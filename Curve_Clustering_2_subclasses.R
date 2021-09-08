# Curve Clustering - 2 subclasses:
# Permet de classer une base de CHCs periveineux en CHCs GOOD ou BAD d'apr�s la m�thode de Curve Clustering

# La fonction prend en entr�e :
#     - Data : Base de donn�es (Format : G�nes x individus)
#     - use.CCR : Si il est d�j� calcul�, l'utilisateur peut renseigner un objet CCR

# La fonction renvoie en sortie :
#     - $Class : Un dataframe contenant les identifiants et la classe obtenue par CurveClustering
#     - $CCR : Un objet CCR � renseigner les prochaines fois pour gagner du temps de calcul

# A titre d'exemple
  # load("D:/home/user/Documents/BDD/Base_ICGC/10-19_Analyse sous groupe MUT vs WT/TCGA_MUT_Survie.RData") # 51923 x 61 CHCs
  # load("D:/home/user/Documents/BDD/11.19 Mise en place Base composite/TCGA Combat N.RData") # TCGA.Combat : 19079 x 370 CHCs
  # Data <- TCGA.Combat[,colnames(TCGA.Combat) %in% colnames(TCGA_Norm_MUT)] # 19079 x 61 CHCs
  # 
  # load("D:/home/user/Documents/BDD/11.19 Reduction de dimensions/CCR_2_Application_TCGA.RData")
  # use.CCR <- CCR

CurvClust_2_subclasses <- function(Data, use.CCR = NULL) {
  
  library("FactoMineR")
  library("factoextra")
  library(curvclust)
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  t.Data <- as.data.frame(t(Data))
  
  # Signal-to-Noise Ratio ---------------------------------------------------
  
  SNR <- function(Base, seuil, plot = T) {
    
    print("Database format: Entities in rows, individuals in columns")
    
    mean.c <- apply(t(Base),2,mean)
    var.c <- apply(t(Base),2,var)
    
    n <- ncol(Base)
    
    v <-  (n * var.c) / (n)
    v[v==0] <- 1
    
    SNR <- abs(mean.c) / sqrt(v)
    print("Summary SNR :")
    print(summary(SNR))
    
    if (plot == TRUE) {
      hist(SNR[SNR < round(summary(SNR)[[4]])], col = "lightblue", main = "SNR")
      abline(v = seuil, col = "red")
    }
    
    keep <- which(SNR >= seuil)
    data_keep <- Base[c(keep),]
    
    print(paste0("Reduced database dimensions : Nrow = ", nrow(data_keep), " x Ncol = ", ncol(data_keep)))
    
    return(data_keep)
    
  }
  
  r.Data <- SNR(Data, 0, plot = T)
  
  print("SNR: completed")
  
  # Variance test for data reduction ----
  
  varianceReduction <- function(Base, seuil) {
    
    print("Database format: Entities in rows, individuals in columns")
    
    median.var <- median(apply(t(Base),2,var))
    Var.accross.samples <- (dim(Base)[2]-1)*apply(t(Base),2,var)/median.var
    
    pval.var.test <- 1 - pchisq(Var.accross.samples, df = dim(Base)[2]-1)
    gene.set <- which(pval.var.test <= seuil)
    
    print(paste0("Number of keeped genes : ", length(gene.set), " (meaning ", round(length(gene.set)*100/nrow(Base)), " %)"))
    gene.set <- as.data.frame(gene.set)
    
    Base.reduite <- Base[rownames(Base) %in% rownames(gene.set),]
    
    print(paste0("Reduced database dimensions : Nrow = ", nrow(Base.reduite), " x Ncol = ", ncol(Base.reduite)))
    
    return(Base.reduite)
    
  }
  
  Data.Pval <- varianceReduction(r.Data, 0.01)
  
  print("Variance test: completed")

  # Graphique d'observation -----
  
  Data <- Data.Pval # 6813 x 61
  
  Data <- log2(Data)
  Data = lapply(Data, function(x) {
    if(any(is.infinite(x))) {
      x[is.infinite(x)] = 0
    }
    return(x)
  })
  Data = as.data.frame(Data)
  
  graphics.off()
  plot(1:nrow(Data),type="n",ylab="Expression (en log)", xaxt = "n", xlab = "", ylim = c(1,max(Data)))
  mtext(1, text = "G�nes", line = 1)
  
  for (i in 1:ncol(Data)) {
    lines(1:nrow(Data),Data[,i], col=c(i))
  }
  
  print("Data transformation: completed")
  
  # CurveClustering ----
  
  Y <- list()
  for (i in 1:ncol(Data)) {
    Y[[i]] <- Data[,i]
  }
  
  n = ncol(Data)
  M = nrow(Data)
  
  CCD <- new("CClustData",Y=Y, filter.number = 1)
  
  CDred <- getUnionCoef(CCD) # dimension reduction
  
  ## options setting
  CCO <- new("CClustO")
  CCO["nbclust"] = 2
  
  ## computing Functional Clustering Mixed Model
  set.seed(111)
  
  if (length(use.CCR) == 0) {CCR <- getFCMM(CDred,CCO)}
  else {CCR <- use.CCR}
  
  print("CCR: completed")
  
  groups <- apply(CCR["Tau"],1,which.max)
  groups <- as.data.frame(groups)
  
  t.Data$groups <- groups$groups
  groups$X <- rownames(t.Data)
  
  Clust<-as.factor(t.Data$groups)
  table(Clust)
  
  classe <- matrix(data = NA, nrow = ncol(Data), ncol = 2)
  classe <- as.data.frame(classe)
  colnames(classe) <- c("id","classe")
  classe$id <- rownames(t.Data)
  classe$classe <- groups$groups
  
  print("Clustering: completed")
  print("Done !")
  
  res <- list("Class" = classe, "CCR" = CCR)
  return(res)
  
}

# Exemple de commande pour lancer la fonction :
  # res <- CurvClust_2_subclasses(Data, use.CCR)