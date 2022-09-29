# Curve Clustering - 4 subclasses : 
# Allows to classify a CHCs base in PP/PV/ECM/STEM according to the Curve Clustering method

# The function takes as input :
#    - Data : Database (Format: Genes x individuals)
#    - use.CCR : If it is already calculated, the user can fill in a CCR file 
#    - couleurs : Vector containing 4 colors

# The function returns as output :
#    - $Class : A dataframe containing the identifiers and the class obtained by Curve Clustering
#    - $CCR : The CCR file to be filled in the next times to save calculation time

CurvClust_4_subclasses <- function(Data, use.CCR = NULL, couleurs = NULL) {
  
  library("FactoMineR")
  library("factoextra")
  library(curvclust)
  
  t.Data <- as.data.frame(t(Data))
  
  ### Data distribution
  
  res <- PCA(t(Data), graph = F)
  plot(fviz_pca_ind(res, col.ind = "#00BFC4"))
  
  print("Data distribution : completed")
  
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
      hist(SNR[SNR < round(summary(SNR)[[4]])], col = "lightblue", main = "SNR", xlim = c(0,10))
      abline(v = seuil, col = "red")
    }
    
    keep <- which(SNR >= seuil)
    data_keep <- Base[c(keep),]
    
    print(paste0("Reduced database dimensions : Nrow = ", nrow(data_keep), " x Ncol = ", ncol(data_keep)))
    
    return(data_keep)
    
  }
  
  r.Data <- SNR(Data, 0.5, plot = T) # 15016 x 61 CHCs
  
  print("SNR: completed")
  
  # Variance test for data reduction ----
  
  varianceReduction <- function(Base, seuil) {
    
    print("Database format: Entities in rows, individuals in columns")
    
    median.var <- median(apply(t(Base),2,var))
    Var.accross.samples <- (dim(Base)[2]-1)*apply(t(Base),2,var)/median.var
    
    pval.var.test <- 1 - pchisq(Var.accross.samples, df = dim(Base)[2]-1)
    gene.set <- which(pval.var.test <= seuil)
    
    print(paste0("Number of keeped genes: ", length(gene.set), " (meaning ", round(length(gene.set)*100/nrow(Base)), " %)"))
    gene.set <- as.data.frame(gene.set)
    
    Base.reduite <- Base[rownames(Base) %in% rownames(gene.set),]
    
    print(paste0("Reduced database dimensions : Nrow = ", nrow(Base.reduite), " x Ncol = ", ncol(Base.reduite)))
    
    return(Base.reduite)
    
  }
  
  Data.Pval <- varianceReduction(r.Data, 0.01)
  
  print("Variance test: completed")
  
  # CurveClustering ----
  
  Data <- Data.Pval # 4202 x 1133
  
  Data = lapply(Data, function(x) {
    if(any(is.infinite(x))) {
      x[is.infinite(x)] = 0
    }
    return(x)
  })
  Data = as.data.frame(Data)
  
  print("Data transformation: completed")
  
  Y <- list()
  for (i in 1:ncol(Data)) {
    Y[[i]] <- Data[,i]
  }
  
  n = ncol(Data)
  M = nrow(Data)
  
  CCD <- new("CClustData",Y=Y, filter.number = 1)
  
  CDred <- getUnionCoef(CCD)
  
  CCO <- new("CClustO")
  CCO["nbclust"] = 4
  
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

  print("Clutering: completed")

  library(kohonen)
  library(pls)
  library(ChemometricsWithR)
  library(colorspace)
  library(ggplot2)
  library(factoextra)
  library(FactoMineR)

  if (length(couleurs) != 4) {couleurs<-hex(HLS(H=seq(1,1000,100), L=0.5, S=0.6))[c(1,2,4,3)]}
  
  col.graph <- couleurs

  Base.PLS<-data.frame(Clust=I(classvec2classmat(Clust)), Data=I(as.matrix(t(Data))))
  detach("package:ChemometricsWithR", unload=TRUE)

  #cppls
  PLS.clust <- cppls(Clust ~ Data, ncomp = 20, data = Base.PLS, validation = "none", scale = T)

  ## Plot 2D

  par(mar = c(5, 5, 2, 2))
  scoreplot(PLS.clust, comps = 1:2, identify = FALSE, type = "n", cex.axis = 1.5, cex.lab = 1.7)
  points(PLS.clust$score[,1][Base.PLS$Clust[,1]==1], PLS.clust$score[,2][Base.PLS$Clust[,1]==1], pch=20, cex=1.2, col=col.graph[1])
  points(PLS.clust$score[,1][Base.PLS$Clust[,2]==1], PLS.clust$score[,2][Base.PLS$Clust[,2]==1], pch=20, cex=1.2, col=col.graph[2])
  points(PLS.clust$score[,1][Base.PLS$Clust[,3]==1], PLS.clust$score[,2][Base.PLS$Clust[,3]==1], pch=20, cex=1.2, col=col.graph[3])
  points(PLS.clust$score[,1][Base.PLS$Clust[,4]==1], PLS.clust$score[,2][Base.PLS$Clust[,4]==1], pch=20, cex=1.2, col=col.graph[4])

  legend("bottomright", legend = unique(Clust),
          pch = rep(20,4), pt.cex=rep(2.1,4), col = c(col.graph[1],col.graph[2], col.graph[3], col.graph[4]), cex = 1.5)
  
  print("Done !")
  
  res <- list("Class" = classe, "CCR" = CCR)
  return(res)
  
}   
