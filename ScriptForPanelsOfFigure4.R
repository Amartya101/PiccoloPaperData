#Script for figure 4 panels

####################################################################
#Functions needed 

#Obtain predicted labels
PredictLabels <- function(PiccoloList,KnownLabels = NULL,k=15,p.thres = 0.5){
  
  if(is.null(KnownLabels)){
    stop("Known labels not provided. Please provide character vector with known cell labels.")
  }
  
  PiccoloList$KnownLabels <- KnownLabels
  
  #Find k nearest neighbors
  NearestNeighbors <- function(PiccoloList,k = 10,query = NULL,sort = TRUE,search = "kdtree",bucketSize = 10,splitRule = "suggest",approx = 0){
    PiccoloList$kNN <- dbscan::kNN(x = PiccoloList$PCA$x,k = k,query = query,sort = sort,search = search,bucketSize = bucketSize,splitRule = splitRule,approx = approx)
    return(PiccoloList)
  }
  
  PiccoloList <- NearestNeighbors(PiccoloList = PiccoloList,k = k)
  
  kNN.Mat <- PiccoloList$kNN$id
  #kNN.Mat <- FNN::get.knn(data = PiccoloList$PCA$x,k = k,algorithm = "kd_tree")
  
  #Find cells in each group
  Labels <- names(table(KnownLabels))
  List.Cluster.Enrichments <- vector(mode = "list",length = length(Labels))
  List.Cluster.Cells <- vector(mode = "list",length = length(Labels))
  for (i in 1:length(Labels))
  {
    List.Cluster.Cells[[i]] <- which(KnownLabels == Labels[i])
    n <- nrow(kNN.Mat)
    p.val.vec <- rep(1,n)
    for (j in 1:nrow(kNN.Mat))
    {
      t <- length(intersect(kNN.Mat[j,],List.Cluster.Cells[[i]]))
      a <- length(List.Cluster.Cells[[i]])
      b <- k
      p.val.vec[j] <- sum(dhyper(t:b,a,n-a,b))
    }
    List.Cluster.Enrichments[[i]] <- p.val.vec
  }
  
  #Predict labels for cells
  Matrix.Cluster.Enrichments <- matrix(unlist(List.Cluster.Enrichments),nrow  = length(List.Cluster.Enrichments[[1]]))
  p.threshold <- p.thres
  Predicted.Cell.Labels <- rep("Ambiguous",nrow(Matrix.Cluster.Enrichments))
  for (i in 1:nrow(Matrix.Cluster.Enrichments))
  {
    Min.p.val <- min(Matrix.Cluster.Enrichments[i,])
    if (Min.p.val < p.threshold){
      Predicted.Cell.Labels[i] <- Labels[which(Matrix.Cluster.Enrichments[i,] == Min.p.val)]
    }
  }
  PiccoloList$PredictedLabels <- Predicted.Cell.Labels
  return(PiccoloList)
}

BasicClassifierMetricsPiccolo <- function(PiccoloList){
  #Confusion matrix
  cm <- as.matrix(table(Actual = PiccoloList$KnownLabels,Predicted = PiccoloList$PredictedLabels)) # create the confusion matrix
  
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class
  rowsums = apply(cm,1,sum) # number of instances per class
  colsums = apply(cm,2,sum) # number of predictions per class
  p = rowsums/n # distribution of instances over the actual classes
  q = colsums/n # distribution of instances over the predicted classes
  
  #Accuracy
  accuracy = sum(diag)/n
  accuracy
  
  #Per-class precision, recall, and F1
  precision = diag/colsums
  recall = diag/rowsums
  f1 = 2*precision*recall/(precision+recall)
  prf1.df <- data.frame(precision,recall,f1)
  
  PiccoloList$PRF1 <- prf1.df
  
  ClassWeights <- as.numeric(table(PiccoloList$KnownLabels)/sum(table(PiccoloList$KnownLabels)))
  
  weightedF1 <- sum(ClassWeights*prf1.df$f1)
  
  #Macro averaged metrics
  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1)
  
  #One vs all
  oneVsAll = lapply(1 : nc,
                    function(i){
                      v = c(cm[i,i],
                            rowsums[i] - cm[i,i],
                            colsums[i] - cm[i,i],
                            n-rowsums[i] - colsums[i] + cm[i,i]);
                      return(matrix(v,nrow = 2,byrow = T))})
  oneVsAll
  
  s = matrix(0,nrow = 2,ncol = 2)
  for(i in 1 : nc){s = s + oneVsAll[[i]]}
  s
  
  #Average accuracy
  avgAccuracy = sum(diag(s))/sum(s)
  avgAccuracy
  
  #Micro-averaged metrics -  Compared to unweighted macro-averaging, micro-averaging favors classes with a larger number of instances
  micro_prf = (diag(s)/apply(s,1,sum))[1];
  micro_prf
  
  #Kappa statistic
  expAccuracy = sum(p*q)
  kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
  kappa
  
  #clustComp
  clustComp.df <- aricode::clustComp(PiccoloList$KnownLabels,PiccoloList$PredictedLabels)
  
  PiccoloList$ClassifierMetrics <- data.frame(Accuracy = accuracy,MacroPrecision = macroPrecision,MacroRecall = macroRecall,MacroF1 = macroF1,WeightedF1 = weightedF1,AvgAccuracy = avgAccuracy,MicroPrecision = micro_prf,MicroRecall = micro_prf,MicroF1 = micro_prf,Kappa = kappa,clustComp.df)
  
  return(PiccoloList)
}

FeatureSelectAndNormScran <- function(PiccoloList,NoOfHVG = 3000){
  PiccoloList1 <- PiccoloList
  PiccoloList1$Counts <- Matrix::t(PiccoloList1$Counts)
  
  MinPercNonZeroCells <- 0.5
  
  No.of.Non.Zero.Per.Feature <- Matrix::diff(PiccoloList1$Counts@p)
  
  Perc.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature/nrow(PiccoloList1$Counts) * 100
  
  Irrelevant.Features <- which(Perc.Non.Zero.Per.Feature <= MinPercNonZeroCells)
  
  PiccoloList1$Counts <- Matrix::t(PiccoloList1$Counts[,-Irrelevant.Features])
  
  if (is.null(dim(PiccoloList1$Genes))){
    PiccoloList1$Genes <- PiccoloList1$Genes[-Irrelevant.Features]
  } else {
    PiccoloList1$Genes <- PiccoloList1$Genes$V1[-Irrelevant.Features]
  }
  
  CountsMat <- as.matrix(PiccoloList1$Counts)
  rownames(CountsMat) <- PiccoloList1$Genes
  colnames(CountsMat) <- PiccoloList1$Barcodes
  
  #Scran
  clusters <- scran::quickCluster(CountsMat)
  SizeFactors <- scran::calculateSumFactors(CountsMat, clusters=clusters)
  
  CountsMat <- t(log1p(t(CountsMat)/SizeFactors))
  
  dec <- scran::modelGeneVar(CountsMat)
  
  PiccoloList$SizeFactors <- SizeFactors
  
  # Get the top 3000 genes.
  Scran.hvg <- scran::getTopHVGs(dec, n=NoOfHVG)
  
  Scran.HVG.SerNos <- vector(mode = "numeric",length = length(Scran.hvg))
  for (i in 1:length(Scran.hvg))
  {
    if (is.null(dim(PiccoloList1$Genes)))
    {
      Scran.HVG.SerNos[i] <- which(PiccoloList1$Genes == Scran.hvg[i])
    } else {
      Scran.HVG.SerNos[i] <- which(PiccoloList1$Genes$V1 == Scran.hvg[i])
    }
  }
  
  PiccoloList$NormCounts <- CountsMat[Scran.HVG.SerNos,]
  
  rm(CountsMat)
  
  #Piccolo top 3000 genes
  #PiccoloList1$Counts <- Matrix::t(PiccoloList1$Counts)
  HVG1 <- data.frame(V1 = Scran.hvg)
  PiccoloList$HVG <- HVG1
  
  Scran.HVG.SerNos <- vector(mode = "numeric",length = length(Scran.hvg))
  for (i in 1:length(Scran.hvg))
  {
    if (is.null(dim(PiccoloList$Genes)))
    {
      Scran.HVG.SerNos[i] <- which(PiccoloList$Genes == Scran.hvg[i])
    } else {
      Scran.HVG.SerNos[i] <- which(PiccoloList$Genes$V1 == Scran.hvg[i])
    }
  }
  
  PiccoloList$HVG.Ser.Nos <- Scran.HVG.SerNos
  
  return(PiccoloList)
}


######################################################################
#For NIH/3T3 simulated counts (Scenario 1) data set

#Load the counts list object
PiccoloList <- Piccolo::CreatePiccoloList(MTX = "ZhengMix8eq_matrix.mtx.gz",Genes = "ZhengMix8eq_features.tsv",Barcodes = "ZhengMix8eq_barcodes.tsv")

PiccoloList$CellLabels <- read.delim("ZhengMix8eq_KnownCellLabels.tsv",header = F)$V1

#Select top 3000 HVGs
PiccoloList <- Piccolo::SelectFeatures(PiccoloList = PiccoloList,NoOfHVG = 3000)

#Normalize
PiccoloList <- Piccolo::Normalize(PiccoloList = PiccoloList)

#Compute PC
PiccoloList <- Piccolo::ComputePC(PiccoloList = PiccoloList)

#Find UMAP coordinates
PiccoloList <- Piccolo::UMAPcoords(PiccoloList = PiccoloList)

#Get UMAP plot
Piccolo::LabelUMAP(PiccoloList = PiccoloList,Labels = PiccoloList$CellLabels,Title = "Piccolo - Zheng Mix 8eq",Size  = 3,BaseSize = 28,LegendPosition = "bottom")

#Predict labels based on kNN approach
PiccoloList <- PredictLabels(PiccoloList = PiccoloList,KnownLabels = PiccoloList$Groups,p.thres = 1)

#Get classification metrics
PiccoloList <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList)

PRF1.Piccolo <- PiccoloList$PRF1
ClassMetrics.Piccolo <- PiccoloList$ClassifierMetrics

write.csv(PRF1.Piccolo,file = "PRF1_Piccolo_ZhengMix8eq.csv",row.names = T)
write.csv(ClassMetrics.Piccolo,file = "ClassMetrics_Piccolo_ZhengMix8eq.csv",row.names = T)

#For Analytic Pearson

#Normalize
PiccoloList <- Piccolo::Normalize(PiccoloList = PiccoloList,Transform = "AnalyticPearson")

#Compute PC
PiccoloList <- Piccolo::ComputePC(PiccoloList = PiccoloList)

#Find UMAP coordinates
PiccoloList <- Piccolo::UMAPcoords(PiccoloList = PiccoloList)

#Get UMAP plot
Piccolo::LabelUMAP(PiccoloList = PiccoloList,Labels = PiccoloList$CellLabels,Title = "Analytic Pearson - Zheng Mix 8eq",Size  = 3,BaseSize = 28,LegendPosition = "bottom")

#Predict labels based on kNN approach
PiccoloList <- PredictLabels(PiccoloList = PiccoloList,KnownLabels = PiccoloList$Groups,p.thres = 1)

#Get classification metrics
PiccoloList <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList)

PRF1.APear <- PiccoloList$PRF1
ClassMetrics.APear <- PiccoloList$ClassifierMetrics

write.csv(PRF1.APear,file = "PRF1_Piccolo_AnalyticPearson_ZhengMix8eq.csv",row.names = T)
write.csv(ClassMetrics.APear,file = "ClassMetrics_Piccolo_AnalyticPearson_ZhengMix8eq.csv",row.names = T)

#Use Seurat - SCTransform v2
rownames(PiccoloList$Counts) <- PiccoloList$Genes$V1
colnames(PiccoloList$Counts) <- PiccoloList$Barcodes

SeuratObj <- Seurat::CreateSeuratObject(counts = PiccoloList$Counts)

SeuratObj <- Seurat::SCTransform(SeuratObj,vst.flavor = "v2")

#Perform normalization using SCTransform v2
SCTransformVariableFeatures <- SeuratObj@assays$SCT@var.features

#Replace scaled matrix in PiccoloList
PiccoloList$NormCounts <- SeuratObj@assays$SCT@scale.data
HVG1 <- data.frame(V1 = SCTransformVariableFeatures)
PiccoloList$HVG <- HVG1

HVG1.Ser.Nos <- vector(mode = "numeric",length = length(PiccoloList$HVG$V1))
for (i in 1:length(PiccoloList$HVG$V1)){
  HVG1.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$HVG$V1[i])
}
PiccoloList$HVG.Ser.Nos <- HVG1.Ser.Nos

#Compute PC
PiccoloList <- Piccolo::ComputePC(PiccoloList = PiccoloList)

#Get UMAP coordinates
PiccoloList <- Piccolo::UMAPcoords(PiccoloList = PiccoloList)

Piccolo::LabelUMAP(PiccoloList = PiccoloList,Labels = PiccoloList$CellLabels,Title = "SCTransform v2 - Zheng Mix 8eq", Size = 3,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList <- PredictLabels(PiccoloList = PiccoloList,KnownLabels = PiccoloList$Groups,p.thres = 1)

#Get classification metrics
PiccoloList <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList)

PRF1.SCT <- PiccoloList$PRF1
ClassMetrics.SCT <- PiccoloList$ClassifierMetrics

write.csv(PRF1.SCT,file = "PRF1_SCTv2_ZhengMix8eq.csv",row.names = T)
write.csv(ClassMetrics.SCT,file = "ClassMetrics_SCTv2_ZhengMix8eq.csv",row.names = T)

#For Scran
PiccoloList <- FeatureSelectAndNormScran(PiccoloList = PiccoloList,NoOfHVG = 3000)

PiccoloList <- Piccolo::ComputePC(PiccoloList = PiccoloList)

PiccoloList <- Piccolo::UMAPcoords(PiccoloList = PiccoloList)

Piccolo::LabelUMAP(PiccoloList = PiccoloList,Labels = PiccoloList$Groups,Title = "Scran logSF - Simulated Counts Scenario 1", Size = 2,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList <- PredictLabels(PiccoloList = PiccoloList,KnownLabels = PiccoloList$Groups,p.thres = 1)

#Get classification metrics
PiccoloList <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList)

PRF1.Scran <- PiccoloList$PRF1
ClassMetrics.Scran <- PiccoloList$ClassifierMetrics

write.csv(PRF1.Scran,file = "PRF1_Scran_ZhengMix8eq.csv",row.names = T)
write.csv(ClassMetrics.Scran,file = "ClassMetrics_Scran_ZhengMix8eq.csv",row.names = T)

