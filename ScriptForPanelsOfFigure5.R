#Script for NIH/3T3 simulated counts

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
PiccoloList.Splatter <- Piccolo::CreatePiccoloList(MTX = "NIH3T3_SimCountsScenario1.mtx.gz",Genes = "NIH3T3_SimCountsScenario1_features.tsv",Barcodes = "NIH3T3_SimCountsScenario1_barcodes.tsv")

PiccoloList.Splatter$Groups <- read.delim("NIH3T3_SimCountsScenario1_groups.tsv",header = F)$V1

#Select top 3000 HVGs
PiccoloList.Splatter <- Piccolo::SelectFeatures(PiccoloList = PiccoloList.Splatter,NoOfHVG = 3000)

#Normalize
PiccoloList.Splatter <- Piccolo::Normalize(PiccoloList = PiccoloList.Splatter)

#Compute PC
PiccoloList.Splatter <- Piccolo::ComputePC(PiccoloList = PiccoloList.Splatter)

#Find UMAP coordinates
PiccoloList.Splatter <- Piccolo::UMAPcoords(PiccoloList = PiccoloList.Splatter)

#Get UMAP plot
Piccolo::LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "Piccolo - Simulated Counts Scenario 1",Size  = 3,BaseSize = 28,LegendPosition = "bottom")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.Piccolo <- PiccoloList.Splatter$PRF1
ClassMetrics.Piccolo <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.Piccolo,file = "PRF1_Piccolo_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.Piccolo,file = "ClassMetrics_Piccolo_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

# write.csv(PRF1.Piccolo,file = "PRF1_Piccolo_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)
# write.csv(ClassMetrics.Piccolo,file = "ClassMetrics_Piccolo_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)

#For Analytic Pearson

#Normalize
PiccoloList.Splatter <- Piccolo::Normalize(PiccoloList = PiccoloList.Splatter,Transform = "AnalyticPearson")

#Compute PC
PiccoloList.Splatter <- Piccolo::ComputePC(PiccoloList = PiccoloList.Splatter)

#Find UMAP coordinates
PiccoloList.Splatter <- Piccolo::UMAPcoords(PiccoloList = PiccoloList.Splatter)

#Get UMAP plot
Piccolo::LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "Analytic Pearson - Simulated Counts Scenario 1",Size  = 3,BaseSize = 28,LegendPosition = "bottom")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.APear <- PiccoloList.Splatter$PRF1
ClassMetrics.APear <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.APear,file = "PRF1_Piccolo_AnalyticPearson_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.APear,file = "ClassMetrics_Piccolo_AnalyticPearson_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

# write.csv(PRF1.APear,file = "PRF1_Piccolo_AnalyticPearson_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)
# write.csv(ClassMetrics.APear,file = "ClassMetrics_Piccolo_AnalyticPearson_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)

#Use Seurat - SCTransform v2
rownames(PiccoloList.Splatter$Counts) <- PiccoloList.Splatter$Genes$V1
colnames(PiccoloList.Splatter$Counts) <- PiccoloList.Splatter$Barcodes

SeuratObj <- Seurat::CreateSeuratObject(counts = PiccoloList.Splatter$Counts)

SeuratObj <- Seurat::SCTransform(SeuratObj,vst.flavor = "v2")

#Perform normalization using SCTransform v2
SCTransformVariableFeatures <- SeuratObj@assays$SCT@var.features

#Replace scaled matrix in PiccoloList
PiccoloList.Splatter$NormCounts <- SeuratObj@assays$SCT@scale.data
HVG1 <- data.frame(V1 = SCTransformVariableFeatures)
PiccoloList.Splatter$HVG <- HVG1

HVG1.Ser.Nos <- vector(mode = "numeric",length = length(PiccoloList.Splatter$HVG$V1))
for (i in 1:length(PiccoloList.Splatter$HVG$V1)){
  HVG1.Ser.Nos[i] <- which(PiccoloList.Splatter$Genes == PiccoloList.Splatter$HVG$V1[i])
}
PiccoloList.Splatter$HVG.Ser.Nos <- HVG1.Ser.Nos

#Compute PC
PiccoloList.Splatter <- Piccolo::ComputePC(PiccoloList = PiccoloList.Splatter)

#Get UMAP coordinates
PiccoloList.Splatter <- Piccolo::UMAPcoords(PiccoloList = PiccoloList.Splatter)

Piccolo::LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "SCTransform v2 - Simulated Counts Scenario 1", Size = 2,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.SCT <- PiccoloList.Splatter$PRF1
ClassMetrics.SCT <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.SCT,file = "PRF1_SCTv2_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.SCT,file = "ClassMetrics_SCTv2_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

# write.csv(PRF1.SCT,file = "PRF1_SCTv2_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)
# write.csv(ClassMetrics.SCT,file = "ClassMetrics_SCTv2_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)

#Piccolo with SCT features
PiccoloList.Splatter <- Piccolo::Normalize(PiccoloList = PiccoloList.Splatter)

PiccoloList.Splatter <- Piccolo::ComputePC(PiccoloList = PiccoloList.Splatter)

PiccoloList.Splatter <- Piccolo::UMAPcoords(PiccoloList = PiccoloList.Splatter)

Piccolo::LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "Piccolo (SCTransform v2) - Sim. Counts Scenario 1", Size = 2,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.PiccoloSCT <- PiccoloList.Splatter$PRF1
ClassMetrics.PiccoloSCT <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.PiccoloSCT,file = "PRF1_PiccoloSCT_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.PiccoloSCT,file = "ClassMetrics_PiccoloSCT_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

#write.csv(PRF1.PiccoloSCT,file = "PRF1_PiccoloSCT_SimCountsNIH3T3_Scenario2_6GroupsSameSize_Scenario2.csv",row.names = T)
#write.csv(ClassMetrics.PiccoloSCT,file = "ClassMetrics_PiccoloSCT_SimCountsNIH3T3_Scenario2_6GroupsSameSize_Scenario2.csv",row.names = T)

#For Scran
PiccoloList.Splatter <- FeatureSelectAndNormScran(PiccoloList = PiccoloList.Splatter,NoOfHVG = 3000)

PiccoloList.Splatter <- Piccolo::ComputePC(PiccoloList = PiccoloList.Splatter)

PiccoloList.Splatter <- Piccolo::UMAPcoords(PiccoloList = PiccoloList.Splatter)

Piccolo::LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "Scran logSF - Simulated Counts Scenario 1", Size = 2,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.Scran <- PiccoloList.Splatter$PRF1
ClassMetrics.Scran <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.Scran,file = "PRF1_Scran_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.Scran,file = "ClassMetrics_Scran_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

# write.csv(PRF1.Scran,file = "PRF1_Scran_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)
# write.csv(ClassMetrics.Scran,file = "ClassMetrics_Scran_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)

#Piccolo with Scran features
PiccoloList.Splatter <- Normalize(PiccoloList = PiccoloList.Splatter)

PiccoloList.Splatter <- ComputePC(PiccoloList = PiccoloList.Splatter)

PiccoloList.Splatter <- UMAPcoords(PiccoloList = PiccoloList.Splatter)

LabelUMAP(PiccoloList = PiccoloList.Splatter,Labels = PiccoloList.Splatter$Groups,Title = "Piccolo (Scran logSF) - Sim. Counts Scenario 2", Size = 2,BaseSize = 28,LegendPosition = "none")

#Predict labels based on kNN approach
PiccoloList.Splatter <- PredictLabels(PiccoloList = PiccoloList.Splatter,KnownLabels = PiccoloList.Splatter$Groups,p.thres = 1)

#Get classification metrics
PiccoloList.Splatter <- BasicClassifierMetricsPiccolo(PiccoloList = PiccoloList.Splatter)

PRF1.PiccoloScran <- PiccoloList.Splatter$PRF1
ClassMetrics.PiccoloScran <- PiccoloList.Splatter$ClassifierMetrics

write.csv(PRF1.PiccoloScran,file = "PRF1_PiccoloScran_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)
write.csv(ClassMetrics.PiccoloScran,file = "ClassMetrics_PiccoloScran_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv",row.names = T)

# write.csv(PRF1.PiccoloScran,file = "PRF1_PiccoloScran_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)
# write.csv(ClassMetrics.PiccoloScran,file = "ClassMetrics_PiccoloScran_SimCountsNIH3T3_6GroupsSameSize_Scenario2.csv",row.names = T)

#Prepare heatmap

#Import result files
PRF1.Piccolo <- data.table::fread("PRF1_Piccolo_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv")
PRF1.APear <- data.table::fread("PRF1_Piccolo_AnalyticPearson_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv")
PRF1.SCT <- data.table::fread("PRF1_SCTv2_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv")
PRF1.Scran <- data.table::fread("PRF1_Scran_SimCountsNIH3T3_6GroupsDiffSize_Scenario1.csv")

F1.Vec <- c(PRF1.Piccolo$f1,PRF1.APear$f1,PRF1.Scran$f1,PRF1.SCT$f1)

Group.Vec <- paste0(substr(PRF1.Piccolo$V1,1,nchar(PRF1.Piccolo$V1)-1)," ",1:length(PRF1.Piccolo$V1))

Method.Vec <- rep(c("Piccolo","Analytic Pearson","Scran logSF","SCTransform v2"),each = 6)

df <- data.frame(Method = Method.Vec,Group = Group.Vec,F1 = F1.Vec)

df$Method <- factor(df$Method, levels=c("SCTransform v2","Scran logSF","Analytic Pearson","Piccolo"))

df1 <- df[df$F1 > 0.75,]
df2 <- df[df$F1 <= 0.75,]

ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Method, fill = F1)) +
  ggplot2::geom_tile() + viridis::scale_fill_viridis() +
  ggplot2::coord_fixed() +
  ggplot2::geom_text(data = df1,ggplot2::aes(label = signif(F1,2)), color = "black", size = 7) +
  ggplot2::geom_text(data = df2,ggplot2::aes(label = signif(F1,2)), color = "white", size = 7) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,hjust = 1)) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::ggtitle("F1 Scores Per Group - Scenario 1") +
  ggplot2::theme(text = ggplot2::element_text(size = 28))

