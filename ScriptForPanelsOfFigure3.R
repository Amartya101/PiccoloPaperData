#Script to prepare figures in Fig. 3 in the Piccolo paper

#Install Piccolo

install.packages("devtools")

devtools::install_github("Amartya101/Piccolo")

########################################################################################
#Functions needed to make the plots

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


PCplot <- function(PiccoloList,NoOfHVG = 3000,DataName,k = 5,Contours = F,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 28,LegendPosition = "right"){
  PiccoloList1 <- PiccoloList
 
  PiccoloList1 <- Piccolo::SelectFeatures(PiccoloList = PiccoloList1,NoOfHVG = 3000)
  
  PiccoloList1$NormCounts <- as.matrix(PiccoloList1$Counts[PiccoloList1$HVG.Ser.Nos,])
  
  PiccoloList1 <- Piccolo::ComputePC(PiccoloList = PiccoloList1)
  
  PiccoloList1$SizeFactors <- log(colSums(PiccoloList1$NormCounts)/mean(colSums(PiccoloList1$NormCounts)))
  
  CorrCoef <- stats::cancor(PiccoloList1$SizeFactors,PiccoloList1$PCA$x[,1:k])$cor
  
  DataName <- paste0(DataName," (",expression(rho)," = ",signif(CorrCoef,digits = 2),")")
  
  df <- data.frame(PC1 = PiccoloList1$PCA$x[,1],PC2 = PiccoloList1$PCA$x[,2],PiccoloList1$SizeFactors)
  colnames(df) <- c("PC1","PC2","log(Size Factors)")
  if (xLabel == T){
    xlabeltext <- "PC1"
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- "PC2"
  } else {
    ylabeltext <- ""
  }
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = PC1,y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(PC1, PC2, color = `log(Size Factors)`),size = 2) + viridis::scale_color_viridis() +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(DataName) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position=LegendPosition) +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
}

PCNormPlot <- function(PiccoloList,Transform = "log",NoOfHVG = 3000,DataName,k = 5,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 28,LegendPosition = "right"){
  PiccoloList1 <- PiccoloList
 
  PiccoloList1 <- Piccolo::SelectFeatures(PiccoloList = PiccoloList1,NoOfHVG = NoOfHVG)
  
  if (Transform == "SCTv2"){
    #Use Seurat - SCTransform v2
    rownames(PiccoloList1$Counts) <- PiccoloList1$Genes
    colnames(PiccoloList1$Counts) <- PiccoloList1$Barcodes
    
    SeuratObj <- Seurat::CreateSeuratObject(counts = PiccoloList1$Counts)
    
    SeuratObj <- Seurat::SCTransform(SeuratObj,vst.flavor = "v2")
    
    #Perform normalization using SCTransform v2
    SCTransformVariableFeatures <- SeuratObj@assays$SCT@var.features
    
    #Replace scaled matrix in PiccoloList
    PiccoloList1$NormCounts <- SeuratObj@assays$SCT@scale.data
    HVG1 <- data.frame(V1 = SCTransformVariableFeatures)
    PiccoloList1$HVG <- HVG1
    
    HVG1.Ser.Nos <- vector(mode = "numeric",length = length(PiccoloList1$HVG$V1))
    for (i in 1:length(PiccoloList1$HVG$V1)){
      HVG1.Ser.Nos[i] <- which(PiccoloList1$Genes == PiccoloList1$HVG$V1[i])
    }
    PiccoloList1$HVG.Ser.Nos <- HVG1.Ser.Nos
    PiccoloList1$SizeFactors <- Matrix::colSums(PiccoloList1$Counts)/mean(Matrix::colSums(PiccoloList1$Counts))
  } else if (Transform == "ScranLogSF"){
    #For Scran
    PiccoloList1 <- FeatureSelectAndNormScran(PiccoloList = PiccoloList1,NoOfHVG = NoOfHVG)
    PiccoloList1$SizeFactors <- Matrix::colSums(PiccoloList1$Counts)/mean(Matrix::colSums(PiccoloList1$Counts))
  } else {
    PiccoloList1 <- Piccolo::Normalize(PiccoloList = PiccoloList1,Transform = Transform)
    PiccoloList1$SizeFactors <- Matrix::colSums(PiccoloList1$Counts)/mean(Matrix::colSums(PiccoloList1$Counts))
  }
  
  PiccoloList1 <- Piccolo::ComputePC(PiccoloList = PiccoloList1)
  
  CorrCoef <- stats::cancor(PiccoloList1$SizeFactors,PiccoloList1$PCA$x[,1:k])$cor
  
  DataName <- paste0(DataName," - Normalized"," (",expression(rho)," = ",signif(CorrCoef,digits = 2),")")
  
  #DataName <- paste0(DataName," - Piccolo Normalized")
  
  df <- data.frame(PC1 = PiccoloList1$PCA$x[,1],PC2 = PiccoloList1$PCA$x[,2],SizeFactor = log(PiccoloList1$SizeFactors))
  colnames(df) <- c("PC1","PC2","log(Size Factors)")
  
  if (xLabel == T){
    xlabeltext <- "PC1"
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- "PC2"
  } else {
    ylabeltext <- ""
  }
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = PC1,y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(PC1, PC2, color = `log(Size Factors)`),size = 2) + viridis::scale_color_viridis() +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(DataName) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position=LegendPosition) + 
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
}

ResidualVarianceMonoPlot <- function(PiccoloList,DataName,Transform = NULL,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 24){
  PiccoloList1 <- PiccoloList
  
  PiccoloList1 <- Piccolo::SelectFeatures(PiccoloList = PiccoloList1,NoOfHVG = 3000)
  
  if(is.null(Transform)){
    Transform <- "log"
  }
  
  PiccoloList1 <- Piccolo::Normalize(PiccoloList = PiccoloList1,Transform = Transform)
  
  PiccoloList1$RowMeans <- Matrix::rowMeans(PiccoloList1$Counts[PiccoloList1$HVG.Ser.Nos,])
  PiccoloList1$RowVar <- matrixStats::rowVars(PiccoloList1$NormCounts)
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- "Variance"
  } else {
    ylabeltext <- ""
  }
  
  if (Transform == "log"){
    Transform <- "Piccolo"
  }
  
  FinalTitle <- paste0(Transform," - ",DataName)
  #FinalTitle <- paste0(DataName," - ","Analytic Pearson")
  
  df <- data.frame(Mean = PiccoloList1$RowMeans,Var = PiccoloList1$RowVar)
  
  ExceedingVar <- which(df$Var > yLim[2])
  if (length(ExceedingVar) != 0){
    df$Var[ExceedingVar] <- yLim[2]
  }
  
  colnames(df) <- c("Mean","Residual Variance")
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = `Residual Variance`)) +
    ggplot2::geom_point(color = "#440154FF",ggplot2::aes(Mean, `Residual Variance`)) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10",scales::math_format(10^.x))) +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(FinalTitle) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position="none")
}

ResidualVarianceMonoPlotSCT <- function(PiccoloList,DataName,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 28){
  
 
  if (is.null(dim(PiccoloList$Genes))){
    rownames(PiccoloList$Counts) <- PiccoloList$Genes
  } else {
    rownames(PiccoloList$Counts) <- PiccoloList$Genes$V1
  }
  
  colnames(PiccoloList$Counts) <- PiccoloList$Barcodes
  
  SeuratObj <- Seurat::CreateSeuratObject(counts = PiccoloList$Counts)
  
  SeuratObj <- Seurat::SCTransform(SeuratObj,vst.flavor = "v2")
  
  SCTransformVariableFeatures <- SeuratObj@assays$SCT@var.features
  
  #Replace scaled matrix in PiccoloList
  PiccoloList$NormCounts <- SeuratObj@assays$SCT@scale.data
  HVG1 <- data.frame(V1 = SCTransformVariableFeatures)
  PiccoloList$HVG <- HVG1
  
  HVG1.Ser.Nos <- vector(mode = "numeric",length = length(PiccoloList$HVG$V1))
  for (i in 1:length(PiccoloList$HVG$V1))
  {
    if (is.null(dim(PiccoloList$Genes))){
      HVG1.Ser.Nos[i] <- which(gsub("_","-",PiccoloList$Genes,fixed = T)  == PiccoloList$HVG$V1[i])
    } else {
      HVG1.Ser.Nos[i] <- which(gsub("_","-",PiccoloList$Genes$V1,fixed = T) == PiccoloList$HVG$V1[i])
    }
  }
  
  PiccoloList$HVG.Ser.Nos <- HVG1.Ser.Nos
  
  PiccoloList$RowMeans <- Matrix::rowMeans(PiccoloList$Counts[PiccoloList$HVG.Ser.Nos,])
  PiccoloList$RowVar <- matrixStats::rowVars(PiccoloList$NormCounts)
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- "Variance"
  } else {
    ylabeltext <- ""
  }
  
  FinalTitle <- paste0("SCTransform v2"," - ",DataName)
  
  df <- data.frame(Mean = PiccoloList$RowMeans,Var = PiccoloList$RowVar)
  
  ExceedingVar <- which(df$Var > yLim[2])
  if (length(ExceedingVar) != 0){
    df$Var[ExceedingVar] <- yLim[2]
  }
  
  colnames(df) <- c("Mean","Residual Variance")
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = `Residual Variance`)) +
    ggplot2::geom_point(color = "#440154FF",ggplot2::aes(Mean, `Residual Variance`)) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10",scales::math_format(10^.x))) +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(FinalTitle) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position="none")
}

ResidualVarianceMonoPlotScran <- function(PiccoloList,DataName,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 28){
  
  #With Scran
  
  PiccoloList <- FeatureSelectAndNormScran(PiccoloList = PiccoloList,NoOfHVG = 3000)
  
  
  PiccoloList$RowMeans <- Matrix::rowMeans(PiccoloList$Counts[PiccoloList$HVG.Ser.Nos,])
  PiccoloList$RowVar <- matrixStats::rowVars(PiccoloList$NormCounts)
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- "Variance"
  } else {
    ylabeltext <- ""
  }
  
  FinalTitle <- paste0("Scran logSF"," - ",DataName)
  
  df <- data.frame(Mean = PiccoloList$RowMeans,Var = PiccoloList$RowVar)
  
  ExceedingVar <- which(df$Var > yLim[2])
  if (length(ExceedingVar) != 0){
    df$Var[ExceedingVar] <- yLim[2]
  }
  
  colnames(df) <- c("Mean","Residual Variance")
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = `Residual Variance`)) +
    ggplot2::geom_point(color = "#440154FF",ggplot2::aes(Mean, `Residual Variance`)) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10",scales::math_format(10^.x))) +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(FinalTitle) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position="none")
}

################################################################
#PC plots for Svensson 1

PiccoloList <- Piccolo::CreatePiccoloList(MTX = "Svensson1_matrix.mtx.gz",Genes = "Svensson1_features.tsv",Barcodes = "Svensson1_barcodes.tsv")

PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 1,MT.Perc = 100,RP.Perc = 100,TotalCountsMADHigh = 3.5,TotalCountsMADLow = 3.5)

PCplot(PiccoloList = PiccoloList,DataName = "Svensson 1")

PCNormPlot(PiccoloList = PiccoloList,DataName = "Svensson 1",Transform = "ScranLogSF")

#PC plots for PBMC 33k

PiccoloList <- Piccolo::CreatePiccoloList(MTX = "PBMC33k.mtx.gz",Genes = "PBMC33k_features.tsv",Barcodes = "PBMC33k_barcodes.tsv")

PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 1,MT.Perc = 100,RP.Perc = 100)

PCplot(PiccoloList = PiccoloList,DataName = "PBMC 33k")

PCNormPlot(PiccoloList = PiccoloList,DataName = "PBMC 33K")

#Residual variance vs mean expression plots

PiccoloList <- Piccolo::CreatePiccoloList(MTX = "PBMC33k.mtx.gz",Genes = "PBMC33k_features.tsv",Barcodes = "PBMC33k_barcodes.tsv")

PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 1,MT.Perc = 100,RP.Perc = 100)

PiccoloList$CountsOriginal <- NULL

#For Piccolo
ResidualVarianceMonoPlot(PiccoloList = PiccoloList,DataName = "PBMC 33k",yLim = c(0,5))

#For Analytic Pearson
ResidualVarianceMonoPlot(PiccoloList = PiccoloList,DataName = "PBMC 33k",Transform = "AnalyticPearson",yLim = c(0,5))

#For SCTransform v2
ResidualVarianceMonoPlotSCT(PiccoloList = PiccoloList,DataName = "PBMC 33k",yLim = c(0,5))

#For Scran logSF
ResidualVarianceMonoPlotScran(PiccoloList = PiccoloList,DataName = "PBMC 33k",yLim = c(0,5))



