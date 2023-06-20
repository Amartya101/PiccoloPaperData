#Script to prepare figures in Fig. 1 (and Supplementary figures 1-4) in the Piccolo paper

#Install Piccolo

install.packages("devtools")

devtools::install_github("Amartya101/Piccolo")

########################################################################################

#Install package MASS

install.packages("MASS")

#Function for plotting alphaQP vs mean log-log plots
AlphaQPMeanLogLogPlot <- function(PiccoloList,DataName,Bins = F,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 32,AxisLabelSize = 32){
  UMI.Mat <- Matrix::t(PiccoloList$Counts)
  
  Gene.IDs <- PiccoloList$Genes
  
  colVarsSPM <- function(X) {
    stopifnot( methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]),function(j) {
      if(X@p[j+1] == X@p[j]) { return(0) } # all entries are 0: var is 0
      mean <- base::sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
      sum((X@x[(X@p[j]+1):X@p[j+1] ] - mean)^2) +
        mean^2 * (X@Dim[1] - (X@p[j+1] - X@p[j]))})/(X@Dim[1] - 1)
    names(ans) <- X@Dimnames[[2]]
    ans
  }
  
  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)
  
  #Filter features
  
  Irrelevant.Features <- which(Var.Arith.Per.Feature == 0)
  
  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (length(dim(Gene.IDs)) > 1){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
  }
  
  colOverdispQPCoef <- function(X,alternative = "greater"){
    stopifnot( methods::is(X,"CsparseMatrix"))
    ans <- sapply( base::seq.int(X@Dim[2]),function(j){
      if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
      #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
      est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
      
      aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                  X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean
      
      mean(aux) + 1})
  }
  
  Alpha.QP <- colOverdispQPCoef(UMI.Mat)
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)
  
  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (length(dim(Gene.IDs)) > 1){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }
  
  Irrelevant.Features <- which(Alpha.QP <= 1)
  
  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (length(dim(Gene.IDs)) > 1){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }
  
  Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature
  
  Theta.NB.Est <- 1/Alpha.NB.Est
  
  df <- data.frame(Mean = Mean.Arith.Per.Feature,AlphaQP = Alpha.QP)
  df$density <- get_density(log10(df$Mean),log10(df$AlphaQP), h = c(1, 1), n = 100)
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- expression("Quasi-Poisson dispersion coefficient" ~ (alpha[QP]))
  } else {
    ylabeltext <- ""
  }
  
  if (Bins == T){
    QuantileIntercepts <- quantile(Mean.Arith.Per.Feature,probs = seq(0,1,0.1))
    
    ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = AlphaQP)) +
      ggplot2::geom_point(ggplot2::aes(Mean, AlphaQP, color = density)) + viridis::scale_color_viridis() +
      ggplot2::geom_vline(size = 0.6,color = "gray40",linetype = "dashed",xintercept = QuantileIntercepts) +
      ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10",scales::math_format(10^.x))) +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10",scales::math_format(10^.x))) +
      ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) + 
      ggplot2::ggtitle(DataName) +
      ggplot2::xlab(xlabeltext) +
      ggplot2::ylab(ylabeltext) +
      ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
      ggplot2::theme(legend.position="none") +
      ggplot2::theme(axis.text = element_text(size = AxisLabelSize)) 
  } else {
    ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = AlphaQP)) +
      ggplot2::geom_point(ggplot2::aes(Mean, AlphaQP, color = density)) + viridis::scale_color_viridis() +
      ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10",scales::math_format(10^.x))) +
      ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10",scales::math_format(10^.x))) +
      ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) + 
      ggplot2::ggtitle(DataName) +
      ggplot2::xlab(xlabeltext) +
      ggplot2::ylab(ylabeltext) +
      ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
      ggplot2::theme(legend.position="none") +
      ggplot2::theme(axis.text = element_text(size = AxisLabelSize)) 
  }
}

#Function for plotting diff alpha vs mean lin-log plots
DiffAlphaLinLogPlot <- function(PiccoloList,DataName,yMax = NULL,yMin = NULL,Reference = 0.1,DiffAlphaThresP = 20,DiffAlphaThresN = -0.5,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 28,AxisLabelSize = 28){
  #For plotting diff alpha
  
  UMI.Mat <- Matrix::t(PiccoloList$Counts)
  
  Gene.IDs <- PiccoloList$Genes
  
  colVarsSPM <- function(X) {
    stopifnot( methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]),function(j) {
      if(X@p[j+1] == X@p[j]) { return(0) } # all entries are 0: var is 0
      mean <- base::sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
      sum((X@x[(X@p[j]+1):X@p[j+1] ] - mean)^2) +
        mean^2 * (X@Dim[1] - (X@p[j+1] - X@p[j]))})/(X@Dim[1] - 1)
    names(ans) <- X@Dimnames[[2]]
    ans
  }
  
  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)
  
  message("Filtering features...")
  
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)
  
  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
  }
  
  message("Estimating overdispersion coefficients...")
  
  colOverdispQPCoef <- function(X,alternative = "greater"){
    stopifnot( methods::is(X,"CsparseMatrix"))
    ans <- sapply( base::seq.int(X@Dim[2]),function(j){
      if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
      est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
      
      aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                  X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean
      
      mean(aux) + 1})
  }
  
  Alpha.QP <- colOverdispQPCoef(UMI.Mat)
  
  Irrelevant.Features <- which(Alpha.QP <= 1)
  
  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }
  
  Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature
  
  #Binning based approach
  
  Mean.Quantiles <- quantile(Mean.Arith.Per.Feature,probs = seq(0.001,1,0.001))
  Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric",length = length(Gene.IDs))
  for (i in 1:length(Mean.Quantiles))
  {
    if (i == 1){
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
    } else {
      Features.In.Bin <- which(Mean.Arith.Per.Feature >= Mean.Quantiles[i-1] & Mean.Arith.Per.Feature <= Mean.Quantiles[i])
    }
    
    Reference.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin],probs = c(Reference))
    Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - Reference.AlphaQP.Bin
  }
  
  df <- data.frame(Mean = Mean.Arith.Per.Feature,DiffAlpha = Diff.AlphaQP.AlphaQPFit)
  if (is.null(dim(Gene.IDs))){
    rownames(df) <- Gene.IDs
  } else {
    if (length(which(duplicated(Gene.IDs$V2))) != 0){
      Gene.IDs$V2[which(duplicated(Gene.IDs$V2))] <- paste0(Gene.IDs$V2[which(duplicated(Gene.IDs$V2))],"_2")
    }
    rownames(df) <- Gene.IDs$V2
  }
  
  if (is.null(yMax) ==  F){
    df$DiffAlpha[df$DiffAlpha >= yMax] <- yMax
  }
  
  if (is.null(yMin) == F){
    df$DiffAlpha[df$DiffAlpha <= yMin] <- yMin
  }
  
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  df$density <- get_density(log10(df$Mean),df$DiffAlpha, h = c(1, 1), n = 100)
  
  df1 <- df[df$DiffAlpha > DiffAlphaThresP,]
  
  df2 <- df[df$DiffAlpha < DiffAlphaThresN,]
 
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- expression(alpha[QP] - alpha["QP(Reference|Bin)"])
  } else {
    ylabeltext <- ""
  }
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = Mean, y = DiffAlpha)) +
    ggplot2::geom_point(ggplot2::aes(Mean, DiffAlpha, color = density)) + viridis::scale_color_viridis() +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10",scales::math_format(10^.x))) +
    ggplot2::geom_hline(size = 1,color = "red",linetype = "dashed",yintercept = DiffAlphaThresP) +
    ggplot2::geom_hline(size = 0.4,color = "black",yintercept = 0) +
    ggplot2::geom_text(size = 6.5,data = df1,ggplot2::aes(label=rownames(df1)),nudge_x = -0.21, nudge_y = 1.6,check_overlap = T) +
    ggplot2::geom_text(size = 6.5,data = df2,ggplot2::aes(label=rownames(df2)),nudge_x = -0.21, nudge_y = -1.6,check_overlap = T) +
    ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
    ggplot2::ggtitle(DataName) +
    ggplot2::ggtitle(DataName) +
    ggplot2::xlab(xlabeltext) +
    ggplot2::ylab(ylabeltext) +
    ggplot2::coord_cartesian(xlim=xLim,ylim=yLim) +
    ggplot2::theme(legend.position="none") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = AxisLabelSize))
}

#####################################################################################
#set working directory to folder containing the counts, gene names, and barcodes files

#As an example we show how to generate the plots for the Svensson 1 technical control data set
setwd("~/Documents/PiccoloPaperData/Svensson1")

#Create Piccolo list object
PiccoloList <- Piccolo::CreatePiccoloList(MTX = "Svensson1_matrix.mtx.gz",Genes = "Svensson1_features.tsv",Barcodes = "Svensson1_barcodes.tsv")

#Filter cells 
PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 100)

#Plot alphaQP vs mean log-log plot with bins
AlphaQPMeanLogLogPlot(PiccoloList = PiccoloList,DataName = "Svensson 1",Bins = T)

#Plot variance vs mean log-log plot
DiffAlphaLinLogPlot(PiccoloList = PiccoloList,DataName = "Svensson 1",DiffAlphaThresP = 20,DiffAlphaThresN = -0.5)


######################################################################
#Script for Figure 2 Panel B 
#CAUTION: Takes a long time to compute

PiccoloList <- Piccolo::CreatePiccoloList(MTX = "PBMC33k.mtx.gz",Genes = "PBMC33k_features.tsv",Barcodes = "PBMC33k_barcodes.tsv")

PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 1,MT.Perc = 100,RP.Perc = 100)

PiccoloList <- Piccolo::SelectFeatures(PiccoloList = PiccoloList,NoOfHVG = 3000)

OriginalHVG <- PiccoloList$HVG
OriginalHVGSerNo <- PiccoloList$HVG.Ser.Nos

Means.Per.Gene.Orig <- Matrix::rowMeans(PiccoloList$Counts[PiccoloList$RelevantGenes.Ser.Nos,]) 

#We need the vegan package for subsampling counts
install.packages("vegan")

#Subsample 30% of the cells to fractions of their total UMI counts

NoOfIterations <- 100
OverlapProp <- vector(mode = "numeric",length = NoOfIterations)
MismatchOrig.Vec <- vector(mode = "numeric",length = NoOfIterations)
SubSampleHVGNo <- vector(mode = "numeric",length = NoOfIterations)
OriginalHVGNo <- rep(length(OriginalHVGSerNo),NoOfIterations)
OrigSet1.Vec <- vector(mode = "numeric",length = NoOfIterations)
OrigSet2.Vec <- vector(mode = "numeric",length = NoOfIterations)
OrigSet3.Vec <- vector(mode = "numeric",length = NoOfIterations)
MismatchedSet1.Vec <- vector(mode = "numeric",length = NoOfIterations)
MismatchedSet2.Vec <- vector(mode = "numeric",length = NoOfIterations)
MismatchedSet3.Vec <- vector(mode = "numeric",length = NoOfIterations)
pval1.Vec <- vector(mode = "numeric",length = NoOfIterations)
pval2.Vec <- vector(mode = "numeric",length = NoOfIterations)
pval3.Vec <- vector(mode = "numeric",length = NoOfIterations)
for (i in 1:NoOfIterations){
  Fractions <- seq(0.3,0.7,0.1)
  
  No.of.Samples.To.Subsample <- ceiling(0.3*ncol(PiccoloList$Counts))
  
  Random.Sample.Nos <- sample(1:ncol(PiccoloList$Counts),size = No.of.Samples.To.Subsample,replace = F)
  
  TempCounts <- matrix(0,nrow = nrow(PiccoloList$Counts),ncol = length(Random.Sample.Nos))
  for (j in 1:length(Random.Sample.Nos)){
    TempJ <- Random.Sample.Nos[j]
    TempJ1 <- sample.int(length(Fractions),1)
    TempCountsVec <- as.vector(PiccoloList$Counts[,TempJ])
    TempCounts[,j] <- vegan::rrarefy(TempCountsVec,ceiling(sum(TempCountsVec)*Fractions[TempJ1]))
  }
  
  TempCounts <- as(TempCounts,"CsparseMatrix")
  
  UnsubsampledCols <- seq(1,ncol(PiccoloList$Counts),1)
  UnsubsampledCols <- UnsubsampledCols[!UnsubsampledCols %in% Random.Sample.Nos]
  SubsampledMat <- cbind(PiccoloList$Counts[,UnsubsampledCols],TempCounts)
  
  PiccoloList1 <- PiccoloList
  PiccoloList1$Counts <- SubsampledMat
  PiccoloList1$Barcodes <- PiccoloList1$Barcodes[c(UnsubsampledCols,Random.Sample.Nos)]
  
  PiccoloList1 <- Piccolo::SelectFeatures(PiccoloList = PiccoloList1,NoOfHVG = 3000)
  
  SubsampleHVG <- PiccoloList1$HVG
  SubsampleHVGSerNo <- PiccoloList1$HVG.Ser.Nos
  
  MisMatchedHVGOrig <- OriginalHVG$V1[!OriginalHVG$V1 %in% SubsampleHVG$V1]
  MisMatchedHVGOrigSerNo <- OriginalHVGSerNo[!OriginalHVGSerNo %in% SubsampleHVGSerNo]
  
  OrigSet1 <- length(which(Means.Per.Gene.Orig < 0.1))
  
  OrigSet2 <- length(which(Means.Per.Gene.Orig >= 0.1 & Means.Per.Gene.Orig < 1))
  
  OrigSet3 <- length(which(Means.Per.Gene.Orig >= 1))
  
  MisMatchedOriginalMeans <- Matrix::rowMeans(PiccoloList$Counts[MisMatchedHVGOrigSerNo,])
  
  MismatchedOrigSet1 <- length(which(MisMatchedOriginalMeans < 0.1))
  
  p.val.Set1 <- sum(dhyper(MismatchedOrigSet1:length(MisMatchedHVGOrigSerNo),OrigSet1,length(PiccoloList$FilteredGenes.Ser.Nos) - OrigSet1,length(MisMatchedHVGOrigSerNo)))
  
  MismatchedOrigSet2 <- length(which(MisMatchedOriginalMeans >= 0.1 & MisMatchedOriginalMeans < 1))
  
  p.val.Set2 <- sum(dhyper(MismatchedOrigSet2:length(MisMatchedHVGOrigSerNo),OrigSet2,length(PiccoloList$FilteredGenes.Ser.Nos) - OrigSet2,length(MisMatchedHVGOrigSerNo)))
  
  MismatchedOrigSet3 <- length(which(MisMatchedOriginalMeans >= 1))
  
  p.val.Set3 <- sum(dhyper(MismatchedOrigSet3:length(MisMatchedHVGOrigSerNo),OrigSet3,length(PiccoloList$FilteredGenes.Ser.Nos) - OrigSet3,length(MisMatchedHVGOrigSerNo)))
  
  OverlapProp[i] <- length(intersect(SubsampleHVGSerNo,OriginalHVGSerNo))/min(c(length(SubsampleHVGSerNo),length(OriginalHVGSerNo)))
  SubSampleHVGNo[i] <- length(PiccoloList1$HVG.Ser.Nos)
  
  OrigSet1.Vec[i] <- OrigSet1
  OrigSet2.Vec[i] <- OrigSet2
  OrigSet3.Vec[i] <- OrigSet3
  pval1.Vec[i] <- p.val.Set1
  pval2.Vec[i] <- p.val.Set2
  pval3.Vec[i] <- p.val.Set3
  MismatchedSet1.Vec[i] <- MismatchedOrigSet1
  MismatchedSet2.Vec[i] <- MismatchedOrigSet2
  MismatchedSet3.Vec[i] <- MismatchedOrigSet3
  
  MismatchOrig.Vec[i] <- length(MisMatchedHVGOrigSerNo)
}

Res.Summ.df <- data.frame(OrigHVGNo = OriginalHVGNo,SubsampleHVGNo = SubSampleHVGNo,MismatchOrig = MismatchOrig.Vec,MatchProp = OverlapProp,Orig.Set1 = OrigSet1.Vec,Mismatched.Set1 = MismatchedSet1.Vec,p.val.1 = pval1.Vec,Orig.Set2 = OrigSet2.Vec,Mismatched.Set2 = MismatchedSet2.Vec,p.val.2 = pval2.Vec,Orig.Set3 = OrigSet3.Vec,Mismatched.Set3 = MismatchedSet3.Vec,p.val.3 = pval3.Vec)

write.csv(Res.Summ.df,file = "PBMC33k_FeatureSelectionSubsampledSummary.csv",row.names = F)

#Prepare violin plot to compare results of consistency

df.3T3 <- data.table::fread("NIH3T310Xv3_FeatureSelectionSubsampledSummary.csv")

df.PBMC3k <- data.table::fread("PBMC33k_FeatureSelectionSubsampledSummary.csv")

df.MouseCortex <- data.table::fread("MouseCortexDorNCseq_FeatureSelectionSubsampledSummary.csv")

df.Angelidis <- data.table::fread("MouseLungAngelidis2019_FeatureSelectionSubsampledSummary.csv")

df.inDrops <- data.table::fread("PBMCinDropsR1_FeatureSelectionSubsampledSummary.csv")

DataLabels <- rep(c("NIH/3T3","PBMC 33k","Mouse Cortex","Mouse Lung","PBMC r1"),each = 100)

PlatformLabels <- rep(c("10X v3","10X v1","DroNc-seq","Drop-seq","inDrops"),each = 100)

df <- data.frame(Data = DataLabels,Platform  = PlatformLabels,MatchProp = 100*c(df.3T3$MatchProp,df.PBMC3k$MatchProp,df.MouseCortex$MatchProp,df.Angelidis$MatchProp,df.inDrops$MatchProp))

df$Data <- factor(df$Data,levels = c("PBMC 33k","NIH/3T3","Mouse Cortex","Mouse Lung","PBMC r1"))

BaseSize <- 28

p <- ggplot2::ggplot(df, ggplot2::aes(Data, MatchProp, fill = Platform)) +
  ggplot2::geom_violin() + 
  ggplot2::geom_boxplot(width = 0.2) +
  ggsci::scale_fill_npg() +
  ggplot2::theme_bw(base_size = BaseSize,base_line_size = 0.4) +
  ggplot2::ylim(c(50,100)) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 45,hjust=1)) +
  ggplot2::ggtitle("HVGs - Subsampled vs Original") +
  ggplot2::xlab("") +
  ggplot2::ylab("Percentage Match")

p

