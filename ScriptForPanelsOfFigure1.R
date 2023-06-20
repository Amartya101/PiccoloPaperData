#Script to prepare figures in Fig. 1 (and Supplementary figures 1-4) in the Piccolo paper

#Install Piccolo

install.packages("devtools")

devtools::install_github("Amartya101/Piccolo")

########################################################################################

#Install package MASS

install.packages("MASS")

#functions needed to make the plots

#Function for plotting variance vs mean log-log plots
VarMeanLogLogPlot <- function(PiccoloList,DataName,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 32,AxisLabelSize = 32){
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
  
  library(MASS)
  library(ggplot2)
  library(viridis)
  #theme_set(theme_bw(base_size = 16))
  
  # Get density of points in 2 dimensions.
  # @param x A numeric vector.
  # @param y A numeric vector.
  # @param n Create a square n by n grid to compute density.
  # @return The density within each square.
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  
  df <- data.frame(Mean = Mean.Arith.Per.Feature,Variance = Var.Arith.Per.Feature)
  #df$density <- get_density(df$Mean, df$Theta, n = 1000)
  df$density <- get_density(log10(df$Mean),log10(df$Variance), h = c(1, 1), n = 100)
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- expression("Variance" ~ (sigma ^ 2))
  } else {
    ylabeltext <- ""
  }
  
  ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = Variance)) +
    #ggplot2::geom_point(size = 0.5,color = "darkblue") +
    geom_point(aes(Mean, Variance, color = density)) + viridis::scale_color_viridis() +
    ggplot2::geom_line(size = 0.85,color = "black",data = df,ggplot2::aes(x = Mean,y = Mean)) +
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

#Function for plotting theta vs mean log-log plots
ThetaMeanLogLogPlot <- function(PiccoloList,DataName,Contours = F,xLabel = T,yLabel = T,xLim = NULL,yLim = NULL,BaseSize = 32,AxisLabelSize = 32){
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
  
  df <- data.frame(Mean = Mean.Arith.Per.Feature,Theta = Theta.NB.Est)
  
  df$density <- get_density(log10(df$Mean),log10(df$Theta), h = c(1, 1), n = 100)
  
  Max.Density.Ser.No <- which(df$density == max(df$density))[1]
  
  Alpha.QP.Est <- Alpha.QP[Max.Density.Ser.No]
  
  if (xLabel == T){
    xlabeltext <- expression("Mean" ~ (mu))
  } else {
    xlabeltext <- ""
  }
  
  if (yLabel == T){
    ylabeltext <- expression("Inverse overdispersion coefficient" ~ (theta))
  } else {
    ylabeltext <- ""
  }
  
  if (Contours == T){
    ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = Theta)) +
      #ggplot2::geom_point(size = 0.5,color = "darkblue") +
      ggplot2::geom_point(ggplot2::aes(Mean, Theta, color = density)) + viridis::scale_color_viridis() +
      ggplot2::stat_density_2d(size = 0.7,ggplot2::aes(fill = density), geom = "polygon", colour="white",alpha  = 0.01)+
      #ggplot2::geom_line(size = 1,color = "red",linetype = "dashed",data = df,ggplot2::aes(x = Mean,y = Mean/(quantile(Alpha.QP[log10(Mean.Arith.Per.Feature) < 0],probs = c(0.5)) - 1))) +
      ggplot2::geom_line(size = 1,color = "red",linetype = "dashed",data = df,ggplot2::aes(x = Mean,y = Mean/(Alpha.QP.Est - 1))) +
      
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
      ggplot2::theme(axis.text = element_text(size = 32)) 
  } else {
    ggplot2::ggplot(data = df,ggplot2::aes(x = Mean,y = Theta)) +
      #ggplot2::geom_point(size = 0.5,color = "darkblue") +
      ggplot2::geom_point(ggplot2::aes(Mean, Theta, color = density)) + viridis::scale_color_viridis() +
      #stat_density_2d(size = 0.25,aes(fill = density), geom = "polygon", colour="white",alpha  = 0.01)+
      #ggplot2::geom_line(size = 1,color = "red",linetype = "dashed",data = df,ggplot2::aes(x = Mean,y = Mean/(quantile(Alpha.QP[log10(Mean.Arith.Per.Feature) < 0],probs = c(0.5)) - 1))) +
      ggplot2::geom_line(size = 1,color = "red",linetype = "dashed",data = df,ggplot2::aes(x = Mean,y = Mean/(Alpha.QP.Est - 1))) +
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

#####################################################################################

#set working directory to folder containing the counts, gene names, and barcodes files

#As an example we show how to generate the plots for the Svensson 1 technical control data set
setwd("~/Documents/PiccoloPaperData/Svensson1")

#Create Piccolo list object
PiccoloList <- Piccolo::CreatePiccoloList(MTX = "Svensson1_matrix.mtx.gz",Genes = "Svensson1_features.tsv",Barcodes = "Svensson1_barcodes.tsv")

#Filter cells 
PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 100)

#Plot variance vs mean log-log plot
VarMeanLogLogPlot(PiccoloList = PiccoloList,DataName = "Svensson 1")

#Plot AlphaQP vs mean log-log plot
AlphaQPMeanLogLogPlot(PiccoloList = PiccoloList,DataName = "Svensson 1")

#Plot Theta vs mean log-log plot
ThetaMeanLogLogPlot(PiccoloList = PiccoloList,DataName = "Svensson 1",Contours = T)

