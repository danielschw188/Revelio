includeIntronicDataAndFilter <- function(dataList,
                                         intronicData){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': copy intronic data: ', sep = ''))

  dataList@DGEs$intronCountData <- intronicData

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': filtering data by intronic data: ', sep = ''))

  filteredCells <- as.character(t(subset(dataList@cellInfo, cellID%in%colnames(dataList@DGEs$intronCountData), select = 'cellID')))
  filteredGenes <- as.character(t(subset(dataList@geneInfo, geneID%in%rownames(dataList@DGEs$intronCountData), select = 'geneID')))
  dataList@DGEs$intronCountData <- dataList@DGEs$intronCountData[filteredGenes, filteredCells]
  dataList@DGEs$countData <- dataList@DGEs$countData[filteredGenes, filteredCells]
  dataList@geneInfo <- subset(dataList@geneInfo, geneID%in%filteredGenes)
  dataList@geneInfo <- subset(dataList@geneInfo, rowSums(dataList@DGEs$countData)>0)
  dataList@DGEs$intronCountData <- dataList@DGEs$intronCountData[rowSums(dataList@DGEs$countData)>0,]
  dataList@DGEs$countData <- dataList@DGEs$countData[rowSums(dataList@DGEs$countData)>0,]
  dataList@cellInfo <- subset(dataList@cellInfo, cellID%in%filteredCells)

  dataList <- getStatisticsOnCellsAndGenes(dataList = dataList)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  return(dataList)
}
#'
#'
#' Get Velocity Matrix.
#'
#' 'getVelocityMatrix' reads the intronic data, filters the existing gene list and calculates the RNA velocity matrix according to Model I in La Manno et al. 2018 ('RNA velocity of single cells').
#'
#' The intronic data should be corresponding to the exonic data previously provided. First, the offset parameter is estimated, then the gamma parameter and finally the velocity matrix from the model I from La Manno et al. 2008 ('RNA velocity of single cells') where the time derivative of spliced counts is assumed constant.
#' @param dataList A Revelio object that contains a raw data matrix with assigned cell cycle phases.
#' @param rawIntronicData A genes-by-cells matrix of raw intronic sequencing data. Rows should contain gene names, columns should contain cell IDs, entries of the matrix are integer numbers. The cell IDs should coincide with the cell IDs provided in the exonic matrix previously.
#' @return Returns the same Revelio object given as input but now with an added intronic data matrix and the velocity matrix.
#'
#' @export
getVelocityMatrix <- function(dataList,
                              rawIntronicData){

  dataList <- includeIntronicDataAndFilter(dataList = dataList,
                                           intronicData = rawIntronicData)

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating velocity offset estimate: ', sep = ''))
  dataList <- offsetEstimation(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating velocity gamma estimate: ', sep = ''))
  dataList <- gammaEstimation(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating RNA velocity matrix: ', sep = ''))
  dataList <- calculateRNAvelocity(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  return(dataList)
}
offsetEstimation <- function(dataList){

  dataUnspliced <- dataList@DGEs$intronCountData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]
  dataSpliced <- dataList@DGEs$countData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]

  normalizedUnspliced <- getNormalizedDataBySpecificVector(data = dataUnspliced,
                                                           divideByCustomVector = colSums(dataUnspliced + dataSpliced))
  normalizedSpliced <- getNormalizedDataBySpecificVector(data = dataSpliced,
                                                         divideByCustomVector = colSums(dataUnspliced + dataSpliced))

  offsetEstimate <- rep(0, dim(normalizedSpliced)[1])
  for (i in 1:dim(normalizedSpliced)[1]){
    if (sum(normalizedSpliced[i,] == 0) == 0){
      offsetEstimate[i] <- 0
    }else{
      offsetEstimate[i] <- mean(as.matrix(normalizedUnspliced[i, normalizedSpliced[i,] == 0]))
    }
  }

  dataList@velocityData$offsetEstimate <- offsetEstimate

  return(dataList)
}
gammaEstimation <- function(dataList){

  dataUnspliced <- dataList@DGEs$intronCountData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]
  dataSpliced <- dataList@DGEs$countData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]

  normalizedUnspliced <- getNormalizedDataBySpecificVector(data = dataUnspliced,
                                                           divideByCustomVector = colSums(dataUnspliced + dataSpliced))
  normalizedSpliced <- getNormalizedDataBySpecificVector(data = dataSpliced,
                                                         divideByCustomVector = colSums(dataUnspliced + dataSpliced))
  offsetEstimate <- dataList@velocityData$offsetEstimate

  gammaEstimate <- rep(0, dim(normalizedSpliced)[1])
  interceptEstimate <- rep(0, dim(normalizedSpliced)[1])

  for (i in 1:dim(normalizedSpliced)[1]){
    #linearRegressionResult <- lsfit(x = t(normalizedSpliced[i,]), y = t(normalizedUnspliced[i,]))
    linearRegressionResult <- lm((t(normalizedUnspliced[i,]) - rep(offsetEstimate[i], dim(normalizedSpliced)[2])) ~ 0 + t(normalizedSpliced[i,]))
    if (is.na(linearRegressionResult[1])){
      gammaEstimate[i] <- 0
    }else{
      gammaEstimate[i] <- linearRegressionResult$coefficients[1]
    }
    #interceptEstimate[i] <- linearRegressionResult$coefficients[1]
  }

  dataList@velocityData$gammaEstimate <- gammaEstimate

  return(dataList)
}
calculateRNAvelocity <- function(dataList){

  dataUnspliced <- dataList@DGEs$intronCountData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]
  dataSpliced <- dataList@DGEs$countData[intersect(rownames(dataList@DGEs$intronCountData), rownames(dataList@DGEs$countData)), intersect(colnames(dataList@DGEs$intronCountData), colnames(dataList@DGEs$countData))]

  normalizedUnspliced <- getNormalizedDataBySpecificVector(data = dataUnspliced,
                                                           divideByCustomVector = colSums(dataUnspliced + dataSpliced))
  normalizedSpliced <- getNormalizedDataBySpecificVector(data = dataSpliced,
                                                         divideByCustomVector = colSums(dataUnspliced + dataSpliced))
  gamma <- dataList@velocityData$gammaEstimate
  offset <- dataList@velocityData$offsetEstimate
  # timeStep <- object@projectParameters$timeStepForVelocity

  velocity <- normalizedUnspliced - replicate(dim(normalizedSpliced)[2], gamma) * normalizedSpliced - replicate(dim(normalizedSpliced)[2], offset)
  # normalizedExtrapolatedState <- normalizedSpliced + velocity*timeStep
  # normalizedExtrapolatedState[normalizedExtrapolatedState < 0] <- 0

  dataList@velocityData$RNAvelocity <- as.data.frame(velocity)
  # object@listExtrapolatedStateModelI@dataNormalized <- normalizedExtrapolatedState

  return(dataList)
}
getNormalizedDataBySpecificVector <- function(data,
                                              divideByCustomVector = NULL){
  return(sweep(data,2,divideByCustomVector,'/'))
}
#'
#'
#' Transform the extrapolated space from RNA velocity into RNA velocity.
#'
#' 'getExtrapolatedStateInDCSpace' calculates the extrapolated state of each cells and transforms it into the coordinate system of the dynamical components.
#'
#' For each cell the extrapolated state is calculated (using the RNA velocity matrix calculated previously) via the model s_ij=s_ij+v. This extrapolated state is then transformed into the dynamical component space using the same transformation matrices generated by only the exonic data.
#' @param dataList A Revelio object that contains a rotation matrix and a velocity matrix.
#' @return Returns the same Revelio object given as input but additionally with the coordinates of the extrapolated state in DC space.
#'
#' @export
getExtrapolatedStateInDCSpace <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating extrapolated state: ', sep = ''))

  extrapolatedStateNormalizedCurrent <- extrapolatedStateFromVelocityModelI(dataList)

  dataList <- transformExtrapolatedStateModelIToRotatedPCA(dataList,
                                                           extrapolatedStateNormalized = extrapolatedStateNormalizedCurrent,
                                                           numberOfDimensionsToGet = sum(dataList@geneInfo$pcaGenes))

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
extrapolatedStateFromVelocityModelI <- function(dataList){

  velocity <- dataList@velocityData$RNAvelocity[dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes], colnames(dataList@DGEs$countData)]
  #    velocity <- dataList@velocityData$RNAvelocity[rownames(dataList@velocityData$RNAvelocity)%in%rownames(object@listDataProcessed@data), colnames(object@listDataProcessed@data)]
  timeStep <- 1#dataList@datasetInfo$velTimeStep

  normalizedExtrapolatedState <- (exp(dataList@DGEs$logOfFractionsData)-1)[dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes],]/dataList@datasetInfo$scalingFactorUMI + velocity*timeStep
  #    normalizedExtrapolatedState <- object@listDataProcessed@dataFractions[rownames(object@listDataProcessed@data)%in%rownames(dataList@velocityData$RNAvelocity),]/object@projectParameters$scalingFactorUMI + velocity*timeStep
  normalizedExtrapolatedState[normalizedExtrapolatedState < 0] <- 0

  return(normalizedExtrapolatedState)
}
transformExtrapolatedStateModelIToRotatedPCA <- function(dataList,
                                                         extrapolatedStateNormalized,
                                                         numberOfDimensionsToGet = 2){

  data <- extrapolatedStateNormalized
  data <- data*dataList@datasetInfo$scalingFactorUMI                                                   #multiplies every entry by 10^4 as a scaling factor
  data <- log(data+1)                                                   #transforms to logarithmic scale (natural logarithm)

  #normalize data by rows
  meanUMIForEachGene <- rowMeans(data)                                   #calculates mean UMI count per gene
  data <- sweep(data,1,meanUMIForEachGene,'-')                           #centralizes gene data
  sdUMIForEachGene <- apply(data,1,sd)                                   #calculates sd of UMI count per gene
  extrapolatedStateScaled <- sweep(data,1,sdUMIForEachGene,'/')          #divides each row by sd of that gene

  extrapolatedStateScaled <- extrapolatedStateScaled[dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes],]
  pcaDataExtrapolatedState <- as.data.frame(t(as.matrix(extrapolatedStateScaled))%*%as.matrix(t(dataList@transformedData$pca$weights)))
  rownames(pcaDataExtrapolatedState) <- colnames(extrapolatedStateScaled)
  rotationMatrixToUse <- as.matrix(dataList@transformedData$dc$rotationMatrix)

  dataList@velocityData$extrapolatedStateInDCSpace <- t(as.data.frame(as.matrix(pcaDataExtrapolatedState) %*% rotationMatrixToUse)[,1:numberOfDimensionsToGet])

  return(dataList)
}
#'
#'
#' Calculate the velocity grid of the cell cycle.
#'
#' 'getVelocityGridCC' places a regular grid on the two-dimensional cell cycle and calculates the average velocity for each grid point as described in La Manno et al. 2018 ('RNA velocity of single cells').
#'
#' A grid is generated overlaying the two-dimensional cell cycle in DC space. For each grid point, the RNA velocity is calculated by using a Gaussian kernel as explained in La Manno et al. 2018 ('RNA velocity of single cells').
#' @param dataList A Revelio object that contains information about the extrapolated state of cells.
#' @param sigma The smoothing parameter for the Gaussian kernel.
#' @param numberOfGridPointsOneSide The number of grid points along one side. The total grid will be a numberOfGridPointsOneSide-by-numberOfGridPointsOneSide grid.
#' @param boolPlotResults TRUE/FALSE if the velocity grid plot should be shown.
#' @return Returns the same Revelio object given as input with the additional information about the RNA velocity at grid points.
#'
#' @export
getVelocityGridCC <- function(dataList,
                              sigma = 1,
                              numberOfGridPointsOneSide = 20,
                              boolPlotResults = FALSE){

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating cell cycle grid for sigma = ', sigma, ': ', sep = ''))

  originalStateInRotatedSpace <- t(dataList@transformedData$dc$data[c(1,2),])
  extrapolatedStateInRotatedStateSpace <- t(dataList@velocityData$extrapolatedStateInDCSpace[c(1,2), rownames(originalStateInRotatedSpace)])
  gridBoundary <- 2.5*max(sd(originalStateInRotatedSpace[,1]), sd(originalStateInRotatedSpace[,2]))
  gridX <- seq(from = -gridBoundary, to = gridBoundary, length.out = numberOfGridPointsOneSide)
  gridY <- seq(from = -gridBoundary, to = gridBoundary, length.out = numberOfGridPointsOneSide)
  gridFinal <- expand.grid(x = gridX, y = gridY)

  euclideanDistance <- function(x,y){
    return(sqrt(sum((x-y)^2)))
  }

  dataDisplacementOfXiDirect <- extrapolatedStateInRotatedStateSpace - originalStateInRotatedSpace
  dataDisplacementGrid <- as.data.frame(matrix(0L, dim(gridFinal)[1], dim(gridFinal)[2]))
  for (i in 1:dim(gridFinal)[1]){
    currentGridPoint <- gridFinal[i,]
    euclideanDistanceVector <- apply(originalStateInRotatedSpace, 1, euclideanDistance, y = currentGridPoint)
    gaussianKernel <- exp(-euclideanDistanceVector^2/(2*sigma^2))
    dataDisplacementGridPC1 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,1])
    dataDisplacementGridPC2 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,2])

    dataDisplacementGrid[i,] <- cbind(dataDisplacementGridPC1, dataDisplacementGridPC2)
  }

  colnames(gridFinal) <- c('gridX', 'gridY')
  colnames(dataDisplacementGrid) <- c('displacementX', 'displacementY')

  dataList@velocityData$grids[[paste('cellCycle2D_sigma_', sigma, sep = '')]] <- cbind(gridFinal, dataDisplacementGrid)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  if (boolPlotResults){
    plotParameters <- list()
    plotParameters$colorPaletteCellCyclePhasesGeneral <- c('#ac4343', '#466caf', '#df8b3f', '#63b558', '#e8d760', '#61c5c7', '#f04ddf', '#a555d4')
    plotParameters$plotLabelTextSize <- 8
    plotParameters$plotDotSize <- 0.1
    plotParameters$plotLineWidth <- 0.2
    plotParameters$fontFamily <- 'Helvetica'
    plotParameters$fontSize <- 8
    plotParameters$colorPaletteCellCyclePhasesGeneral <- plotParameters$colorPaletteCellCyclePhasesGeneral[1:length(levels(dataList@cellInfo$ccPhase))]
    names(plotParameters$colorPaletteCellCyclePhasesGeneral) <- levels(dataList@cellInfo$ccPhase)

    plotBoundary <- max(dataList@transformedData$dc$data[c('DC1', 'DC2'),])*0.78

    ccPhaseBorderPathPolar <- data.frame(angle = rep(0,7), radius = rep(0,7))
    ccPhaseBorderPathPolar[c(1,3,5,7),1] <- 2*pi*(1-cumsum(c(dataList@datasetInfo$ccDurationG1, dataList@datasetInfo$ccDurationS, dataList@datasetInfo$ccDurationG2, dataList@datasetInfo$ccDurationM))/dataList@datasetInfo$ccDurationTotal)
    ccPhaseBorderPathPolar[c(1,3,5,7),2] <- 50
    ccPhaseBorderPathCartesian <- ccPhaseBorderPathPolar
    ccPhaseBorderPathCartesian[,1] <- ccPhaseBorderPathPolar[,2]*cos(ccPhaseBorderPathPolar[,1])
    ccPhaseBorderPathCartesian[,2] <- ccPhaseBorderPathPolar[,2]*sin(ccPhaseBorderPathPolar[,1])
    colnames(ccPhaseBorderPathCartesian) <- c('xValue', 'yValue')

    labelPositionHelp <- append(2*pi, ccPhaseBorderPathPolar[c(1,3,5,7),1])
    labelPositionHelp <- (labelPositionHelp[1:(length(labelPositionHelp)-1)]-labelPositionHelp[2:length(labelPositionHelp)])/2+labelPositionHelp[2:length(labelPositionHelp)]
    labelRadiusHelp <- labelPositionHelp
    labelRadiusHelp[(labelRadiusHelp>pi/4)&(labelRadiusHelp<3*pi/4)] <- labelRadiusHelp[(labelRadiusHelp>pi/4)&(labelRadiusHelp<3*pi/4)]-pi/2
    labelRadiusHelp[(labelRadiusHelp>5*pi/4)&(labelRadiusHelp<7*pi/4)] <- labelRadiusHelp[(labelRadiusHelp>5*pi/4)&(labelRadiusHelp<7*pi/4)]-pi/2

    labelPositionPolar <- data.frame(angle = labelPositionHelp, radius = sqrt((plotBoundary*tan(labelRadiusHelp))^2+plotBoundary^2), label = c('G1', 'S', 'G2', 'M'))
    labelPositionCartesian <- labelPositionPolar
    labelPositionCartesian[,1] <- labelPositionPolar[,2]*cos(labelPositionPolar[,1])
    labelPositionCartesian[,2] <- labelPositionPolar[,2]*sin(labelPositionPolar[,1])
    colnames(labelPositionCartesian) <- c('xValue', 'yValue', 'label')

    plotDC1DC2 <- ggplot(data = cbind(as.data.frame(t(dataList@transformedData$dc$data)), ccPhase = dataList@cellInfo$ccPhase))+
      theme_gray(base_size = plotParameters$plotLabelTextSize)+
      theme(text=element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
            axis.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
            axis.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
            legend.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
            legend.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize))+
      geom_point(aes(x = DC1, y = DC2, color = ccPhase), size = plotParameters$plotDotSize*5)+
      coord_cartesian(xlim = c(-plotBoundary, plotBoundary),ylim = c(-plotBoundary, plotBoundary))+
      scale_color_manual(values = plotParameters$colorPaletteCellCyclePhasesGeneral,
                         labels = levels(dataList@cellInfo$ccPhase))+
      theme(legend.position = 'right',
            axis.text.y=element_text(angle=90, hjust=0.5),
            legend.title = element_blank(),
            aspect.ratio=1)+
      guides(color = guide_legend(override.aes = list(size = 1)))

    scaleFactorArrow <- 1/max(sqrt(dataList@velocityData$grids[[paste('cellCycle2D_sigma_', sigma, sep = '')]][,'displacementX']^2+dataList@velocityData$grids[[paste('cellCycle2D_sigma_', sigma, sep = '')]][,'displacementY']^2))*3#2.5
    arrowLengthThreshold <- 0.2
    dataToUse <- filter(dataList@velocityData$grids[[paste('cellCycle2D_sigma_', sigma, sep = '')]], sqrt(displacementX^2+displacementY^2)*scaleFactorArrow>arrowLengthThreshold)
    dataToUse[,c('displacementX', 'displacementY')] <- dataToUse[,c('displacementX', 'displacementY')]*scaleFactorArrow
    plotDC1DC2_WithVelocityGrid <- plotDC1DC2+
      geom_point(data = dataToUse,
                 aes(x = gridX, y = gridY),
                 color = 'black',
                 size = plotParameters$plotDotSize/20)+
      geom_segment(data = dataToUse,
                   aes(x = gridX, xend = gridX+displacementX,y = gridY, yend = gridY+displacementY),
                   color = 'black',
                   size = plotParameters$plotLineWidth*1.35,
                   arrow=arrow(length = unit(plotParameters$plotLineWidth/3*1.5, "cm")))

    grid.arrange(plotDC1DC2_WithVelocityGrid)

  }
  return(dataList)
}
getCartesianCoordinates <- function(dataPolar){
  data <- as.data.frame(cbind(dataPolar[,'radius']*cos(dataPolar[,'angle']), dataPolar[,'radius']*sin(dataPolar[,'angle'])))
  rownames(data) <- rownames(dataPolar)
  return(data)
}
#'
#'
#' Calculate the velocity grid along a third dimension.
#'
#' 'getVelocityGridAlongThirdDimensions' places a regular grid on the two-dimensional plot consisting of cell cycle pseudotime and a higher order DC and calculates the average velocity for each grid point as described in La Manno et al. 2018 ('RNA velocity of single cells').
#'
#' A grid is generated overlaying the two-dimensional plot of cell cycle pseudotime and a higher order DC. For each grid point, the RNA velocity is calculated by using a Gaussian kernel as explained in La Manno et al. 2018 ('RNA velocity of single cells').
#' @param dataList A Revelio object that contains information about the extrapolated state of cells.
#' @param whichThirdDimension The number of the dimension to be investigated. Needs to be an integer.
#' @param sigma The smoothing parameter for the Gaussian kernel.
#' @param numberOfGridPointsOneSide The number of grid points along one side. The total grid will be a numberOfGridPointsOneSide-by-numberOfGridPointsOneSide grid.
#' @param boolPlotResults TRUE/FALSE if the velocity grid plot should be shown.
#' @return Returns the same Revelio object given as input with the additional information about the RNA velocity at grid points.
#'
#' @export
getVelocityGridAlongThirdDimensions <- function(dataList,
                                                whichThirdDimension = 3,
                                                sigma = 1,
                                                numberOfGridPointsOneSide = 20,
                                                boolPlotResults = FALSE){

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating velocity grid along DC', whichThirdDimension, ': ', sep = ''))

  meanRadiusOfCells <- mean(dataList@cellInfo$ccRadius)
  sdRadiusOfCells <- sd(dataList@cellInfo$ccRadius)
  minimalRadius <- meanRadiusOfCells-2*sdRadiusOfCells
  if (minimalRadius<0){
    minimalRadius <- meanRadiusOfCells-1*sdRadiusOfCells
  }
  if (minimalRadius<0){
    minimalRadius <- 0
  }
  cellsOutsideEye <- (dataList@cellInfo$cellID[dataList@cellInfo$ccPositionIndex])[dataList@cellInfo$ccRadius[dataList@cellInfo$ccPositionIndex]>=minimalRadius]

  originalStateInRotatedSpace <- t(dataList@transformedData$dc$data[c(1,2), cellsOutsideEye])
  extrapolatedStateInRotatedStateSpace <- t(dataList@velocityData$extrapolatedStateInDCSpace[c(1,2,whichThirdDimension), rownames(originalStateInRotatedSpace)[rownames(originalStateInRotatedSpace)%in%cellsOutsideEye]])
  meanRadiusOfCells <- median(t(subset(dataList@cellInfo, dataList@cellInfo$cellID%in%cellsOutsideEye, select = ccRadius)))

  #new approach: project cells onto mean circle and infer change in cell cycle angle from there
  velocityVector <- extrapolatedStateInRotatedStateSpace[,c(1,2)]-originalStateInRotatedSpace
  originalStateProjectedToCirclePolar <- getPolarCoordinates(originalStateInRotatedSpace)
  originalStateProjectedToCirclePolar[,'radius'] <- rep(meanRadiusOfCells, dim(originalStateProjectedToCirclePolar)[1])
  originalStateProjectedToCircle <- getCartesianCoordinates(originalStateProjectedToCirclePolar)
  extrapolatedStateFromProjectedCircle <- originalStateProjectedToCircle+velocityVector

  currentStateAngle <- t(subset(dataList@cellInfo[dataList@cellInfo$ccPositionIndex,], dataList@cellInfo$cellID%in%cellsOutsideEye, select = ccAngle))
  currentStateAngle <- 2*pi-currentStateAngle
  polarDataExtrapolatedState <- getPolarCoordinates(extrapolatedStateFromProjectedCircle[,c(1,2)])
  polarDataExtrapolatedState <- polarDataExtrapolatedState[,'angle']
  polarDataExtrapolatedState <- 2*pi - polarDataExtrapolatedState

  #get any cells that cross x>0,y=0 during extrapolation
  if (sum((polarDataExtrapolatedState-currentStateAngle)<(-pi))>0){
    polarDataExtrapolatedState[(polarDataExtrapolatedState-currentStateAngle)<(-pi)] <- polarDataExtrapolatedState[(polarDataExtrapolatedState-currentStateAngle)<(-pi)]+2*pi
  }
  if (sum((polarDataExtrapolatedState-currentStateAngle)>pi)>0){
    polarDataExtrapolatedState[(polarDataExtrapolatedState-currentStateAngle)>pi] <- polarDataExtrapolatedState[(polarDataExtrapolatedState-currentStateAngle)>pi]-2*pi
  }
  polarDataExtrapolatedState <- polarDataExtrapolatedState*meanRadiusOfCells
  currentStateAngle <- currentStateAngle*meanRadiusOfCells

  dataCurrentStateCCAndRotatedThirdDimension <- as.data.frame(t(rbind(currentStateAngle, dataList@transformedData$dc$data[whichThirdDimension, cellsOutsideEye])))
  colnames(dataCurrentStateCCAndRotatedThirdDimension) <- c('V1', 'V2')

  dataExtrapolatedStateCCAndRotatedThirdDimension <- as.data.frame(cbind(polarDataExtrapolatedState, extrapolatedStateInRotatedStateSpace[,3]))
  colnames(dataExtrapolatedStateCCAndRotatedThirdDimension) <- c('V1', 'V2')

  arrowLengths <- sqrt(rowSums((dataExtrapolatedStateCCAndRotatedThirdDimension-dataCurrentStateCCAndRotatedThirdDimension)^2))
  arrowsProperLength <- arrowLengths<(mean(arrowLengths)+4*sd(arrowLengths))
  dataCurrentStateCCAndRotatedThirdDimension <- dataCurrentStateCCAndRotatedThirdDimension[arrowsProperLength,]
  dataExtrapolatedStateCCAndRotatedThirdDimension <- dataExtrapolatedStateCCAndRotatedThirdDimension[arrowsProperLength,]

  gridX <- seq(from = 0, to = 2*pi*meanRadiusOfCells, length.out = numberOfGridPointsOneSide)
  gridY <- seq(from = -pi*meanRadiusOfCells, to = pi*meanRadiusOfCells, length.out = numberOfGridPointsOneSide)
  gridFinal <- expand.grid(x = gridX, y = gridY)

  euclideanDistance <- function(x,y){
    return(sqrt(sum((x-y)^2)))
  }

  dataCurrentStateCCAndRotatedThirdDimensionThreePeriods <- rbind(dataCurrentStateCCAndRotatedThirdDimension,dataCurrentStateCCAndRotatedThirdDimension,dataCurrentStateCCAndRotatedThirdDimension)
  dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[1:dim(dataCurrentStateCCAndRotatedThirdDimension)[1],1] <- dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[1:dim(dataCurrentStateCCAndRotatedThirdDimension)[1],1] - 2*pi*meanRadiusOfCells
  dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[(2*dim(dataCurrentStateCCAndRotatedThirdDimension)[1]+1):dim(dataCurrentStateCCAndRotatedThirdDimensionThreePeriods)[1],1] <- dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[(2*dim(dataCurrentStateCCAndRotatedThirdDimension)[1]+1):dim(dataCurrentStateCCAndRotatedThirdDimensionThreePeriods)[1],1] + 2*pi*meanRadiusOfCells
  orderSaveThreePeriods <- order(dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[,1])
  dataCurrentStateCCAndRotatedThirdDimensionThreePeriods <- dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[orderSaveThreePeriods,]
  dataCurrentStateCCAndRotatedThirdDimensionThreePeriods <- dataCurrentStateCCAndRotatedThirdDimensionThreePeriods[floor((dim(dataCurrentStateCCAndRotatedThirdDimension)[1])/4*3):floor((dim(dataCurrentStateCCAndRotatedThirdDimensionThreePeriods)[1])/12*9),]
  dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods <- rbind(dataExtrapolatedStateCCAndRotatedThirdDimension,dataExtrapolatedStateCCAndRotatedThirdDimension,dataExtrapolatedStateCCAndRotatedThirdDimension)
  dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[1:dim(dataExtrapolatedStateCCAndRotatedThirdDimension)[1],1] <- dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[1:dim(dataExtrapolatedStateCCAndRotatedThirdDimension)[1],1] - 2*pi*meanRadiusOfCells
  dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[(2*dim(dataExtrapolatedStateCCAndRotatedThirdDimension)[1]+1):dim(dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods)[1],1] <- dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[(2*dim(dataExtrapolatedStateCCAndRotatedThirdDimension)[1]+1):dim(dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods)[1],1] + 2*pi*meanRadiusOfCells
  dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods <- dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[orderSaveThreePeriods,]
  dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods <- dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods[floor((dim(dataExtrapolatedStateCCAndRotatedThirdDimension)[1])/4*3):floor((dim(dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods)[1])/12*9),]

  dataDisplacementOfXiDirect <- dataExtrapolatedStateCCAndRotatedThirdDimensionThreePeriods - dataCurrentStateCCAndRotatedThirdDimensionThreePeriods
  dataDisplacementGrid <- as.data.frame(matrix(0L, dim(gridFinal)[1], dim(gridFinal)[2]))
  for (i in 1:dim(gridFinal)[1]){
    currentGridPoint <- gridFinal[i,]
    euclideanDistanceVector <- apply(dataCurrentStateCCAndRotatedThirdDimensionThreePeriods, 1, euclideanDistance, y = currentGridPoint)
    gaussianKernel <- exp(-euclideanDistanceVector^2/(2*sigma^2))
    dataDisplacementGridPC1 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,1])
    dataDisplacementGridPC2 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,2])

    dataDisplacementGrid[i,] <- cbind(dataDisplacementGridPC1, dataDisplacementGridPC2)
  }

  colnames(gridFinal) <- c('gridX', 'gridY')
  colnames(dataDisplacementGrid) <- c('displacementX', 'displacementY')

  dataList@velocityData$grids[[paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = '')]] <- cbind(gridFinal, dataDisplacementGrid)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  if (boolPlotResults){
    plotParameters <- list()
    plotParameters$colorPaletteCellCyclePhasesGeneral <- c('#ac4343', '#466caf', '#df8b3f', '#63b558', '#e8d760', '#61c5c7', '#f04ddf', '#a555d4')
    plotParameters$plotLabelTextSize <- 8
    plotParameters$plotDotSize <- 0.1
    plotParameters$plotLineWidth <- 0.2
    plotParameters$fontFamily <- 'Helvetica'
    plotParameters$fontSize <- 8
    plotParameters$colorPaletteCellCyclePhasesGeneral <- plotParameters$colorPaletteCellCyclePhasesGeneral[1:length(levels(dataList@cellInfo$ccPhase))]
    names(plotParameters$colorPaletteCellCyclePhasesGeneral) <- levels(dataList@cellInfo$ccPhase)

    ccPhaseTranscitionPercentages <- cumsum(c(0, dataList@datasetInfo$ccDurationG1, dataList@datasetInfo$ccDurationS, dataList@datasetInfo$ccDurationG2, dataList@datasetInfo$ccDurationM))/dataList@datasetInfo$ccDurationTotal
    tickBreaksHelp <- (ccPhaseTranscitionPercentages[1:(length(ccPhaseTranscitionPercentages)-1)]+ccPhaseTranscitionPercentages[2:length(ccPhaseTranscitionPercentages)])/2
    tickBreaks <- rbind(ccPhaseTranscitionPercentages[1:(length(ccPhaseTranscitionPercentages)-1)], tickBreaksHelp)
    tickBreaks <- c(tickBreaks)
    tickBreaks <- append(tickBreaks, ccPhaseTranscitionPercentages[length(ccPhaseTranscitionPercentages)])
    tickColour <- rbind(rep('#000000', length(ccPhaseTranscitionPercentages)-1), rep('#ffffff', length(ccPhaseTranscitionPercentages)-1))
    tickColour <- c(tickColour)
    tickColour <- as.character(append(tickColour, tickColour[1]))
    tickLabels <- rbind(rep('', length(ccPhaseTranscitionPercentages)-1), c('G1','S','G2','M'))
    tickLabels <- c(tickLabels)
    tickLabels <- append(tickLabels, tickLabels[1])

    tickInfo <- data.frame(position = tickBreaks, color = tickColour, label = tickLabels)

    meanRadiusScaleFactor <- max(dataList@velocityData$grids[[paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = '')]][,'gridX'])/(2*pi)

    thirdDimensionToDisplay <- as.numeric(substr(paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = ''), start = 11, stop = 11))
    dataToUse <- data.frame(ccCircleDistance = dataList@cellInfo$ccPercentage*2*pi*meanRadiusScaleFactor, thirdDimension = dataList@transformedData$dc$data[thirdDimensionToDisplay,], ccPhase = dataList@cellInfo$ccPhase)

    scaleFactorArrow <- 1/max(sqrt(dataList@velocityData$grids[[paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = '')]][,'displacementX']^2+dataList@velocityData$grids[[paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = '')]][,'displacementY']^2))*5
    arrowLengthThreshold <- 0.2
    dataToUseGrid <- filter(dataList@velocityData$grids[[paste('velAlongDC', whichThirdDimension, '_sigma_', sigma, sep = '')]], sqrt(displacementX^2+displacementY^2)*scaleFactorArrow>arrowLengthThreshold)
    dataToUseGrid[,c('displacementX', 'displacementY')] <- dataToUseGrid[,c('displacementX', 'displacementY')]*scaleFactorArrow

    plotDC3Velocity <- ggplot(data = dataToUse)+
      theme_gray(base_size = plotParameters$plotLabelTextSize)+
      theme(text=element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
            axis.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
            axis.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
            legend.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
            legend.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
            plot.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize))+
      geom_point(aes(x = ccCircleDistance, y = thirdDimension, color = ccPhase), size = plotParameters$plotDotSize*5)+
      ylim(c(-10,10))+
      scale_color_manual(values = plotParameters$colorPaletteCellCyclePhasesGeneral,
                         labels = levels(dataList@cellInfo$ccPhase))+
      xlab('cell cycle progression')+
      ylab(paste('DC', thirdDimensionToDisplay, sep = ''))+
      guides(color = guide_legend(override.aes = list(size = 1)))+
      geom_vline(xintercept = cumsum(c(0, dataList@datasetInfo$ccDurationG1, dataList@datasetInfo$ccDurationS, dataList@datasetInfo$ccDurationG2, dataList@datasetInfo$ccDurationM))/dataList@datasetInfo$ccDurationTotal*2*pi*meanRadiusScaleFactor, color = 'black', size = plotParameters$plotLineWidth*3/2, linetype = 'dotted')+
      theme(legend.position = 'top',
            legend.title = element_blank(),
            legend.box.margin = margin(0,0,-5,0),
            legend.margin = margin(0,0,0,0),
            legend.key.size = unit(0.65,"line"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_line(colour = as.character(tickInfo$color)),
            panel.grid.major.x = element_blank(),
            axis.text.y=element_text(angle=90, hjust=0.5))+
      guides(color = guide_legend(override.aes = list(size = 0.8)))+
      scale_x_continuous(breaks = tickInfo$position*2*pi*meanRadiusScaleFactor,
                         labels = tickInfo$label,
                         minor_breaks = seq(from = (0+10^(-10))*2*pi*meanRadiusScaleFactor, to = (1+10^(-10))*2*pi*meanRadiusScaleFactor, length.out = 11))+
      geom_point(data = dataToUseGrid,
                 aes(x = gridX, y = gridY),
                 color = 'black',
                 size = plotParameters$plotDotSize/20)+
      geom_segment(data = dataToUseGrid,
                   aes(x = gridX, xend = gridX+displacementX, y = gridY, yend = gridY+displacementY),
                   color = 'black',
                   size = plotParameters$plotLineWidth*1.35,
                   arrow=arrow(length = unit(plotParameters$plotLineWidth/3*1.35, "cm")))

    grid.arrange(plotDC3Velocity)

  }

  return(dataList)
}
getIdealizedCell <- function(dataList,
                             numberOfAnchorsForIdealizedCell = 10,
                             boolReturnAnchorPointPolarCoordinates = FALSE,
                             numberOfPointsInSpline = 10^4,
                             filterVector = NULL){
  if (is.null(filterVector)){
    filterVector <- colnames(dataList@transformedData$dc$data)
  }
  data <- t(dataList@transformedData$dc$data[c(1,2),filterVector])
  indexOfCells <- which(colnames(dataList@transformedData$dc$data)%in%filterVector)
  data <- data[order(dataList@cellInfo$ccAngle[indexOfCells], decreasing = TRUE),]
  dataPolar <- getPolarCoordinates(data = data)
  dataAngle <- dataPolar[,'angle']
  resultList <- getBordersAndMiddlesAlongAngleVectorForSpecificNumberOfIntervals(dataAngle,
                                                                                 numberOfIntervals = numberOfAnchorsForIdealizedCell)
  {
    angleOfBorders <- resultList[[1]]
    angleOfBorders <- append(angleOfBorders,angleOfBorders[1]-2*pi)
    angleOfMiddles <- resultList[[2]]
  }

  radiusValues <- rep(0,numberOfAnchorsForIdealizedCell)
  for (i in 1:numberOfAnchorsForIdealizedCell){
    radiusValues[i] <- mean(dataPolar[(dataAngle<angleOfBorders[i])&(dataAngle>angleOfBorders[i+1]),'radius'])
  }

  angleOfMiddles <- as.numeric(append(angleOfMiddles, angleOfMiddles[1]-2*pi))
  radiusValues <- as.numeric(append(radiusValues, radiusValues[1]))
  splineIdealizedCell <- spline(x = angleOfMiddles, y =  radiusValues, n = numberOfPointsInSpline, method = 'periodic')
  splineXValues <- splineIdealizedCell[[1]]
  splineIdealizedCell <- splineIdealizedCell[[2]]

  numberOfPointsToDistort <- sum(splineXValues<0)
  splineXValues <- c(splineXValues[(numberOfPointsToDistort+1):(length(splineXValues)-1)], splineXValues[1:numberOfPointsToDistort]+2*pi)
  splineXValues <- rev(splineXValues)
  splineIdealizedCell <- append(splineIdealizedCell[(numberOfPointsToDistort+1):(length(splineIdealizedCell)-1)], splineIdealizedCell[1:numberOfPointsToDistort])
  splineIdealizedCell <- rev(splineIdealizedCell)
  idealizedCellPolar <- as.data.frame(cbind(splineXValues, splineIdealizedCell))
  colnames(idealizedCellPolar) <- c('angle', 'radius')

  idealizedCellCart <- getCartesianCoordinates(dataPolar = idealizedCellPolar)
  colnames(idealizedCellCart) <- c('x', 'y')

  idealizedCell <- as.data.frame(cbind(idealizedCellCart, idealizedCellPolar))

  if (boolReturnAnchorPointPolarCoordinates){
    return(list(idealizedCell = idealizedCell, anchorPointsPolar = data.frame(angle = angleOfMiddles, radius = radiusValues)))
  }else{
    return(idealizedCell)
  }
}
getBordersAndMiddlesAlongAngleVectorForSpecificNumberOfIntervals <- function(dataAngle,
                                                                             numberOfIntervals){

  numberOfCellsPerInterval <- rep(floor(length(dataAngle)/numberOfIntervals),numberOfIntervals)
  if ((length(dataAngle)%%numberOfIntervals)>0){
    additionalIntervals <- length(dataAngle)%%numberOfIntervals
    numberOfCellsPerInterval[1:additionalIntervals] <- numberOfCellsPerInterval[1:additionalIntervals]+1
  }
  angleOfIntervalBorder <- cumsum(numberOfCellsPerInterval[1:(length(numberOfCellsPerInterval)-1)])
  angleOfIntervalBorder <- (dataAngle[angleOfIntervalBorder] + dataAngle[angleOfIntervalBorder+1])/2
  angleOfIntervalBorder <- as.numeric(append((dataAngle[1]+dataAngle[length(dataAngle)]+2*pi)/2, angleOfIntervalBorder))

  angleOfIntervalMiddle <- (angleOfIntervalBorder[1:(length(angleOfIntervalBorder)-1)] + angleOfIntervalBorder[2:length(angleOfIntervalBorder)])/2
  angleOfIntervalMiddle <- append(angleOfIntervalMiddle, (angleOfIntervalBorder[length(angleOfIntervalBorder)]-2*pi+angleOfIntervalBorder[1])/2)

  return(list(angleOfIntervalBorder, angleOfIntervalMiddle))
}
getVelocityGridForIdealizedCell <- function(dataList,
                                            numberOfAnchorsForIdealizedCell = 10,
                                            factorForSigmaDecrease = 1.5){

  resultList <- getIdealizedCell(dataList = dataList,
                                 numberOfAnchorsForIdealizedCell = numberOfAnchorsForIdealizedCell,
                                 boolReturnAnchorPointPolarCoordinates = TRUE)
  idealizedCell <- resultList[['idealizedCell']]
  anchorPointsPolar <- resultList[['anchorPointsPolar']]
  anchorPointsPolar <- anchorPointsPolar[1:(dim(anchorPointsPolar)[1]-1),]
  rm(resultList)
  anchorPointsCartesian <- getCartesianCoordinates(anchorPointsPolar)
  originalStateInRotatedSpace <- t(dataList@transformedData$dc$data[c(1,2),dataList@cellInfo$ccPositionIndex])
  extrapolatedStateInRotatedStateSpace <- t(dataList@velocityData$extrapolatedStateInDCSpace[c(1,2), rownames(originalStateInRotatedSpace)])
  gridFinal <- anchorPointsCartesian
  sigma <- (2*pi*median(dataList@cellInfo$ccRadius))/(numberOfAnchorsForIdealizedCell*factorForSigmaDecrease)

  euclideanDistance <- function(x,y){
    return(sqrt(sum((x-y)^2)))
  }

  dataDisplacementOfXiDirect <- extrapolatedStateInRotatedStateSpace - originalStateInRotatedSpace
  dataDisplacementGrid <- as.data.frame(matrix(0L, dim(gridFinal)[1], dim(gridFinal)[2]))
  for (i in 1:dim(gridFinal)[1]){
    currentGridPoint <- gridFinal[i,]
    euclideanDistanceVector <- apply(originalStateInRotatedSpace, 1, euclideanDistance, y = currentGridPoint)
    gaussianKernel <- exp(-euclideanDistanceVector^2/(2*sigma^2))
    dataDisplacementGridPC1 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,1])
    dataDisplacementGridPC2 <- sum(gaussianKernel*dataDisplacementOfXiDirect[,2])

    dataDisplacementGrid[i,] <- cbind(dataDisplacementGridPC1, dataDisplacementGridPC2)
  }

  colnames(gridFinal) <- c('gridX', 'gridY')
  colnames(dataDisplacementGrid) <- c('displacementX', 'displacementY')

  return(cbind(gridFinal, dataDisplacementGrid))
}
#'
#'
#' Calculate the velocity along an averaged trajectory.
#'
#' 'getVelocityOnAveragedTrajectory' places an average trajectory through the data cloud and calculates the average velocity for each anchor point as described in La Manno et al. 2018 ('RNA velocity of single cells').
#'
#' A trajectory is placed through the data cloud as a spline in polar coordinates with a number of anchor points specified by the user. For each anchor point, the RNA velocity is calculated by using a Gaussian kernel as explained in La Manno et al. 2018 ('RNA velocity of single cells').
#' @param dataList A Revelio object that contains information about the extrapolated state of cells.
#' @param numberOfAnchors The number of anchor points used for the spline of the averaged trajectory. This also defines the points at which RNA velocity will be calculated.
#' @return Returns the same Revelio object given as input.
#'
#' @export
plotVelocityOnAveragedTrajectory <- function(dataList,
                                             numberOfAnchors = 10){
  idealizedCellAnchor <- getIdealizedCell(dataList = dataList,
                                          numberOfAnchorsForIdealizedCell = numberOfAnchors)

  plotParameters <- list()
  plotParameters$colorPaletteCellCyclePhasesGeneral <- c('#ac4343', '#466caf', '#df8b3f', '#63b558', '#e8d760', '#61c5c7', '#f04ddf', '#a555d4')
  plotParameters$plotLabelTextSize <- 8
  plotParameters$plotDotSize <- 0.1
  plotParameters$plotLineWidth <- 0.2
  plotParameters$fontFamily <- 'Helvetica'
  plotParameters$fontSize <- 8
  plotParameters$colorPaletteCellCyclePhasesGeneral <- plotParameters$colorPaletteCellCyclePhasesGeneral[1:length(levels(dataList@cellInfo$ccPhase))]
  names(plotParameters$colorPaletteCellCyclePhasesGeneral) <- levels(dataList@cellInfo$ccPhase)
  plotParameters$colorOptimalTraj <- c('#9428bf')

  plotBoundary <- max(dataList@transformedData$dc$data[c('DC1', 'DC2'),])*0.78

  ccPhaseBorderPathPolar <- data.frame(angle = rep(0,7), radius = rep(0,7))
  ccPhaseBorderPathPolar[c(1,3,5,7),1] <- 2*pi*(1-cumsum(c(dataList@datasetInfo$ccDurationG1, dataList@datasetInfo$ccDurationS, dataList@datasetInfo$ccDurationG2, dataList@datasetInfo$ccDurationM))/dataList@datasetInfo$ccDurationTotal)
  ccPhaseBorderPathPolar[c(1,3,5,7),2] <- 50
  ccPhaseBorderPathCartesian <- ccPhaseBorderPathPolar
  ccPhaseBorderPathCartesian[,1] <- ccPhaseBorderPathPolar[,2]*cos(ccPhaseBorderPathPolar[,1])
  ccPhaseBorderPathCartesian[,2] <- ccPhaseBorderPathPolar[,2]*sin(ccPhaseBorderPathPolar[,1])
  colnames(ccPhaseBorderPathCartesian) <- c('xValue', 'yValue')

  labelPositionHelp <- append(2*pi, ccPhaseBorderPathPolar[c(1,3,5,7),1])
  labelPositionHelp <- (labelPositionHelp[1:(length(labelPositionHelp)-1)]-labelPositionHelp[2:length(labelPositionHelp)])/2+labelPositionHelp[2:length(labelPositionHelp)]
  labelRadiusHelp <- labelPositionHelp
  labelRadiusHelp[(labelRadiusHelp>pi/4)&(labelRadiusHelp<3*pi/4)] <- labelRadiusHelp[(labelRadiusHelp>pi/4)&(labelRadiusHelp<3*pi/4)]-pi/2
  labelRadiusHelp[(labelRadiusHelp>5*pi/4)&(labelRadiusHelp<7*pi/4)] <- labelRadiusHelp[(labelRadiusHelp>5*pi/4)&(labelRadiusHelp<7*pi/4)]-pi/2

  labelPositionPolar <- data.frame(angle = labelPositionHelp, radius = sqrt((plotBoundary*tan(labelRadiusHelp))^2+plotBoundary^2), label = c('G1', 'S', 'G2', 'M'))
  labelPositionCartesian <- labelPositionPolar
  labelPositionCartesian[,1] <- labelPositionPolar[,2]*cos(labelPositionPolar[,1])
  labelPositionCartesian[,2] <- labelPositionPolar[,2]*sin(labelPositionPolar[,1])
  colnames(labelPositionCartesian) <- c('xValue', 'yValue', 'label')

  plotDC1DC2 <- ggplot(data = cbind(as.data.frame(t(dataList@transformedData$dc$data)), ccPhase = dataList@cellInfo$ccPhase))+
    theme_gray(base_size = plotParameters$plotLabelTextSize)+
    theme(text=element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
          axis.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
          axis.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
          legend.text = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize-1),
          legend.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize),
          plot.title = element_text(family=plotParameters$fontFamily, size=plotParameters$fontSize))+
    geom_point(aes(x = DC1, y = DC2, color = ccPhase), size = plotParameters$plotDotSize*2)+
    geom_path(data = ccPhaseBorderPathCartesian, aes(x = xValue, y = yValue), color = 'black', size = plotParameters$plotLineWidth*3, linetype = 'dotted')+
    coord_cartesian(xlim = c(-plotBoundary, plotBoundary),ylim = c(-plotBoundary, plotBoundary))+
    scale_color_manual(values = plotParameters$colorPaletteCellCyclePhasesGeneral,
                       labels = levels(dataList@cellInfo$ccPhase))+
    theme(legend.position = 'right',
          legend.title = element_blank(),
          aspect.ratio=1,
          legend.box.margin = margin(0,0,0,0),
          legend.margin = margin(0,0,0,0),
          legend.key.size = unit(0.35,"line"),
          axis.text.y=element_text(angle=90, hjust=0.5))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    geom_text(data = labelPositionCartesian, aes (x = xValue, y = yValue, label = label), size = plotParameters$plotLabelTextSize/2)


  plotIdealizedCell <- plotDC1DC2+
    geom_path(data = idealizedCellAnchor, aes (x = x, y = y), size = plotParameters$plotLineWidth*4, color = plotParameters$colorOptimalTraj)

  velocityGridIdealizedCell <- getVelocityGridForIdealizedCell(dataList = dataList,
                                                               numberOfAnchorsForIdealizedCell = as.numeric(numberOfAnchors),
                                                               factorForSigmaDecrease = 1.5)
  scaleFactorArrow <- 1/max(sqrt(velocityGridIdealizedCell[,'displacementX']^2+velocityGridIdealizedCell[,'displacementY']^2))*3.5
  dataToUse <- velocityGridIdealizedCell
  dataToUse[,c('displacementX', 'displacementY')] <- dataToUse[,c('displacementX', 'displacementY')]*scaleFactorArrow
  plotIdealizedCellWithVel <- plotIdealizedCell+
    geom_segment(data = dataToUse,
                 aes(x = gridX, y = gridY, xend = gridX+displacementX, yend = gridY+displacementY),color = 'black',
                 size = plotParameters$plotLineWidth*3,
                 arrow=arrow(length = unit(plotParameters$plotLineWidth*1, "cm")))
  # theme(legend.margin = margin(0, 10, 0, -10),
  #       legend.box.margin = margin(0, 10, 0, -10))

  grid.arrange(plotIdealizedCellWithVel)

  return(dataList)
}
