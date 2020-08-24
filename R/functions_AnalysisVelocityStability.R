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
  timeStep <- dataList@datasetInfo$velTimeStep
  
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
getVelocityGridsCC <- function(dataList){
  
  for (i in dataList@datasetInfo$velCCGridSigma){
    startTime <- Sys.time()
    cat(paste(Sys.time(), ': calculating cell cycle grid for sigma = ', i, ': ', sep = ''))
    dataList <- calculateVelocityCCGridDisplacement(dataList = dataList,
                                                    sigma = i)
    cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  }
  
  return(dataList)
}
calculateVelocityCCGridDisplacement <- function(dataList,
                                                sigma = 1){
  originalStateInRotatedSpace <- t(dataList@transformedData$dc$data[c(1,2),])
  extrapolatedStateInRotatedStateSpace <- t(dataList@velocityData$extrapolatedStateInDCSpace[c(1,2), rownames(originalStateInRotatedSpace)])
  gridBoundary <- 2.5*max(sd(originalStateInRotatedSpace[,1]), sd(originalStateInRotatedSpace[,2]))
  gridX <- seq(from = -gridBoundary, to = gridBoundary, length.out = dataList@datasetInfo$velCCGridAmountOfPoints)
  gridY <- seq(from = -gridBoundary, to = gridBoundary, length.out = dataList@datasetInfo$velCCGridAmountOfPoints)
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
  
  return(dataList)
}
getVelocityGridsAlongThirdDimensions <- function(dataList){
  
  for (i in dataList@datasetInfo$velWhichThirdDimensionsToGetVelocityFor){
    startTime <- Sys.time()
    cat(paste(Sys.time(), ': calculating velocity grid along DC', i, ': ', sep = ''))
    dataList <- calculateVelocityGridAlongThirdDimension(dataList = dataList,
                                                         whichThirdDimension = i,
                                                         sigma = dataList@datasetInfo$velCCGridSigma[length(dataList@datasetInfo$velCCGridSigma)])
    cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  }
  
  return(dataList)
}
calculateVelocityGridAlongThirdDimension <- function(dataList,
                                                     whichThirdDimension = 3,
                                                     sigma = 1){
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
  
  gridX <- seq(from = 0, to = 2*pi*meanRadiusOfCells, length.out = dataList@datasetInfo$velCCGridAmountOfPoints)
  gridY <- seq(from = -pi*meanRadiusOfCells, to = pi*meanRadiusOfCells, length.out = dataList@datasetInfo$velCCGridAmountOfPoints)
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
  
  return(dataList)
}
getStabilityIndexAlongCC <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': doing stability analysis: ', sep = ''))
  
  numberOfIntervals <- dataList@datasetInfo$stabilityIndexNumberOfIntervals
  correlationType <- dataList@datasetInfo$stabilityIndexCorrelationType
  pValueForCorrelation <- dataList@datasetInfo$stabilityIndexpValueForCorrelation
  
  rawData <- dataList@DGEs$countData
  dataToDownsample <- as.matrix(rawData)
  baseCountToDownsampleTo <- (quantile(colSums(dataToDownsample)))[2]
  downsampleFractions <- baseCountToDownsampleTo/colSums(dataToDownsample)
  downsampleFractions <- pmin(downsampleFractions, 1)
  dataDownsampled <- downsampleMatrix(dataToDownsample, prop = downsampleFractions, bycol = TRUE)
  countUMIPerCell <- colSums(dataDownsampled)                                   #gives total UMI count per cell
  dataDownsampledLogOfFractions <- sweep(dataDownsampled,2,countUMIPerCell,'/')                          #divides each column by total UMI count for that cell
  dataDownsampledLogOfFractions <- dataDownsampledLogOfFractions*baseCountToDownsampleTo                                                     #multiplies every entry by 10^4 as a scaling factor
  dataDownsampledLogOfFractions <- log(dataDownsampledLogOfFractions+1)                                                   #transforms to logarithmic scale (natural logarithm)
  
  #cell-cell and gene-gene correlations
  #dataRawMatrix <- as.matrix(dataList@DGEs$logOfFractionsData[dataList@geneInfo$pcaGenes, dataList@cellInfo$ccPositionIndex])
  dataRawMatrix <- as.matrix(dataDownsampledLogOfFractions[dataList@geneInfo$pcaGenes, dataList@cellInfo$ccPositionIndex])
  #dataRawMatrix <- as.matrix(dataList@DGEs$logOfFractionsData[genesCCDC12, dataList@cellInfo$ccPositionIndex])
  
  cellCorrelationMatrix <- rep(0, numberOfIntervals)
  geneCorrelationMatrix <- rep(0, numberOfIntervals)
  numberOfZeroEntries <- rep(0, numberOfIntervals)
  
  numberOfCellsPerInterval <- rep(floor(dim(dataRawMatrix)[2]/numberOfIntervals),numberOfIntervals)
  if ((dim(dataRawMatrix)[2]%%numberOfIntervals)>0){
    additionalIntervals <- dim(dataRawMatrix)[2]%%numberOfIntervals
    numberOfCellsPerInterval[1:additionalIntervals] <- numberOfCellsPerInterval[1:additionalIntervals]+1
  }
  indexOfIntervalBorder <- append(0,cumsum(numberOfCellsPerInterval))
  dataList@stabilityIndex[[paste(correlationType, 'CorrelationPValue', pValueForCorrelation, 'And', numberOfIntervals, 'Intervals_DataMatrices', sep = '')]] <- list()
  
  for (i in 1:numberOfIntervals){
    currentCells <- (dataList@cellInfo$cellID[dataList@cellInfo$ccPositionIndex])[(indexOfIntervalBorder[i]+1):indexOfIntervalBorder[i+1]]
    
    if (length(currentCells)>2){
      # startTimeNew <- Sys.time()
      # 
      # cellCorVectorCurrent <- rep(0, length(currentCells)*(length(currentCells)-1)/2)
      # counterCellVector <- 0
      # if (pValueForCorrelation>0){
      #   for (j in 1:(length(currentCells)-2)){
      #     counterCellVector <- counterCellVector + 1
      #     corSave <- apply(dataRawMatrix[,currentCells[(j+1):length(currentCells)]], 2, cor.test, y=as.vector(dataRawMatrix[,currentCells[j]]), method = correlationType)
      #     corSaveCor <- as.numeric(lapply(corSave, '[[', 4))
      #     corSavePValue <- as.numeric(lapply(corSave, '[[', 3))
      #     corSaveCor[is.na(corSaveCor)] <- 0
      #     corSavePValue[is.na(corSavePValue)] <- 1
      #     corSave <- corSaveCor[corSavePValue<pValueForCorrelation]
      #     if (length(corSave)>0){
      #       cellCorVectorCurrent[counterCellVector:(counterCellVector+length(corSave)-1)] <- corSave
      #       counterCellVector <- counterCellVector+length(corSave)-1
      #     }else{
      #       counterCellVector <- counterCellVector-1
      #     }
      #   }
      #   counterCellVector <- counterCellVector + 1
      #   corSave <- cor.test(dataRawMatrix[,currentCells[length(currentCells)]], as.vector(dataRawMatrix[,currentCells[length(currentCells)-1]]), method = correlationType)
      #   corSaveCor <- corSave[[4]]
      #   corSavePValue <- corSave[[3]]
      #   corSaveCor[is.na(corSaveCor)] <- 0
      #   corSavePValue[is.na(corSavePValue)] <- 1
      #   corSave <- corSaveCor[corSavePValue<pValueForCorrelation]
      #   if (length(corSave)>0){
      #     cellCorVectorCurrent[counterCellVector:(counterCellVector+length(corSave)-1)] <- corSave
      #     counterCellVector <- counterCellVector+length(corSave)-1
      #   }else{
      #     counterCellVector <- counterCellVector-1
      #   }
      # }else{
      #   for (j in 1:(length(currentCells)-2)){
      #     counterCellVector <- counterCellVector + 1
      #     corSave <- apply(dataRawMatrix[,currentCells[(j+1):length(currentCells)]], 2, cor, y=as.vector(dataRawMatrix[,currentCells[j]]), method = correlationType)
      #     corSave[is.na(corSave)] <- 0
      #     cellCorVectorCurrent[counterCellVector:(counterCellVector+length(corSave)-1)] <- corSave
      #     counterCellVector <- counterCellVector+length(corSave)-1
      #   }
      #   counterCellVector <- counterCellVector + 1
      #   corSave <- cor(dataRawMatrix[,currentCells[length(currentCells)]], as.vector(dataRawMatrix[,currentCells[length(currentCells)-1]]), method = correlationType)
      #   corSave[is.na(corSave)] <- 0
      #   cellCorVectorCurrent[counterCellVector:(counterCellVector+length(corSave)-1)] <- corSave
      #   counterCellVector <- counterCellVector+length(corSave)-1
      # }
      # 
      # geneCorVectorCurrent <- rep(0, dim(dataRawMatrix)[1]*(dim(dataRawMatrix)[1]-1)/2)
      # counterGeneVector <- 0
      # if (pValueForCorrelation>0){
      #   for (j in 1:(dim(dataRawMatrix)[1]-2)){
      #     counterGeneVector <- counterGeneVector + 1
      #     corSave <- apply(dataRawMatrix[(j+1):dim(dataRawMatrix)[1],currentCells], 1, cor.test, y=as.vector(dataRawMatrix[j,currentCells]), method = correlationType)
      #     corSaveCor <- as.numeric(lapply(corSave, '[[', 4))
      #     corSavePValue <- as.numeric(lapply(corSave, '[[', 3))
      #     corSaveCor[is.na(corSaveCor)] <- 0
      #     corSavePValue[is.na(corSavePValue)] <- 1
      #     corSave <- corSaveCor[corSavePValue<pValueForCorrelation]
      #     if (length(corSave)>0){
      #       geneCorVectorCurrent[counterGeneVector:(counterGeneVector+length(corSave)-1)] <- corSave
      #       counterGeneVector <- counterGeneVector+length(corSave)-1
      #     }else{
      #       counterGeneVector <- counterGeneVector-1
      #     }
      #   }
      #   counterGeneVector <- counterGeneVector + 1
      #   corSave <- cor.test(dataRawMatrix[dim(dataRawMatrix)[1],currentCells], as.vector(dataRawMatrix[dim(dataRawMatrix)[1]-1,currentCells]), method = correlationType)
      #   corSaveCor <- corSave[[4]]
      #   corSavePValue <- corSave[[3]]
      #   corSaveCor[is.na(corSaveCor)] <- 0
      #   corSavePValue[is.na(corSavePValue)] <- 1
      #   corSave <- corSaveCor[corSavePValue<pValueForCorrelation]
      #   if (length(corSave)>0){
      #     geneCorVectorCurrent[counterGeneVector:(counterGeneVector+length(corSave)-1)] <- corSave
      #     counterGeneVector <- counterGeneVector+length(corSave)-1
      #   }else{
      #     counterGeneVector <- counterGeneVector-1
      #   }
      # }else{
      #   for (j in 1:(dim(dataRawMatrix)[1]-2)){
      #     counterGeneVector <- counterGeneVector + 1
      #     corSave <- apply(dataRawMatrix[(j+1):dim(dataRawMatrix)[1],currentCells], 1, cor, y=as.vector(dataRawMatrix[j,currentCells]), method = correlationType)
      #     corSave[is.na(corSave)] <- 0
      #     geneCorVectorCurrent[counterGeneVector:(counterGeneVector+length(corSave)-1)] <- corSave
      #     counterGeneVector <- counterGeneVector+length(corSave)-1
      #   }
      #   counterGeneVector <- counterGeneVector + 1
      #   corSave <- cor(dataRawMatrix[dim(dataRawMatrix)[1],currentCells], as.vector(dataRawMatrix[dim(dataRawMatrix)[1]-1,currentCells]), method = correlationType)
      #   corSave[is.na(corSave)] <- 0
      #   geneCorVectorCurrent[counterGeneVector:(counterGeneVector+length(corSave)-1)] <- corSave
      #   counterGeneVector <- counterGeneVector+length(corSave)-1
      # }
      
      corTestValuesCell <- corr.test(dataRawMatrix[,currentCells], adjust = "none", ci = FALSE, method = correlationType)
      corTestValuesGene <- corr.test(t(dataRawMatrix[,currentCells]), adjust = "none", ci = FALSE, method = correlationType)
      
      cellCorVectorCurrent <- corTestValuesCell$r
      cellPValueVectorCurrent <- corTestValuesCell$p
      geneCorVectorCurrent <- corTestValuesGene$r
      genePValueVectorCurrent <- corTestValuesGene$p
      
      cellCorVectorCurrent[is.na(cellCorVectorCurrent)] <- 0
      cellPValueVectorCurrent[is.na(cellPValueVectorCurrent)] <- 1
      geneCorVectorCurrent[is.na(geneCorVectorCurrent)] <- 0
      genePValueVectorCurrent[is.na(genePValueVectorCurrent)] <- 1
      
      cellCorVectorCurrent[cellPValueVectorCurrent>=pValueForCorrelation] <- 0
      geneCorVectorCurrent[genePValueVectorCurrent>=pValueForCorrelation] <- 0
      
      cellCorrelationValues <- cellCorVectorCurrent[upper.tri(cellCorVectorCurrent, diag = FALSE)]
      geneCorrelationValues <- geneCorVectorCurrent[upper.tri(geneCorVectorCurrent, diag = FALSE)]
      
      cellCorrelationMatrix[i] <- mean(cellCorrelationValues[cellCorrelationValues!=0])
      geneCorrelationMatrix[i] <- mean(abs(geneCorrelationValues[geneCorrelationValues!=0]))
      numberOfZeroEntries[i] <- sum(geneCorrelationValues==0)
      
      if (!dataList@datasetInfo$pcaGenes=='allGenes'){
        dataList@stabilityIndex[[paste(correlationType, 'CorrelationPValue', pValueForCorrelation, 'And', numberOfIntervals, 'Intervals_DataMatrices', sep = '')]][[as.character(i)]] <- list(cellCorMatrix = cellCorVectorCurrent,
                                                                                                                                                                                              cellPValueMatrix = cellPValueVectorCurrent,
                                                                                                                                                                                              geneCorMatrix = geneCorVectorCurrent,
                                                                                                                                                                                              genePValueMatrix = genePValueVectorCurrent)
      }
    }else{
      cellCorrelationMatrix[i] <- 0
      geneCorrelationMatrix[i] <- 0
      numberOfZeroEntries[i] <- 0
    }
  }
  
  cellCorrelationMatrix <- c((cellCorrelationMatrix[1]+cellCorrelationMatrix[length(cellCorrelationMatrix)])/2, cellCorrelationMatrix, (cellCorrelationMatrix[1]+cellCorrelationMatrix[length(cellCorrelationMatrix)])/2)
  geneCorrelationMatrix <- c((geneCorrelationMatrix[1]+geneCorrelationMatrix[length(geneCorrelationMatrix)])/2, geneCorrelationMatrix, (geneCorrelationMatrix[1]+geneCorrelationMatrix[length(geneCorrelationMatrix)])/2)
  dataList@stabilityIndex[[paste(correlationType, 'CorrelationPValue', pValueForCorrelation, 'And', numberOfIntervals, 'Intervals', sep = '')]] <- data.frame(ccPercentage = c(0, (1:numberOfIntervals)/numberOfIntervals-1/numberOfIntervals/2, 1),
                                                                                                                                                              cellcell = cellCorrelationMatrix,
                                                                                                                                                              genegene = geneCorrelationMatrix,
                                                                                                                                                              stabilityIndex = geneCorrelationMatrix/cellCorrelationMatrix)
  
  
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
