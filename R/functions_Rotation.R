#' @export
getOptimalRotations <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating optimal rotations: ', sep = ''))

  dataList@transformedData$dc <- calculateOptimalRotation(pcaData = dataList@transformedData$pca$data,
                                                          pcaWeights = dataList@transformedData$pca$weights,
                                                          pcaIsPCAssociatedWithCC = dataList@transformedData$pca$pcProperties$isComponentAssociatedWithCC,
                                                          ccPhaseInformation = dataList@cellInfo$ccPhase)

  cellCycleScore <- getCellCycleScoreForPCA(dataList@transformedData$dc$data, dataList@cellInfo[,'ccPhase'])
  boolOutliers <- rep(FALSE, length(cellCycleScore))
  boolOutliers[1:(min(length(cellCycleScore),100))] <- calculateOutliersInVector(data = as.vector(cellCycleScore[1:(min(length(cellCycleScore),100))]),
                                                                                 outlierThreshold = 2)
  dataList@transformedData$dc$dcProperties <- cbind(dataList@transformedData$dc$dcProperties,
                                                    ccScore = cellCycleScore,
                                                    isComponentAssociatedWithCC = boolOutliers)

  thresholdForPCWeightToBeSignificant <- dataList@datasetInfo$thresholdsForPCWeightSignificance
  ccGenesHelp <- rep(FALSE, dim(dataList@geneInfo)[1])
  names(ccGenesHelp) <- dataList@geneInfo[,'geneID']
  ccGenesHelp[dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes][which((abs(dataList@transformedData$dc$weights[1,])>thresholdForPCWeightToBeSignificant)|(abs(dataList@transformedData$dc$weights[2,])>thresholdForPCWeightToBeSignificant))]] <- TRUE
  dataList@geneInfo <- cbind(dataList@geneInfo, ccGenes = ccGenesHelp)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
calculateOptimalRotation <- function(pcaData,
                                     pcaWeights,
                                     pcaIsPCAssociatedWithCC,
                                     ccPhaseInformation){
  projVectorDimensionsToBeRotated <- which(pcaIsPCAssociatedWithCC)
  projVectorDimensionsToBeRotated <- projVectorDimensionsToBeRotated[projVectorDimensionsToBeRotated>2]
  whichDimensionsWereRotatedBefore <- NULL
  if (length(projVectorDimensionsToBeRotated) > 0){
    sequenceOfViewingAxis <- as.data.frame(matrix(0L, nrow = 3, ncol = max(projVectorDimensionsToBeRotated)))
    rotationMatrix <- diag(1, dim(pcaData)[1])
    rotatedData <- t(pcaData)
    rotatedWeightMatrix <- t(pcaWeights)

    gridOnSphere <- getGoldenSpiral(numberOfPoints = 10^4)
    gridOnSphere <- gridOnSphere[gridOnSphere[,3]>=0,]
    for (j in projVectorDimensionsToBeRotated){
      # print(paste('PC to rotate: ', j, sep = ''))
      # start.Time <- Sys.time()
      firstViewingAxis <- getBestViewingAxisMinimizingDispersionInThirdDimension(data = rotatedData[,c(1,2,j)],
                                                                                 gridToUse = gridOnSphere,
                                                                                 ccPhaseInformation = ccPhaseInformation,
                                                                                 boolInitialAxis = FALSE)
      radiusForAxis <- 3*sort(apply(gridOnSphere,1,euclideanDistance,y=firstViewingAxis),partial=2)[2]
      gridAroundAxis <- getGridAroundAxis(viewingAxisForGrid = firstViewingAxis,
                                          numberOfPoints = 10^4,
                                          maxRadius = radiusForAxis)
      # allGridPoints <- rbind(gridOnSphere, gridAroundAxis)
      # plot3DTest <- plot_ly(allGridPoints, x = ~V1, y = ~V2, z = ~V3, marker = list(size = 1))
      # add_markers(p = plot3DTest)
      newViewingAxis <- getBestViewingAxisMinimizingDispersionInThirdDimension(data = rotatedData[,c(1,2,j)],
                                                                               gridToUse = gridAroundAxis,
                                                                               ccPhaseInformation = ccPhaseInformation,
                                                                               boolInitialAxis = TRUE,
                                                                               initialAxis = firstViewingAxis)
      counterLoops <- 0

      while ((euclideanDistance(firstViewingAxis, newViewingAxis)>(radiusForAxis/4*3))&(counterLoops <= 5)){
        counterLoops <- counterLoops+1
        firstViewingAxis <- newViewingAxis
        radiusForAxis <- 3*sort(apply(gridOnSphere,1,euclideanDistance,y=firstViewingAxis),partial=2)[2]
        gridAroundAxis <- getGridAroundAxis(viewingAxisForGrid = firstViewingAxis,
                                            numberOfPoints = 10^4,
                                            maxRadius = radiusForAxis)
        newViewingAxis <- getBestViewingAxisMinimizingDispersionInThirdDimension(data = rotatedData[,c(1,2,j)],
                                                                                 gridToUse = gridAroundAxis,
                                                                                 ccPhaseInformation = ccPhaseInformation,
                                                                                 boolInitialAxis = TRUE,
                                                                                 initialAxis = firstViewingAxis)
      }

      # end.Time <- Sys.time()
      # print(end.Time - start.Time)
      # rm(start.Time, end.Time)

      newRotationMatrix3D <- getRotationMatrix3D(viewingAxis = newViewingAxis)
      rotatedData[,c(1,2,j)] <- rotatedData[,c(1,2,j)] %*% newRotationMatrix3D
      rotatedWeightMatrix[,c(1,2,j)] <- rotatedWeightMatrix[,c(1,2,j)] %*% newRotationMatrix3D
      rotationMatrix[,c(1,2,j)] <- rotationMatrix[,c(1,2,j)] %*% newRotationMatrix3D

      sequenceOfViewingAxis[,j] <- newViewingAxis
    }
  }else{
    rotationMatrix <- diag(1, dim(pcaData)[1])
    rotatedData <- t(pcaData)
    rotatedWeightMatrix <- t(pcaWeights)
    sequenceOfViewingAxis <- as.data.frame(matrix(c(0,0,0,0,0,0,0,0,1), nrow = 3, ncol = 3))
  }
  colnames(rotatedData) <- paste('DC', 1:dim(rotatedData)[2], sep = '')
  colnames(rotatedWeightMatrix) <- paste('DC', 1:dim(rotatedData)[2], sep = '')
  colnames(sequenceOfViewingAxis) <- paste('PC', 1:dim(sequenceOfViewingAxis)[2], sep = '')

  dataAngle <- getPolarCoordinates(rotatedData[,c(1,2)])[,'angle']
  intervalBorders <- seq(from = 0, to = 2*pi, length.out = 11)
  cellsInInterval <- which((dataAngle<intervalBorders[7])&(dataAngle>intervalBorders[6]))
  phaseToConsider <- as.numeric(which.max(table(ccPhaseInformation[cellsInInterval])))
  medianCurrentPhase <- median(dataAngle[ccPhaseInformation == levels(ccPhaseInformation)[phaseToConsider]])
  if (phaseToConsider<length(levels(ccPhaseInformation))){
    nextPhase <- phaseToConsider+1
  }else{
    nextPhase <- 1
  }
  medianNextPhase <- median(dataAngle[ccPhaseInformation == levels(ccPhaseInformation)[nextPhase]])
  if (medianNextPhase>medianCurrentPhase){
    rotatedData[,2] <- (-rotatedData[,2])
    rotatedWeightMatrix[,2] <- (-rotatedWeightMatrix[,2])
    rotationMatrix[,2] <- (-rotationMatrix[,2])
  }

  return(list(data = t(rotatedData), weights = t(rotatedWeightMatrix), rotationMatrix = rotationMatrix, dcProperties = data.frame(dcID = paste('DC', 1:dim(rotatedData)[2], sep = ''), diagCovOfData = diag(cov(rotatedData))), sequenceOfViewingAxes = sequenceOfViewingAxis))
}
#' @export
getGoldenSpiral <- function(numberOfPoints = 10^4){

  golden_angle <- pi * (1 + sqrt(5))

  theta <- golden_angle * (0:(numberOfPoints-1))
  z <- seq(from = 1 - 1.0 / numberOfPoints, to = 1.0 / numberOfPoints - 1, length.out = numberOfPoints)
  radius <- sqrt(1 - z * z)

  points <- as.data.frame(matrix(0L, numberOfPoints, 3))
  points[,1] <- radius * cos(theta)
  points[,2] <- radius * sin(theta)
  points[,3] <- z

  # plot3DTest <- plot_ly(points, x = ~V1, y = ~V2, z = ~V3, marker = list(size = 1))
  # add_markers(p = plot3DTest)

  return (points)
}
#' @export
getGridAroundAxis <- function(viewingAxisForGrid,
                              numberOfPoints = 10^4,
                              maxRadius = 0.1){

  radius <- sqrt((0:(numberOfPoints-1)) / numberOfPoints)*maxRadius

  golden_angle <- pi * (1 + sqrt(5))
  theta <- golden_angle * (0:(numberOfPoints-1))

  points <- as.data.frame(matrix(0L, numberOfPoints, 3))
  points[,1] <- radius*cos(theta)
  points[,2] <- radius*sin(theta)
  points[,3] <- rep(0, numberOfPoints)

  rotationMatrix <- getRotationMatrix3D(viewingAxisForGrid)
  rotatedPoints <- as.matrix(points) %*% t(rotationMatrix)

  viewingAxisForGridExtended <- matrix(viewingAxisForGrid,nrow=numberOfPoints,ncol=length(viewingAxisForGrid), byrow = TRUE)
  rotatedPoints <- rotatedPoints + viewingAxisForGridExtended

  rotatedPointsOnSphere <- rotatedPoints
  vectorLengths <- apply(rotatedPointsOnSphere,1,euclideanDistance,y=0)
  rotatedPointsOnSphere <- replicate(3, 1/vectorLengths)*rotatedPointsOnSphere

  # testData <-as.data.frame(rbind(pointsForSphere, rotatedPointsOnSphere))
  # plot3DTest <- plot_ly(testData, x=~V1, y=~V2, z=~V3, marker = list(size = 1))
  # add_markers(p = plot3DTest)

  return(rotatedPointsOnSphere)
}
#' @export
getBestViewingAxisMinimizingDispersionInThirdDimension <- function(data,
                                                                   gridToUse,
                                                                   ccPhaseInformation,
                                                                   boolInitialAxis = FALSE,
                                                                   initialAxis = NULL){

  getOnlyThirdDimensionIfRotateDataToNewViewingAxis <- function(dataToRotate,
                                                                rotationMatrixToUse){
    return((as.data.frame(as.matrix(dataToRotate) %*% rotationMatrixToUse))[,3])
  }
  numberOfRuns <- dim(gridToUse)[1]
  dataPhaseNames <- levels(ccPhaseInformation)

  if (!boolInitialAxis){
    newRotationMatrix <- getRotationMatrix3D(viewingAxis = c(0,0,1))
    rotatedDataCurrent <- getOnlyThirdDimensionIfRotateDataToNewViewingAxis(dataToRotate = data,
                                                                            rotationMatrixToUse = newRotationMatrix)
    #first angle
    meanThirdDimensionPerPhase <- matrix(0L, length(dataPhaseNames),1)
    for (i in 1:length(dataPhaseNames)){
      currentDataSet <- rotatedDataCurrent[ccPhaseInformation == dataPhaseNames[i]]
      meanThirdDimensionPerPhase[i] <- median(currentDataSet)
    }
    phaseDispersionThirdDimension <- sd(meanThirdDimensionPerPhase)
    bestViewingAxis <- c(0,0,1)
    smallestDispersion <- phaseDispersionThirdDimension
    bestMeanThirdDimensionPerPhase <- meanThirdDimensionPerPhase
  }else{
    newRotationMatrix <- getRotationMatrix3D(viewingAxis = initialAxis)
    rotatedDataCurrent <- getOnlyThirdDimensionIfRotateDataToNewViewingAxis(dataToRotate = data,
                                                                            rotationMatrixToUse = newRotationMatrix)
    #first angle
    meanThirdDimensionPerPhase <- matrix(0L, length(dataPhaseNames),1)
    for (i in 1:length(dataPhaseNames)){
      currentDataSet <- rotatedDataCurrent[ccPhaseInformation == dataPhaseNames[i]]
      meanThirdDimensionPerPhase[i] <- median(currentDataSet)
    }
    phaseDispersionThirdDimension <- sd(meanThirdDimensionPerPhase)
    bestViewingAxis <- initialAxis
    smallestDispersion <- phaseDispersionThirdDimension
    bestMeanThirdDimensionPerPhase <- meanThirdDimensionPerPhase
  }


  #all angles defined by grid
  for (k in 1:numberOfRuns){
    newRotationMatrix <- getRotationMatrix3D(viewingAxis = as.vector(t(gridToUse[k,])))
    rotatedDataCurrent <- getOnlyThirdDimensionIfRotateDataToNewViewingAxis(dataToRotate = data,
                                                                            rotationMatrixToUse = newRotationMatrix)

    meanThirdDimensionPerPhase <- matrix(0L, length(dataPhaseNames),1)
    for (i in 1:length(dataPhaseNames)){
      currentDataSet <- rotatedDataCurrent[ccPhaseInformation == dataPhaseNames[i]]
      meanThirdDimensionPerPhase[i] <- median(currentDataSet)
    }
    phaseDispersionThirdDimension <- sd(meanThirdDimensionPerPhase)

    if (phaseDispersionThirdDimension < smallestDispersion){
      bestViewingAxis <- as.vector(t(gridToUse[k,]))
      smallestDispersion <- phaseDispersionThirdDimension
      bestMeanThirdDimensionPerPhase <- meanThirdDimensionPerPhase
    }
  }

  return(bestViewingAxis)
}
#' @export
getRotationMatrix3D <- function(viewingAxis){

  zRotation <- 0
  z1 <- viewingAxis[1]/norm(viewingAxis,"2")
  z2 <- viewingAxis[2]/norm(viewingAxis,"2")
  z3 <- viewingAxis[3]/norm(viewingAxis,"2")
  sinRotAngleX <- z2/(z2^2+z3^2)*sqrt(1-z1^2)
  cosRotAngleX <- z3/(z2^2+z3^2)*sqrt(1-z1^2)
  sinRotAngleY <- (-z1)
  cosRotAngleY <- sqrt(1-z1^2)
  sinRotAngleZ <- sin(zRotation)
  cosRotAngleZ <- cos(zRotation)
  rotMatrixX <- matrix(c(1,0,0,0,cosRotAngleX,sinRotAngleX,0,-sinRotAngleX,cosRotAngleX), nrow=3)
  rotMatrixY <- matrix(c(cosRotAngleY,0,-sinRotAngleY,0,1,0,sinRotAngleY,0,cosRotAngleY), nrow=3)
  rotMatrixZ <- matrix(c(cosRotAngleZ,sinRotAngleZ,0,-sinRotAngleZ,cosRotAngleZ,0,0,0,1), nrow=3)

  return(t(rotMatrixZ %*% rotMatrixY %*% rotMatrixX))
}
#' @export
euclideanDistance <- function(x,y){
  return(sqrt(sum((x-y)^2)))
}
#' @export
getPolarCoordinates <- function(data,
                                boolGetCCTime = FALSE){
  locDataX <- data[,1]
  locDataY <- data[,2]

  locDataRadius <- sqrt(locDataX^2+locDataY^2)
  locDataAngle <- 0*locDataRadius

  #calculate angle in interval (-pi,pi]
  locDataAngle[locDataX>0] <- atan(locDataY[locDataX>0]/locDataX[locDataX>0])
  locDataAngle[(locDataX<0)&(locDataY>=0)] <- atan(locDataY[(locDataX<0)&(locDataY>=0)]/locDataX[(locDataX<0)&(locDataY>=0)])+pi
  locDataAngle[(locDataX<0)&(locDataY<0)] <- atan(locDataY[(locDataX<0)&(locDataY<0)]/locDataX[(locDataX<0)&(locDataY<0)])-pi
  locDataAngle[(locDataX=0)&(locDataY>0)] <- pi/2
  locDataAngle[(locDataX=0)&(locDataY<0)] <- -pi/2
  locDataAngle[(locDataX=0)&(locDataY=0)] <- 0

  #transform angle to interval [0,2pi)
  locDataAngle[locDataAngle<0] <- locDataAngle[locDataAngle<0]+2*pi

  if (boolGetCCTime){
    locCCProgression <- (2*pi-locDataAngle)/(2*pi)
    dataNew <- as.data.frame(cbind(locDataAngle, locDataRadius, locCCProgression))
    colnames(dataNew) <- c('angle', 'radius', 'cctime')
  }else{
    dataNew <- as.data.frame(cbind(locDataAngle, locDataRadius))
    colnames(dataNew) <- c('angle', 'radius')
  }
  rownames(dataNew) <- rownames(data)


  return(dataNew)
}
#' @export
getRotation2D <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': finding optimal 2D alignment: ', sep = ''))

  dataList <- alignCellDivisionToXAxis(dataList)

  # if (!(any(dataList@datasetInfo$goldStandardWeights == '')|is.null(dataList@datasetInfo$goldStandardWeights))){
  #   dataList <- alignCellCycleToGoldStandard(dataList,
  #                                            weightsGoldStandard = dataList@datasetInfo$goldStandardWeights)
  # }

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
alignCellDivisionToXAxis <- function(dataList){
  source('~/Studium Doktorarbeit/Projects/CellCycleSorting/Scripts/getPolarCoordinates.R')
  source('~/Studium Doktorarbeit/Projects/CellCycleSorting/Scripts/calculateAverageUMIPerIntervalBasedOnIndex.R')

  dataCurrent <- t(dataList@transformedData$dc$data[c(1,2),])
  dataPolar <- getPolarCoordinates(dataCurrent)
  dataAngle <- dataPolar[, 'angle']
  #names(dataAngle) <- rownames(dataPolar)
  cellNamesOrdered <- rev(rownames(dataPolar)[order(dataAngle)])
  countUMIperInterval <- as.vector(t(calculateAverageUMIPerIntervalBasedOnIndex(countData = dataList@DGEs$countData,
                                                                                inputPolarMatrix = dataPolar)))
  countUMIperIntervalSorted <- countUMIperInterval[order(countUMIperInterval)]
  minCountUMI <- mean(countUMIperIntervalSorted[1:max(2,floor(0.2*length(countUMIperInterval)))])
  maxCountUMI <- mean(countUMIperIntervalSorted[min(length(countUMIperInterval)-1,ceiling(0.9*length(countUMIperInterval))):length(countUMIperInterval)])
  linearFunction <- seq(from = minCountUMI, to = maxCountUMI, length.out = length(countUMIperInterval))
  residualsHelp <- rep(0, length(countUMIperInterval))
  residualsHelp[1] <- sum((countUMIperInterval-linearFunction)^2)
  for (i in 2:length(countUMIperInterval)){
    residualsHelp[i] <- sum((append(countUMIperInterval[i:length(countUMIperInterval)], countUMIperInterval[1:(i-1)])-linearFunction)^2)
  }



  divisionIndex <- which.min(residualsHelp)
  if (divisionIndex>1){
    if (countUMIperInterval[divisionIndex-1]<(minCountUMI+(maxCountUMI-minCountUMI)/4)){
      divisionIndex <- divisionIndex-1
    }else{
      if (countUMIperInterval[divisionIndex]>(minCountUMI+(maxCountUMI-minCountUMI)/3)){
        if (divisionIndex<length(countUMIperInterval)){
          divisionIndex <- divisionIndex+1
        }else{
          divisionIndex <- 1
        }
      }
      if (countUMIperInterval[divisionIndex]>(minCountUMI+(maxCountUMI-minCountUMI)/3)){
        if (divisionIndex<length(countUMIperInterval)){
          divisionIndex <- divisionIndex+1
        }else{
          divisionIndex <- 1
        }
      }
    }
  }else{
    if (countUMIperInterval[length(countUMIperInterval)]<(minCountUMI+(maxCountUMI-minCountUMI)/4)){
      divisionIndex <- length(countUMIperInterval)
    }else{
      if (countUMIperInterval[divisionIndex]>(minCountUMI+(maxCountUMI-minCountUMI)/3)){
        divisionIndex <- 2
      }
      if (countUMIperInterval[divisionIndex]>(minCountUMI+(maxCountUMI-minCountUMI)/3)){
        divisionIndex <- 3
      }
    }
  }
  divisionCellIndex <- floor(floor(dim(dataCurrent)[1]/(length(countUMIperInterval)))*(divisionIndex-1)+min(divisionIndex-1,dim(dataCurrent)[1]%%(length(countUMIperInterval))))+1
  divisionCell <- cellNamesOrdered[divisionCellIndex]
  if (divisionCellIndex>1){
    divisionAngle <- (dataPolar[divisionCell, 'angle']+dataPolar[cellNamesOrdered[divisionCellIndex-1], 'angle'])/2
  }else{
    divisionAngle <- (dataPolar[divisionCell, 'angle']+dataPolar[cellNamesOrdered[dim(dataCurrent)[1]], 'angle']-2*pi)/2
  }

  angleToRotate <- divisionAngle
  rotationMatrix2D <- matrix(c(cos(angleToRotate), sin(angleToRotate), -sin(angleToRotate), cos(angleToRotate)), nrow = 2)
  dataList@transformedData$dc$data[c(1,2),] <- t(t(dataList@transformedData$dc$data[c(1,2),])%*%rotationMatrix2D)
  dataList@transformedData$dc$weights[c(1,2),] <- t(t(dataList@transformedData$dc$weights[c(1,2),])%*%rotationMatrix2D)
  dataList@transformedData$dc$rotationMatrix[,c(1,2)] <- dataList@transformedData$dc$rotationMatrix[,c(1,2)]%*%rotationMatrix2D
  return(dataList)
}
alignCellCycleToGoldStandard <- function(dataList,
                                         weightsGoldStandard){
  numberOfIntervals <- 1000
  weightsDC1GoldStandard <- weightsGoldStandard[,1]
  names(weightsDC1GoldStandard) <- rownames(weightsGoldStandard)
  weightsDC2GoldStandard <- weightsGoldStandard[,2]
  names(weightsDC2GoldStandard) <- rownames(weightsGoldStandard)

  dataToUse <- t(as.matrix(dataList@transformedData$dc$data[c(1,2),]))
  dataWeightsToUse <- t(as.matrix(dataList@transformedData$dc$weights[c(1,2),]))
  if (dataList@datasetInfo$cellType == '3T3'){
    convertMouseGeneList <- function(x){

      #require("biomaRt")
      human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

      #genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human,  uniqueRows=FALSE)
      #humanx <- genesV2[, 2]

      mouseGenes <- mygene::queryMany(x, scopes="symbol", fields="ensembl.gene", species="mouse")
      mouseGenesID <- lapply(mouseGenes@listData$ensembl,'[[',1)
      for (i in 1:length(mouseGenesID)){
        if (length(mouseGenesID[[i]])>1){
          mouseGenesID[[i]] <- mouseGenesID[[i]][1]
        }
      }
      mouseGenesID <- as.character(mouseGenesID)
      mouseGenesNames <- as.character(lapply(mouseGenes@listData$query,'[[',1))
      mouseGenesID <- mouseGenesID[match(unique(mouseGenesNames), mouseGenesNames)]
      mouseGenesNames <- mouseGenesNames[match(unique(mouseGenesNames), mouseGenesNames)]

      #indexWithNullEntry <- mouseGenesID=='NULL'
      #xShortened <- x[!(mouseGenesID=='NULL')]
      mouseGenesIDShortened <- mouseGenesID[!(mouseGenesID=='NULL')]
      #mouseGenesNamesShortened <- mouseGenesNames[!(mouseGenesID=='NULL')]

      genesV2 = biomaRt::getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), filters = "ensembl_gene_id", values = mouseGenesIDShortened , mart = mouse, attributesL = c("hgnc_symbol"), martL = human,  uniqueRows=FALSE)

      indexOfGenesWhereNoAssociatedHumanGeneIsFound <- !(mouseGenesID%in%genesV2[,1])
      matchOfIndexes <- match(mouseGenesID[(mouseGenesID%in%genesV2[,1])], genesV2[,1])

      humanGeneNames <- as.data.frame(matrix(0L, length(x), 1))
      humanGeneNames[which(!(indexOfGenesWhereNoAssociatedHumanGeneIsFound)),] <- genesV2[matchOfIndexes,3]
      humanGeneNames <- as.character(t(humanGeneNames))
      helpVar <- table(humanGeneNames[!(humanGeneNames=='0')])
      helpVar <- names(helpVar[helpVar>1])
      for (i in 1:length(helpVar)){
        indexesOfDoubleName <- which(humanGeneNames == helpVar[i])
        bestFit <- 0
        for (j in 1:length(indexesOfDoubleName)){
          if (tolower(x[indexesOfDoubleName[j]])==tolower((humanGeneNames[indexesOfDoubleName])[1])){
            bestFit <- j
          }
        }
        if (bestFit>0){
          humanGeneNames[indexesOfDoubleName[-bestFit]] <- '0'
        }else{
          humanGeneNames[indexesOfDoubleName[-1]] <- '0'
        }
      }

      return(humanGeneNames)
    }
    mouseGenesUsedDuringPCAConverted <- as.character(t(convertMouseGeneList(dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes])))
    mouseGenesUsedDuringPCAConvertedWithoutZeros <- mouseGenesUsedDuringPCAConverted[!(mouseGenesUsedDuringPCAConverted=='0')]
    #A <- convertMouseGeneList(rownames(dataWeightsToUse))
    dataWeightsToUse <- dataWeightsToUse[!(mouseGenesUsedDuringPCAConverted=='0'),]
    rownames(dataWeightsToUse) <- mouseGenesUsedDuringPCAConvertedWithoutZeros

    genesUsedHere <- mouseGenesUsedDuringPCAConvertedWithoutZeros
  }else{
    genesUsedHere <- dataList@geneInfo$geneID[dataList@geneInfo$pcaGenes]
  }
  maximalCorrelation <- -1
  for (i in 1:numberOfIntervals){
    angleToRotate <- (i-1)*2*pi/numberOfIntervals
    rotationByAngle <- matrix(c(cos(angleToRotate), sin(angleToRotate), -sin(angleToRotate), cos(angleToRotate)), nrow = 2)
    dataRotated <- dataToUse%*%rotationByAngle
    dataWeightsRotated <- dataWeightsToUse%*%rotationByAngle

    weightsDC1DataRotated <- dataWeightsRotated[,1]
    weightsDC2DataRotated <- dataWeightsRotated[,2]

    jointGenesCurrent <- intersect(rownames(weightsGoldStandard), genesUsedHere)
    reducedDC1GoldStandard <- weightsDC1GoldStandard[jointGenesCurrent]
    reducedDC2GoldStandard <- weightsDC2GoldStandard[jointGenesCurrent]
    reducedDC1DataRotated <- weightsDC1DataRotated[names(weightsDC1DataRotated)%in%jointGenesCurrent]
    reducedDC2DataRotated <- weightsDC2DataRotated[names(weightsDC2DataRotated)%in%jointGenesCurrent]

    currentCorrelationDC1 <- cor(reducedDC1GoldStandard, reducedDC1DataRotated)
    currentCorrelationDC2 <- cor(reducedDC2GoldStandard, reducedDC2DataRotated)
    currentCorrelationMean <- mean(c(currentCorrelationDC1, currentCorrelationDC2))

    if (currentCorrelationMean>maximalCorrelation){
      optimalRotationAngle <- angleToRotate
      optimalInterval <- i
      optimalGeneSet <- jointGenesCurrent
      maximalCorrelation <- currentCorrelationMean
      maximalCorrelationDC1 <- currentCorrelationDC1
      maximalCorrelationDC2 <- currentCorrelationDC2
    }
  }

  dataList@transformedData$dc$angleRotatedToGoldStandard <- optimalRotationAngle
  rotationByAngle <- matrix(c(cos(optimalRotationAngle), sin(optimalRotationAngle), -sin(optimalRotationAngle), cos(optimalRotationAngle)), nrow = 2)
  dataList@transformedData$dc$data[c(1,2),] <- t(t(dataList@transformedData$dc$data[c(1,2),])%*%rotationByAngle)
  dataList@transformedData$dc$weights[c(1,2),] <- t(t(dataList@transformedData$dc$weights[c(1,2),])%*%rotationByAngle)
  dataList@transformedData$dc$rotationMatrix[,c(1,2)] <- dataList@transformedData$dc$rotationMatrix[,c(1,2)]%*%rotationByAngle

  return(dataList)
}
#' @export
getCCSorting <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': getting cell cycle sorting info: ', sep = ''))

  dcDataPolar <- getPolarCoordinates(t(dataList@transformedData$dc$data[c(1,2),]))
  newInfo <- c('ccAngle', 'ccRadius', 'ccTime', 'ccPercentage', 'ccPercentageUniformlySpaced', 'ccPositionIndex')
  oldInfo <- newInfo[newInfo%in%colnames(dataList@cellInfo)]
  newInfo <- newInfo[!(newInfo%in%colnames(dataList@cellInfo))]
  newData <- cbind(ccAngle = dcDataPolar[,'angle'],
                   ccRadius = dcDataPolar[,'radius'],
                   ccTime = (2*pi-dcDataPolar[,'angle'])/(2*pi)*dataList@datasetInfo$ccDurationTotal,
                   ccPercentage = (2*pi-dcDataPolar[,'angle'])/(2*pi),
                   ccPercentageUniformlySpaced = match(1:dim(dcDataPolar)[1], rev(order(dcDataPolar[,'angle'])))/dim(dcDataPolar)[1],
                   ccPositionIndex = rev(order(dcDataPolar[,'angle'])))
  if (length(oldInfo)>0){
    dataList@cellInfo[,oldInfo] <- newData[,oldInfo]
  }
  if (length(newInfo)>0){
    dataList@cellInfo <- cbind(dataList@cellInfo, newData[,newInfo])
  }

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
