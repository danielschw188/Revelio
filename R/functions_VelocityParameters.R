#' @export
readIntronicDataAndFilter <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': reading intronic data: ', sep = ''))

  dataList@DGEs$intronCountData <- readRDS(dataList@datasetInfo$intronicDataPath)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': filtering data by intronic data: ', sep = ''))

  filteredCells <- as.character(t(subset(dataList@cellInfo, cellID%in%colnames(dataList@DGEs$intronCountData), select = 'cellID')))
  filteredGenes <- as.character(t(subset(dataList@geneInfo, geneID%in%rownames(dataList@DGEs$intronCountData), select = 'geneID')))
  dataList@DGEs$intronCountData <- dataList@DGEs$intronCountData[filteredGenes, filteredCells]
  dataList@DGEs$countData <- dataList@DGEs$countData[filteredGenes, filteredCells]
  dataList@geneInfo <- dplyr::filter(dataList@geneInfo, geneID%in%filteredGenes)
  dataList@geneInfo <- dplyr::filter(dataList@geneInfo, rowSums(dataList@DGEs$countData)>0)
  dataList@DGEs$intronCountData <- dataList@DGEs$intronCountData[rowSums(dataList@DGEs$countData)>0,]
  dataList@DGEs$countData <- dataList@DGEs$countData[rowSums(dataList@DGEs$countData)>0,]
  dataList@cellInfo <- dplyr::filter(dataList@cellInfo, cellID%in%filteredCells)

  dataList <- getStatisticsOnCellsAndGenes(dataList = dataList)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
getVelocityMatrix <- function(dataList){

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating velocity offset estimate: ', sep = ''))
  dataList <- offsetEstimation(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating velocity gamma estimate: ', sep = ''))
  dataList <- gammaEstimation(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  # startTime <- Sys.time()
  # cat(paste(Sys.time(), ': filtering genes according to velocity paper: ', sep = ''))
  # dataList <- geneFilteringVelocityPaper(dataList = dataList)
  # cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating RNA velocity matrix: ', sep = ''))
  dataList <- calculateRNAvelocity(dataList = dataList)
  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))

  return(dataList)
}
#' @export
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
#' @export
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
#' @export
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
#' @export
getNormalizedDataBySpecificVector <- function(data,
                                              divideByCustomVector = NULL){
  return(sweep(data,2,divideByCustomVector,'/'))
}
