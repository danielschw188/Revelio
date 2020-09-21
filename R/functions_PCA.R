#' @export
getNormalizedData <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating z-score data: ', sep = ''))

  dataList@datasetInfo$scalingFactorUMI <- median(dataList@cellInfo$nUMI)
  dataList@DGEs$logOfFractionsData <- getLogOfFractionsData(data = dataList@DGEs$countData,
                                                            scalingFactorUMI = dataList@datasetInfo$scalingFactorUMI)

  dataList@DGEs$scaledData <- getScaledDataFromLogOfFractions(dataList@DGEs$logOfFractionsData)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
getScaledDataFromLogOfFractions <- function(data){

  #normalize data by rows
  meanUMIForEachGene <- rowMeans(data)                                   #calculates mean UMI count per gene
  data <- sweep(data,1,meanUMIForEachGene,'-')                           #centralizes gene data
  sdUMIForEachGene <- apply(data,1,sd)                                   #calculates sd of UMI count per gene
  data <- sweep(data,1,sdUMIForEachGene,'/')                             #divides each row by sd of that gene

  #  return(list(data, dataLogOfFractions, dataFractions, meanUMIForEachGene, sdUMIForEachGene))
  return(data)
}
#' @export
getVariableGenes <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating variable genes: ', sep = ''))

  variableGenes <- rep(FALSE, dim(dataList@geneInfo)[1])
  names(variableGenes) <- dataList@geneInfo[,'geneID']
  resultList <- calculateVariableGenes(data = dataList@DGEs$logOfFractionsData,
                                       batchID = dataList@cellInfo[,'batchID'])
  variableGenesHelp <- resultList[['varGenes']]
  meanGeneExpression <- resultList[['meanGeneExpr']]
  normalizedDispersionPerGene <- resultList[['normDispPerGene']]
  variableGenes[variableGenesHelp] <- TRUE
  dataList@geneInfo <- cbind(dataList@geneInfo, variableGenes, meanGeneExpression, normalizedDispersionPerGene)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
calculateVariableGenes <- function(data,
                                   batchID,
                                   numberOfBins = 20,
                                   minGeneExpression = 0.2,
                                   maxGeneExpression = 4,
                                   minNormDispersion = 0.5,
                                   maxNormDispersion = 10){

  exponentialMean <- function(x) {
    return(log(mean(exp(x))))
  }
  logVarDividedMean <- function(x) {
    return(log(var(exp(x)-1)/mean(exp(x)-1)))
  }

  # batchNames <- unique(batchID)
  # if (length(batchNames)>1){
  #   variableGenesAll <- as.data.frame(matrix(NA, nrow = dim(data)[1], ncol = length(batchNames)))
  #   meanGeneExpression <- rep(0, dim(data)[1])
  #   normalizedDispersionPerGene <- rep(0, dim(data)[1])
  #   for (i in 1:length(batchNames)){
  #     dataCurrent <- data[,batchID==batchNames[i]]
  #
  #     meanGeneExpressionSingle <- apply(dataCurrent,1,exponentialMean)
  #     meanGeneExpressionSingle[is.na(meanGeneExpressionSingle)] <- 0
  #     dispersionPerGene <- apply(dataCurrent,1,logVarDividedMean)
  #     dispersionPerGene[is.na(dispersionPerGene)] <- 0
  #
  #     binsMeanGeneExpression <- cut(meanGeneExpressionSingle,numberOfBins)
  #     binnedMeanForDispersion <- tapply(dispersionPerGene,binsMeanGeneExpression,mean)
  #     binnedMeanForDispersion[is.na(binnedMeanForDispersion)] <- 0
  #     binnedSDForDispersion <- tapply(dispersionPerGene,binsMeanGeneExpression,sd)
  #     binnedSDForDispersion[is.na(binnedSDForDispersion)] <- 0
  #     normalizedDispersionPerGeneSingle <- (dispersionPerGene-binnedMeanForDispersion[as.numeric(binsMeanGeneExpression)])/binnedSDForDispersion[as.numeric(binsMeanGeneExpression)]
  #     normalizedDispersionPerGeneSingle[is.na(normalizedDispersionPerGeneSingle)] <- 0
  #
  #     #select variable genes
  #     variableGenesSingle <- rownames(dataCurrent)[(meanGeneExpressionSingle>minGeneExpression)&(meanGeneExpressionSingle<maxGeneExpression)&(normalizedDispersionPerGeneSingle>minNormDispersion)&(normalizedDispersionPerGeneSingle<maxNormDispersion)]
  #
  #     variableGenesAll[1:length(variableGenesSingle),i] <- variableGenesSingle
  #     meanGeneExpression <- meanGeneExpression + meanGeneExpressionSingle*dim(dataCurrent)[2]/dim(data)[2]
  #     normalizedDispersionPerGene <- normalizedDispersionPerGene + normalizedDispersionPerGeneSingle*dim(dataCurrent)[2]/dim(data)[2]
  #   }
  #   variableGenesAll <- variableGenesAll[rowSums(is.na(variableGenesAll))<length(batchNames),]
  #   variableGenesVector <- as.vector(variableGenesAll[,1])
  #   for (i in 2:length(batchNames)){
  #     variableGenesVector <- append(variableGenesVector, as.vector(variableGenesAll[,i]))
  #   }
  #   variableGenesVector <- variableGenesVector[!is.na(variableGenesVector)]
  #   amountOfDataSetAVariableGeneIsFoundIn<-table(variableGenesVector)
  #   if (length(batchNames)==2){
  #     numberOfDataSetsAVariableGeneShouldBeIn <- 1  #used to be 2
  #   }else{
  #     numberOfDataSetsAVariableGeneShouldBeIn <- length(batchNames)-ceiling((length(batchNames)-1)/3)
  #   }
  #   variableGenes <- names(amountOfDataSetAVariableGeneIsFoundIn)[amountOfDataSetAVariableGeneIsFoundIn>=numberOfDataSetsAVariableGeneShouldBeIn]
  #
  # }else{
  meanGeneExpression <- apply(data,1,exponentialMean)
  meanGeneExpression[is.na(meanGeneExpression)] <- 0
  dispersionPerGene <- apply(data,1,logVarDividedMean)
  dispersionPerGene[is.na(dispersionPerGene)] <- 0

  binsMeanGeneExpression <- cut(meanGeneExpression,numberOfBins)
  binnedMeanForDispersion <- tapply(dispersionPerGene,binsMeanGeneExpression,mean)
  binnedMeanForDispersion[is.na(binnedMeanForDispersion)] <- 0
  binnedSDForDispersion <- tapply(dispersionPerGene,binsMeanGeneExpression,sd)
  binnedSDForDispersion[is.na(binnedSDForDispersion)] <- 0
  normalizedDispersionPerGene <- (dispersionPerGene-binnedMeanForDispersion[as.numeric(binsMeanGeneExpression)])/binnedSDForDispersion[as.numeric(binsMeanGeneExpression)]
  normalizedDispersionPerGene[is.na(normalizedDispersionPerGene)] <- 0

  #select variable genes
  variableGenes <- rownames(data[(meanGeneExpression>minGeneExpression)&(meanGeneExpression<maxGeneExpression)&(normalizedDispersionPerGene>minNormDispersion)&(normalizedDispersionPerGene<maxNormDispersion),])
  # }

  # object@listVariableGenes$variableGenes <- variableGenes
  # object@listVariableGenes$meanGeneExpression <- meanGeneExpression
  # object@listVariableGenes$normalizedDispersionForEachGene <- normalizedDispersionPerGene

  return(list(varGenes = variableGenes, meanGeneExpr = meanGeneExpression, normDispPerGene = normalizedDispersionPerGene))
}
#' @export
getPCAGenes <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': determining PCA genes: ', sep = ''))

  if (length(dataList@datasetInfo$pcaGenes) == 1){
    if (dataList@datasetInfo$pcaGenes == 'variableGenes'){
      # if (dataList@datasetInfo$intronDataExists){
      #   dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = dataList@geneInfo[,'variableGenes']&(rowSums(dataList@DGEs$intronCountData[dataList@geneInfo$geneID,])>0))
      # }else{
        dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = dataList@geneInfo[,'variableGenes'])
      # }
    }else{
      if (dataList@datasetInfo$pcaGenes == 'allGenes'){
        # if (dataList@datasetInfo$intronDataExists){
        #   dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = rep(TRUE, dim(dataList@geneInfo)[1])&(rowSums(dataList@DGEs$intronCountData[dataList@geneInfo$geneID,])>0))
        # }else{
          dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = rep(TRUE, dim(dataList@geneInfo)[1]))
        # }
      }
    }
  }else{
    pcaGenesHelp <- rep(FALSE, dim(dataList@geneInfo)[1])
    names(pcaGenesHelp) <- dataList@geneInfo[,'geneID']
    pcaGenesHelp[dataList@datasetInfo$pcaGenes[dataList@datasetInfo$pcaGenes%in%names(pcaGenesHelp)]] <- TRUE
    # if (dataList@datasetInfo$intronDataExists){
    #   dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = pcaGenesHelp&(rowSums(dataList@DGEs$intronCountData[dataList@geneInfo$geneID,])>0))
    # }else{
      dataList@geneInfo <- cbind(dataList@geneInfo, pcaGenes = pcaGenesHelp)
    # }
  }

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
getPCAData <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': calculating PCA: ', sep = ''))

  pcaHelp <- dplyr::filter(dataList@DGEs$scaledData, dataList@geneInfo[,'pcaGenes'])
  rownames(pcaHelp) <- t(subset(dataList@geneInfo, dataList@geneInfo[,'pcaGenes'], select = geneID))
  dataList@transformedData$pca <- doPCA(pcaHelp)

  cellCycleScore <- getCellCycleScoreForPCA(dataList@transformedData$pca$data, dataList@cellInfo[,'ccPhase'])
  boolOutliers <- rep(FALSE, length(cellCycleScore))
  boolOutliers[1:(min(length(cellCycleScore),100))] <- calculateOutliersInVector(data = as.vector(cellCycleScore[1:(min(length(cellCycleScore),100))]),
                                                                                 outlierThreshold = 2)
  dataList@transformedData$pca$pcProperties <- cbind(dataList@transformedData$pca$pcProperties,
                                                     ccScore = cellCycleScore,
                                                     isComponentAssociatedWithCC = boolOutliers)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
#' @export
doPCA <- function(data){

  #calculate covariance matrix
  covData <- cov(t(data))
  eigenProperties <- eigen(covData)
  eigenValuesCov <- eigenProperties$values
  eigenVectorsCov <- eigenProperties$vectors

  #pcadata after rotation for 20 dimensions
  weightMatrix <- eigenVectorsCov
  rownames(weightMatrix) <- rownames(data)
  pcaData <- as.data.frame(t(data)%*%weightMatrix)
  rownames(pcaData) <- colnames(data)
  colnames(pcaData) <- paste('PC', 1:dim(pcaData)[2], sep = '')
  colnames(weightMatrix) <- paste('PC', 1:dim(pcaData)[2], sep = '')
  rownames(eigenVectorsCov) <- rownames(data)
  colnames(eigenVectorsCov) <- paste('PC', 1:dim(pcaData)[2], sep = '')
  names(eigenValuesCov) <- paste('PC', 1:dim(pcaData)[2], sep = '')

  return(list(data = t(pcaData), weights = as.data.frame(t(weightMatrix)), pcProperties = data.frame(pcID = names(eigenValuesCov), eigenValue = eigenValuesCov, row.names = NULL), eigenVectors = as.data.frame(t(eigenVectorsCov))))

}
#' @export
getCellCycleScoreForPCA <- function(pcaData,
                                    ccPhaseInformation){
  numberOfPCs <- dim(pcaData)[1]
  scoreCorrelationOfCellCycleToEachPC <- matrix(0L, numberOfPCs, 1)
  for (i in 1:numberOfPCs){
    scoreCorrelationOfCellCycleToEachPC[i] <- calculateCCScoreSingleDimension(pcaData[i,],
                                                                              ccPhaseInformation)
  }
  return(scoreCorrelationOfCellCycleToEachPC)
}
#' @export
calculateOutliersInVector <- function(data,
                                      outlierThreshold = 2){
  return((abs(data - mean(data)) / sd(data))>outlierThreshold)
  #return((abs(data - median(data)) / mad(data))>outlierThreshold)
}
#' @export
calculateCCScoreSingleDimension <- function(data,
                                            ccPhaseInformation){
  ccPhaseNames <- levels(ccPhaseInformation)

  meanOfCellCycleClusters <- matrix(0L, length(ccPhaseNames),1)
  for (i in 1:length(ccPhaseNames)){
    currentDataSet <- data[ccPhaseInformation == ccPhaseNames[i]]
    meanOfCellCycleClusters[i] <- mean(currentDataSet)
  }
  return (sd(meanOfCellCycleClusters))
}
