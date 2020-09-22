#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom methods new
#'
setClass('Revelio',
         slots = list(datasetInfo = 'list',
                      DGEs = 'list',
                      cellInfo = 'list',
                      geneInfo = 'list',
                      transformedData = 'list',
                      velocityData = 'list',
                      stabilityIndex = 'list'))
#'
#'
#'
#'Create Revelio Object.
#'
#' 'createRevelioObject' returns a Revelio object according to the specified parameters.
#'
#'testtest
#'
#' @param datasetID ID should correspond to ID provided in supplementary file.
#' @return Returns Revelio object.
#'
#' @export
createRevelioObject <- function(rawData = NULL,
                                cyclicGenes = NULL,
                                datasetID = 'data1',
                                cellType = 'unknown',
                                sequTechnique = 'unknown',
                                lowernGeneCutoff = 1200,
                                uppernUMICutoff = 10^7,
                                ccPhaseAssignBasedOnIndividualBatches = TRUE,
                                ccPhaseAssignThresholdForCorSingleGeneToAvgExpression = 0.2,
                                ccPhaseAssignMaxDiffOfCCPhasesBetweenHighestScores = 1,
                                ccPhaseAssignThresholdHighestPhaseScore = 0.75,
                                ccPhaseAssignThresholdSecondHighestPhaseScore = 0.5,
                                pcaGenes = 'variableGenes',
                                thresholdsForPCWeightSignificance = 0.05,
                                ccDurationG1 = 0.487*19.33,
                                ccDurationS = 0.392*19.33,
                                ccDurationG2 = 0.093*19.33,
                                ccDurationM = 0.028*19.33,
                                ccDurationTotal = 19.33){

  startTime <- Sys.time()
  cat(paste(Sys.time(), ': reading data: ', sep = ''))

  dataList <- new('Revelio')

  dataList@datasetInfo <- list(datasetID = datasetID,
                               cellType = cellType,
                               sequTechnique = sequTechnique,
                               lowernGeneCutoff = lowernGeneCutoff,
                               uppernUMICutoff = uppernUMICutoff,
                               ccPhaseAssignBasedOnIndividualBatches = ccPhaseAssignBasedOnIndividualBatches,
                               ccPhaseAssignThresholdForCorSingleGeneToAvgExpression = ccPhaseAssignThresholdForCorSingleGeneToAvgExpression,
                               ccPhaseAssignMaxDiffOfCCPhasesBetweenHighestScores = ccPhaseAssignMaxDiffOfCCPhasesBetweenHighestScores,
                               ccPhaseAssignThresholdHighestPhaseScore = ccPhaseAssignThresholdHighestPhaseScore,
                               ccPhaseAssignThresholdSecondHighestPhaseScore = ccPhaseAssignThresholdSecondHighestPhaseScore,
                               pcaGenes = pcaGenes,
                               thresholdsForPCWeightSignificance = thresholdsForPCWeightSignificance,
                               ccDurationG1 = ccDurationG1,
                               ccDurationG2 = ccDurationG2,
                               ccDurationM = ccDurationM,
                               ccDurationTotal = ccDurationTotal)

  if (!is.null(rawData)){
    dataList@DGEs$countData <- as.data.frame(rawData)
  }
  if (!is.null(cyclicGenes)){
    dataList@datasetInfo$cyclicGenes <- as.data.frame(cyclicGenes)
  }

  #filter cells by nGene and nUMI
  dataList <- getFilteredDataBynUMIAndnGeneThresholds(dataList = dataList)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
getFilteredDataBynUMIAndnGeneThresholds <- function(dataList){
  # startTime <- Sys.time()
  # cat(paste(Sys.time(), ': filtering data by nGene and nUMI thresholds: ', sep = ''))

  dataList <- getStatisticsOnCellsAndGenes(dataList = dataList)

  filteredCells <- as.character(t(subset(dataList@cellInfo, nUMI<dataList@datasetInfo$uppernUMICutoff & nGene>dataList@datasetInfo$lowernGeneCutoff, select = 'cellID')))
  dataList@DGEs$countData <- dataList@DGEs$countData[,filteredCells]
  dataList@geneInfo <- dplyr::filter(dataList@geneInfo, rowSums(dataList@DGEs$countData)>0)
  dataList@DGEs$countData <- dataList@DGEs$countData[rowSums(dataList@DGEs$countData)>0,]
  dataList@cellInfo <- dplyr::filter(dataList@cellInfo, cellID%in%filteredCells)

  dataList <- getStatisticsOnCellsAndGenes(dataList = dataList)

  # cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
getStatisticsOnCellsAndGenes <- function(dataList){
  #raw data for current list element
  currentData <- dataList@DGEs$countData
  datasetName <- dataList@datasetInfo$datasetID
  #get batch information from column names (format expected: 'batchID_barcode')
  batchInfo <- colnames(currentData)%>%
    {do.call(rbind, strsplit(., '_'))}
  if (dim(batchInfo)[2]>1){
    batchInfo <- batchInfo[,1]
  }else{
    batchInfo <- rep(paste(datasetName, 'Batch1', sep = ''), dim(currentData)[2])
  }
  #calculate statistics on cells and write them to a data frame
  statisticsOnCells <- data.frame('datasetID' = rep(datasetName, dim(currentData)[2]), 'batchID' = batchInfo,  'cellID' = colnames(currentData), 'nUMI' = colSums(currentData), 'nGene' = colSums(currentData>0))
  statisticsOnCells[,'datasetID'] <- factor(statisticsOnCells[,'datasetID'], levels = datasetName)
  statisticsOnCells[,'batchID'] <- factor(statisticsOnCells[,'batchID'], levels = unique(batchInfo))
  statisticsOnCells[,'cellID'] <- as.character(statisticsOnCells[,'cellID'])

  #calculate statistics on genes and write them to a data frame
  statisticsOnGenes <- data.frame('datasetID' = rep(datasetName, dim(currentData)[1]), 'geneID' = rownames(currentData), 'nUMI' = rowSums(currentData), 'nCell' = rowSums(currentData>0))
  statisticsOnGenes[,'datasetID'] <- factor(statisticsOnGenes[,'datasetID'], levels = datasetName)
  statisticsOnGenes[,'geneID'] <- as.character(statisticsOnGenes[,'geneID'])

  # if (boolIncludeStatisticsOnBatches){
  #   for (j in levels(statisticsOnCells[,'batchID'])){
  #     currentDataBatch <- currentData[, statisticsOnCells[,'batchID']==j]
  #     currentDataBatch <- currentDataBatch[rowSums(currentDataBatch)>0,]
  #     listOfStatisticsOnCells[[j]] <- data.frame('datasetID' = rep(datasetName, dim(currentDataBatch)[2]), 'batchID' = rep(j, dim(currentDataBatch)[2]),  'cellID' = colnames(currentDataBatch), 'nUMI' = colSums(currentDataBatch), 'nGene' = colSums(currentDataBatch>0))
  #     listOfStatisticsOnCells[[j]][,'datasetID'] <- factor(listOfStatisticsOnCells[[j]][,'datasetID'], levels = datasetName)
  #     listOfStatisticsOnCells[[j]][,'batchID'] <- factor(listOfStatisticsOnCells[[j]][,'batchID'], levels = j)
  #     listOfStatisticsOnCells[[j]][,'cellID'] <- as.character(listOfStatisticsOnCells[[j]][,'cellID'])
  #
  #     #calculate statistics on genes and write them to a data frame
  #     listOfStatisticsOnGenes[[j]] <- data.frame('datasetID' = rep(datasetName, dim(currentDataBatch)[1]), 'geneID' = rownames(currentDataBatch), 'nUMI' = rowSums(currentDataBatch), 'nCell' = rowSums(currentDataBatch>0))
  #     listOfStatisticsOnGenes[[j]][,'datasetID'] <- factor(listOfStatisticsOnGenes[[j]][,'datasetID'], levels = datasetName)
  #     listOfStatisticsOnGenes[[j]][,'geneID'] <- as.character(listOfStatisticsOnGenes[[j]][,'geneID'])
  #   }
  # }
  if (length(dataList@cellInfo)==0){
    dataList@cellInfo <- statisticsOnCells
  }else{
    dataList@cellInfo[,colnames(statisticsOnCells)] <- statisticsOnCells
  }
  if (length(dataList@geneInfo)==0){
    dataList@geneInfo <- statisticsOnGenes
  }else{
    dataList@geneInfo[,colnames(statisticsOnGenes)] <- statisticsOnGenes
  }
  #rm(currentData, batchInfo)
  return(dataList)
}
#' @export
getCellCyclePhaseAssignInformation <- function(dataList){
  startTime <- Sys.time()
  cat(paste(Sys.time(), ': assigning cell cycle phases: ', sep = ''))

  dataList@datasetInfo$scalingFactorUMI <- median(dataList@cellInfo$nUMI)
  dataList@DGEs$logOfFractionsData <- getLogOfFractionsData(data = dataList@DGEs$countData,
                                                            scalingFactorUMI = dataList@datasetInfo$scalingFactorUMI)

  resultList <- cellCyclePhaseAssign(data = dataList@DGEs$logOfFractionsData,
                                     batchIDInformation = dataList@cellInfo[,'batchID'],
                                     boolBaseColouringOnIndividualReplicates = dataList@datasetInfo$ccPhaseAssignBasedOnIndividualBatches,
                                     geneBucketListFileLocation = dataList@datasetInfo$ccPhaseAssignGeneBucketListFileLocation,
                                     geneBucketListSheetIndex = dataList@datasetInfo$ccPhaseAssignGeneBucketListSheetIndex,
                                     corScoreGeneToAvgExprThreshold = dataList@datasetInfo$ccPhaseAssignThresholdForCorSingleGeneToAvgExpression,
                                     cyclicGenes = dataList@datasetInfo$cyclicGenes)
  phaseScore <- resultList[['phaseScore']]
  phaseScoreNormalized <- resultList[['phaseScoreNormalized']]
  colnames(phaseScoreNormalized) <- paste(colnames(phaseScoreNormalized), '_zScore', sep = '')
  phaseAssignStatistics <- data.frame(ccPhase = resultList[['ccPhase']], meanOfPhaseScore = apply(resultList[['phaseScore']], 1, mean), sdOfPhaseScore = apply(resultList[['phaseScore']], 1, sd), highestPhaseScore = resultList[['highestPhaseScore']], secondHighestPhaseScore = resultList[['secondHighestPhaseScore']])
  phaseAssignStatistics <- cbind(phaseAssignStatistics, phaseScore, phaseScoreNormalized)
  dataList@geneInfo <- cbind(dataList@geneInfo, ccPhase = resultList[['geneCCPhase']])

  #get outliers
  phaseDifferenceBetweenHighestAndSecondHighest <- t(subset(phaseAssignStatistics, select = colnames(phaseScoreNormalized)))
  highestPhaseScoreIndex <- apply(phaseDifferenceBetweenHighestAndSecondHighest, 2, which.max)
  highestPhaseScoreIndexOfMatrix <- (1:length(highestPhaseScoreIndex))*dim(phaseDifferenceBetweenHighestAndSecondHighest)[1]-dim(phaseDifferenceBetweenHighestAndSecondHighest)[1]+highestPhaseScoreIndex
  phaseDifferenceBetweenHighestAndSecondHighest[highestPhaseScoreIndexOfMatrix] <- -5
  secondHighestPhaseScoreIndex <- apply(phaseDifferenceBetweenHighestAndSecondHighest, 2, which.max)
  differenceHighestIndexToSecondIndex <- apply(cbind(abs(highestPhaseScoreIndex - secondHighestPhaseScoreIndex), abs(highestPhaseScoreIndex + dim(phaseDifferenceBetweenHighestAndSecondHighest)[1] - secondHighestPhaseScoreIndex), abs(highestPhaseScoreIndex - dim(phaseDifferenceBetweenHighestAndSecondHighest)[1] - secondHighestPhaseScoreIndex)), 1, min)
  isOutlierSuspectedDoublets <- as.logical(rep(FALSE, dim(phaseAssignStatistics)[1]))
  isOutlierNoConfidenceInPhaseScore <- as.logical(rep(FALSE, dim(phaseAssignStatistics)[1]))
  isOutlierSuspectedDoublets[which((differenceHighestIndexToSecondIndex>dataList@datasetInfo$ccPhaseAssignMaxDiffOfCCPhasesBetweenHighestScores)&(phaseAssignStatistics[,'secondHighestPhaseScore']>dataList@datasetInfo$ccPhaseAssignThresholdSecondHighestPhaseScore))] <- TRUE
  isOutlierNoConfidenceInPhaseScore[which(phaseAssignStatistics[,'highestPhaseScore']<dataList@datasetInfo$ccPhaseAssignThresholdHighestPhaseScore)] <- TRUE
  phaseAssignStatistics <- cbind(phaseAssignStatistics, isOutlierSuspectedDoublets, isOutlierNoConfidenceInPhaseScore)

  #combine data
  dataList@cellInfo <- cbind(dataList@cellInfo, phaseAssignStatistics)

  #filter outliers
  dataList <- getFilteredDataByCCPhaseAssignOutliers(dataList = dataList)

  cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
cellCyclePhaseAssign <- function(data,
                                 batchIDInformation,
                                 boolBaseColouringOnIndividualReplicates,
                                 nameSeparationSymbol = '_',
                                 geneBucketListFileLocation,
                                 geneBucketListSheetIndex,
                                 corScoreGeneToAvgExprThreshold,
                                 cyclicGenes){

  #read excel file containing cell cycle genes for different phases
  geneCategoriesExternal <- cyclicGenes

  #get category names
  ccPhaseNames <- colnames(geneCategoriesExternal)

  dataAll <- data
  scoreForEachCCPhase <- as.data.frame(matrix(0L, dim(dataAll)[2], length(ccPhaseNames)))
  rownames(scoreForEachCCPhase) <- colnames(dataAll)
  colnames(scoreForEachCCPhase) <- ccPhaseNames
  ccPhaseInformationHelpAll <- NULL
  namesCCPhaseInformation <- NULL
  phaseScoreCellsBeforeNormalizationAll <- NULL
  phaseScoreCellsAll <- NULL
  highestPhaseScoreAll <- NULL
  secondHighestPhaseScoreAll <- NULL
  if (!boolBaseColouringOnIndividualReplicates){
    replicateNamesSave <- sort(unique(batchIDInformation))
    replicateNames <- 1
    replicateIDInformationSave <- batchIDInformation
    batchIDInformation <- rep(1, length(batchIDInformation))
  }else{
    replicateNames <- sort(unique(batchIDInformation))
  }
  for (k in replicateNames){
    data <- dataAll[, batchIDInformation==k]
    #initiate data frame
    geneNames <- t(data.frame(matrix(NA, nrow = length(ccPhaseNames), ncol = dim(geneCategoriesExternal)[1]), row.names = ccPhaseNames))

    #get all genes from excel file which are actually present in data
    for (i in 1:length(ccPhaseNames)){
      geneNames[1:(length(geneCategoriesExternal[,i])), i] <- gsub(" ", "", geneCategoriesExternal[,i])
      number.Of.Genes.In.Data <- sum(geneNames[,i] %in% rownames(data))
      geneNames[1:number.Of.Genes.In.Data, i] <- geneNames[(geneNames[,i] %in% rownames(data)), i]
      geneNames[(number.Of.Genes.In.Data+1):dim(geneCategoriesExternal)[1], i] <- NA
    }

    #delete zero genes
    for (i in 1:length(ccPhaseNames)){
      notNAGeneNames <- geneNames[!is.na(geneNames[, i]), i]
      indexOfZeroGenes <- which(rowSums(data[notNAGeneNames,]) == 0)
      if (length(indexOfZeroGenes)>0){
        geneNames[1:(length(notNAGeneNames)-length(indexOfZeroGenes)), i] <- notNAGeneNames[-indexOfZeroGenes]
        geneNames[(length(notNAGeneNames)-length(indexOfZeroGenes)+1):dim(geneCategoriesExternal)[1], i] <- NA
      }
    }

    #delete empty rows
    geneNames <- geneNames[rowSums(is.na(geneNames))<length(ccPhaseNames),]
    geneNamesSave <- geneNames

    #average gene expression pattern per phase
    avgExprPattern <- data.frame(matrix(0L, dim(data)[2], length(ccPhaseNames)))
    colnames(avgExprPattern) <- colnames(geneNames)
    for (i in 1:length(ccPhaseNames)){
      avgExprPattern[,i] <- colMeans(data[geneNames[!is.na(geneNames[,i]), i],])
    }

    #score for each gene compared to genes in the same category
    corScore <- matrix(0L, dim(geneNames)[1], length(ccPhaseNames))
    for (i in 1:length(ccPhaseNames)){
      for (j in 1:sum(!is.na(geneNames[,i]))){
        corScore[j,i] <- cor(as.numeric(t(data[geneNames[j,i],])), as.numeric(avgExprPattern[,i]))
      }
    }

    #dplyr::filter out genes behaving differently in our data
    for (i in 1:length(ccPhaseNames)){
      numberOfValidGenes <- sum(corScore[,i] > corScoreGeneToAvgExprThreshold)
      geneNames[1:numberOfValidGenes, i] <- geneNames[corScore[,i] > corScoreGeneToAvgExprThreshold, i]
      if (numberOfValidGenes<(dim(geneNames)[1])){
        geneNames[(numberOfValidGenes+1):dim(geneNames)[1], i] <- NA
      }
    }

    #delete empty rows
    geneNames <- geneNames[rowSums(is.na(geneNames))<length(ccPhaseNames),]

    #score for each cell in each category
    phaseScoreCells <- data.frame(matrix(0L, dim(data)[2], length(ccPhaseNames)))
    rownames(phaseScoreCells) <- colnames(data)
    colnames(phaseScoreCells) <- ccPhaseNames
    for (j in 1:length(ccPhaseNames)){
      phaseScoreCells[,j] <- colMeans(data[geneNames[!is.na(geneNames[,j]), j], ], na.rm=TRUE)
    }
    phaseScoreCellsBeforeNormalization <- phaseScoreCells

    #scale the scores
    phaseScoreCells <- as.data.frame(apply(phaseScoreCells, 2, scale))
    phaseScoreCells <- as.data.frame(t(apply(phaseScoreCells, 1, scale)))
    rownames(phaseScoreCells) <- colnames(data)
    colnames(phaseScoreCells) <- ccPhaseNames

    #assign highest score
    phaseAssign <- apply(phaseScoreCells, 1, which.max)

    phaseScoreCellsSecond <- t(phaseScoreCells)
    maxIndexPerCell <- (1:length(phaseAssign))*length(ccPhaseNames)-length(ccPhaseNames)+phaseAssign
    phaseScoreCellsSecond[maxIndexPerCell] <- -9
    phaseScoreCellsSecond <- t(phaseScoreCellsSecond)

    highestPhaseScore <- apply(phaseScoreCells, 1, max)
    secondHighestPhaseScore <- apply(phaseScoreCellsSecond, 1, max)

    ccPhaseNames <- colnames(phaseScoreCells)
    ccPhaseInformation <- ccPhaseNames[phaseAssign]
    ccPhaseNames <- ccPhaseNames[as.numeric(names(table(phaseAssign)))]

    orderNecessary <- match(ccPhaseNames,ccPhaseNames[order(ccPhaseNames)])

    ccPhaseInformationHelpAll <- append(ccPhaseInformationHelpAll,as.character(ccPhaseInformation))
    namesCCPhaseInformation <- append(namesCCPhaseInformation,colnames(data))

    phaseScoreCellsBeforeNormalizationAll <- rbind(phaseScoreCellsBeforeNormalizationAll, phaseScoreCellsBeforeNormalization)
    phaseScoreCellsAll <- rbind(phaseScoreCellsAll, phaseScoreCells)
    highestPhaseScoreAll <- append(highestPhaseScoreAll, highestPhaseScore)
    secondHighestPhaseScoreAll <- append(secondHighestPhaseScoreAll, secondHighestPhaseScore)
  }
  if (!boolBaseColouringOnIndividualReplicates){
    batchIDInformation <- replicateIDInformationSave
    replicateNames <- sort(unique(batchIDInformation))
  }
  data <- dataAll
  matchHelp <- match(colnames(data),namesCCPhaseInformation)
  ccPhaseInformationHelpAll <- ccPhaseInformationHelpAll[matchHelp]
  namesCCPhaseInformation <- namesCCPhaseInformation[matchHelp]
  ccPhaseInformationHelpAll <- as.factor(ccPhaseInformationHelpAll)
  ccPhaseInformationAlgorithm <- factor(ccPhaseInformationHelpAll, levels = levels(ccPhaseInformationHelpAll)[orderNecessary])
  names(ccPhaseInformationAlgorithm) <- namesCCPhaseInformation

  ccPhaseAssignOfGenes <- rep('-', dim(data)[1])
  names(ccPhaseAssignOfGenes) <- rownames(data)
  for (i in levels(ccPhaseInformationAlgorithm)){
    ccPhaseAssignOfGenes[geneNamesSave[!is.na(geneNamesSave[,i]),i]] <- i
  }
  ccPhaseAssignOfGenes <- factor(ccPhaseAssignOfGenes, levels = c(levels(ccPhaseInformationAlgorithm), '-'))

  return(list(ccPhase = ccPhaseInformationAlgorithm, phaseScore = phaseScoreCellsBeforeNormalizationAll, phaseScoreNormalized = phaseScoreCellsAll, highestPhaseScore = highestPhaseScoreAll, secondHighestPhaseScore = secondHighestPhaseScoreAll, geneCCPhase = ccPhaseAssignOfGenes))
}
getLogOfFractionsData <- function(data,
                                  scalingFactorUMI){

  countUMIPerCell <- colSums(data)                                   #gives total UMI count per cell
  data <- sweep(data,2,countUMIPerCell,'/')                          #divides each column by total UMI count for that cell
  data <- data*scalingFactorUMI                                                     #multiplies every entry by 10^4 as a scaling factor
  data <- log(data+1)                                                   #transforms to logarithmic scale (natural logarithm)

  return(data)
}
getFilteredDataByCCPhaseAssignOutliers <- function(dataList){
  # startTime <- Sys.time()
  # cat(paste(Sys.time(), ': filtering data by cc phase outliers: ', sep = ''))

  filteredCells <- as.character(t(subset(dataList@cellInfo, !(isOutlierSuspectedDoublets|isOutlierNoConfidenceInPhaseScore), select = 'cellID')))
  dataList@DGEs$countData <- dataList@DGEs$countData[,filteredCells]
  dataList@geneInfo <- dplyr::filter(dataList@geneInfo, rowSums(dataList@DGEs$countData)>0)
  dataList@DGEs$countData <- dataList@DGEs$countData[rowSums(dataList@DGEs$countData)>0,]
  dataList@cellInfo <- dplyr::filter(dataList@cellInfo, cellID%in%filteredCells)

  dataList <- getStatisticsOnCellsAndGenes(data = dataList)

  # cat(paste(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units'), '\n', sep = ''))
  return(dataList)
}
