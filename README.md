# Revelio

## Installation
```
devtools::install_github('danielschw188/Revelio')
```

## How to Use
```
library(Revelio)
myData <- createRevelioObject(rawData = dataSource,
                              cyclicGenes = cyclicGenesSource)
myData <- getCellCyclePhaseAssignInformation(dataList = myData)
myData <- getPCAData(dataList = myData)
myData <- getOptimalRotation(dataList = myData)
```

## References
