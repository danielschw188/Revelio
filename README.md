# Revelio

The Revelio package is designed to extract cell cycle information from scRNA-seq experiments of cycling cells. What sets Revelio apart is that the transformation from normalized data to final outcome is fully linear. The advantage is that any conclusion drawn on the transformed data can be transferred back to the original data since the algorithm is easily inverted. One simple application of this principle is the removal of cell cycle effects from the normalized data matrix.

## Installation
```
devtools::install_github('danielschw188/Revelio')
```

## How to Use

The algorithm needs two inputs from the user:
1. 'rawData': A gene-by-cell data matrix that contains the UMI counts as entries. Comments about formatting:
    - gene names are placed in the row names
    - cell IDs are placed in the column names
    - entries of the matrix are typically positive integers
    - NA entries should be replaced by 0
2. 'cyclicGenes': A matrix of gene names that will be used for marker gene lists. Comments about formatting:
    - the entries of each column of the matrix represents a separate marker gene list for a specific phase or phase transition
    - the name of the phase or phase transition is placed in the column names
    - the entries of the matrix are the gene names
    - if one of the marker gene lists is shorter than the longest marker gene list, all additionally entries of the shorter list should be set to NA
    
We have included sample data for these two inputs into the package. The sample raw data matrix is called 'revelioTestData_rawDataMatrix' and the table of marker gene lists is called 'revelioTestData_cyclicGenes'. This is the data illustrated in Figure 1 of the publication. During the following example, we will utilize these two inputs. When analyzing your own data, the variable names of your own data have to be exchanged.

Now the Revelio algorithm can be executed. First, a Revelio object is created:
```
library(Revelio)
myData <- createRevelioObject(rawData = revelioTestData_rawDataMatrix,
                              cyclicGenes = revelioTestData_cyclicGenes)
```
Next, the provided table of marker gene lists is used to computationally infer cell cycle phases for each cell:
```
myData <- getCellCyclePhaseAssignInformation(dataList = myData)
```
The data is now prepared for analysis which is done in two steps. First, a standard PCA is performed and afterwards an additional linear rotation is applied to the data:
```
myData <- getPCAData(dataList = myData)
myData <- getOptimalRotation(dataList = myData)
```
If results should be plotted right away, the parameter 'boolPlotResults' has to be set TRUE:
```
myData <- getPCAData(dataList = myData, boolPlotResults = TRUE)
myData <- getOptimalRotation(dataList = myData, boolPlotResults = TRUE)
```

Due to the linearity of our algorithm and the fact that cell cycle effects are functionally isolated into two dimensions, we can now proceed to remove cell cycle effects from the normalized data matrix:
```
normalizedDataWithoutCCEffects <- removeCCEffects(dataList = myData)
```
The output is a gene-by-cell data matrix of normalized counts where cell cycle effects have been removed. This data can now be utilized for further analysis (such as doing another PCA or differential gene expression analysis).

## References

Schwabe D, Formichetti S, Junker JP, Falcke M, Rajewsky N. The transcriptome dynamics of single cells during the cell cycle. Mol Syst Biol. 2020 Nov;16(11):e9946. doi: 10.15252/msb.20209946. PMID: 33205894.
