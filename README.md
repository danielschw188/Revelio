# Revelio

The Revelio package is designed to extract cell cycle information from scRNA-seq experiments of cycling cells. What sets Revelio apart is that the transformation from normalized data to final outcome is fully linear. The advantage is that any conclusion drawn on the transformed data can be transferred back to the original data since the algorithm is easily inverted. One simple application of this principle is the removal of cell cycle effects from the normalized data matrix.

## Installation
```
devtools::install_github('danielschw188/Revelio')
```

## How to Use
### Cell Cycle only, no intronic data required

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

### Velocity analysis, intronic data required

For the velocity analysis, the algorithm requires another input from the user:  
3. 'rawIntronicData': A gene-by-cell data matrix that contains the UMI counts as entries. Formatting should be the same as 'rawData'. Additionally, the cell IDs should match the cell IDs found in 'rawData'.

We have included a sample data into the package as well, called 'revelioTestData_rawIntronicMatrix'. Together with the data from the previous section, this data is illustrated in Figure 2 of the publication. Again, please adjust your variable name when using your own data.

Ideally, you have already run the previous section. Now, we include the intronic data into the data object and generate an RNA velocity matrix as described in model I in La Manno et al. 2018 'RNA velocity of single cells':
```
myData <- getVelocityMatrix(dataList = myData,
                            rawIntronicData = revelioTestData_rawIntronicMatrix)
```
This calculation can take up to 30 minutes depending on the size of your data. We now rerun the PCA and rotation steps from the previous sections:
```
myData <- getPCAData(dataList = myData)
myData <- getOptimalRotation(dataList = myData, boolPlotResults = TRUE)
```
We then calculate the extrapolated state within the dynamical component space by utilizing model I from La Manno et al. 2018 and transforming the resulting extrapolated state into the dynamical component space with the previously calculated rotation matrix:
```
myData <- getExtrapolatedStateInDCSpace(dataList = myData)
```
Next, we can overlay the 2D cell cycle with a grid and calculate average velocities at each grid point:
```
myData <- getVelocityGridCC(dataList = myData,
                            sigma = 0.6,
                            numberOfGridPointsOneSide = 20,
                            boolPlotResults = TRUE)
```
This is also computationally expensive. In the publication we used a 50x50 grid ('numberOfGridPointsOneSide' has to be set to 50) which might require an hour of computation. We can also look at the averaged velocity along an averaged trajectory:
```
myData <- plotVelocityOnAveragedTrajectory(dataList = myData,
                                           numberOfAnchors = 10)
```
Lastly, we compute the RNA velocity along the cylinder axis:
```
myData <- getVelocityGridAlongThirdDimensions(dataList = myData,
                                              whichThirdDimension = 3,
                                              sigma = 1.5,
                                              numberOfGridPointsOneSide = 20,
                                              boolPlotResults = TRUE)
```
Again, we used a 50x50 grid for the publication which can take up to 1.5 hours.


## References

Schwabe D, Formichetti S, Junker JP, Falcke M, Rajewsky N. The transcriptome dynamics of single cells during the cell cycle. Mol Syst Biol. 2020 Nov;16(11):e9946. doi: 10.15252/msb.20209946. PMID: 33205894.
