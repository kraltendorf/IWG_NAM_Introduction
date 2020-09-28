# IWG_NAM_Introduction

**Introduction to the intermediate wheatgrass (IWG) nested association mapping (NAM) population. Submitted for publication to G3.**

Each analysis folder includes three sub-directories with "scripts", "data", and "output" and includes all files necessary to reproduce the analyses.


### Step 1: Growing Degree Calculation
Using data from NOAA and The Land Institute Weather Station, calculate and plot the accumulation of growing degree days for each of the environments.

### Step 2: Phenotypic Data Analysis
Using phenotypic emmeans derived from a previous analysis (see: IWG_Yield_Components), I calculate parental means, determine whether variable parents predict variable progeny, run a t-test between progeny derived between the common and donor parents, generate boxplots of variation within families, and model the relationship between the two flowering time traits. I also look at the difference between selfed and outcrossed progeny, and prepare .qua files for MapQTL. 

### Step 3: Variant Calling
Here I follow the GBarleyS pipeline developed by @neyhartj to call SNPs from GBS sequence data using GATK. The pipeline was adapted for IWG where two barcodes are used in GBS library development.

### Step 4: Variant Filtering and Imputation
SNPs are filtered using VCFTools and imputed using LinkImpute.

### Step 5: Additional Filtering and JoinMap File Prep
Here I filter for "impossible" allele calls based on parental genotypes and allele frequencies and prepare files in the JoinMap format. 

### Step 6: Population Structure Analysis
Here I create a PCA of the population's relationship matrix and include scripts and output from the Structure analysis.

### Step 7: Linkage Disequillibrium
First I calculate pairwise LD within families and model decay using a spline approach over cM.

### Step 8: JoinMap
After linkage map creation, files are appended together and converted from the cross pollinated format to the two-way psuedo testrcoss, or doubled haploid model for use in MapQTL.

### Step 9: Marker Summary
Determine how many markers were in common with the previous map, determine how many were significantly distorted, and how many ended up on the final map. 

### Step 10: GWAS
Here I run GAPIT for each unique trait and environment combination, and narrow QTL windows using the Sommer package in R.

### Step 11: QTL Linkage Mapping
Here I calculate LOD thresholds from output in MapQTL and calcualte 2-LOD dropoff intervals for QTL mapping results.

### Step 12: Visualization of QTL Mapping Results
Finally, I plot the results from both analyses using the R program LinkageMapView. 
