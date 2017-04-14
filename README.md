# CAGI-p16-assessment

## Overview

This repository contains all the R scrips used for the assessment of the CAGI-3 ["p16 challenge"](https://genomeinterpretation.org/content/predict-how-variants-p16-tumor-suppressor-protein-affect-cell-proliferation).

All the scripts can be reused for performance assessment of bioinformatics tools to predict phenotypic effects of genetic variants of unknown significance (VUS).

## Dependencies
The analysis requires the following R packages to be installed

* [ROCR](https://cran.r-project.org/web/packages/ROCR/index.html)
* [plotrix](https://cran.r-project.org/web/packages/plotrix/index.html)

## Usage

Each script can be run from a terminal as the example below

``` 
  Rscript 1_main_numerical_assessment.R 
```

## Details

In this section a brief description of each script is given

Script | Description
------------ | -------------
1 | calculates the main numerical measures (i.e. correlations, AUC, RMSD).In additions it produces the tables needed by other scripts to generate assessment figures and tables.
2 | calculates correlation measure among all predictions and produce an heatmap figure to visualize results
3 | calculates correlation measure among performance indices and produce an heatmap figure to visualize results
4 | calculates the pairwise significance of challenge evaluations and produce an heatmap figure
5 | draws experimental values versus predicted values graph
6 | calculates only PSWD10 and produce a table to identify difficut targets
