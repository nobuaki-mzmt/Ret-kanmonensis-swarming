# README
## Article Information
This repository provides access to the data and source code used for the manuscript   
### **The significance of social interactions in synchronized swarming flight in a termite**    
Nobuaki Mizumoto, Tomonari Nozaki  
<!-- Author names are commented out for DBR.

Preprint will be available at bioRxiv. [![DOI:XXX](http://img.shields.io/badge/DOI-10.1101/XXX.svg)]  
The all data will be uploaded in Zenodo upon acceptance: [![DOI](https://zenodo.org/badge/DOI/XXXDOIXXX.svg)](https://doi.org/XXXDOIXXX) -->

This study proposes that synchronized termite swarming results from collective decision-making within a group of alates. We observed the swarming behavior in _Reticulitermes kanmonensis_ under both semi-natural and laboratory conditions to confirm 1) termites suppress minor dispersal flights under lower temperatures, 2) they can synchronize flight even without environmental cues, 3) group size facilitates swarming.
This includes data obtained from empirical observations and R scripts to analyze them.

## Table of Contents
* [README](./README.md)
* [scripts](./scripts)
  * [output.R](./scripts/output.R) - output all plots and statistical analysis of empirical works
  <!-- * [cdm_model.R](./scripts/cdm_model.R) - simple simulation model for collective decision making. ourput all plots. -->
* [output](./output) - all outputs are stored
* [data](./data)
  * [raw](./data/raw) - raw data in .csv    

## Session information
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] plyr_1.8.8          viridis_0.6.3       viridisLite_0.4.2  
 [4] ggplot2_3.4.2       fitdistrplus_1.1-11 multcomp_1.4-25    
 [7] TH.data_1.1-2       survival_3.5-5      mvtnorm_1.2-3      
[10] lme4_1.1-34         Matrix_1.6-1        car_3.1-2          
[13] carData_3.0-5       MASS_7.3-60        

