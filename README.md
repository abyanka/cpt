# Replication Code

This repository contains code and data to replicate the analyses in the manuscript.  More specifically, the contents of the repository are:
* The cpt R package, which implements the CPT
* The cpt.paper R package, which contains the datasets and also some helper functions
* The Replication folder, which contains the R scripts that perform the analyses

To replicate the analyses:
1. Install the cpt and cpt.paper R packages
2. Ensure that all other required R packages are installed (see list below)
3. Change into the Replication folder and run script_run.sh

On a 6-core Ryzen CPU, the code takes about 1 week to run.

The required R packages are:

install.packages("rgenoud","Matching","Hotelling","nnet","randomForest",
"glmnet","tidyverse","energy","crossmatch","MASS","reshape","ggplot2","ggthemes",
"rdd","xtable","stargazer","BalanceCheck","asbio","rpart","ggdendro")
