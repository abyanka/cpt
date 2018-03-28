#!/bin/bash

mkdir output/figures

### Replicate results/analysis of Lyall (2009) application
R CMD BATCH --vanilla '--args 1000' lyall_2009_analysis.R 

### MPs for Sale application
# performs preliminary computationally intensive estimation
# First argument is number of random permutations, second argument is the number of different window sizes to evaluate at.
R CMD BATCH --vanilla '--args 1000 80' MPs.R 

# less computationally intensive estimation and figure/tables
R CMD BATCH --vanilla mps_analysis.R 

### Judges, Green and Winik (2010) 
R CMD BATCH --vanilla '--args 1000' judges_analysis.R 

### Community College, Rouse (2010) 
R CMD BATCH --vanilla '--args 1000' community_college_analysis.R

### Sims
R CMD BATCH --vanilla '--args 500 400' sims.R
R CMD BATCH --vanilla simplots.R

cp *.Rout output
cp *.rda output



