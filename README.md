# BayesianACMimputation
Bayesian imputation project for autoregression correlation models

# scripts here are adapted from the original http://www.spatial-econometrics.com/ Toolbox for Matlab into R code.

# R
gibbs_zsamplingfunctions.r = support functions for the Gibbs sampling procedure, 
  includes lots of functions that are defined separately in the spatial econometrics toolbox. 
  Function probit_g() is for Gibbs sampling (equal to probit_g.m) and
  function probit_g_mp is for Gibbs sampling with imputation (equal to probit_g_mp.r
 
gibbs_logistic.r = R script to conduct Gibbs sampling 
  This script includes the data of n 100 and k 3 generated in matlab (with seed as: rng(2017)) to make sure 
  that the R results are identical to the matlab results (small deviations due to random Gibbs sampling aside).

# Matlab 
probit_g.m = support function probit_g() for the Gibbs sampling procedure
probit_g_mp.m = support function probit_g_mp() for the Gibbs sampling with multiple imputation procedure
probit_gd.m = demo script to create dataset, define priors, call the function and describe output
