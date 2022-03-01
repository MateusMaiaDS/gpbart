# Gaussian Process BART.
Implementation of GP-BART algorithm in R language, from the article entitled "GP-BART: a novel Bayesian Additive Regression Trees approach using Gaussian Processes" by Mateus Maia, Keefe Murphy, and Andrew C. Parnell.

To run all functions necessary, it's just necessary to source the codes and require the libraries from the code snippet below:

```{r}
# Importing libraries and function
rm(list=ls())
library(mlbench)
library(Rcpp)

source("gpbart.R")
source("bart.R")
source("tree_manipulation_objects.R")
source("common_help_functions.R")
source("fast_gp_single_tau.R")
sourceCpp("dist_matrix.cpp")
```

To see a demonstration of how you could use the data just check the files `unit_test.R` and `unit_test_gpbart.R` to check the models of BART and GP-BART respectively.
