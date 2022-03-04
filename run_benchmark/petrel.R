# Cleaning the workspace
rm(list=ls())
# Importing all the packages and functions
library(tidyverse)
library(doParallel)
library(foreach)
library(spData)

# Sourcing the kfold_file
source("run_benchmark/spatial_cross_validation.R") # Spatial K-Fold value

# petrel
load("run_benchmark/data/petrel.Rdata")
database<-petrel


# Setting the cross-validation parameterisation
N <-  1000
n_iter <- 2000
K <- 20
beta <- 2
rotation_boolean <- TRUE
phi_sample <- TRUE
kappa_value <- 0.5

cross_validation_object <- spatial_k_fold(data = database,
                                          dependent_variable = "y",
                                          k_partitions = 10,seed = 42)


# Separating the clusters
number_cores <- 1 # Take care with this line to not break your computer
cl <- makeCluster(number_cores)
doParallel::registerDoParallel(cl)


# Testing the simple K
simple_K_example <- foreach(i = 1:10, .packages = c("BART","SoftBart","tgp")) %dopar%{
  
  # Loading the functions
  source("run_benchmark/fast_cross_validation_saving_models.R") # Getting the function that it will apply the cv
  source("gpbart.R") # New simulated scenarios
  source("fast_gp_single_tau.R") # GP-BART function
  
  aux <- gpbart_comparison(cross_validation_data = cross_validation_object,
                           n_rep = i,K = K, beta = beta, 
                           phi_sample = phi_sample,rotation_boolean = rotation_boolean,
                           n_iter_value = n_iter, kappa_value = kappa_value)
  # BART
  # GP-BART
  # soft-BART
  # tGP
}

stopCluster(cl)

saveRDS(object = simple_K_example,
        file = paste0("run_benchmark/CORRECTED_new_rotation_FINAL_scale_gp_spatial_parallel_march_petrel_gpbart_mod_result_",N,"_K_",K,"_beta_",round(beta,0),"_phi_sample_",
                      phi_sample,
                      "rotation_boolean_",rotation_boolean,"_n_rep_",n_iter,".Rds"))
