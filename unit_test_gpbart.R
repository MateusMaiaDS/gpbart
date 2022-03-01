# This file is just to measure if the main code for GP-BART is running well, isn't made for 
# evaluation purposes.

# Setting a seed
set.seed(42)

# Importing data and library
rm(list=ls())
library(mlbench)
library(Rcpp)

source("gpbart.R")
source("bart.R")
source("tree_manipulation_objects.R")
source("common_help_functions.R")
source("fast_gp_single_tau.R")
sourceCpp("dist_matrix.cpp")

# Getting a simple train and test dataset
train_data <- mlbench.friedman1(n = 50, sd = 1)
x_train <- train_data$x
y_train <- train_data$y

test_data <- mlbench.friedman1(n = 10, sd = 1)
x_test <- test_data$x
y_test <- test_data$y
colnames(x_train) <- colnames(x_test) <- paste0("x.", seq_len(ncol(x_train)))

# Running the model
gpbart_mod <- my_rBart_gp(x = x_train,
                          y = y_train,
                          number_trees = 20, node_min_size = 5, kappa = 0.5,
                          n_iter = 2000, alpha = 0.5, beta = 20, scale_boolean = TRUE)

# Running the prediction
pred_gpbart <- my_predict_rBART(rBart_model = gpbart_mod, x_test = x_test, type = "mean")

# Running the standard from dBARTs
dbart_mod <- dbarts::bart(x.train = x_train,
                          y.train = y_train,
                          ntree = 20, keeptrees = TRUE)

crossprod(pred_gpbart$out$pred - y_test)
crossprod(y_test - predict(dbart_mod,newdata = x_test) %>% colMeans)

pred_gpbart$out$pred
predict(dbart_mod,newdata = x_test) %>% colMeans
y_test