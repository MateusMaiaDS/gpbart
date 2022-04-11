# GP-BART: Gaussian Processses Bayesian Additive Regression Trees
Implementation of GP-BART algorithm.

## Installation

To install the **GP-BART** package `GP-BART` in the development version from GitHub, just use the command:

```r
# install.packages("devtools")
devtools::install_github("MateusMaiaDS/gpbart")
```


To run an example of the GP-BART algorithm just ran the following code

```r
# This file is just to measure if the main code for GP-BART is running well, isn't made for 
# evaluation purposes.
library(gpbart)

# Importing data and library
rm(list=ls())

# Setting a seed
set.seed(42)

# Getting a simple train and test dataset
train_data <- mlbench::mlbench.friedman1(n = 50, sd = 1)
x_train <- train_data$x
y_train <- train_data$y

test_data <- mlbench::mlbench.friedman1(n = 10, sd = 1)
x_test <- test_data$x
y_test <- test_data$y
colnames(x_train) <- colnames(x_test) <- paste0("x.", seq_len(ncol(x_train)))

# Running the model
gpbart_mod <- gp_bart(x = x_train,
                          y = y_train,
                          number_trees = 20, node_min_size = 5, kappa = 0.5,
                          n_iter = 2000, alpha = 0.5, beta = 20, scale_boolean = TRUE)

# Running the prediction
pred_gpbart <- predict(rBart_model = gpbart_mod, x_test = x_test, type = "mean")


```
