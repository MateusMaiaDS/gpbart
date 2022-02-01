# Importing data and library
rm(list=ls())
library(mlbench)
source("bart.R")
source("tree_manipulation_objects.R")
# Getting a simple train and test dataset
train_data <- mlbench.friedman1(n = 50,sd = 1)
x_train <- train_data$x
y_train <- train_data$y

test_data <- mlbench.friedman1(n = 10, sd = 1)
x_test <- test_data$x
y_test <- test_data$y


bart_mod <- bart(x = x_train,y = y_train,number_trees = 20)
bart_pred <- predict_bart(bart_mod = bart_mod, newdata = x_test,type = "mean")

crossprod( (bart_pred-y_test))

x_train <- cross_validation_object[[1]]$x_train
y_train <- cross_validation_object[[1]]$y_train
x_test <- cross_validation_object[[1]]$x_test

dbarts::bart(x.train = x_train,
             y.train = y_train,
             ntree = 20)