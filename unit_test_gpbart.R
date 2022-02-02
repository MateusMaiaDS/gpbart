# Importing data and library
rm(list=ls())
library(mlbench)
source("gpbart.R")
source("bart.R")
source("tree_manipulation_objects.R")

# Getting a simple train and test dataset
train_data <- mlbench.friedman1(n = 50,sd = 1)
x_train <- train_data$x
y_train <- train_data$y

test_data <- mlbench.friedman1(n = 10, sd = 1)
x_test <- test_data$x
y_test <- test_data$y
colnames(x_train) <- colnames(x_test) <- paste0("x.",1:ncol(x_train))
gpbart_mod <- my_rBart_gp(x = x_train,
                          y = y_train,
                          number_trees = 20,node_min_size = 5,kappa = 1,
                          n_iter = 2000,bart_number_iter = 2000)

gpbart_mod$nu_store<- matrix(1e-10, nrow = nrow(gpbart_mod$nu_store),ncol = ncol(gpbart_mod$nu_store))

pred_gpbart <- my_predict_rBART(rBart_model = gpbart_mod, x_test = x_test,type  = c("mean"),pred_bart_only = TRUE)
crossprod( (pred_gpbart$out$pred-y_test))

# x_train <- cross_validation_object[[1]]$x_train
# y_train <- cross_validation_object[[1]]$y_train
# x_test <- cross_validation_object[[1]]$x_test

dbart_mod <- dbarts::bart(x.train = x_train,
                          y.train = y_train,
                          ntree = 20,keeptrees = TRUE)

crossprod( (y_test- predict(dbart_mod,newdata = x_test) %>% apply(2,mean)))

pred_gpbart$out$pred
predict(dbart_mod,newdata = x_test) %>% apply(2,mean)
y_test


points(x_test[,1],pred_gpbart$list_matrix_pred[[1]][1,], pch = 20, col = "blue")
plot(x_train[,1],gpbart_mod$current_partial_residuals_list[[1]][1,],pch = 20)
