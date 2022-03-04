# Creating the function that it will return a matrix with all the metrics
source("bart.R")
source("gpbart.R")
source("tree_manipulation_objects.R")


gpbart_comparison <- function(cross_validation_data, # A list of cross-validation objects
                              n_rep = 1,K,beta,
                              phi_sample,rotation_boolean,
                              n_iter_value, kappa_value
){
  
  # Creating the object that it will be return after this iteration
  cross_validation_matrix <- data.frame(n_rep = n_rep,
                                        model = rep(c("gpbart","bart","softbart","tgp")),
                                        rmse = NA,
                                        pi = NA,
                                        pred = NA)
  
  # Associating the training and test data
  x_train <- cross_validation_data[[ n_rep ]]$x_train
  x_test <- cross_validation_data[[ n_rep ]]$x_test
  y_train <- cross_validation_data[[ n_rep ]]$y_train
  y_test <- cross_validation_data[[ n_rep ]]$y_test
  
  
  # Calculating the predictions
  gpbart_pred_mean <- numeric()
  bart_pred_mean <- numeric()
  softbart_pred_mean <- numeric()
  tgp_pred_mean <- numeric()
  
  
  # Setting the models parameters
  beta_value <- beta

  # # Generating the GP_BART model
  gpbart_mod <- my_rBart_gp(x = x_train,
                            y = y_train,
                            number_trees = K, p = 1,
                            rotation = rotation_boolean,node_min_size = 15,
                            beta = beta_value,alpha = 0.5,
                            discrete_phi_boolean = FALSE,
                            mu = 1,
                            tau = 1,
                            a_tau = 5,
                            d_tau = rate_tau(x = x_train,
                                             y = y_train,prob = 0.9,shape = 3),
                            phi_vector = rep(1,K),
                            predictions = matrix(rep((y_train/K),K),nrow = K,byrow = TRUE),
                            n_iter = n_iter_value,x_scale = TRUE,
                            scale_boolean = TRUE,nu_multiplier = 1,
                            nu_update = FALSE, phi_update = TRUE,
                            K_bart = 2,
                            kappa = kappa_value,gp_variables = c("lat","lon"),
                            bart_boolean = TRUE,bart_number_iter = 250,
                            burn = 500)
  
  
  # Predicting the GP-BART model
  pred_gpbart <- predict(rBart_model = gpbart_mod,x_test = x_test,
                         type = c("all"),pred_bart_only = FALSE)
  
  
  # For BART case only
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(rmse = ifelse((model == "gpbart" & n_rep == n_rep),
                                rmse(obs = y_test, pred = colMeans(pred_gpbart$out$pred)),rmse))
  
  # Storing predictions
  gpbart_pred_mean <- colMeans(pred_gpbart$out$pred)
  
  # Hyperparameters calibration
  N_test <- length(y_test) # Getting the number of observations on the test set
  
  low_ci <- pred_gpbart$out$pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[1,]
  
  up_ci <- pred_gpbart$out$pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[2,]
  
  
  # -- Calculating the callibration measure
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(pi = ifelse((model == "gpbart" & n_rep == n_rep),
                              sum(  (low_ci <= y_test)  &
                                      (y_test <= up_ci))/N_test,pi))
  
  # ==  BART-model ===
  bart_model <- dbarts::bart(x.train = x_train,y.train = y_train,ntree = 200,keeptrees = TRUE)
  
  
  bart_pred <- predict(bart_model,newdata = x_test)
  
  bart_pred_mean <- colMeans(bart_pred)
  
  
  # Predciting the RMSE for BART
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(rmse = ifelse((model == "bart" & n_rep == n_rep),
                                rmse(obs = y_test,pred = bart_pred_mean),
                                rmse))
  
  # Calculating the PI for bart
  low_ci <- bart_pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[1,]
  
  # -- Getting the up CI
  up_ci <- bart_pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[2,]
  
  # -- Calculating the callibration measure
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(pi = ifelse((model == "bart" & n_rep == n_rep),
                              sum(low_ci <= y_test &
                                    y_test <= up_ci)/N_test,pi))
  
  # =====
  # Getting the SoftBART package
  # =====
  library(SoftBart)
  if(dim(x_train)[2]==1){
    softbart_mod <- SoftBart::softbart(X = cbind(x_train,x_train), Y = y_train,X_test = cbind(x_test,x_test))
  } else {
    softbart_mod <- SoftBart::softbart(X = x_train, Y = y_train,X_test = x_test )
  }
  
  softbart_pred <- softbart_mod$y_hat_test
  
  softbart_pred_mean <- colMeans(softbart_pred)
  
  # Predciting the RMSE for softBART
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(rmse = ifelse((model == "softbart" & n_rep == n_rep),
                                rmse(obs = y_test,pred = softbart_pred_mean),
                                rmse))
  
  # cross_validation_matrix
  
  # Preciting the calibration for the softBART
  low_ci <- softbart_pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[1,]
  
  # -- Getting the up CI
  up_ci <-softbart_pred %>% apply(2, function(x){ quantile(x, probs = c(0.25,0.75))}) %>% .[2,]
  
  # -- Calculating the callibration measure
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(pi = ifelse((model == "softbart" & n_rep == n_rep),
                              sum(low_ci <= y_test &
                                    y_test <= up_ci)/N_test,pi))
  
  
  # Using Treed Gaussian Processes patterns
  mod_tgp <- tgp::btgp(X = x_train,Z = y_train,XX = x_test)
  
  
  # Storing the tgp
  tgp_pred_mean <- mod_tgp$ZZ.mean
  
  # Predciting the RMSE for tGP
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(rmse = ifelse((model == "tgp" & n_rep == n_rep),
                                rmse(obs = y_test,pred = mod_tgp$ZZ.mean),
                                rmse))
  
  # Preciting the calibration for the softBART
  low_ci <- mod_tgp$ZZ.mean-
    0.67*sqrt(mean(mod_tgp$ZZ.s2))
  
  # -- Getting the up CI
  up_ci <- mod_tgp$ZZ.mean+
    0.67*sqrt(mean(mod_tgp$ZZ.s2))
  
  # -- Calculating the callibration measure
  cross_validation_matrix <- cross_validation_matrix %>% 
    dplyr::mutate(pi = ifelse((model == "tgp" & n_rep == n_rep),
                              sum(low_ci <= y_test &
                                    y_test <= up_ci)/N_test,pi))
  
  
  # cross_validation_matrix
  
  # return(as.matrix(cross_validation_matrix))
  return(list(cross_validation_matrix = cross_validation_matrix,
              gpbart_pred = gpbart_pred_mean,
              bart_pred = bart_pred_mean,
              softbart_pred = softbart_pred_mean,
              tgp_pred = tgp_pred_mean,
              cross_validation_data = cross_validation_data,
              gpbart_mod = gpbart_mod))
  
}

