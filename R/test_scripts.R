# rm(list=ls())
# # source(file = "R/gpbart.R")
# library(gpbart)
# set.seed(42)
# # # Creating a simple example
# x <- sort(runif(n = 100,min = -pi,max = pi))
# x <- as.matrix(x)
# colnames(x) <- "x"
# y <- sin(x)+rnorm(n = length(x),sd = 0.1)
# 
# # Running the model
# gpbart_mod <- gpbart::gp_bart(x = x,y = y,number_trees = 2,kappa = 0.5,
#                               beta = 20,alpha = 0.9,rotation = FALSE,scale_boolean = TRUE)
# rBart_model <- gpbart_mod
# x_test <- x
# # x_t
# pred_gpbart <- predict(gpbart_mod,x_test = x_test,type = "mean")
# 
# # Comparing the up sd from the quantile with the from \tau
# # up_pi <- pred_gpbart$out$pred %>% apply(2,function(x)quantile(x,probs = c(0.75)))
# # mean_pi <- pred_gpbart$out$pred %>% colMeans()
# 
# # SD from pi
# (up_pi-mean_pi)/qnorm(0.75) %>% mean
# 
# # SD from \tau
# # pred_gpbart$out$sd %>% mean
# 
# plot(x,y,pch=20)
# lines(x,rBart_model$y_hat %>% colMeans())
# lines(x_test,pred_gpbart$out$pred,col = "blue")
