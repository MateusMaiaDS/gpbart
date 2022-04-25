# # # Creating a simple example
# x <- sort(runif(n = 100,min = -pi,max = pi))
# x <- as.matrix(x)
# colnames(x) <- "x"
# y <- sin(x)+rnorm(n = length(x),sd = 0.1)
#
# # Running the model
# gpbart_mod <- gpbart::gp_bart(x = x,y = y,number_trees = 2)
# rBart_model <- gpbart_mod
# x_test <- x[1,,drop= FALSE]
# pred_gpbart <- predict(gpbart_mod,x_test = x_test,type = "mean")
#
# # Comparing the up sd from the quantile with the from \tai
# up_pi <- pred_gpbart$out$pred %>% apply(2,function(x)quantile(x,probs = c(0.75)))
# mean_pi <- pred_gpbart$out$pred %>% colMeans()
#
# # SD from pi
# (up_pi-mean_pi)/qnorm(0.75) %>% mean
#
# # SD from \tau
# pred_gpbart$out$sd %>% mean

