rm(list=ls())
library(gpbart)
set.seed(42)
# # Creating a simple example
x <- sort(runif(n = 100,min = -pi,max = pi))
x <- as.matrix(x)
colnames(x) <- "x"
y <- sin(x)+rnorm(n = length(x),sd = 0.01)

# Running the model
gpbart_mod <- gpbart::gp_bart(x = x,y = y,number_trees = 1,kappa = 0.5,node_min_size = 5,
                              beta = 20,alpha = 0.9,
                              rotation = FALSE,scale_boolean = TRUE,x_scale = TRUE)
rBart_model <- gpbart_mod
x_test <- sort(runif(n = 100,min = -pi,max = pi))
x_test <- as.matrix(x_test)
colnames(x_test) <- colnames(x)
pred_gpbart <- predict(gpbart_mod,x_test = x,type = "all")

# Comparing the up sd from the quantile with the from \tau
up_pi <- pred_gpbart$out$pred %>% apply(2,function(x)quantile(x,probs = c(0.75)))
low_pi <- pred_gpbart$out$pred %>% apply(2,function(x)quantile(x,probs = c(0.25)))

# mean_pi <- pred_gpbart$out$pred %>% colMeans()

# SD from pi
(up_pi-mean_pi)/qnorm(0.75) %>% mean

# SD from \tau
# pred_gpbart$out$sd %>% mean

plot(x,y,pch=20)
# lines(x,colMeans(rBart_model$y_hat) )
# lines(x,pred_gpbart$out$pred[120,],col = "blue")
lines(x,colMeans(pred_gpbart$out$pred), col = "red")
lines(x,low_pi, col = "red", lty = "dashed")
lines(x,up_pi, col = "red", lty = "dashed")
lines(x,colMeans(pred_gpbart$out$pred)+qnorm(0.25)*mean(pred_gpbart$out$sd), col = "red", lty = "dashed")
# lines(x,colMeans(pred_gpbart$out$pred)+qnorm(0.75)*mean(pred_gpbart$out$sd), col = "red", lty = "dashed")

