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
rm(list=ls())
# Loading package
library(gpbart)

# Setting a seed
set.seed(42)

# Creating a simple example
x <- sort(runif(n = 100,min = -pi,max = pi))
x <- as.matrix(x)
colnames(x) <- "x"
y <- sin(x)+rnorm(n = length(x),sd = 0.1)

# Running the model
gpbart_mod <- gpbart::gp_bart(x = x,y = y,number_trees = 5,kappa = 0.5,
                              beta = 20,alpha = 0.9)
rBart_model <- gpbart_mod
x_test <- sort(runif(n = 100,min = -pi,max = pi))
x_test <- as.matrix(x_test)
colnames(x_test) <- colnames(x)
pred_gpbart <- predict(gpbart_mod,x_test = x,type = "all")

# Comparing the up sd from the quantile with the from \tau
up_pi <- apply(pred_gpbart$out$pred ,2,function(x)quantile(x,probs = c(0.75)))
low_pi <-  apply(pred_gpbart$out$pred,2,function(x)quantile(x,probs = c(0.25)))

# Plotting the result
plot(x,y,pch=20)
lines(x,colMeans(pred_gpbart$out$pred), col = "red")
lines(x,low_pi, col = "red", lty = "dashed")
lines(x,up_pi, col = "red", lty = "dashed")



```

## References
For more details: [**Maia, Mateus, Keefe Murphy, and Andrew C. Parnell. "GP-BART: a novel Bayesian additive regression trees approach using Gaussian processes." arXiv preprint arXiv:2204.02112 (2022)**](https://doi.org/10.48550/arXiv.2204.02112).
