### All helper functions which are common to bart.R & gpbart.R & some other helper functions
# other common functions related to tree structures are found in
source("tree_manipulation_objects.R")

tree_prior <- function(tree, alpha, beta) {
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
  
  # Selecting internal nodes names
  names_internal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 0))
  
  # Selecting the depth of the terminal nodes
  depth_terminal <- vapply(tree[names_terminal_nodes], "[[", numeric(1), "depth_node")
  
  # Selecting the depth of the internal nodes
  depth_internal <- vapply(tree[names_internal_nodes], "[[", numeric(1), "depth_node")
  
  # Case for stump (No internal node)
  if(length(depth_internal) == 0) {
    log_p <- log(1 - alpha)
  } else {
    # Calculating the log-likelihood
    log_p <- sum(log(1 - alpha * (1 + depth_terminal)^(-beta))) + sum(log(alpha) - beta * log1p(depth_internal))
  }
    return(log_p)
}

# Update tau_j values
update_tau <- function(x, 
                       y,
                       a_tau,
                       d_tau, 
                       predictions) {
  
  # Calculating the values of a and d
  n <- nrow(x)
  
  # Getting the shape parameter from the posterior
  shape_tau_post <- 0.5 * n + a_tau
  
  # Getting the ratio parameter
  rate_tau_post <- 0.5 * crossprod(y - predictions) + d_tau
  
  # Updating the \tau
  tau_sample <- rgamma(n = 1, shape = shape_tau_post, rate = rate_tau_post)
    return(tau_sample)
}

# Functions to find the zero for tau
zero_tau_prob <- function(x, naive_tau_value, prob, shape) {
  
  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined   
  return(pgamma(naive_tau_value,
                shape = shape,
                rate = x) - (1 - prob))
}

zero_tau_prob_squared <- function(rate, naive_tau_value, prob, shape) {
  
  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined   
  return((pgamma(naive_tau_value,
                 shape = shape,
                 rate = rate) - (1 - prob))^2)
}

# Naive tau_estimation
naive_tau <- function(x, y) {
  
  # Getting the valus from n and p
  n <- length(y)
  
  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))
  
  # Naive lm_mod 
  lm_mod <- lm(formula = y ~ ., data =  data.frame(y, x))
  
  sigma <- sqrt(sum((lm_mod$residuals)^2)/(n - p))
  
  tau <- 1/sigma^2
    return(tau)
}

# Naive sigma_estimation
naive_sigma <- function(x,y){
  
  # Getting the valus from n and p
  n <- length(y)
  
  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))
  
  # Naive lm_mod 
  lm_mod <- lm(formula = y ~ ., data =  data.frame(y,x))
  
  sigma <- sqrt(sum((lm_mod$residuals)^2)/(n - p))
    return(sigma)
}

# Return rate parameter from the tau prior
rate_tau <- function(x, # X value
                     y, # Y value
                     prob = 0.9,
                     shape) {
  # Find the tau_ols
  tau_ols <- naive_tau(x = x,
                       y = y)
  
  # Getting the root
  min_root <-  try(uniroot(f = zero_tau_prob, interval = c(1e-2, 100),
                           naive_tau_value = tau_ols,
                           prob = prob, shape = shape)$root, silent = TRUE)
  
  if(inherits(min_root, "try-error")) {
    # Verifying the squared version
    min_root <- optim(par = runif(1), fn = zero_tau_prob_squared,
                      method = "L-BFGS-B", lower = 0,
                      naive_tau_value = tau_ols,
                      prob = prob, shape = shape)$par
  }
    return(min_root)
}

# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {
  
  # Defining the a and b
  a <- min(y)
  b <- max(y)
  
  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
    return(y) 
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
    return(y)
}

rmse <- function(obs, pred) {
  return(sqrt(mean((obs - pred)^2)))
}

rMVN_var <- function(mean, Sigma) {
  if(is.matrix(Sigma)) {
    drop(mean + crossprod(PD_chol(Sigma), rnorm(length(mean))))
  } else {
    mean + sqrt(Sigma) * rnorm(length(mean))
  }
}

is_diag_matrix <- function(m) all(m[!diag(nrow(m))] == 0)

PD_chol  <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
    d    <- nrow(x)
    eigs <- eigen(x, symmetric = TRUE)
    eval <- eigs$values
    evec <- eigs$vectors
      return(chol(x + evec %*% tcrossprod(diag(pmax.int(1e-8, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec), ...))
  }
)
