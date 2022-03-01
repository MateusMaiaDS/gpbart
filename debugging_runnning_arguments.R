# Getting a generic dataset
library(mlbench)
friedman_example <- mlbench.friedman1(n = 50)
x <- friedman_example$x # Covariate matrix
y <- friedman_example$y # Target variable
number_trees = 2 # Number of trees
control = list(node_min_size = 5,
               scale_boolean = TRUE)

prior_parameter = list(a_tau = 2, # Parameters from the prior
                       d_tau = 3,
                       k_bart = 2,
                       alpha = 0.95,
                       beta = 2,
                       prob_tau = 0.9)

mcmc_parameter = list(n_iter = 2000, # Parameters from MCMC
                      burn = 1000,
                      thin = 1)

init = list(tau = 1, # Initial values
            mu = 0)

# list(x = x, # Covariate matrix
#     y = y_scale, # Target variable
#     tree_store = tree_store, # Tree store
#     residuals_store = residuals_store,
#     predictions_store = predictions_store,
#     tau_store = tau_store,
#     y_hat_store = y_hat_store,
#     number_trees = number_trees, # Number of trees
#     control = list(node_min_size = node_min_size,
#                    scale_boolean = scale_boolean),
#     prior_parameter = list(a_tau = a_tau, # Parameters from the prior
#                            d_tau = d_tau,
#                            k_bart = k_bart,
#                            alpha = alpha,
#                            beta = beta,
#                            prob_tau = prob_tau),
#     mcmc_parameter = list(n_iter = n_iter, # Parameters from MCMC
#                           burn = burn,
#                           thin = thin))

