## GP-Bart
#' @useDynLib gpbart
#' @importFrom Rcpp sourceCpp
# ==================================#
# Objects to test the tree_complete_conditional function
# ==================================#

tree_complete_conditional_gpbart <- function(tree, x, residuals, nu = 1, phi = 1,
                                             tau_mu,
                                             number_trees = number_trees) {

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Calculating value S
  S <- unlist(mapply(terminal_nodes, FUN=function(z, x=z$Omega_plus_I_inv) {
    if(z$is_Omega_diag) sum(diag(x)) else sum(x)
  }, SIMPLIFY = TRUE)) + tau_mu

  # Log Omega
  log_det_Omega_plus_I_tau <- vapply(terminal_nodes, function(z, x=z$Omega_plus_I_tau) {
    if(z$is_Omega_diag) sum(log(diag(x))) else determinant(x, logarithm = TRUE)$modulus
  }, numeric(1))

  # Defining RT_Omega_I_R
  RTR <- unlist(mapply(terminal_nodes, residuals_terminal_nodes,
                       FUN = function(nodes, resid, x=nodes$Omega_plus_I_inv) {
    if(nodes$is_Omega_diag) sum(resid^2 * diag(x)) else crossprod(resid, crossprod(x, resid))
  }, SIMPLIFY = FALSE))

  # The term R^{T} solve(Omega + I ) 1
  R_Omega_I_one <- unlist(mapply(terminal_nodes, residuals_terminal_nodes,
                                 FUN = function(nodes, residuals, x=nodes$Omega_plus_I_inv) {
    if(nodes$is_Omega_diag) sum(residuals * diag(x)) else rowSums(crossprod(residuals, x))
  }, SIMPLIFY = FALSE))

  log_posterior <- -0.5 * sum(log_det_Omega_plus_I_tau) - 0.5 * sum(log(S) - log(tau_mu)) - 0.5 * sum(RTR) + 0.5 * sum(R_Omega_I_one^2 / S)

    return(list(log_posterior = log_posterior,
                S = S,
                RTR = RTR,
                R_Omega_I_one = R_Omega_I_one))
}

# tree_complete_conditional_gpbart_grow <- function(old_tree,new_tree,
#                                                   x, residuals, nu = 1, phi = 1,
#                                              tau_mu,
#                                              number_trees = number_trees) {
#   
#   # Selecting the terminal nodes
#   old_terminal_nodes <- old_tree[names(which(vapply(old_tree, "[[", numeric(1), "terminal") == 1))]
#   new_terminal_nodes <- new_tree[names(which(vapply(new_tree, "[[", numeric(1), "terminal") == 1))]
#   
#   # New nodes
#   nodes_new_grow <- new_terminal_nodes[!(names(new_terminal_nodes) %in% names(old_terminal_nodes))]
#   parent_old_node <- old_terminal_nodes[paste0("node_",nodes_which_grow[[1]]$parent_node)]
#   
#   # Number of nodes
#   n_node_new <- length(nodes_new_grow)
#   
#   # Picking each node size
#   nodes_size_new <- vapply(nodes_new_grow, function(x) {
#     length(x$observations_index)
#   }, numeric(1))
#   
#   # Residuals terminal nodes
#   residuals_terminal_nodes_new_grow <- lapply(nodes_new_grow, function(x) {
#     residuals[x$observations_index]
#   })
#   
#   residuals_terminal_nodes_old_grow <- lapply(parent_old_node, function(x) {
#     residuals[x$observations_index]
#   })
#   
#   # Calculating value S
#   S_new <- unlist(mapply(nodes_new_grow, FUN=function(z, x=z$Omega_plus_I_inv) {
#     if(z$is_Omega_diag) sum(diag(x)) else sum(x)
#   }, SIMPLIFY = TRUE)) + tau_mu
#   
#   S_old <- unlist(mapply(parent_old_node, FUN=function(z, x=z$Omega_plus_I_inv) {
#     if(z$is_Omega_diag) sum(diag(x)) else sum(x)
#   }, SIMPLIFY = TRUE)) + tau_mu
#   
#   
#   # Log Omega
#   log_det_Omega_plus_I_tau_new <- vapply(nodes_new_grow, function(z, x=z$Omega_plus_I_tau) {
#     if(z$is_Omega_diag) sum(log(diag(x))) else determinant(x, logarithm = TRUE)$modulus
#   }, numeric(1))
#   
#   log_det_Omega_plus_I_tau_old <- vapply(parent_old_node, function(z, x=z$Omega_plus_I_tau) {
#     if(z$is_Omega_diag) sum(log(diag(x))) else determinant(x, logarithm = TRUE)$modulus
#   }, numeric(1))
#   
#   
#   # Defining RT_Omega_I_R
#   RTR_new <- unlist(mapply(nodes_new_grow, residuals_terminal_nodes_new_grow,
#                        FUN = function(nodes, resid, x=nodes$Omega_plus_I_inv) {
#                          if(nodes$is_Omega_diag) sum(resid^2 * diag(x)) else crossprod(resid, crossprod(x, resid))
#                        }, SIMPLIFY = FALSE))
#   
#   RTR_old <- unlist(mapply(parent_old_node, residuals_terminal_nodes_old_grow,
#                            FUN = function(nodes, resid, x=nodes$Omega_plus_I_inv) {
#                              if(nodes$is_Omega_diag) sum(resid^2 * diag(x)) else crossprod(resid, crossprod(x, resid))
#                            }, SIMPLIFY = FALSE))
#   
#   
#   
#   # The term R^{T} solve(Omega + I ) 1
#   R_Omega_I_one_new <- unlist(mapply(nodes_new_grow, residuals_terminal_nodes_new_grow,
#                                  FUN = function(nodes, residuals, x=nodes$Omega_plus_I_inv) {
#                                    if(nodes$is_Omega_diag) sum(residuals * diag(x)) else rowSums(crossprod(residuals, x))
#                                  }, SIMPLIFY = FALSE))
#   
#   R_Omega_I_one_old <- unlist(mapply(parent_old_node, residuals_terminal_nodes_old_grow,
#                                      FUN = function(nodes, residuals, x=nodes$Omega_plus_I_inv) {
#                                        if(nodes$is_Omega_diag) sum(residuals * diag(x)) else rowSums(crossprod(residuals, x))
#                                      }, SIMPLIFY = FALSE))
#   
#   log_posterior_new <- -0.5 * sum(log_det_Omega_plus_I_tau_new) - 0.5 * sum(log(S_new) - log(tau_mu)) - 0.5 * sum(RTR_new) + 0.5 * sum(R_Omega_I_one_new^2 / S_new)
#   log_posterior_old <- -0.5 * sum(log_det_Omega_plus_I_tau_old) - 0.5 * sum(log(S_old) - log(tau_mu)) - 0.5 * sum(RTR_old) + 0.5 * sum(R_Omega_I_one_old^2 / S_old)
#   
#   return(list(log_posterior_new = log_posterior_new,
#               S_new = S_new,
#               RTR_new = RTR_new,
#               R_Omega_I_one_new = R_Omega_I_one_new,
#               log_posterior_old = log_posterior_old,
#               S_old = S_old,
#               RTR_old = RTR_old,
#               R_Omega_I_one_old = R_Omega_I_one_old))
# }


# Generate mu_j values
update_mu <- function(tree,
                      x,
                      residuals,
                      likelihood_object,
                      seed = NULL) {

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Mu mean value
  mu_mean <- likelihood_object$R_Omega_I_one / likelihood_object$S

  # Calculating mu values
  mu <- mapply(mu_mean, sqrt(likelihood_object$S),
               FUN = function(x, y) {
                 stats::rnorm(
                   n = 1,
                   mean = x,
                   sd = 1/y
                 )
               }, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- mu[[names_terminal_nodes[i]]]
  }
    return(tree)
}

update_residuals <- function(tree, x, nu, phi, residuals, tau, seed = NULL) {
  # set.seed(seed)

  # New g (new vector prediction for g)
  residuals_new <- rep(NA, length(residuals))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  
  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")

  # Calculating Omega matrix
  Omega_matrix <- lapply(terminal_nodes, "[[", "Omega_matrix")
  
  # Getting Omega_matrix plust tau
  Omega_matrix_plus_tau <- lapply(Omega_matrix, function(x) { x + diag(tau^-1,nrow = nrow(x))})
  
  # Checking if diagonal
  is_Omega_diag <- lapply(terminal_nodes, "[[", "is_Omega_diag")

  # Calculating Omega matrix plus I INVERSE
  Omega_matrix_inverse_plus_I  <- lapply(terminal_nodes, "[[", "Omega_plus_I_inv")

  # Calculating g_mean posterior
  residuals_mean <- mapply(Omega_matrix,
                           Omega_matrix_inverse_plus_I,
                           residuals_terminal_nodes,
                           mu_values,
                           is_Omega_diag,
                           FUN = function(omega, omega_i_inv, residuals, mu, omega_is_diag) {
                             # Getting the values
                             if(omega_is_diag) {
                               mu + diag(omega) * diag(omega_i_inv) * (residuals - mu)
                             } else {
                               mu + crossprod(omega, crossprod(omega_i_inv, residuals - mu))
                             }
                           }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  residuals_variance <- mapply(Omega_matrix,
                               Omega_matrix_inverse_plus_I,
                               residuals_terminal_nodes,
                               mu_values,
                               is_Omega_diag,
                               FUN = function(omega, omega_i_inv, residuals, mu, omega_is_diag) {
                                 # Getting the Omega value
                                 if(omega_is_diag) {
                                   diag(omega) - diag(omega)^2 * diag(omega_i_inv) ### Old version

                                 } else {
                                   omega - crossprod(omega, crossprod(omega_i_inv, omega))  ### Old version

                                 }
                               }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  residuals_sample <- mapply(FUN=rMVN_var, residuals_mean, residuals_variance, SIMPLIFY = FALSE)
  # residuals_sample <- mapply(FUN=rMVN_var, residuals_terminal_nodes, Omega_matrix_plus_tau, SIMPLIFY = FALSE)
  
  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    residuals_new[terminal_nodes[[i]]$observations_index] <- residuals_sample[[i]]
    # residuals_new[terminal_nodes[[i]]$observations_index] <- residuals_mean[[i]]
    
    
  }
    return(residuals_new)
}

# Prediction function values
get_prediction <- function(trees, x, single_tree = FALSE) {

  # Verying single tree
  if(length(trees) == 1) {
    single_tree <- TRUE
  }

  # Count the number of trees
  n_trees <- length(trees)

  # Stop nesting problems in case of multiple trees
  if(!is.na(pmatch("tree", names(trees))) && (length(trees) == 1)) trees <- trees[[1]]

  # Selecting just one tree
  if(single_tree) {

    # Creating the predictions vector
    predictions <- rep(0, nrow(x))

    # Getting terminal node indexes
    terminal_nodes <- trees[names(which(vapply(trees, "[[", numeric(1), "terminal") == 1))]

    for(i in seq_along(terminal_nodes)) {

      # Retrieving mean values for each obsevartion index
      predictions[terminal_nodes[[i]]$observations_index] <- terminal_nodes[[i]]$mu
      # Round lowest values of predictions
      predictions <- ifelse(abs(predictions) < 1e-15, 0, predictions)
    }
  } else { # Recursive call to multiple trees to calculate the sum of means
    frac_tree <- trees
    frac_tree[[1]] <- NULL # Removing one of trees (the first one)
    predictions <- get_prediction(trees = trees[[1]], x = x, single_tree = TRUE) + # Getting the value of the tree that was blanked
      get_prediction(
        trees = frac_tree, x = x,
        single_tree = (length(frac_tree) == 1)
      ) # Getting the other trees until left just one
  }
    return(predictions)
}

# ==============#
# rBart-GP FUNCTION
#' Fit a GP-BART model
#'
#' @param x A named matrix x of the covariates from the training data.
#' @param y A named matriy y of the response from the training data
#' @param number_trees Number of trees used in the GP-BART model.
#' @param node_min_size Minimum number of observations in a terminal node
#' @param mu Initial values for the \eqn{\mu} parameter
#' @param alpha Parameter \eqn{\alpha} from tree prior
#' @param beta Parameter \eqn{\beta} from tree prior
#' @param tau Initial value for \eqn{\tau} parameter
#' @param n_iter Total number of iterations in GP-BART
#' @param burn Number of burn MCMC iterations
#' @param thin Number of thinning from MCMC
#' @param rotation Boolean to decide if rotated splits will be used
#' @param theta Fixed value for \eqn{\theta}, if NULL a random value for \eqn{\theta} will be used.
#' @param seed Setting a seed for the model
#' @param scale_boolean Boolean to scale or not the response \eqn{y}
#' @param a_tau Scale parameter from \eqn{\tau} prior
#' @param d_tau Rate parameter from \eqn{\tau} prior
#' @param discrete_phi_boolean Boolean to decide if it will be used a discrete grid for \eqn{\phi} proposals
#' @param x_scale Boolean to scale x or not
#' @param nu_vector A constant vector of length of number of trees with the \eqn{\nu} constant value
#' @param gp_variables Covariates used to build the covariance matrix from GP
#' @param K_bart Prior parameter from \eqn{\tau_{\mu}}
#' @param prob_tau Prior parameter from \eqn{\tau}
#' @param kappa Weight parameter \eqn{\kappa} from the proposed mixture of BART and GP-BART
#' @param bart_boolean Boolean to decide or not to warmup the model with BART samples
#' @param bart_number_iter Number of BART initial iterations to be used
#'
#' @return A 'gpbart_GPBART' model object,
#' @export
#'
gp_bart <- function(x, y,
                        number_trees = 2, # Setting the number of trees
                        node_min_size = 15, # Min node size,
                        mu = 0,
                        alpha = 0.5, # Alpha from prior
                        beta = 5, # Beta from prior
                        tau = 1, # Tau from prior,
                        n_iter = 2000, # Number of iterations
                        burn = 500, # Number of burn
                        thin = 1, # Number of thin
                        rotation = TRUE, # If rotated lon and lat will be used in tree building
                        theta = NULL, # If theta is NULL, then the rotation angle will be randomly selected
                        seed = NULL, # Alpha vector values from the Dirichlet prior
                        scale_boolean = TRUE,
                        # This will be defining the nu the default value
                        nu_vector = NULL,
                        a_tau = 3, # Prior from a_v_ratio gamma
                        d_tau = 1, # Prior from d_v_ratio gamma,
                        discrete_phi_boolean = FALSE,
                        x_scale =  TRUE,
                        gp_variables = colnames(x),   # Selecting the GP-Variables
                        K_bart = 2,
                        prob_tau = 0.9,
                        kappa = 0.5, bart_boolean = TRUE, bart_number_iter = 250) {

  # Changing the node_min_size
  if(node_min_size>=nrow(x)){
    stop("Node min size is greater than the number of observation")
  }

  # This parameter is a "scale paramter" to the GP
  phi_vector = rep(0.1 / (sqrt(number_trees)), number_trees)


  # Creating the prediction elements to be stored
  predictions = matrix(0, nrow = number_trees, ncol = nrow(x))
  predictions_list =  NULL

  # If there's only one covariate
  rotation <- ncol(x) != 1

  mean_x <- NULL
  sd_x <- NULL
  
  # Scaling the x
  if(x_scale){
    x_original <- x
    # Scaled version
    xscale <- scale(x)
    mean_x <- attr(xscale,"scaled:center")
    sd_x <- attr(xscale,"scaled:scale")
    x <- as.matrix(xscale)
  } else{ 
    x_original <- x
  }

  # Adjusting the kappa (avoiding the Infinity error)
  if(kappa == 1 ){
    kappa <- kappa - 2*.Machine$double.eps
  }

  if(kappa == 0 ){
    kappa <- kappa + 2*.Machine$double.eps
  }

  # Getting the maximum and minimum values from a distance matrix
  distance_matrix_x <- symm_distance_matrix(m1 = x[,gp_variables, drop = FALSE])
  distance_range <- range(distance_matrix_x[upper.tri(distance_matrix_x)])
  distance_min <- sqrt(distance_range[1])
  distance_max <- sqrt(distance_range[2])

  # Setting seed
  set.seed(seed)
  acc_ratio <- 0
  acc_ratio_phi <- 0
  acc_ratio_nu <- 0

  # Saving a_min and b_max
  a_min <- NULL
  b_max <- NULL

  a_min <- min(y)
  b_max <- max(y)
  
  # Scale values
  if(scale_boolean) {
    # Normalizing y
    y_scale <- normalize_bart(y = y)

    # Defing the nu vector if not in default
    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4 * number_trees * K_bart^2)/(1 - kappa), number_trees)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees * K_bart^2)
    tau_mu_gpbart <- tau_mu_bart/kappa

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)

  } else {

    # Not scaling the y
    y_scale <- y

    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4 * number_trees * K_bart^2)/((1 - kappa) * (max(y_scale) - min(y_scale))^2), number_trees)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees * K_bart^2)/((max(y_scale) - min(y_scale))^2)
    tau_mu_gpbart <- tau_mu_bart/kappa

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
  }

   
  # Creating the likelihood object list
  likelihood_object <- list()
  
  # Getting the number of observations
  n <- length(y)

  # Creating the stor of accepted tree verbs and which split it was
  verb_store_list <- list()

  # Getting the current trees
  current_trees <- list()
  current_trees_proposal <- list()

  # Fixed trees with true parameters
  fixed_trees <- list()

  # Creating the current_partial_residuals
  current_partial_residuals_matrix <-
  current_predictions_matrix <- matrix(NA, nrow = number_trees, ncol = nrow(x))

  # Creating the predictions saving
  current_partial_residuals_list <- list()
  current_predictions_list <- list()

  # Error of the matrix
  if(is.null(colnames(x))) {
    stop("Insert a valid NAMED matrix")
  }

  if(node_min_size == 0) {
    stop("Node Minimum Size need to be greater than 0")
  }

  if(length(nu_vector) != number_trees) {
    stop("Insert a valid \\nu vector for the number of trees")
  }

  if(length(phi_vector) != number_trees) {
    stop("Insert a valid \\phi vector for the number of trees")
  }

  # Recommendation about the min_node_size
  if(node_min_size < 15) {
    warning("\n It is recommended that the node_min_size should be of at least 15 observations.", immediate.=TRUE)
  }

  # Storage containers
  store_size <- (n_iter - burn) 
  tree_store <- vector("list", store_size)
  tau_store <- c()
  signal_pc_store <- matrix(NA, ncol = number_trees, nrow = store_size)

  y_hat_store <-
  y_hat_store_proposal <- matrix(NA, ncol = length(y), nrow = store_size)

  # Storing the likelihoods
  log_lik_store <-
  log_lik_store_fixed_tree <- rep(NA,store_size)

  # Getting the numbers
  loglike_fixed_tree_residuals <- numeric()
  loglike_tree_residuals <- numeric()
  loglike_fixed_tree_residuals_matrix <-
  loglike_tree_residuals_matrix <- matrix(NA, nrow = store_size, ncol = number_trees)

  full_cond_store <-
  phi_store <-
  phi_proposal_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  phi_vector_proposal <- rep(0.1, number_trees)

  # Creating the list of trees stumps
  for(i in seq_len(number_trees)) {
    # Creating the fixed two split trees
    current_trees[[i]] <- stump(x = x, tau = tau, mu = mu)
    current_trees_proposal[[i]] <- stump(x = x, tau = tau, mu =  mu)
  }

  names(current_trees) <-
  names(current_trees_proposal) <- vapply(seq_len(number_trees), function(x) paste0("tree_", x), character(1)) # Naming each tree

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )

  # Getting the parameters to unormalize the data
  a <- min(y)
  b <- max(y)

  # Setting initial values for phi vector
  for(i in seq_len(n_iter)) {

    utils::setTxtProgressBar(progress_bar, i)

    # Changing the bart boolean, when reach the maximum
    if(i >= bart_number_iter){
      bart_boolean <- FALSE
      tau_mu <- tau_mu_gpbart
    } else tau_mu <- tau_mu_bart

    if((i > burn) && ((i %% thin) == 0)) {

      # Saving the store of the other ones
      curr <- (i - burn) / thin
      tree_store[[curr]] <- lapply(current_trees, remove_omega_plus_I_inv)
      tau_store[curr] <- tau

      y_hat_store[curr, ] <- if(scale_boolean){
        
        # Getting the unnormalized version from tau
        unnormalize_bart(colSums(predictions),a = a_min,b = b_max)
      } else {
        
        colSums(predictions)
      }
      # Saving the current partial
      current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix

      # Saving the predictions
      current_predictions_list[[curr]] <- if(scale_boolean){
        unnormalize_bart(predictions, a = a_min, b = b_max)
      } else {
        predictions
      }

      phi_store[curr, ] <- phi_vector
      phi_proposal_store[curr, ] <- phi_vector_proposal
      verb_store_list[[curr]] <- verb_store
    }

    # Creating a boolean to create the first trees only using BART model
    if(bart_boolean) {

      # Verb Store
      verb_store <- data.frame(verb = rep(NA, number_trees),
                               accepted = rep(NA, number_trees),
                               identical = rep(NA, number_trees))

      for(j in seq_len(number_trees)) {

        # Getting the verb list

        # Calculating the residuals for each tree
        if(number_trees > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change", "swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- sample(c("grow", "grow_projection"),
                                                                                               size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate

        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x = x,
          gp_variables = gp_variables,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )
        
        # Checking if the update tree generated a valid tree, if not skip likelihood calculations
        if( !identical(current_trees[[j]],new_trees[[j]]) ){


          # Calculating the likelihood of the new tree
          likelihood_new <- tree_complete_conditional_bart(
            tree = new_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x = x, tau_mu = tau_mu, tau = tau
          )
  
          # Calculating the likelihood of the old tree
          likelihood_old <- tree_complete_conditional_bart(
            tree = current_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x = x, tau_mu = tau_mu, tau = tau
          )
  
          # Extracting only the likelihood
          l_new <- likelihood_new +
            tree_prior(
              tree = new_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )
  
          # Extracting only the likelihood
          l_old <- likelihood_old +
            tree_prior(
              tree = current_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )
  
          # Getting the log of transitin prob
          log_transition <- log_transition_prob(current_tree = current_trees[[j]],
                                                new_tree = new_trees[[j]],verb = verb)
          
          # (log) Probability of accept the new proposed tree
          acceptance <- (l_new - l_old + log_transition)
  
          # If Storage or not based on thin and burn parameters
          if((i > burn) && ((i %% thin) == 0)) {
            full_cond_store[curr, j] <- l_old
          }
          
        } else { # Accepting exact same trees
          # Calculating the likelihood of the new and old tree
          likelihood_new <- likelihood_old <- tree_complete_conditional_bart(
            tree = new_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x = x, tau_mu = tau_mu, tau = tau
          )
          acceptance <- 0
        }
        # In case of acceptance

        if(acceptance > 0 || acceptance > -stats::rexp(1)) {

          # Counting acc ratio
          acc_ratio <- acc_ratio + 1

          # Make changes if accept
          current_trees <- new_trees

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- TRUE


          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_new

        } else {

          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_old

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- FALSE

        } # End of accept if statement

        # To update the mu values
        current_trees[[j]] <- update_mu_bart(
          tree = current_trees[[j]],
          x = x,
          residuals = current_partial_residuals,
          tau = tau,
          tau_mu = tau_mu)

        # Updating the BART predictions
        predictions[j, ] <- update_predictions_bart(
          tree = current_trees[[j]], x = x
        )

        # Updating the current residuals
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]

      } # End of Loop through the trees

      # =================
      # ATTENTION HERE!!!
      # =================
    } else { # Going over the case where the BART-boolean is no more valid

      # Verb Store
      verb_store <- data.frame(verb = rep(NA,number_trees),
                               accepted = rep(NA,number_trees),
                               identical = rep(NA,number_trees))

      for(j in seq_len(number_trees)) {

        # Calculating the residuals for each tree
        if(number_trees > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change","swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- sample(c("grow","grow_projection"),
                                                                                              size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate

        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x = x,
          gp_variables = gp_variables,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )
        
        
        # ==================== #
        # Getting the Omega Inverse the current and the future tree
        # ==================== #

        # Checking if the previous tree were identitcal or not
        if(!identical(current_trees[[j]],new_trees[[j]])){
        
            # Getting the inverse for the current terminal nodes
            current_trees[[j]] <- inverse_omega_plus_I(tree = current_trees[[j]],
                                                       x = x, tau = tau,
                                                       nu = nu_vector[j],
                                                       phi = phi_vector[j])
    
            # Getting the inverse for the new tree terminal nodes
            new_trees[[j]] <- inverse_omega_plus_I(tree = new_trees[[j]],
                                                   x = x,tau = tau,
                                                   nu = nu_vector[j],
                                                   phi = phi_vector[j])
    
            # Calculating the likelihood of the new tree
            likelihood_new <- tree_complete_conditional_gpbart(
              tree = new_trees[[j]],  # Calculate the full conditional
              residuals = current_partial_residuals,
              x = x, tau_mu = tau_mu,
              nu = nu_vector[j], phi = phi_vector[j],
              number_trees = number_trees
            )
    
            # Calculating the likelihood of the old tree
            likelihood_old <- tree_complete_conditional_gpbart(
              tree = current_trees[[j]], # Calculate the full conditional
              residuals = current_partial_residuals,
              x = x, tau_mu = tau_mu,
              nu = nu_vector[j], phi = phi_vector[j],
              number_trees = number_trees
            )
    
            # Extracting only the likelihood
            l_new <- likelihood_new$log_posterior +
              tree_prior(
                tree = new_trees[[j]], # Calculate the tree prior
                alpha = alpha,
                beta = beta
              )
    
            # Extracting only the likelihood
            l_old <- likelihood_old$log_posterior +
              tree_prior(
                tree = current_trees[[j]], # Calculate the tree prior
                alpha = alpha,
                beta = beta
              )
    
            # Getting the log of transitin prob
            log_transition <- log_transition_prob(current_tree = current_trees[[j]],
                                                  new_tree = new_trees[[j]],verb = verb)
            
            # (log) Probability of accept the new proposed tree
            acceptance <- (l_new - l_old + log_transition)
    
            # If Storage or not based on thin and burn parameters
            if((i > burn) && ((i %% thin) == 0)) {
              full_cond_store[curr, j] <- l_old
            }
            # In case of them being identical
        } else {
          
          acceptance <- 0.1
          # Checking if none Omega matrix element were calculated
          if(is.null(current_trees[[j]][unlist(lapply(current_trees[[j]],
                                                      function(node) node$terminal==1))][[1]]$Omega_matrix)){
              # Creating the current tree object
              # Getting the inverse for the current terminal nodes
              current_trees[[j]] <- new_trees[[j]] <- inverse_omega_plus_I(tree = current_trees[[j]],
                                                         x = x, tau = tau,
                                                         nu = nu_vector[j],
                                                         phi = phi_vector[j])
              
              # Creating the likelihood object
              likelihood_new <- likelihood_old <- 
                tree_complete_conditional_gpbart(
                tree = current_trees[[j]], # Calculate the full conditional
                residuals = current_partial_residuals,
                x = x, tau_mu = tau_mu,
                nu = nu_vector[j], phi = phi_vector[j],
                number_trees = number_trees
              )
            } else{ # Replacing the likelihood_object
              likelihood_new <- likelihood_old <- likelihood_object[[j]]
            } # end of the if for null likelihood and current tree 
          
          }

        if(acceptance > 0 || acceptance > -stats::rexp(1)) { #
          acc_ratio <- acc_ratio + 1

          # Checking whether the trees are identical
          if(identical(current_trees[[j]], new_trees[[j]])){
            verb_store[j,"identical"] <- TRUE
          } else {
            verb_store[j,"identical"] <- FALSE
          }

          # Make changes if accept
          current_trees <- new_trees

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- TRUE

          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_new

        } else {
          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_old

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- FALSE
          verb_store[j,"identical"] <- FALSE

        } # End of accept if statement

        # # # To update the mu values
        current_trees[[j]] <- update_mu(
          tree = current_trees[[j]],
          x = x,
          residuals = current_partial_residuals,
          likelihood_object = likelihood_object[[j]])

        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j, ] <- update_residuals(
          tree = current_trees[[j]], x = x,
          residuals = current_partial_residuals,
          phi = phi_vector[j], nu = nu_vector[j], tau = tau
        )

        # To update phi
        mh_update_phi <- update_phi_marginal(current_tree_iter = current_trees[[j]],
                                               residuals = current_partial_residuals,
                                               x = x,nu = nu_vector[j],phi = phi_vector[j],
                                               gp_variables = gp_variables,
                                               likelihood_object = likelihood_object[[j]],
                                               number_trees = number_trees,
                                               discrete_phi = discrete_phi_boolean,
                                               tau = tau,
                                               tau_mu = tau_mu,
                                               distance_min = distance_min,
                                               distance_max = distance_max)

          # In case of accept the update over \phi update everything
          if(mh_update_phi$phi_boolean) {

            # Updating the tree and the \phi object from the tree
            current_trees[[j]] <- mh_update_phi$tree

            # Updating the likelihood objects
            likelihood_object[[j]] <- mh_update_phi$likelihood_object

            # Updating the phi value
            phi_vector[j] <- mh_update_phi$phi_proposal

          } # If doesn't accept, nothing changes.



        # current_partial_residuals_matrix<-
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]

      } # End of Loop through the trees
    }

    tau <- update_tau(x = x,
                      y = y_scale,
                      a_tau = a_tau,
                      d_tau = d_tau,
                      predictions = colSums(predictions))

  } # End of Loop through the n_inter
  cat("\n")

  # Returning X to its original scale
  if(x_scale) {
    x <- x_original
  }

  results <- list(trees = tree_store,
                  tau_store = tau_store,
                  y_hat = y_hat_store,
                  log_lik = log_lik_store,
                  log_lik_fixed_tree = log_lik_store_fixed_tree,
                  loglike_fixed_tree_residuals_matrix = loglike_fixed_tree_residuals_matrix,
                  loglike_tree_residuals_matrix = loglike_tree_residuals_matrix,
                  full_cond = full_cond_store,
                  phi_store = phi_store,
                  phi_proposal_store = phi_proposal_store,
                  nu_vector = nu_vector,
                  y = y_scale,
                  X = x_original,
                  x_scale = x_scale,
                  mean_x = mean_x,
                  sd_x = sd_x,
                  scale_boolean = scale_boolean,
                  acc_ratio = acc_ratio,
                  acc_ratio_phi = acc_ratio_phi,
                  iter = n_iter,
                  burn = burn,
                  thin = thin,
                  store_size = store_size,
                  number_trees = number_trees,
                  node_min_size = node_min_size,
                  a_min = a_min,
                  b_max = b_max,
                  a_tau = a_tau,
                  d_tau = d_tau,
                  current_partial_residuals_list = current_partial_residuals_list,
                  beta = beta,
                  current_predictions_list = current_predictions_list,
                  tau_mu = tau_mu, kappa = kappa,
                  verb_store_list = verb_store_list)
  class(results) <- "gpbart_GPBART"
    return(results)
}

# #Do a MH for PHI
update_phi_marginal <- function(x, current_tree_iter,residuals,
                                seed = NULL,
                                tau,
                                tau_mu,
                                phi, nu,number_trees,
                                likelihood_object, p, gp_variables,
                                discrete_phi = TRUE,
                                distance_min,
                                distance_max) {

  # Increased the range of tree proposal
  if(discrete_phi){
    phi_proposal <- sample(c(0.1,0.5,1,5,10), size = 1)
  } else {
    phi_proposal <- stats::runif(1, min = distance_min, max = distance_max)
  }

  # Calculating the likelihood from the new step
  tree_from_phi_proposal <- inverse_omega_plus_I(tree = current_tree_iter,
                                                 x = x,nu = nu, tau = tau,
                                                 phi = phi_proposal, gp_variables = gp_variables)

  likelihood_phi_proposal <- tree_complete_conditional_gpbart(tree = tree_from_phi_proposal,
                                                              x = x,
                                                              residuals = residuals,
                                                              nu = nu, tau_mu = tau_mu,
                                                              phi = phi_proposal,
                                                              number_trees = number_trees)

  # Old phi likelhood
  l_old_phi <- likelihood_object$log_posterior

  # Proposal likelihood
  l_proposal_phi <- likelihood_phi_proposal$log_posterior

  # (log) Probability of accept the new proposed tree
  acceptance_phi <- l_proposal_phi - l_old_phi

  # If storage for phi

  if(acceptance_phi > 0 || acceptance_phi > -stats::rexp(1)) { #

    # Nu boolean to see if was accepted or not
    phi_boolean <- TRUE
    return(list(phi_boolean = phi_boolean,
                likelihood_object = likelihood_phi_proposal,
                tree = tree_from_phi_proposal,
                phi_proposal = phi_proposal)) # Returning the proposal value for phi
  } else {
    # Case of not accepting
    phi_boolean <- FALSE
    return(list(phi_boolean = phi_boolean)) # Returning the old value for phi
  } #
}


# Function to return the depth trees

tree_depth_hist <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- max(vapply(tree, "[[", numeric(1), "depth_node"))
    }
  }
    return(tree_depth)
}

# Get the covariate splits
tree_var_hist <- function(gpbart_model) {
  tree_depth <- c()

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for(i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth <- c(tree_depth, sapply(tree, "[[", "node_var"))
    }
  }
    return(tree_depth)
}

# Function to count th enumber of terminal nodes in a tree
tree_count_terminals <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for(i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- sum(vapply(tree, "[[", numeric(1), "terminal"))
    }
  }
    return(tree_depth)
}

# Tau values
get_tau_values <- function(gpbart_model) {

  # n_iter
  n_iter <- length(gpbart_model$trees)

  tau_iters <- rep(list(NA), gpbart_model$number_trees)

  num_trees <- gpbart_model$number_trees

  for(j in seq_len(n_iter)) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for tau terminals
    tau_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(tau_terminals) <- paste0("node_", 0:49)

    for(k in seq_len(num_trees)) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)

      for(nodes in seq_along(all_nodes)) {
        if(tree[[all_nodes[nodes]]]$terminal == 1) {
          # tree[[nodes]]$tau <- sample(c(0.001,seq(0.005,10,by = 0.005),10,50,100), 1)
          tau_terminals[k, as.numeric(stringr::str_extract(all_nodes[nodes], pattern = "\\(?[0-9,.]+\\)?")) + 1] <- tree[[all_nodes[nodes]]]$tau
        }
      }
    }

    # Saving all
    tau_iters[[j]] <- tau_terminals
  }
    return(tau_iters)
}

# mu values
get_mu_values <- function(gpbart_model) {

  # n_iter
  n_iter <- length(gpbart_model$trees)

  mu_iters <- rep(list(NA), gpbart_model$number_trees)

  num_trees <- gpbart_model$number_trees

  for(j in seq_len(n_iter)) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for mu terminals
    mu_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(mu_terminals) <- paste0("node_", 0:49)

    for(k in seq_len(num_trees)) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)

      for(nodes in seq_len(all_nodes)) {
        if (tree[[all_nodes[nodes]]]$terminal == 1) {
          # tree[[nodes]]$mu <- sample(c(0.001,seq(0.005,10,by = 0.005),10,50,100), 1)
          mu_terminals[k, as.numeric(stringr::str_extract(all_nodes[nodes], pattern = "\\(?[0-9,.]+\\)?")) + 1] <- tree[[all_nodes[nodes]]]$mu
        }
      }
    }

    # Saving all
    mu_iters[[j]] <- mu_terminals
  }
    return(mu_iters)
}

# ORIGINAL PREDICT GAUSSIAN FROM MULTIPLE TREES
predict_gaussian_from_multiple_trees <- function(multiple_trees, # A list of trees
                                                 phi_vector, # A vector of phi values
                                                 nu_vector, # The nu value to be used
                                                 x_train, # The x of the training model
                                                 x_new, # The x that will be predicted
                                                 partial_residuals, # The partial_residual values
                                                 tau,
                                                 pred_bart_only # Boolean argument to predict a BART object
                                                 ) {
  # Defining objects
  y_pred_final <-

  # Calculating the sd from the prediction interval
  y_pred_pi_final <- matrix(0, nrow = length(multiple_trees), ncol = nrow(x_new))

  # # Covariance matrix
  # cov_pred_final <- list()
  # if(isTRUE(get_cov_star)) { variance <- matrix(0, nrow = nrow(x_new), ncol = nrow(x_new)) }

  # Iterating over all trees
  for(m in seq_along(multiple_trees)) {

    # Creating the list to be predicted (The if is just in case of of just one tree)
    new_tree <- multiple_trees[[m]]
    phi <- phi_vector[m]
    nu <- nu_vector[m]

    # Creating the pred  vector
    y_pred <- numeric(nrow(x_new))
    y_pred_pi <- numeric(nrow(x_new))

    # if(isTRUE(get_cov_star)) { variance <- matrix(0, nrow = nrow(x_new), ncol = nrow(x_new)) }

    # Setting the root node with the new observations
    new_tree[["node_0"]]$test_index <- seq_len(nrow(x_new))

    # Creating the list of nodes
    list_nodes <- names(new_tree)[-1]

    # IN CASE OF JUST ROOT NODE
    if(length(new_tree) == 1) {
      list_nodes <- "node_0"
    }

    # Updating all nodes
    for(i in seq_along(list_nodes)) {

      current_node_aux <- new_tree[[list_nodes[i]]]

      # In case of more than one node
      if(length(list_nodes) > 1) {
        # Veryfing the type of the current node
        if (is.list(current_node_aux$node_var)) {

          # Rotated Lon and Lat

          rotated_x <- tcrossprod((A(current_node_aux$theta)), x_new[,current_node_aux$node_var$node_var_pair,drop = FALSE])
          rownames(rotated_x) <- current_node_aux$node_var$node_var_pair

          # Updating observations from the left node
          if(current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index]  < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] >= current_node_aux$node_var_split)]
          }
        } else { # Evaluating the case where is not a rotated lat/lon

          # To continous covariates
          if(is.numeric(x_new[, current_node_aux$node_var, drop = FALSE])) {

            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var]  < current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
            }

            # To categorical covariates
          } else {
            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
            }
          }
        }
      }
      # Here I will calculate the predicted values based on the terminal nodes AND only on terminal nodes which have test

      if(new_tree[[list_nodes[[i]]]]$terminal == 1 && length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {

        # Selecting the observations from the current node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED
        x_current_node <- matrix(x_train[new_tree[[list_nodes[[i]]]]$observations_index, ],
                                 nrow = length(new_tree[[list_nodes[[i]]]]$observations_index))

        # Selecting the observations from the test node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED

        x_star <- matrix(x_new[new_tree[[list_nodes[[i]]]]$test_index, ],
                         nrow = length(new_tree[[list_nodes[[i]]]]$test_index))

        # Calcualting the distance matrix
        distance_matrix_current_node <- symm_distance_matrix(m1 = x_current_node)

        # Generating the scenario where there are only BART predictions
        if(isFALSE(pred_bart_only)) {

          # Getting the GP from a terminal node
          gp_process <- gp_main_sample(
            x_train = x_current_node, distance_matrix_train = distance_matrix_current_node,
            y_train = matrix((partial_residuals[m,new_tree[[list_nodes[[i]]]]$observations_index]) - new_tree[[list_nodes[[i]]]]$mu,
                              nrow = nrow(x_current_node)),
            x_star = x_star, tau = tau,
            nu = nu, phi = phi,get_cov_star = FALSE
          )

          # Creating the mu vector
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- gp_process$mu_pred + new_tree[[list_nodes[[i]]]]$mu

        } else {

          # Attributing the predicted value as the sampled \mu from the terminal node
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- new_tree[[list_nodes[[i]]]]$mu
        }

        # Creating the sd of the mu vector
        y_pred_pi[new_tree[[list_nodes[i]]]$test_index] <- (1/tau)

      }
    }

    y_pred_pi_final[m,] <- y_pred_pi
    y_pred_final[m, ] <- y_pred
    # if(isTRUE(get_cov_star)) { cov_pred_final[[m]] <- variance }
  }
  # Return the new tree
  return(list(y_pred = colSums(y_pred_final),
              all_tree_pred = y_pred_final
              ))
}

# Function to count the number of terminal nodes
count_terminal_nodes <- function(tree) {
  return(sum(vapply(tree, "[[", numeric(1), "terminal") == 1))
}

# Predicting a gpbart
#' @method predict "gpbart_GPBART"
#' @rdname gpbart_GPBART The fitted gpBART model
#' @param x_test the test set
#' @param pred_bart_only boolean if there are only bart predictions
#' @param type select the prediction outputs among 'c("all", "mean","median"))'
#' @param ... other parameters
#' @usage
#' \method{predict}{gpbart_GPBART}(object,
#'         x_test,
#'         pred_bart_only = FALSE,
#'         type = c('all', 'median', 'mean'),
#'         ...)
#' @export
predict.gpbart_GPBART <- function(object, x_test, 
                                  pred_bart_only = FALSE,
                                  type = c("all","mean","median"),...) { # type argument Can be "all", "mean" or "meadian"

  # Adjusting the type
  type <- match.arg(type)
  
  # Loading x_train
  if(object$x_scale){
    x_train <- as.matrix(scale(object$X,center = object$mean_x,scale = object$sd_x))
    x_test  <- as.matrix(scale(x_test,center = object$mean_x,scale = object$sd_x))

    # Retrieving the col.names
    colnames(x_train) <- colnames(x_test)  <- colnames(object$X)
  } else {
    x_train <- object$X
  }

  # Number of iters of bayesian simulation
  n_iter <- length(object$tau_store)

  # The number of columns is the number of test observations and the rows are the iterations
  y_hat_matrix <- matrix(NA, ncol = nrow(x_test),nrow = n_iter)

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running rBART..."
  )
  y_list_matrix <- list()

  # Creating the final vector
  y_pred_final <- matrix(0, nrow = object$number_trees, ncol = nrow(x_test))


  # Looping around the trees
  for(i in seq_len(n_iter)) {

    utils::setTxtProgressBar(progress_bar, i)

    # Selecting one tree from BART model
    current_tree <- object$trees[[i]]

    # Getting the predictions from the test observations
    y_pred_aux <- predict_gaussian_from_multiple_trees(
      multiple_trees = current_tree, x_train = x_train,
      x_new = x_test, partial_residuals = object$current_partial_residuals_list[[i]],
      phi_vector = object$phi_store[i, ],
      nu_vector = object$nu_vector,
      tau = object$tau_store[[i]], pred_bart_only = pred_bart_only
    )

    # Iterating over all trees (test)
    y_pred_final <- y_pred_aux$all_tree_pred

    if(object$scale_boolean) {
      # Recovering the prediction interval from test
      y_hat_matrix[i, ] <- unnormalize_bart(colSums(y_pred_final), a = object$a_min, b = object$b_max)
    } else {
      y_hat_matrix[i, ] <- colSums(y_pred_aux$all_tree_pred)
    }

    y_list_matrix[[i]] <- y_pred_final

  }
  # Chaging the value of \tau in case of scaling
  object$tau_store <- if(object$scale_boolean){
    (object$tau_store/((object$b_max-object$a_min)^2))
  } else {
    object$tau_store
  }
  
  out <- list(
    pred = switch(type,
                  all = y_hat_matrix,
                  mean = colMeans(y_hat_matrix),
                  median = apply(y_hat_matrix, 2, "median"),
    ),
    sd = switch(type,
                all = object$tau_store^(-1/2),
                mean = mean(object$tau_store^(-1/2)),
                median = stats::median(object$tau_store^(-1/2))
    )
  )

  return(list(out = out, list_matrix_pred = y_list_matrix))
}

# A function to count the mean values of observations in terminal nodes
gpbart_count_terminal_nodes <- function(mod_gpbart){

  # Auxiliar matrix
  all_tree_terminal_nodes <- array(NA,dim = c(mod_gpbart$number_trees, 50, mod_gpbart$store_size),
                                   dimnames = list(paste0("tree_", seq_len(mod_gpbart$number_trees)),
                                                   paste0("node_", seq_len(50)),
                                                   paste0("iter_", seq_len(mod_gpbart$store_size))))

  # Iterating over MH
  for(k in seq_len(mod_gpbart$store_size)) {
    tree_iter <- mod_gpbart$trees[[k]]

    # Iterating with for over the terminal nodes
    for(i in seq_len(mod_gpbart$number_trees)) {

      # Gathering the node number
      node_number <- unlist(lapply(tree_iter[[i]], function(x) { x[x$terminal == 1]$node_number}))
      all_tree_terminal_nodes[i,node_number,k] <- unlist(lapply(tree_iter[[i]][names(node_number)],function(x) {
        length(x[x$terminal == 1]$observations_index)
      }))

    }
  }
  # Three splits
  split_nodes <- unique(colMeans(all_tree_terminal_nodes, na.rm = TRUE))
    return(split_nodes[!is.na(split_nodes)])
}



# Retrieve new points on the cubic scale
unit_cube_scale_new_points_uni <- function(x, x_new) {

  # Scaling to -1 to 1 function
  scaling <-
    (2 * x_new - (max(x) + min(x))) / (max(x) - min(x))

  # Applying on all covariates
    return(scaling)
}

# Getting the likelihood from the Gaussian Processes
neg_loglike <- function(prediction_object,
                        y_test) {

  # Getting the prediction
  y_pred <- colMeans(prediction_object$mcmc_pi_mean)

  # Getting the variance
  sd_pred <- colMeans(prediction_object$mcmc_pi_sd)

  return(-stats::dnorm(x = y_test,
                mean = y_pred,
                sd = sqrt(sd_pred), log = TRUE))
}

# sum(neg_loglike(prediction_object = pred_gpbart,y_test = y_test))

# Get tau from terminal nodes
get_tau_values_from_single_tree <- function(gpbart_mod, tree_number, mh_iter = 100){
  # Getting the vector of tau values from the terminal nodes
  return(unlist(lapply(gpbart_mod$trees[[mh_iter]][[tree_number]],function(x) x[x$terminal == 1]$tau)))
}

# Getting Omega Inverse argument list
# get_omega_inverse_terminal_nodes <- function(tree,
#                                              x = x,
#                                              nu = nu_vector[j], phi = phi_vector[j]){
#
#   # Selecting terminal nodes names
#   names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
#
#   # Selecting the terminal nodes
#   terminal_nodes <- tree[names_terminal_nodes]
#
#   # Number of nodes
#   n_node <- length(terminal_nodes)
#
#   # Picking each node size
#   nodes_size <- sapply(terminal_nodes, function(x) {
#     length(x$observations_index)
#   })
#
#   # Calculating Omega matrix INVERSE
#   Omega_matrix_INV <- mapply(terminal_nodes, FUN = function(y) {
#     chol2inv(  chol(kernel_function(
#       x = matrix(x[y$observations_index, ], nrow = length(y$observations_index)),
#       nu = nu, phi = phi
#     )) )
#   }, SIMPLIFY = FALSE)
#
#   # Adding the mu values calculated
#   for(i in seq_along(names_terminal_nodes)) {
#     tree[[names_terminal_nodes[i]]]$Omega_inv <- Omega_matrix_INV[[names_terminal_nodes[i]]]
#   }
#     return(tree)
# }

# Getting Omega Inverse + Diag Inverse
inverse_omega_plus_I <- function(tree,
                                 x = x,
                                 nu, phi,
                                 tau,
                                 number_trees = number_trees,
                                 gp_variables = colnames(x)  # Selecting which gp-variables to use

) {
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })

  # Calculating Omega matrix INVERSE
  distance_matrices <- mapply(terminal_nodes, FUN = function(y) {
    symm_distance_matrix(matrix(x[y$observations_index, gp_variables], nrow = length(y$observations_index)))
  }, SIMPLIFY = FALSE)

  # Calculating Omega
  Omega_matrix <- mapply(distance_matrices, FUN = function(dist_m) {
    kernel_function(
    squared_distance_matrix = dist_m,
    nu = nu, phi = phi)
  }, SIMPLIFY = FALSE)

  # Checking if diagonal
  is_Omega_diag <- lapply(Omega_matrix, is_diag_matrix)

  # Calculating Omega_plus_I*tau^(-1)
  Omega_matrix_plus_I <- mapply(Omega_matrix, FUN = function(omega) {
    omega + diag(1/(tau), nrow = nrow(omega))
  }, SIMPLIFY = FALSE)

  # Calculating Omega matrix plus I INVERSE
  Omega_matrix_plus_I_INV <- mapply(Omega_matrix_plus_I, FUN = function(omega_plus_I_tau) { # p is the shrinkage factor
    chol2inv(PD_chol(omega_plus_I_tau))
  }, SIMPLIFY = FALSE)

  # Adding the Omega_matrix_plus_I_Inv
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$Omega_plus_I_tau <- Omega_matrix_plus_I[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$Omega_plus_I_inv <- Omega_matrix_plus_I_INV[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$distance_matrix <- distance_matrices[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$Omega_matrix <- Omega_matrix[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$is_Omega_diag <- is_Omega_diag[[names_terminal_nodes[i]]]
  }
    return(tree)
}

# # Removing the Omega_plus_I_inv object
remove_omega_plus_I_inv <- function(current_tree_iter) {

  # Selecting terminal nodes names
  # names_terminal_nodes <- names(which(vapply(current_tree_iter, "[[", numeric(1), "terminal") == 1))
  names_terminal_nodes <- names(current_tree_iter)
  
  for(i in names_terminal_nodes) {
    current_tree_iter[[i]]$Omega_plus_I_inv <-
    current_tree_iter[[i]]$distance_matrix <-
    current_tree_iter[[i]]$Omega_plus_I_tau <-
    current_tree_iter[[i]]$Omega_matrix <-
    current_tree_iter[[i]]$is_Omega_diag <- NULL
  }
    return(current_tree_iter)
}


# Getting the average values from the training predictions
gpbart_train_mean <- function(gpbart_mod) {
  colSums(Reduce("+", gpbart_mod$current_predictions_list)/length(gpbart_mod$current_predictions_list))
}

# Calculating the variance from the training set.
gpbart_training_var <-  function(gpbart_mod) {

  # Number of MCMC samples
  n_mcmc <- length(gpbart_mod$trees)
  # Creating the var vector
  var_train <- matrix(0, nrow = n_mcmc, ncol = length(gpbart_mod$y))

  # Iterating over all trees and getting the terminal nodes
  for(i in seq_len(n_mcmc)) {

    for(m in seq_len(gpbart_mod$number_trees)) {
      # Selecting the current tree
      tree <- gpbart_mod$trees[[i]][[m]]

      # Selecting the terminal nodes
      terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

      # Iterating over the terminal nodes
      for(node in terminal_nodes){

        # In this line I adding up the quantity of 1/tau for each tree
        var_train[i,node$observations_index] <- var_train[i,node$observations_index] + 1/node$tau
      }
    }
  }

  # Returning the var_train vector
    return(var_train)
}


# Calculating the residuals log-likelihood
loglike_residuals <- function(tree,
                              x,
                              current_partial_residuals,
                              phi,
                              nu,p){

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    current_partial_residuals[x$observations_index]
  })

  # Getting the \tau values
  tau_terminal <- lapply(terminal_nodes, "[[", "tau")

  # Getting the \mu values
  mu_terminal <- lapply(terminal_nodes, "[[", "mu")

  Omega_matrix_plus_I <- lapply(terminal_nodes, function(nodes){
    kernel_function(
      squared_distance_matrix = nodes$distance_matrix,
      nu = nu, phi = phi) + diag(p, nrow = nrow(nodes$distance_matrix))
  })

  # Creating the loglikeresiduals_vec
  loglike_residuals_vec <- numeric()
  # Doing a quick for
  for(i in seq_along(residuals_terminal_nodes)){
    loglike_residuals_vec[i] <-  mvtnorm::dmvnorm(x = residuals_terminal_nodes[[i]],
                                                  mean = rep(mu_terminal[[i]], length(residuals_terminal_nodes[[i]])),
                                                  sigma = (1/tau_terminal[[i]]) * Omega_matrix_plus_I[[i]], log = TRUE)
  }
  # Returning the sum of the log_like_residual_vec
    return(sum(loglike_residuals_vec))
}


# NEW UPDATE G
# Update tau_j values
update_g <- function(tree, x, nu, phi, residuals, seed = NULL, p) {
  # set.seed(seed)

  # New g (new vector prediction for g)
  g_new <- rep(NA, length(residuals))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal")))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")

  # Selecting the tau values
  tau_j <- vapply(terminal_nodes, "[[", numeric(1), "tau")

  # Calculating Omega matrix
  Omega_matrix_inverse <- mapply(terminal_nodes, FUN = function(nodes) {
    chol2inv(PD_chol(
      kernel_function(squared_distance_matrix = nodes$distance_matrix,
                      nu = nu,
                      phi = phi)
    ))
  }, SIMPLIFY = FALSE)

  # Getting the A matrix inverse
  A_matrix_inv <- mapply(Omega_matrix_inverse, FUN = function(x) {
    chol2inv(PD_chol(diag(p, nrow = nrow(x)) + x))
  }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  g_mean <- mapply(A_matrix_inv,
                   residuals_terminal_nodes,
                   mu_values,Omega_matrix_inverse, FUN = function(A_inv, res, mu, omg) {
                     crossprod(A_inv, p * res + mu * rowSums(omg))
                   }, SIMPLIFY = FALSE)

  # Putting in the Keefe's speed order
  g_sd <- mapply(tau_j, Omega_matrix_inverse, FUN = function(tau, omg) {
    tau * (omg + diag(p, nrow=nrow(omg)))
  }, SIMPLIFY = FALSE)

  g_sample <- mapply(FUN=rMVN_var, g_mean, g_sd, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    g_new[terminal_nodes[[i]]$observations_index] <- g_sample[[i]]
  }
    return(g_new)
}

# Some tests over nu parameter
calculate_nu <- function(nu) {
  (nu + 1)/nu
}

calculate_p_nu <- function(nu, p) {
  (p * nu + 1)/nu
}

# Get train predictions
get_train_predictions <- function(gpbart_mod) {

  # Getting the quantile
  gpbart_sum_pred <- do.call(rbind, lapply(gpbart_mod$current_predictions_list, colSums))

  # Returning the matrix of final predictions
    return(gpbart_sum_pred)
}

# Calculating a PI coverage
#' @export
pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){
  
  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)
  
  full_post_draw <- list()
  
  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )
  
  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)
      
      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }
  
  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }
  
  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})
  
  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)
  
  return(pi_cov)
}
