# GP-Bart

# Loading test data
library(magrittr)
library(Rcpp)

# Load in functions which are common to BART
source("common_help_functions.R")

# ==================================#
# Objects to test the tree_complete_conditional function
# ==================================#

tree_complete_conditional_gpbart <- function(tree, x, residuals, nu = 1, phi = 1, 
                                             tau_mu,
                                             number_trees= number_trees,
                                             tau_multiplier) {

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  
  # Calculating value S
  S <- unlist(mapply( terminal_nodes, FUN=function(z) {
    (sum(z$Omega_plus_I_inv) + c(tau_mu ))
  },SIMPLIFY = TRUE))
  
  
  # Log Omega
  log_det_Omega_plus_I_tau <- unlist(lapply(terminal_nodes, function(z) {
    determinant(z$Omega_plus_I_tau, logarithm = TRUE)$modulus
  }))
  
  log_det_Omega_plus_I_tau
  # Defining RT_Omega_I_R
  RTR <- unlist(mapply(terminal_nodes,residuals_terminal_nodes, FUN = function(nodes,resid) {
    crossprod(resid,crossprod(nodes$Omega_plus_I_inv,resid))
  },SIMPLIFY = FALSE))
  
  
  # The term R^{T} solve(Omega + I ) 1 
  R_Omega_I_one <- unlist(mapply(terminal_nodes, residuals_terminal_nodes, FUN = function(nodes, residuals) {
    rowSums(crossprod(residuals, nodes$Omega_plus_I_inv))
  },SIMPLIFY = FALSE))
  
  
  log_posterior <- -0.5*sum(log_det_Omega_plus_I_tau) - 0.5 * sum(log(S/tau_mu)) - 0.5*sum(RTR) + 0.5 * sum(S^(-1) *R_Omega_I_one^2)  
  
  
  return(list(log_posterior = log_posterior,
              S = S,
              RTR = RTR,
              R_Omega_I_one =  R_Omega_I_one))
}


# Generate mu_j values
update_mu <- function(tree,
                      x,
                      likelihood_object,
                      seed = NULL) {
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })
  
  
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  
  # Mu mean value
  mu_mean <- likelihood_object$R_Omega_I_one*likelihood_object$S^(-1)
  
  
  # Mu mean SD
  mu_var <- (likelihood_object$S)^-1
  
  # cat("Mu mean equals to: ",unlist(mu_mean),"\n")
  # cat("Mu sd equals to: ",unlist(mu_var),"\n")
  
  # Calculating mu values
  mu <- mapply(mu_mean, mu_var,
               FUN = function(x, y) {
                 rnorm(
                   n = 1,
                   mean = x,
                   sd = sqrt(y)
                 )
               }, SIMPLIFY = FALSE
  )
  
  # Adding the mu values calculated
  for (i in 1:length(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- c(mu[[names_terminal_nodes[i]]])
  }
  
  return(tree)
}

update_residuals <- function(tree, x, nu, phi, residuals, tau, error_handling_residuals, seed = NULL) {
  # set.seed(seed)
  
  # New g (new vector prediction for g)
  residuals_new <- rep(NA, length(residuals))
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  
  
  # Getting the \mu_{j} vector
  mu_values <- unlist(
    lapply(terminal_nodes, function(x) {
      x$mu
    })
  )
  
  
  # Calculating Omega matrix plus I INVERSE
  Omega_matrix <- lapply(terminal_nodes, FUN = function(nodes) {
    kernel_function(
      squared_distance_matrix =  nodes$distance_matrix,
      nu = nu, phi = phi)
  })
  
  # Calculating Omega matrix
  Omega_matrix_inverse_plus_I  <- lapply(terminal_nodes, FUN = function(nodes) {
    nodes$Omega_plus_I_inv
  })
  
  # print(diag((Omega_matrix_inverse_plus_I[[1]]%*%Omega_matrix[[1]]))[1:5])
  
  # Calculating g_mean posterior
  residuals_mean <- mapply(Omega_matrix,
                           Omega_matrix_inverse_plus_I,
                           residuals_terminal_nodes,
                           mu_values,
                           FUN = function(omega, omega_i_inv, residuals, mu) {
                             # Getting the values
                             rep(mu,dim(omega)[1]) + crossprod(omega,
                                                               crossprod(omega_i_inv,
                                                                         (residuals - rep(mu,dim(omega)[1]))
                                                               )
                             )
                           }, SIMPLIFY = FALSE)
  
  
  # Calculating g_mean posterior
  residuals_variance <- mapply(Omega_matrix,
                               Omega_matrix_inverse_plus_I,
                               residuals_terminal_nodes,
                               mu_values,
                               FUN = function(omega, omega_i_inv, residuals, mu) {
                                 # Getting the Omega value
                                 (omega - crossprod(omega,crossprod(omega_i_inv,omega)))  #+ diag(1/tau,nrow = dim(omega)[1])
                               }, SIMPLIFY = FALSE)
  
  
  # # Calculating g_mean posterior
  residuals_sample <- mapply(residuals_mean,residuals_variance,
                             FUN =  function(mean,var){
                               
                               if(error_handling_residuals){
                                 
                                 # Using tryCatch to handle the error
                                 r_sample <- try(rMVN2(b = mean,Q = var),silent = TRUE)
                                 
                                 # Veryfing if it's an error
                                 if(class(r_sample)=="try-error"){
                                   error_counter_rmvn <<- error_counter_rmvn + 1 # Global variable so ugly!!!
                                   return(MASS::mvrnorm(n = 1,mu = mean,Sigma = var))
                                 }else{
                                   return(r_sample)
                                 }
                               } else{
                                 # Run without error handling
                                 return(MASS::mvrnorm(n = 1,mu = mean,Sigma = var))
                                 # rMVN2(b = mean,Q = (var+diag(1e-8, nrow = size(var)[1])))
                               }
                             }, SIMPLIFY = FALSE)
  
  
  # Adding the mu values calculated
  for (i in 1:length(terminal_nodes)) {
    # Saving g
    residuals_new[terminal_nodes[[i]]$observations_index] <- residuals_sample[[i]]
    
  }
  
  # cat("Residuals_sample",residuals_sample[[1]][1:5],"\n\n")
  # cat("Residuals_mean",residuals_mean[[1]][1:5],"\n\n")
  
  return(residuals_new)
}

# Prediction function values
get_prediction <- function(trees, x, single_tree = FALSE) {
  
  # Verying single tree
  if (length(trees) == 1) {
    single_tree <- TRUE
  }
  
  # Count the number of trees
  n_trees <- length(trees)
  
  # Stop nesting problems in case of multiple trees
  if (!is.na(pmatch("tree", names(trees))) & (length(trees) == 1)) trees <- trees[[1]]
  
  
  # Selecting just one tree
  if (single_tree) {
    
    # Creating the predictions vector
    predictions <- rep(0, nrow(x))
    
    # Getting terminal node indexes]
    terminal_nodes <- trees[names(which(sapply(trees, function(x) {
      x$terminal == 1
    })))]
    
    for (i in 1:length(terminal_nodes)) {
      
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
my_rBart_gp <- function(x, y,
                        number_trees = 2, # Setting the number of trees
                        node_min_size = round(0.05 * length(y)), # Min node size,
                        mu = 0,
                        alpha = 0.95, # Alpha from prior
                        beta = 2, # Beta from prior
                        tau = 1, # Tau from prior,
                        nu_vector = NULL, # This will be defining the nu the default value
                        phi_vector = rep(0.1 / (sqrt(number_trees)), number_trees), # This parameter is a "scale paramter" to the GP
                        n_iter = 2000, # Number of iterations
                        burn = 500, # Number of burn
                        thin = 1, # Number of thin
                        rotation = TRUE, # If rotated lon and lat will be used in tree building
                        theta = NULL, # If theta is NULL, then the rotation angle will be randomly selected
                        seed = NULL, # Alpha vector values from the Dirichlet prior
                        scale_boolean = TRUE,
                        a_tau = 1, # Prior from a_v_ratio gamma
                        d_tau = 3, # Prior from d_v_ratio gamma,
                        predictions=matrix(0, nrow = number_trees, ncol = nrow(x)),
                        predictions_list =  NULL,
                        discrete_phi_boolean = FALSE,
                        nu_multiplier = 1,tau_multiplier = 1, scale_multiplier = 1,
                        x_scale =  TRUE,
                        nu_update = TRUE,
                        phi_update = TRUE,
                        gp_variables = colnames(x),   # Selecting the GP-Variables
                        p = 1, # Shrink main parameter from GP-BART
                        K_bart = 2,
                        prob_tau = 0.9,
                        error_handling_residuals = FALSE,
                        kappa = 0.5, bart_boolean = TRUE, bart_number_iter = 250
                        
) {
  
  # If there's only one covariate
  if(dim(x)[2]==1){
    rotation <- FALSE
  }
  
  # Saving mean values from x
  mean_x <- colMeans(x)
  sd_x <- apply(x,2,sd)
  
  
  # Scaling the x
  if(x_scale){
    x <- t(apply(x, 1, function(x) {(x - mean_x)/sd_x}))
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
  
  # Scale values
  if(scale_boolean){
    # Normalizing y
    a_min <- min(y)
    b_max <- max(y)
    
    y_scale <- normalize_bart(y = y)
    
    
    # Defing the nu vector if not in default
    if( is.null(nu_vector) ){
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4*number_trees*K_bart^2)/((max(y_scale)-min(y_scale))^2),number_trees)/(1-kappa)
    }
    
    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- ((4*(number_trees)*K_bart^2)/((max(y_scale)-min(y_scale))^2))/kappa
    
    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
    
  } else {
    
    # Not scaling the y
    y_scale <- y
    
    if( is.null(nu_vector) ){
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4*number_trees*K_bart^2)/((max(y_scale)-min(y_scale))^2),number_trees)
    }
    
    
    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- (4*(number_trees)*K_bart^2)/((1-kappa)*(max(y_scale)-min(y_scale))^2)
    
    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
    
  }
  
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
  current_partial_residuals_matrix <- matrix(NA, nrow = number_trees, ncol = nrow(x))
  current_predictions_matrix <- matrix(NA, nrow = number_trees, ncol = nrow(x))
  
  
  # Creating the predictions saving
  current_partial_residuals_list <- list()
  current_predictions_list <- list()
  
  # Error of the matrix
  if (is.null(colnames(x))) {
    stop("Insert a valid NAMED matrix")
  }
  
  if (node_min_size == 0) {
    stop("Node Minimum Size need to be greater than 0")
  }
  
  if (length(nu_vector) != number_trees) {
    stop("Insert a valid \\nu vector for the number of trees")
  }
  
  if (length(phi_vector) != number_trees) {
    stop("Insert a valid \\phi vector for the number of trees")
  }
  
  
  # Recommendation about the min_node_size
  if (node_min_size < 15) {
    warning("\n It is recommended that the min_node_size should be of at least 15 observations.")
  }
  
  # Storage containers
  store_size <- (n_iter - burn) / thin
  tree_store <- vector("list", store_size)
  tau_store <- c()
  signal_pc_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  
  y_hat_store <- matrix(NA, ncol = length(y), nrow = store_size)
  y_hat_store_proposal <- matrix(NA, ncol = length(y), nrow = store_size)
  
  # Storing the likelihoods
  log_lik_store <- rep(NA, store_size)
  log_lik_store_fixed_tree <- rep(NA,store_size)
  
  
  # Getting the numbers
  loglike_fixed_tree_residuals <- numeric()
  loglike_tree_residuals <- numeric()
  loglike_fixed_tree_residuals_matrix <- matrix(NA, nrow = store_size, ncol = number_trees)
  loglike_tree_residuals_matrix <- matrix(NA, nrow = store_size, ncol = number_trees)
  
  full_cond_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  phi_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  phi_proposal_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  phi_vector_proposal <- rep(0.1, number_trees)
  
  # Creating the list of trees stumps
  for (i in 1:number_trees) {
    # Creating the fixed two split trees
    current_trees[[i]] <- stump(x = x, tau = tau, mu = mu)
    current_trees_proposal[[i]] <- stump(x = x, tau = tau, mu =  mu)
  }
  
  
  names(current_trees) <- (sapply(1:number_trees, function(x) paste0("tree_", x))) # Naming each tree
  names(current_trees_proposal) <- (sapply(1:number_trees, function(x) paste0("tree_", x))) # Naming each tree
  
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
  for (i in 1:n_iter) {
    
    utils::setTxtProgressBar(progress_bar, i)
    
    # Changing the bart boolean, when reach the maximum
    if(i >= bart_number_iter){
      bart_boolean <- FALSE
    }
    
    if ((i > burn) & ((i %% thin) == 0)) {
      
      # Saving the store of the other ones
      curr <- (i - burn) / thin
      tree_store[[curr]] <- current_trees
      tau_store[[curr]] <- tau
      
      y_hat_store[curr, ] <-  colSums(predictions)
      
      # Saving the current partial
      current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix
      
      
      # Saving the predictions
      current_predictions_list[[curr]] <- predictions
      
      phi_store[curr, ] <- phi_vector
      phi_proposal_store[curr, ] <- phi_vector_proposal
      verb_store_list[[curr]] <- verb_store
      
      
    }
    
    
    
    # Creating a boolean to create the first trees only using BART model
    if (bart_boolean){
      
      # Verb Store
      verb_store <- data.frame(verb = rep(NA,number_trees),
                               accepted = rep(NA,number_trees),
                               identical = rep(NA,number_trees))
      
      for (j in 1:number_trees) {
        
        # Getting the verb list
        
        # Calculating the residuals for each tree
        if (number_trees > 1) {
          
          # Calculating what Chipman called as R(j) = y - g_others_trees
          if (number_trees > 2) {
            # current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }
        
        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow","grow_projection", "prune", "change","change_projection","swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change","swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }
        
        # Case of rotation
        if(rotation){
          if (i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- sample(c("grow","grow_projection"),
                                                                                               size = 1) # Grow the tree for the first few iterations
        } else {
          if (i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }  
        
        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate
        
        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x = x,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )
        
        
        
        # Calculating the likelihood of the new tree
        likelihood_new <- tree_complete_conditional_bart(
          tree = new_trees[[j]], # Calculate the full conditional
          residuals_values = current_partial_residuals,
          x = x,tau_mu = tau_mu, tau = tau
        )
        
        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional_bart(
          tree = current_trees[[j]], # Calculate the full conditional
          residuals_values = current_partial_residuals,
          x = x,tau_mu = tau_mu, tau = tau
        )
        
        # Extracting only the likelihood
        l_new <- likelihood_new+
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
        
        
        
        # Probability of accept the new proposed tree
        acceptance <- exp(l_new - l_old)
        
        # If Storage or not based on thin and burn parameters
        if ((i > burn) & ((i %% thin) == 0)) {
          full_cond_store[curr, j] <- l_old
        }
        
        # In case of acceptance
        
        if (runif(1) < acceptance) {
          
          # Counting acc ratio
          acc_ratio <- acc_ratio + 1
          
          # Make changes if accept
          current_trees <- new_trees
          
          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- TRUE
          
          
          # Storing likelihood matrix objects
          likelihood_object <- likelihood_new
          
        } else {
          
          
          # Storing likelihood matrix objects
          likelihood_object <- likelihood_old
          
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
      
      for (j in 1:number_trees) {
        
        # Calculating the residuals for each tree
        if (number_trees > 1) {
          
          # Calculating what Chipman called as R(j) = y - g_others_trees
          if (number_trees > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }
        
        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow","grow_projection", "prune", "change","change_projection","swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change","swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }
        
        # Case of rotation
        if(rotation){
          if (i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- sample(c("grow","grow_projection"),
                                                                                               size = 1) # Grow the tree for the first few iterations
        } else {
          if (i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }  
        
        
        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate
        
        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x = x,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )
        
        # ==================== #
        # Getting the Omega Inverse the current and the future tree
        # ==================== #
        
        # Getting the inverse for the current terminal nodes
        current_trees[[j]] <- inverse_omega_plus_I(tree = current_trees[[j]],
                                                   x = x,tau = tau,
                                                   nu = nu_vector[j],
                                                   phi = phi_vector[j])
        
        # Getting the inverse for the new tree terminal nodes
        new_trees[[j]] <- inverse_omega_plus_I(tree = new_trees[[j]],
                                               x = x,tau = tau,
                                               nu = nu_vector[j],
                                               phi = phi_vector[j])
        
        
        
        # Calculating the likelihood of the new tree
        likelihood_new <- tree_complete_conditional_gpbart(
          tree = new_trees[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          nu = nu_vector[j], phi = phi_vector[j], 
          number_trees = number_trees,
          tau_multiplier = tau_multiplier
        )
        
        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional_gpbart(
          tree = current_trees[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          nu = nu_vector[j], phi = phi_vector[j], 
          number_trees = number_trees,
          tau_multiplier = tau_multiplier
        )
        
        # Extracting only the likelihood
        l_new <- likelihood_new$log_posterior+
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
        
        
        
        # Probability of accept the new proposed tree
        acceptance <- exp(l_new - l_old)
        
        # If Storage or not based on thin and burn parameters
        if ((i > burn) & ((i %% thin) == 0)) {
          full_cond_store[curr, j] <- l_old
        }
        
        if (runif(1) < acceptance) { #
          acc_ratio <- acc_ratio + 1
          
          
          # Checking whether the trees are identical
          if(identical(current_trees[[j]],new_trees[[j]])){
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
          likelihood_object <- likelihood_new
          
        } else {
          # Storing likelihood matrix objects
          likelihood_object <- likelihood_old
          
          
          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- FALSE
          verb_store[j,"identical"] <- FALSE
          
        } # End of accept if statement
        
        # # # To update the mu values
        current_trees[[j]] <- update_mu(
          tree = current_trees[[j]],
          x = x,
          likelihood_object = likelihood_object)
        
        
        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j, ] <- update_residuals(
          tree = current_trees[[j]], x = x,
          residuals = current_partial_residuals,
          phi = phi_vector[j], nu = nu_vector[j],tau = tau,
          error_handling_residuals = error_handling_residuals # Boolean to check if handle the residuals or not
        )
        
        
        # To update phi
        if(phi_update){
          mh_update_phi<- update_phi_marginal(current_tree_iter = current_trees[[j]],
                                              residuals = current_partial_residuals,
                                              x = x,nu = nu_vector[j],phi = phi_vector[j],
                                              gp_variables = gp_variables,
                                              likelihood_object = likelihood_object,
                                              number_trees = number_trees,
                                              discrete_phi = discrete_phi_boolean,
                                              tau = tau,
                                              tau_mu = tau_mu,
                                              distance_min = distance_min,
                                              distance_max = distance_max
          )
          
          #In case of accept the update over \phi update everything
          if(mh_update_phi$phi_boolean){
            
            # Updating the tree and the \phi object from the tree
            current_trees[[j]] <- mh_update_phi$tree
            
            # Updating the likelihood objects
            likelihood_object <- mh_update_phi$likelihood_object
            
            # Updating the phi value
            phi_vector[j] <- mh_update_phi$phi_proposal
            
          } # If doesn't accept, nothing changes.
          
        }
        
        # Getting the residuals likelihood
        # loglike_tree_residuals[j] <- loglike_residuals(tree = current_trees[[j]],
        #                                                x = x,current_partial_residuals = current_partial_residuals,
        #                                                phi = phi_vector[j],
        #                                                nu = nu_vector[j], p = p )
        
        
        # current_partial_residuals_matrix<-
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]
        
        current_trees[[j]] <- remove_omega_plus_I_inv(current_tree_iter = current_trees[[j]])
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
    x <- t(apply(x, 1, function(y) {y * sd_x + mean_x}))
  }
  return(list(
    trees = tree_store,
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
    X = x,
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
    number_trees = number_trees, node_min_size = node_min_size, 
    a_min = a_min,
    b_max = b_max,
    a_tau = a_tau,
    d_tau = d_tau,
    current_partial_residuals_list = current_partial_residuals_list,
    beta = beta, 
    current_predictions_list = current_predictions_list,
    nu_multiplier = nu_multiplier , scale_multiplier = scale_multiplier,
    p = p, tau_mu = tau_mu, kappa = kappa,
    verb_store_list = verb_store_list
  ))
}


# #Do a MH for PHI
update_phi_marginal <- function(x, current_tree_iter,residuals,
                                seed = NULL,
                                tau, 
                                tau_mu,
                                phi, nu,number_trees,
                                likelihood_object,p,gp_variables,
                                discrete_phi = TRUE,
                                distance_min,
                                distance_max
){
  # Increased the range of tree proposal
  if(discrete_phi){
    phi_proposal <- sample(c(0.1,0.5,1,5,10),size = 1)
  } else {
    phi_proposal <- runif(1,min = distance_min,max = distance_max)
  }
  
  # Calculating the likelihood from the new step
  tree_from_phi_proposal <- inverse_omega_plus_I(tree = current_tree_iter,
                                                 x = x,nu = nu, tau = tau,
                                                 phi = phi_proposal,gp_variables = gp_variables)
  
  
  likelihood_phi_proposal <- tree_complete_conditional_gpbart(tree = tree_from_phi_proposal,
                                                              x = x,
                                                              residuals = residuals,
                                                              nu = nu, tau_mu = tau_mu,
                                                              phi = phi_proposal,
                                                              number_trees = number_trees,tau_multiplier = tau_multiplier)
  
  # Old phi likelhood
  l_old_phi <- likelihood_object$log_posterior
  
  # Proposal likelihood
  l_proposal_phi <- likelihood_phi_proposal$log_posterior
  
  
  # Rejecting the new value
  if(l_proposal_phi == -Inf){
    # Case of not accepting
    phi_boolean <- FALSE
    return( list(phi_boolean = phi_boolean)) # Returning the old value for phi
    
  }
  
  # Probability of accept the new proposed tree
  acceptance_phi <- exp(l_proposal_phi- l_old_phi)
  
  # If storage for phi
  
  if (runif(1) < acceptance_phi) { #
    
    # Nu boolean to see if was accepted or not
    phi_boolean <- TRUE
    
    return(list(phi_boolean = phi_boolean,
                likelihood_object = likelihood_phi_proposal,
                tree = tree_from_phi_proposal,
                phi_proposal = phi_proposal)) # Returning the proposal value for phi
  }else{
    # Case of not accepting
    phi_boolean <- FALSE
    return( list(phi_boolean = phi_boolean)) # Returning the old value for phi
  } #
}


# #Do a MH for PHI
update_phi <- function(x, current_tree_iter,residuals,
                       seed = NULL,
                       tau, 
                       tau_mu,
                       phi, nu,number_trees,
                       likelihood_object,p,gp_variables,
                       discrete_phi = TRUE,
                       distance_min,
                       distance_max
){
  
  # Increased the range of tree proposal
  if(discrete_phi){
    phi_proposal <- sample(c(0.1,0.5,1,5,10),size = 1)
  } else {
    phi_proposal <- runif(1,min = distance_min,max = distance_max)
  }
  
  # Calculating the likelihood from the new step
  tree_from_phi_proposal <- inverse_omega_plus_I(tree = current_tree_iter,
                                                 x = x,nu = nu, tau = tau,
                                                 phi = phi_proposal,gp_variables = gp_variables)
  
  
  likelihood_phi_proposal <- tree_complete_conditional_gpbart(tree = tree_from_phi_proposal,
                                                              x = x,
                                                              residuals = residuals,
                                                              nu = nu, tau_mu = tau_mu,
                                                              phi = phi_proposal, 
                                                              number_trees = number_trees,tau_multiplier = tau_multiplier)
  
  # Selecting the terminal nodes
  terminal_nodes <- current_tree_iter[names(which(sapply(current_tree_iter, function(x) {
    x$terminal == 1
  })))]
  
  # Selecting the new terminal nodes
  terminal_nodes_phi_proposal <- tree_from_phi_proposal[names(which(sapply(tree_from_phi_proposal, function(x) {
    x$terminal == 1
  })))]
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  # Calculating the old likelihood
  g_mean_old <- mapply(terminal_nodes,
                       residuals_terminal_nodes,
                       FUN = function(nodes,residuals_nodes){
                         nodes$mu + crossprod(nodes$Omega_matrix, crossprod(nodes$Omega_plus_I_inv, (residuals_nodes - nodes$mu)))
                       },SIMPLIFY = FALSE)
  
  g_variance_old <- mapply(terminal_nodes,
                           FUN = function(nodes){
                             (nodes$Omega_matrix) - crossprod(nodes$Omega_matrix,
                                                              crossprod(nodes$Omega_plus_I_inv, nodes$Omega_matrix))
                           },SIMPLIFY = FALSE)
  
  
  # Calculating the new likelihood
  g_mean_new <- mapply(terminal_nodes_phi_proposal,
                       residuals_terminal_nodes,
                       FUN = function(nodes,residuals){
                         nodes$mu + crossprod(nodes$Omega_matrix, crossprod(nodes$Omega_plus_I_inv, (residuals - nodes$mu)))
                       },SIMPLIFY = FALSE)
  
  g_variance_new <- mapply(terminal_nodes_phi_proposal,
                           FUN = function(nodes){
                             (nodes$Omega_matrix) - crossprod(nodes$Omega_matrix,
                                                              crossprod(nodes$Omega_plus_I_inv, nodes$Omega_matrix))
                           },SIMPLIFY = FALSE)
  
  log_likelihood_old <- mapply(terminal_nodes, residuals_terminal_nodes,
                               FUN = function(nodes,residuals){
                                 mvtnorm::dmvnorm(x = residuals,
                                                  mean = rep(nodes$mu,dim(nodes$Omega_matrix)[1]),
                                                  sigma = nodes$Omega_matrix + diag(nodes$tau^-1, 
                                                                                    nrow = dim(nodes$Omega_matrix)[1]),log = TRUE)
                               })
  
  log_likelihood_new <- mapply(terminal_nodes_phi_proposal, residuals_terminal_nodes,
                               FUN = function(nodes,residuals){
                                 mvtnorm::dmvnorm(x = residuals,
                                                  mean = rep(nodes$mu,dim(nodes$Omega_matrix)[1]),
                                                  sigma = nodes$Omega_matrix + diag(nodes$tau^-1, 
                                                                                    nrow = dim(nodes$Omega_matrix)[1]),log = TRUE)
                               })
  
  
  # Old phi likelhood
  l_old_phi <- sum(log_likelihood_old)
  
  # Proposal likelihood
  l_proposal_phi <- sum(log_likelihood_new)
  
  
  # Rejecting the new value
  if(l_proposal_phi == -Inf){
    # Case of not accepting
    phi_boolean <- FALSE
    return( list(phi_boolean = phi_boolean)) # Returning the old value for phi
    
  }
  
  # Probability of accept the new proposed tree
  acceptance_phi <- exp(l_proposal_phi- l_old_phi)
  
  # If storage for phi
  
  if (runif(1) < acceptance_phi) { #
    
    # Nu boolean to see if was accepted or not
    phi_boolean <- TRUE
    
    return(list(phi_boolean = phi_boolean,
                likelihood_object = likelihood_phi_proposal,
                tree = tree_from_phi_proposal,
                phi_proposal = phi_proposal)) # Returning the proposal value for phi
  }else{
    # Case of not accepting
    phi_boolean <- FALSE
    return( list(phi_boolean = phi_boolean)) # Returning the old value for phi
  } #
}

# Function to return the depth trees

tree_depth_hist <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)
  
  for (k in 1:length(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in 1:gpbart_model$number_trees) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- lapply(tree, function(x) {
        x$depth
      }) %>%
        unlist() %>%
        max()
    }
  }
  return(tree_depth)
}

# Get the covariate splits
tree_var_hist <- function(gpbart_model) {
  tree_depth <- c()
  
  for (k in 1:length(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in 1:gpbart_model$number_trees) {
      tree <- tree_iter[[i]]
      tree_depth <- c(tree_depth,lapply(tree, function(x) {
        x$node_var
      }) %>%
        unlist() )
    }
  }
  return(tree_depth)
}

# Function to count th enumber of terminal nodes in a tree
tree_count_terminals <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)
  
  for (k in 1:length(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in 1:gpbart_model$number_trees) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- lapply(tree, function(x) {
        x$termina
      }) %>%
        unlist() %>%
        sum()
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
  
  for (j in 1:n_iter) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for tau terminals
    tau_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(tau_terminals) <- paste0("node_", 0:49)
    
    for (k in 1:num_trees) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)
      
      
      for (nodes in 1:length(all_nodes)) {
        if (tree[[all_nodes[nodes]]]$terminal == 1) {
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
  
  for (j in 1:n_iter) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for mu terminals
    mu_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(mu_terminals) <- paste0("node_", 0:49)
    
    for (k in 1:num_trees) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)
      
      
      for (nodes in 1:length(all_nodes)) {
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
  y_pred_final <- matrix(0, nrow = length(multiple_trees), ncol = nrow(x_new))
  
  # Calculating the sd from the prediction interval
  y_pred_pi_final <- matrix(0, nrow = length(multiple_trees), ncol = nrow(x_new))
  
  
  # Covariance matrix
  cov_pred_final <- list()
  variance <- matrix(0, nrow = nrow(x_new), ncol = nrow(x_new))
  
  
  # Iterating over all trees
  for (m in 1:length(multiple_trees)) {
    
    # Creating the list to be predicted (The if is just in case of of just one tree)
    new_tree <- multiple_trees[[m]]
    phi <- phi_vector[m]
    nu <- nu_vector[m]
    
    # Creating the pred  vector
    y_pred <- numeric(nrow(x_new))
    y_pred_pi <- numeric(nrow(x_new))
    
    variance <- matrix(0, nrow = nrow(x_new), ncol = nrow(x_new))
    
    # Setting the root node with the new observations
    new_tree[["node_0"]]$test_index <- 1:nrow(x_new)
    
    # Creating the list of nodes
    list_nodes <- names(new_tree)[-1]
    
    # IN CASE OF JUST ROOT NODE
    if (length(new_tree) == 1) {
      list_nodes <- "node_0"
    }
    
    # Updating all nodes
    for (i in 1:length(list_nodes)) {
      
      current_node_aux <- new_tree[[list_nodes[i]]]
      
      # In case of more than one node
      if (length(list_nodes) > 1) {
        # Veryfing the type of the current node
        if (is.list(current_node_aux$node_var)) {
          
          # Rotated Lon and Lat
          
          rotated_x <- tcrossprod((A(current_node_aux$theta)), x_new[,current_node_aux$node_var$node_var_pair])
          rownames(rotated_x) <- current_node_aux$node_var$node_var_pair
          
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] >= current_node_aux$node_var_split)]
          }
        } else { # Evaluating the case where is not a rotated lat/lon
          
          # To continous covariates
          if (is.numeric(x_new[, current_node_aux$node_var])) {
            
            # Updating observations from the left node
            if (current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] < current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
            }
            
            
            # To categorical covariates
          } else {
            # Updating observations from the left node
            if (current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
            }
          }
        }
      }
      # Here I will calculate the predicted values based on the terminal nodes AND only on terminal nodes which have test
      
      if (new_tree[[list_nodes[[i]]]]$terminal == 1 & length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {
        
        # Selecting the observations from the current node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED
        x_current_node <- matrix(x_train[new_tree[[list_nodes[[i]]]]$observations_index, ],
                                 nrow = length(new_tree[[list_nodes[[i]]]]$observations_index)
        )
        
        # Selecting the observations from the test node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED
        
        x_star <- matrix(x_new[new_tree[[list_nodes[[i]]]]$test_index, ], nrow = length(new_tree[[list_nodes[[i]]]]$test_index))
        
        
        
        # Calcualting the distance matrix
        distance_matrix_current_node <- symm_distance_matrix(m1 = x_current_node)
        
        # Generating the scenario where there are only BART predictions
        if(!pred_bart_only){
          
          # Getting the GP from a terminal node
          gp_process <- gp_main(
            x_train = x_current_node,distance_matrix_train = distance_matrix_current_node,
            y_train = (matrix((partial_residuals[m,new_tree[[list_nodes[[i]]]]$observations_index]) - new_tree[[list_nodes[[i]]]]$mu,
                              nrow = nrow(x_current_node)
            )), x_star = x_star, tau = tau,
            nu = nu, phi = phi
          ) 
          
          # Creating the mu vector
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- gp_process$mu_pred + new_tree[[list_nodes[[i]]]]$mu
          
        } else {
          
          # Attributing the predicted value as the sampled \mu from the terminal node
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- new_tree[[list_nodes[[i]]]]$mu
          
        }
        
        # Creating the sd of the mu vector
        y_pred_pi[new_tree[[list_nodes[i]]]$test_index] <- (1/tau)
        
        if(!pred_bart_only){
          # Variance is a matrix n_testXn_test
          variance[new_tree[[list_nodes[i]]]$test_index, new_tree[[list_nodes[i]]]$test_index] <- gp_process$cov
        }
      }
    }
    # points(x_new,y_pred,col='red')
    # points(x_current_node,g_values,col='orange',pch=20)
    
    y_pred_pi_final[m,] <- y_pred_pi
    y_pred_final[m, ] <- y_pred
    cov_pred_final[[m]] <- variance
  }
  # Return the new tree
  return(list(y_pred = colSums(y_pred_final),
              y_pred_pi_sd = y_pred_pi_final, new_tree = new_tree,
              all_tree_pred = y_pred_final, variance_matrix = cov_pred_final))
  # return(y_pred)
}


# Function to count the number of terminal nodes
count_terminal_nodes <- function(tree) {
  return(sum(unlist(lapply(tree, function(x) {
    x$terminal == 1
  }))))
}

# Predict rBART model

my_predict_rBART <- function(rBart_model, x_test, type = c("all", "mean", "median"),
                             pred_bart_only = FALSE) {
  
  # Loading x_train
  if(rBart_model$x_scale){
    x_train <- t(apply(rBart_model$X, 1, function(z){ (z - rBart_model$mean_x)/rBart_model$sd_x}))
    x_test  <- t(apply(x_test, 1, function(z){ (z - rBart_model$mean_x)/rBart_model$sd_x}))
    
    # Retrieving the col.names
    colnames(x_train) <- 
    colnames(x_test)  <- colnames(rBart_model$X)
  } else {
    x_train <- rBart_model$X
  }

  # Number of iters of bayesian simulation
  n_iter <- length(rBart_model$tau_store)
  
  # The number of columns is the number of test observations and the rows are the iterations
  y_hat_matrix <- matrix(0, nrow = n_iter, ncol = nrow(x_test))
  y_sd_pi_matrix <- matrix(0, nrow = n_iter, ncol = nrow(x_test))
  
  # Getting the training objects
  y_hat_matrix_train <- matrix(0, nrow = n_iter, ncol = nrow(x_train))
  y_sd_pi_matrix_train <- matrix(0, nrow = n_iter, ncol = nrow(x_train))
  
  cov_hat_matrix <- matrix(0, nrow = n_iter, ncol = nrow(x_test))
  cov_hat_matrix_list <- list()
  cov_hat_matrix_sqrt_sum <- list()
  
  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running rBART..."
  )
  y_list_matrix <- list()
  
  # Creating the final vector
  y_pred_final <- matrix(0, nrow = rBart_model$number_trees, ncol = nrow(x_test))
  cov_pred_final <- list()
  variance <- matrix(0, nrow = nrow(x_test), ncol = nrow(x_test))
  
  # Looping around the trees
  for (i in 1:n_iter) {
    
    utils::setTxtProgressBar(progress_bar, i)
    
    
    
    # Selecting one tree from BART model
    current_tree <- rBart_model$trees[[i]]
    
    # Getting the predictions from the test observations
    y_pred_aux <- predict_gaussian_from_multiple_trees(
      multiple_trees = current_tree, x_train = x_train,
      x_new = x_test, partial_residuals = rBart_model$current_partial_residuals_list[[i]],
      phi_vector = rBart_model$phi_store[i, ], nu_vector = rBart_model$nu_vector,
      tau = rBart_model$tau_store[[i]],pred_bart_only = pred_bart_only
    )
    
    # Iterating over all trees (test)
    y_pred_final <- y_pred_aux$all_tree_pred
    y_pred_final_pi <- y_pred_aux$y_pred_pi_sd
    cov_pred_final <- y_pred_aux$variance_matrix[[1]]
    
    
    if(rBart_model$scale_boolean){
      # Recovering the prediction interval from test
      y_hat_matrix[i, ] <- unnormalize_bart(colSums(y_pred_final), a = rBart_model$a_min, b = rBart_model$b_max)
      y_sd_pi_matrix[i, ] <- colSums(y_pred_final_pi)*sqrt((rBart_model$b_max-rBart_model$a_min))
    } else {
      y_hat_matrix[i, ] <- colSums(y_pred_aux$all_tree_pred)
      y_sd_pi_matrix[i, ] <- colSums(y_pred_final_pi)
    }
    
    y_list_matrix[[i]] <- y_pred_final
    
    cov_hat_matrix_list[[i]] <- Reduce("+", cov_pred_final)
    cov_hat_matrix_sqrt_sum[[i]] <- Reduce("+", lapply(cov_pred_final, function(x){sqrt(abs(x))}))
    
  }
  
  # Getting the parameters a and b
  a <- min(rBart_model$y)
  b <- max(rBart_model$y)
  
  
  out <- list(
    pred = switch(type,
                  all = y_hat_matrix,
                  mean = colMeans(y_hat_matrix),
                  median = apply(y_hat_matrix, 2, "median"),
                  
    ),
    sd = switch(type,
                all = cov_hat_matrix,
                mean = mean(unlist(rBart_model$tau_store)^(-1/2)),
                median = median(unlist(rBart_model$tau_store)^(-1/2))
    )
  )
  
  
  
  return(list(out = out, list_matrix_pred = y_list_matrix,
              list_matrix_cov = cov_hat_matrix_list,
              list_matrix_cov_sqrt_sum = cov_hat_matrix_sqrt_sum,
              mcmc_pi_mean =  y_hat_matrix,
              mcmc_pi_sd = y_sd_pi_matrix))
  
}

# A function to count the mean values of observations in terminal nodes
gpbart_count_terminal_nodes <- function(mod_gpbart){
  
  # Auxiliar matrix 
  all_tree_terminal_nodes <- array(NA,dim = c(mod_gpbart$number_trees,50,mod_gpbart$store_size),
                                   dimnames = list(paste0("tree_",1:mod_gpbart$number_trees),
                                                   paste0("node_",1:50),
                                                   paste0("iter_",1:mod_gpbart$store_size)))
  
  
  # Iterating over MH 
  for(k in 1:mod_gpbart$store_size){
    tree_iter <- mod_gpbart$trees[[k]]
    
    # Iterating with for over the terminal nodes
    for(i in 1:mod_gpbart$number_trees){
      
      # Gathering the node number
      node_number <- unlist(tree_iter[[i]] %>% lapply(function(x) { x[x$terminal == 1]$node_number})) 
      all_tree_terminal_nodes[i,node_number,k] <- unlist(tree_iter[[i]][names(node_number)] %>% lapply(function(x) {
        length(x[x$terminal == 1]$observations_index)
      }))
      
    }
  }
  # Three splits
  split_nodes <- unique(apply(all_tree_terminal_nodes, 2, mean,na.rm = TRUE))
  return(split_nodes[!is.na(split_nodes)])
  
}


# A function to return the split values in a tree
gpbart_split <- function(mod_gpbart, k){
  
  # Auxiliar matrix 
  all_tree_splits <- array(NA,dim = c(mod_gpbart$number_trees,50,mod_gpbart$store_size),
                           dimnames = list(paste0("tree_",1:mod_gpbart$number_trees),
                                           paste0("node_",1:50),
                                           paste0("iter_",1:mod_gpbart$store_size)))
  
  
  # Iterating over MH 
  # for(k in 1:gpbart_mod$store_size){
  tree_iter <- mod_gpbart$trees[[k]]
  
  # Iterating with for over the terminal nodes
  for(i in 1:mod_gpbart$number_trees){
    
    # Gathering the node number
    node_number <- unlist(tree_iter[[i]] %>% lapply(function(x) { x[x$terminal == 1]$node_number })) 
    all_tree_splits[i,node_number,k] <- unlist(tree_iter[[i]][names(node_number)] %>%
                                                 lapply(function(x) { x[x$terminal == 1]$node_var_split}))
    
  }
  # }
  # Three splits
  split_nodes <- apply(all_tree_splits, c(1,2), median,na.rm = TRUE)[,1]
  # split_nodes <- split_nodes[!is.na(split_nodes)]
  # split_nodes <- c(all_tree_splits)
  
  return(split_nodes[!is.na(split_nodes)])
  
}

gpbart_split_mean <- function(mod_gpbart){
  
  # Auxiliar matrix 
  all_tree_splits <- array(NA,dim = c(mod_gpbart$number_trees,50,mod_gpbart$store_size),
                           dimnames = list(paste0("tree_",1:mod_gpbart$number_trees),
                                           paste0("node_",1:50),
                                           paste0("iter_",1:mod_gpbart$store_size)))
  
  
  # Iterating over MH 
  for(k in 1:gpbart_mod$store_size){
    tree_iter <- mod_gpbart$trees[[k]]
    
    # Iterating with for over the terminal nodes
    for(i in 1:mod_gpbart$number_trees){
      
      # Gathering the node number
      node_number <- unlist(tree_iter[[i]] %>% lapply(function(x) { x[x$terminal == 1]$node_number })) 
      all_tree_splits[i,node_number,k] <- unlist(tree_iter[[i]][names(node_number)] %>%
                                                   lapply(function(x) { x[x$terminal == 1]$node_var_split}))
      
    }
  }
  # Three splits
  split_nodes <- apply(all_tree_splits, c(1,2), median,na.rm = TRUE)[,1]
  # split_nodes <- split_nodes[!is.na(split_nodes)]
  # split_nodes <- c(all_tree_splits)
  
  return(split_nodes[!is.na(split_nodes)])
  
}



# Retrieve new points on the cubic scale
unit_cube_scale_new_points_uni <- function(x,x_new){
  
  
  # Scaling to -1 to 1 function
  scaling <- 
    (2 * x_new - (max(x) + min(x))) / (max(x) - min(x))
  
  
  # Applying on all covariates
  
  return(scaling)
  
}


# Getting the likelihood from the Gaussian Processes
neg_loglike<-function(prediction_object,
                      y_test){
  
  # Getting the prediction
  y_pred <- colMeans(prediction_object$mcmc_pi_mean)
  
  # Getting the variance
  sd_pred <- colMeans(prediction_object$mcmc_pi_sd)
  
  return(-dnorm(x = y_test,
                mean = y_pred,
                sd = sqrt(sd_pred),log = TRUE))
}

# sum(neg_loglike(prediction_object = pred_gpbart,y_test = y_test))


# Get tau from terminal nodes
get_tau_values_from_single_tree <- function(gpbart_mod,tree_number, mh_iter = 100){
  # Getting the vector of tau values from the terminal nodes
  return(unlist(gpbart_mod$trees[[mh_iter]][[tree_number]] %>%
                  lapply(function(x) x[x$terminal==1]$tau)))
}

# Getting Omega Inverse argument lsit
# get_omega_inverse_terminal_nodes <- function(tree,
#                                              x = x,
#                                              nu = nu_vector[j], phi = phi_vector[j]){
#   
#   
#   # Selecting terminal nodes names
#   names_terminal_nodes <- names(which(sapply(tree, function(x) {
#     x$terminal == 1
#   })))
#   
#   # Selecting the terminal nodes
#   terminal_nodes <- tree[names(which(sapply(tree, function(x) {
#     x$terminal == 1
#   })))]
#   
#   # Number of nodes
#   n_node <- length(terminal_nodes)
#   
#   # Picking each node size
#   nodes_size <- sapply(terminal_nodes, function(x) {
#     length(x$observations_index)
#   })
#   
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
#   for (i in 1:length(names_terminal_nodes)) {
#     tree[[names_terminal_nodes[i]]]$Omega_inv <- Omega_matrix_INV[[names_terminal_nodes[i]]]
#   }
#   
#   return(tree)
#   
# }

# Getting Omega Inverse + Diag Inverse
inverse_omega_plus_I <- function(tree,
                                 x = x,
                                 nu, phi,
                                 tau,
                                 gp_variables = colnames(x)  # Selecting which gp-variables to use
                                 
){
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })
  
  # Calculating Omega matrix INVERSE
  distance_matrix <- mapply(terminal_nodes, FUN = function(y) {
    symm_distance_matrix(matrix(x[y$observations_index, gp_variables], nrow = length(y$observations_index)))
  }, SIMPLIFY = FALSE)
  
  # Calculating Omega
  Omega_matrix <- mapply(distance_matrix, FUN = function(dist_m){
    ( kernel_function(
      squared_distance_matrix =  dist_m,
      nu = nu, phi = phi)
    )
  },SIMPLIFY = FALSE)
  
  # Calculating Omega_plus_I*tau^(-1)
  Omega_matrix_plus_I <- mapply(Omega_matrix, FUN = function(omega){
    ( omega + diag((1/tau),nrow = dim(omega)[1])
    )
  },SIMPLIFY = FALSE)
  
  # print(Omega_matrix_plus_I[[1]][1:5,1:5])
  
  # Calculating Omega matrix plus I INVERSE
  Omega_matrix_plus_I_INV <- mapply(Omega_matrix_plus_I,  FUN = function(omega_plus_I_tau) { # p is the shrinkage factor
    chol2inv(  
      chol(
        omega_plus_I_tau
      )
    )
  }, SIMPLIFY = FALSE)
  
  # Adding the Omega_matrix_plus_I_Inv
  for (i in 1:length(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$Omega_plus_I_tau <- Omega_matrix_plus_I[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$Omega_plus_I_inv <- Omega_matrix_plus_I_INV[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$distance_matrix <- distance_matrix[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$Omega_matrix <- Omega_matrix[[names_terminal_nodes[i]]]
    
  }
  
  return(tree)
}


# # Removing the Omega_plus_I_inv object
remove_omega_plus_I_inv <- function(current_tree_iter){
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(current_tree_iter, function(x) {
    x$terminal == 1
  })))
  
  for(i in names_terminal_nodes){
    current_tree_iter[[i]]$Omega_plus_I_inv <- NULL
    current_tree_iter[[i]]$distance_matrix <- NULL
    current_tree_iter[[i]]$Omega_plus_I_tau <- NULL
    current_tree_iter[[i]]$Omega_matrix <- NULL
    
  }
  
  return(current_tree_iter)
}


# Read one of the trees as example
rMVN2 <- function(b, Q) {
  Q_inv <- chol2inv(chol(Q))
  p   <- ncol(Q)
  U   <- chol(Q_inv)
  z   <- rnorm(p)
  backsolve(U,matrix(z),transpose = FALSE, k=p)+b
}





# Creating one tree with 2 terminal nodes only
tree_two_terminal <- function(x,
                              node_split, # Setting the value of the node split
                              tau_left, # Setting the left value for tau
                              tau_right, # Settting the right value for tau
                              mu_left, # Set the left. for mu
                              mu_right, # Setting the right for mu
                              node_min_size = 5){
  # Creating a stump object from the tree
  stump_object <- stump(x = x,tau = 1)
  
  # Growing to the tree with a fixed split
  two_terminal_tree <- grow_tree_fixed_split(tree = stump_object,
                                             x = x,
                                             node_min_size = node_min_size,
                                             rotation = FALSE,
                                             theta = NULL,
                                             fixed_split = node_split)
  
  # Defining the fixed values for the parameters
  if(two_terminal_tree[["node_1"]]$left == 1){
    
    # Defining the fixed parameters for the left node
    two_terminal_tree[["node_1"]]$tau <- tau_left
    two_terminal_tree[["node_1"]]$mu <- mu_left
    
  } else {
    
    # Defining the fixed parameters for the left node
    two_terminal_tree[["node_1"]]$tau <- tau_right
    two_terminal_tree[["node_1"]]$mu <- mu_right
    
  }
  
  # Defining the fixed parameters for the second node
  
  # Defining the fixed values for the parameters
  if(two_terminal_tree[["node_2"]]$left == 1){
    
    # Defining the fixed parameters for the left node
    two_terminal_tree[["node_2"]]$tau <- tau_left
    two_terminal_tree[["node_2"]]$mu <- mu_left
    
  } else {
    
    # Defining the fixed parameters for the left node
    two_terminal_tree[["node_2"]]$tau <- tau_right
    two_terminal_tree[["node_2"]]$mu <- mu_right
    
  }
  
  # Returing the terminal tree
  return(two_terminal_tree)
}

# Creating the tree one
# tree_two_terminal(x = x,node_split = 0.3,
#                   tau_left = 1,
#                   tau_right = 5,
#                   mu_left = -10,
#                   mu_right = 5)

# Getting the average values from the training predictions
gpbart_train_mean <- function(gpbart_mod){
  colSums(Reduce("+",gpbart_mod$current_predictions_list)/length(gpbart_mod$current_predictions_list))
}


# Calculating the variance from the training set.
gpbart_training_var <-  function(gpbart_mod){
  
  # Number of MCMC samples
  n_mcmc <- length(gpbart_mod$trees)
  # Creating the var vector
  var_train <- matrix(0,
                      nrow = n_mcmc,
                      ncol = length(gpbart_mod$y))
  
  # Iterating over all trees and getting the terminal nodes
  for(i in 1:n_mcmc){
    
    for(m in 1:gpbart_mod$number_trees){
      # Selecting the current tree 
      tree <- gpbart_mod$trees[[i]][[m]]
      
      # Selecting the terminal nodes
      terminal_nodes <- tree[names(which(sapply(tree, function(x) {
        x$terminal == 1
      })))]
      
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



# # Getting the PI estimation
# pi_calibration <- function(gpbart_mod,pi_perc){
#   
#   # Calculating alpha as 1-pi_perc
#   alpha <- 1-pi_perc
#   
#   # Calculating the mean pred vaue and the var
#   train_mean_pred <- gpbart_train_mean(gpbart_mod = gpbart_mod)
#   train_var <- apply(new_gpbart_training_var(gpbart_mod),2,mean)
#   
#   # Creating the alpha% PI interval
#   pi_up <- train_mean_pred + qnorm(1-alpha/2) * sqrt(train_var)
#   pi_low <- train_mean_pred + qnorm(alpha/2) * sqrt(train_var)
#   
#   # Returning the proportion
#   prop_alpha <- sum( (pi_low <= gpbart_mod$y) & (gpbart_mod$y <= pi_up) )/length(gpbart_mod$y)
#   
#   # Returning the proportion within the alpha% PI interval
#   return(prop_alpha)
# }



# Calculating the residuals log-likelihood
loglike_residuals <- function(tree,
                              x,
                              current_partial_residuals,
                              phi, 
                              nu,p){
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    current_partial_residuals[x$observations_index]
  })
  
  # Getting the \tau values
  tau_terminal <- lapply(terminal_nodes, function(nodes){
    nodes$tau
  })
  
  # Getting the \mu values
  mu_terminal <- lapply(terminal_nodes, function(nodes){
    nodes$mu
  })
  
  Omega_matrix_plus_I <- lapply(terminal_nodes, function(nodes){
    kernel_function(
      squared_distance_matrix = nodes$distance_matrix,
      nu = nu, phi = phi) + diag(p,nrow = (dim(nodes$distance_matrix)[1]))
  })
  
  # Creating the loglikeresiduals_vec
  loglike_residuals_vec <- numeric()
  # Doing a quick for
  for( i in 1:length(residuals_terminal_nodes)){
    loglike_residuals_vec[i] <-  mvtnorm::dmvnorm(x = residuals_terminal_nodes[[i]],
                                                  mean = rep(mu_terminal[[i]],length(residuals_terminal_nodes[[i]])),
                                                  sigma = (1/tau_terminal[[i]])*Omega_matrix_plus_I[[i]],log = TRUE)
  }
  # Returning the sum of the log_like_residual_vec
  return(sum(loglike_residuals_vec))
}

# # Plot the fixed structure from a tree
# plot_tree_nodes <- function(gpbart_mod,
#                             tree_number){
#   
# }


# NEW UPDATE G
# Update tau_j values
update_g <- function(tree, x, nu, phi, residuals, seed = NULL,p) {
  # set.seed(seed)
  
  # New g (new vector prediction for g)
  g_new <- rep(NA, length(residuals))
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  
  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  
  
  # Getting the \mu_{j} vector
  mu_values <- unlist(
    lapply(terminal_nodes, function(x) {
      x$mu
    })
  )
  
  # Seleciting the tau values
  tau_j <- sapply(terminal_nodes, function(x) {
    x$tau
  })
  
  
  # Calculating Omega matrix
  Omega_matrix_inverse <- mapply(terminal_nodes, FUN = function(nodes) {
    chol2inv(chol(
      (kernel_function(squared_distance_matrix = nodes$distance_matrix,
                       nu = nu,
                       phi = phi)+diag(1e-4,nrow = nrow(nodes$distance_matrix))) # To be sure that is invertible
    ))
  }, SIMPLIFY = FALSE)
  
  # Getting the A matrix inverse
  A_matrix_inv <- mapply(Omega_matrix_inverse, FUN = function(x) {
    chol2inv(chol(diag(p,nrow = nrow(x)) + x))
  }, SIMPLIFY = FALSE)  
  
  # Calculating g_mean posterior
  g_mean <- mapply(A_matrix_inv,
                   residuals_terminal_nodes,
                   mu_values,Omega_matrix_inverse, FUN = function(A_inv, res, mu, omg) {
                     crossprod(
                       A_inv,
                       (p*res + mu*rowSums(omg))
                     )
                   }, SIMPLIFY = FALSE)
  
  
  
  # Putting in the Keefe's speed order
  g_sd <- mapply(tau_j, Omega_matrix_inverse, FUN = function(tau, omg) {
    (tau) * (omg + diag(p,nrow=nrow(omg)))
  }, SIMPLIFY = FALSE)
  
  g_sample <- mapply(g_mean, g_sd,
                     FUN = function(x, y) {
                       rMVN2(b = x,
                             Q = y)
                     },
                     SIMPLIFY = FALSE
  )
  
  
  # Adding the mu values calculated
  for (i in 1:length(terminal_nodes)) {
    # Saving g
    g_new[terminal_nodes[[i]]$observations_index] <- g_sample[[i]]
    
  }
  
  return(g_new)
}

# Some tests over nu parameter
calculate_nu <- function(nu){
  (nu+1)/nu
}

calculate_p_nu <- function(nu,p){
  (p*nu+1)/nu
}

# New PI calibration dataset
pi_calibration <- function(gpbart_mod,prob = 0.5, quantile = FALSE){
  
  # Decoding of I will get it from the quantile or not
  if(quantile){
    # Getting the quantile
    gpbart_sum_pred <- do.call(rbind,lapply(gpbart_mod$current_predictions_list, colSums))
    get_quantiles <- apply(gpbart_sum_pred,2, function(x) quantile(x, prob = c(prob/2,1-prob/2)))
    gpbart_up_ci <- get_quantiles[2,]
    gpbart_low_ci <- get_quantiles[1,]
    
    return( sum(((Ey <= gpbart_up_ci) & (Ey  >= (gpbart_low_ci)) )) / length(gpbart_mod$y))
    
  } else {
    # Getting the intervals from the training set
    gpbart_mod_train_sum_samples <- do.call(rbind,lapply(gpbart_mod$current_predictions_list, function(x){ apply(x,2,sum)}))
    gpbart_mod_mean <- gpbart_mod_train_sum_samples %>% colMeans
    gpbart_mod_train_low_ci <- gpbart_mod_mean + qnorm(prob/2)*sqrt(1/mean(gpbart_mod$tau_store))
    gpbart_mod_train_up_ci <- gpbart_mod_mean + qnorm(1-prob/2)*sqrt(1/mean(gpbart_mod$tau_store))
    
    
    return( sum((gpbart_mod$y  <= gpbart_mod_train_up_ci) & (gpbart_mod$y  >= gpbart_mod_train_low_ci)) / length(gpbart_mod$y)) 
  }
  
}

# Get train predictions
get_train_predictions <- function(gpbart_mod){
  
  # Getting the quantile
  gpbart_sum_pred <- do.call(rbind,lapply(gpbart_mod$current_predictions_list, function(x) {apply(x,2,sum)}))
  
  # Returning the matrix of final predictions
  return(gpbart_sum_pred)
}