# A file for a simple version implementation of BART algorithm ---
# with all functions needed to it
# Calculate Tree Prior
tree_prior <- function(tree, alpha, beta) {
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  # Selecting internal nodes names
  names_internal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 0
  })))
  
  # Selecting the depth of the terminal nodes
  depth_terminal <- sapply(tree[names_terminal_nodes], function(x) {
    x$depth_node
  })
  
  # Selecting the depth of the internal nodes
  depth_internal <- sapply(tree[names_internal_nodes], function(x) {
    x$depth_node
  })
  
  # Case for stump (No internal node)
  if (length(depth_internal) == 0) {
    log_p <- log(1 - alpha)
  } else {
    # Calculating the log-likelihood
    log_p <- sum(log(1 - alpha * (1 + depth_terminal)^(-beta))) + sum(log(alpha) - beta * log(1 + depth_internal))
  }
  return(log_p)
}

# Update tau_j values
update_tau <- function(x, 
                       y,
                       a_tau,
                       d_tau, 
                       predictions
) {
  
  # Calculating the values of a and d
  n <- nrow(x)
  
  # Getting the shape parameter from the posterior
  shape_tau_post <- 0.5*n+a_tau
  
  # Getting the ratio parameter
  rate_tau_post <- 0.5*crossprod( (y-predictions) ) + d_tau
  
  # Updating the \tau
  tau_sample <- rgamma(n = 1,shape = shape_tau_post,rate = rate_tau_post)
  
  return(tau_sample)
}




# Return rate parameter from the tau prior
rate_tau <- function(x, # X value
                     y, # Y value
                     prob = 0.9,
                     shape){
  # Find the tau_ols
  tau_ols <- naive_tau(x = x,
                       y = y)
  
  # Function to find the zero
  zero_tau_prob <- function(x,naive_tau_value, prob , shape){
    
    # Find the zero to the function P(tau < tau_ols ) = 0.1, for a defined   
    return (pgamma(naive_tau_value,
                   shape = shape,
                   rate = x)-(1-prob) )
    
  }
  
  # Function to find the zero
  zero_tau_prob_squared <- function(rate,naive_tau_value, prob , shape){
    
    # Find the zero to the function P(tau < tau_ols ) = 0.1, for a defined   
    return ((pgamma(naive_tau_value,
                    shape = shape,
                    rate = rate)-(1-prob))^2 )
    
  }
  
  
  # Getting the root
  min_root <-  try(uniroot(f = zero_tau_prob,interval = c(1e-2,100),
                           naive_tau_value = tau_ols,
                           prob = prob, shape = shape)$root,silent = TRUE)
  
  if( class(min_root) =="try-error"){
    # Verifying the squared version
    min_root <- optim(par = runif(1),fn = zero_tau_prob_squared,
                      method = "L-BFGS-B",lower = c(0),
                      naive_tau_value = tau_ols,
                      prob = prob, shape = shape)$par
    
  }
  
  
  return(min_root)
  
}


# Naive tau_estimation
naive_tau <- function(x,y){
  
  # Getting the valus from n and p
  n <- length(y)
  
  # Getting the value from p
  if(is.null(ncol(x))){
    p <- 1
  } else { 
    p <- ncol(x)
  }
  
  # Naive lm_mod 
  lm_mod <- lm(formula = y ~ ., data =  data.frame(y,x))
  
  sigma <- sqrt( sum((lm_mod$residuals)^2)/(n-p) )
  
  tau <- sigma^(-2)
  return(tau)
  
}


# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y){
  
  # Defining the a and b
  a <- min(y)
  b <- max(y)
  
  # This will normalize y between -0.5 and 0.5
  y  <- (y-a)/(b-a)-0.5
  
  return(y) 
}


# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  
  # Just getting back to the regular BART
  y <- (b-a)*(z+0.5) + a
  
  return(y)
}

# Function to calculate the tree complete conditional using BART model
tree_complete_conditional_bart <- function(x,
                                           residuals_values,
                                           tree,
                                           tau_mu,
                                           tau
){
  
  # Getting the number of observations of the data
  n <- nrow(x)
  
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
  
  # Retrieving Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals_values[x$observations_index]
  })
  
  
  # Defining RT_Omega_I_R
  RTR <- unlist(mapply(terminal_nodes,residuals_terminal_nodes, FUN = function(nodes,resid) {
    crossprod(resid)
  },SIMPLIFY = FALSE))
  
  
  # The term R^{T} solve(Omega + I ) 1 
  R_Omega_I_one <- unlist(mapply(terminal_nodes, residuals_terminal_nodes, FUN = function(nodes, residuals) {
    sum(residuals)
  },SIMPLIFY = FALSE))
  
  # Retrieve all nodes values and calculate all of them
  log_posterior <- 0.5*n*log(tau) + 0.5 * sum( log( tau_mu/(nodes_size*tau + tau_mu)) ) -
    0.5*tau*sum(RTR) + 0.5 * (tau^2) * sum((R_Omega_I_one^2)/(nodes_size*tau + tau_mu) )  
  
  
  
  return(log_posterior)
  
}


# Update \mu using the BART simple version
update_mu_bart <- function(tree,
                           x,
                           tau,
                           tau_mu,
                           residuals,
                           seed = NULL) {
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })
  
  
  # Residuals terminal nodes
  residuals_terminal_nodes_sum <- unlist(lapply(terminal_nodes, function(node) {
    sum(residuals[node$observations_index])
  }))
  
  # Remember that S = n_node*tau
  S <- (nodes_size*tau+tau_mu)
  
  
  # Mu mean value
  mu_mean <- (S^-1)*(tau)*residuals_terminal_nodes_sum
  
  
  # Mu mean SD
  mu_var <- (S^-1)
  
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



update_predictions_bart <- function(tree, x) {
  
  # New g (new vector prediction for g)
  predictions_new <- rep(NA, length(residuals))
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))
  
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(sapply(tree, function(x) {
    x$terminal == 1
  })))]
  
  
  # Getting the \mu_{j} vector
  mu_values <- unlist(
    lapply(terminal_nodes, function(x) {
      x$mu
    })
  )
  
  
  # Adding the mu values calculated
  for (i in 1:length(terminal_nodes)) {
    # Saving g
    predictions_new[terminal_nodes[[i]]$observations_index] <- mu_values[[i]]
    
  }
  
  
  return(predictions_new)
}


# Calculating the BART object model


bart <- function(x, # Covarariate matrix
                 y, # Target variable
                 number_trees = 2, # Number of trees
                 control = list(node_min_size = 5,
                                scale_boolean = TRUE,
                                rotation_boolean = FALSE),
                 prior_parameter = list(a_tau = 10, # Parameters from the prior
                                        d_tau = 3,
                                        k_bart = 2,
                                        alpha = 0.95,
                                        beta = 2,
                                        prob_tau = 0.9),
                 mcmc_parameter = list(n_iter = 2000, # Parameters from MCMC
                                       burn = 1000,
                                       thin = 1),
                 init = list(tau = 1, # Initial values
                             mu = 0)){
  
  # Setting the initial values for each one of them, 
  node_min_size <- control[["node_min_size"]]
  scale_boolean <- control[["scale_boolean"]]
  rotation_boolean <- control[["rotation_boolean"]]
  
  a_tau <- prior_parameter[["a_tau"]]
  d_tau <- prior_parameter[["d_tau"]]
  k_bart <- prior_parameter[["k_bart"]]
  alpha <- prior_parameter[["alpha"]]
  beta <- prior_parameter[["beta"]]
  prob_tau <- prior_parameter[["prob_tau"]]
  acc_ratio <- 0
  
  n_iter <- mcmc_parameter[["n_iter"]]
  burn <- mcmc_parameter[["burn"]]
  thin <- mcmc_parameter[["thin"]]
  
  tau <- init$tau
  mu <- init$mu
  
  # Giving names for col of "x" if they don't have it
  if(is.null(colnames(x))){
    colnames(x) <- paste0("x.",1:ncol(x))
  }
  
  # Create the storage for objects to be exported later
  store_size <- (n_iter-burn)/thin
  tree_store <- list()
  tau_store <- list()
  residuals_store <- list()
  y_hat_store <- matrix(NA, ncol = nrow(x),nrow = store_size)
  predictions_store <- list()
  
  # Creating a object to store the current trees
  current_trees <- list()
  
  # Creating the list of initial stumps and predictions
  for(j in 1:number_trees){
    current_trees[[j]] <- stump(x = x,
                                tau = tau,
                                mu = mu)
  }
  
  # Naming the trees to get the objects name
  names(current_trees) <- (sapply(1:number_trees, function(x) paste0("tree_", x))) # Naming each tree
  
  
  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )
  
  # Getting the maxium and the minimum from the data
  a_min <- min(y)
  b_max <- max(y)
  
  # Scaling the data or not (and chaning the hyperparameters with respect to it)
  if(scale_boolean){
    
    # Scaling the y between -0.5 and 0.5
    y_scale <- normalize_bart(y = y)
    
    # \tau_mu depends over the scale of y
    tau_mu <-   tau_mu <- 4*(k_bart^2)*number_trees/((max(y_scale)-min(y_scale))^2)
    
    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
    
  } else {
    
    y_scale <- y
    
    # \tau_mu depends over the scale of y
    tau_mu <-   tau_mu <- 4*(k_bart^2)*number_trees/((max(y_scale)-min(y_scale))^2)
    
    
    # Getting the optimal tau value
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
    
  }
  
  # Creating the init values of predictions
  predictions <- matrix(0,
                        ncol = nrow(x),
                        nrow = number_trees)
  
  current_partial_residuals_matrix <- matrix(0,
                                             ncol = nrow(x),
                                             nrow = number_trees)
  
  # Iterating over the MCMC sampling 
  for (i in 1:n_iter) {
    
    utils::setTxtProgressBar(progress_bar, i)
    
    
    
    # Saving after the warmup step
    if ((i > burn) & ((i %% thin) == 0)) {
      
      curr <- (i - burn) / thin
      tree_store[[curr]] <- current_trees
      tau_store[[curr]] <- tau
      y_hat_store[curr,] <- colSums(predictions)
      residuals_store[[curr]] <- current_partial_residuals_matrix
      predictions_store[[curr]] <- predictions
      
    }
    
    
    # Verb Store
    verb_store <- data.frame(verb = rep(NA,number_trees),
                             accepted = rep(NA,number_trees),
                             identical = rep(NA,number_trees))
    
    
    # Iterating over the trees 
    for(j in 1:number_trees){
      
      # Making a exception for the case of only one tree
      if( number_trees == 1){
        
        # Getting the current partial values
        current_partial_residuals <- matrix(y_scale,ncol = length(y_scale))
        
      } else {
        
        # Getting the current partial values
        current_partial_residuals <- y_scale - colSums(predictions[-j,,drop = FALSE])
      }
      
      # Propose a new tree based on the verbs: grow/prune/change/swap
      if(rotation_boolean){
        verb <- sample(c("grow","grow_projection", "prune", "change","change_projection","swap"),
                       prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
      } else{
        verb <- sample(c("grow", "prune", "change","swap"),
                       prob = c(0.25,0.25,0.4,0.1), size = 1)
      }
      
      
      
      # Case of rotation
      if(rotation_boolean){
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
        verb = verb
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
      
      # In case of acceptance
      if (runif(1) < acceptance) {
        
        # Counting acc ratio
        acc_ratio <- acc_ratio + 1
        
        # Make changes if accept
        current_trees <- new_trees
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store[j,"verb"] <- verb
        verb_store[j,"accepted"] <- TRUE
        
        
        
      } else {
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store[j,"verb"] <- verb
        verb_store[j,"accepted"] <- FALSE
        
      } # End of accept for MH sample
      
      # # # To update the mu values
      current_trees[[j]] <- update_mu_bart(
        tree = current_trees[[j]],
        x = x,
        residuals = current_partial_residuals, 
        tau = tau,
        tau_mu = tau_mu)
      
      
      # EQUATION FROM SECTION 4
      # ==== Using the prediction from R_star_bar
      predictions[j,] <- update_predictions_bart(
        tree = current_trees[[j]], x = x
      )
      
      current_partial_residuals_matrix[j,] <- current_partial_residuals
      
    } # End of iterations over trees
    
    
    # Calling the function to update tau
    tau <- update_tau(x = x,
                      y = y_scale,
                      a_tau = a_tau,
                      d_tau = d_tau,
                      predictions = colSums(predictions))
    
  } # End of the iteration over MCMC sampling
  
  
  return(list(x = x, # Covarariate matrix
              y = y_scale, # Target variable
              tree_store = tree_store, # Tree store
              residuals_store = residuals_store,
              predictions_store = predictions_store,
              tau_store = tau_store,
              y_hat_store = y_hat_store,
              number_trees = number_trees, # Number of trees
              control = list(node_min_size = node_min_size,
                             scale_boolean = scale_boolean,
                             a_min = a_min,
                             b_max = b_max,
                             rotation_boolean = rotation_boolean),
              prior_parameter = list(a_tau = a_tau, # Parameters from the prior
                                     d_tau = d_tau,
                                     k_bart = k_bart,
                                     alpha = alpha,
                                     beta = beta,
                                     prob_tau = prob_tau),
              mcmc_parameter = list(n_iter = n_iter, # Parameters from MCMC
                                    burn = burn,
                                    thin = thin),
              verb_store = verb_store)
  )
  
  
}

# Getting the predictions for each tree
get_predictions_tree <- function(tree,
                                 x){
  
  new_tree <- tree
  pred_new <- numeric(nrow(x))
  # Setting the root node with the new observations
  new_tree[["node_0"]]$test_index <- 1:nrow(x)
  
  # Creating the list of nodes
  list_nodes <- names(new_tree)[-1]
  
  # IN CASE OF JUST ROOT NODE
  if (length(new_tree) == 1) {
    list_nodes <- "node_0"
  }
  
  # Iterating over the list of nodes
  for(i in 1:length(list_nodes)){
    
    # Selecting the current node
    current_node_aux <- new_tree[[list_nodes[[i]]]]
    
    
    # Veryfing the type of the current node
    if (is.list(current_node_aux$node_var)){
      
      
      # Rotation function
      A <- function(theta) {
        matrix(c(
          cos(theta), -sin(theta),
          sin(theta), cos(theta)
        ),
        ncol = 2, nrow = 2, byrow = TRUE
        )
      }
      
      # Rotated Lon and Lat
      rotated_x <- tcrossprod(A(current_node_aux$theta), x[, current_node_aux$node_var$node_var_pair])
      rownames(rotated_x) <- current_node_aux$node_var$node_var_pair
      
      # Updating observations from the left node
      if (current_node_aux$left == 1) {
        new_tree[[list_nodes[i]]]$test_index<- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] < current_node_aux$node_var_split)] # Updating the left node
      } else {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] >= current_node_aux$node_var_split)]
      }
      
    } else { # Checking the case where there is no rotated variable
      
      # To continous covariates
      if (is.numeric(x[, current_node_aux$node_var])) {
        
        # Updating observations from the left node
        if (current_node_aux$left == 1) {
          new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] < current_node_aux$node_var_split)] # Updating the left node
        } else {
          new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
        }
        
        
        # To categorical covariates
      } else {
        # Updating observations from the left node
        if (current_node_aux$left == 1) {
          new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
        } else {
          new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
        }
      } # Ending the else for type of variable
      
      # Checking if it is a terminal node or not
      if (new_tree[[list_nodes[[i]]]]$terminal == 1 & length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {
        
        # Assign the mu value for that terminal node
        pred_new[new_tree[[list_nodes[[i]]]]$test_index] <- new_tree[[list_nodes[[i]]]]$mu
        
      } # End of check terminal
    } # End of the rotation check
  }# End Terminal node iterations
  
  return(pred_new)
  
  
}


# Predicting for new observation
predict_bart <- function(bart_mod, newdata, type = c("all")){
  
  # Getting the MCMC prediction values
  mcmc_post_pred <- matrix(NA,
                           nrow = length(bart_mod$tree_store),
                           ncol = nrow(newdata))
  
  # Setting up the colnames
  colnames(newdata) <- colnames(bart_mod$x)
  
  # Retrieving all information necessary to predict for new observations
  for(i in 1:length(bart_mod$tree_store)){
    
    current_trees <- bart_mod$tree_store[[i]]
    
    # Creating a pred_aux matrix to store the predictions for new obs
    pred_aux <- matrix(0,
                       nrow = bart_mod$number_trees,
                       ncol = nrow(newdata))
    
    for(j in 1:bart_mod$number_trees){
      
      pred_aux[j,] <- get_predictions_tree(tree = current_trees[[j]],x = newdata)
      
    }
    
    # Adding up all trees
    mcmc_post_pred[i,] <- colSums(pred_aux)
  }
  
  # Getting the values
  out <- switch(type,
                all = unnormalize_bart(z = mcmc_post_pred, a = bart_mod$control$a_min, b = bart_mod$control$b_max),
                mean = unnormalize_bart(z = colMeans(mcmc_post_pred), a = bart_mod$control$a_min, b = bart_mod$control$b_max),
                median = unnormalize_bart(z = colMeans(mcmc_post_pred), a = bart_mod$control$a_min, b = bart_mod$control$b_max))
  
  return(out)
}




