# GP-Bart

# Loading test data
library(magrittr)
library(Rcpp)
source("gpbart/seed_bart/seed_bart_simple.R")
source("gpbart/fast_gp_multiple_tau.R")

error_counter_rmvn <<- 0

# set.seed(42)
# Doing quantile normalization and linear transformation to uniform distribution
normalizing_and_scaling_x <- function(data) {
  
  # Creating index to keep the initial order
  data <- cbind(matrix(1:nrow(data),
                       nrow = nrow(data)), data)
  
  # Ordering the data
  data <- data[order(data[, 1]), , drop = FALSE]
  
  # Saving only the covariates
  x <- data[, -ncol(data), drop = FALSE]
  
  # Saving the target variable
  y <- data[, ncol(data), drop = FALSE]
  
  # Getting ordering from the minimum to  the maximum
  data_rank <- apply(x, 2, rank, ties.method = "min")
  
  data_sort <- apply(x, 2, sort)
  data_mean <- apply(x, 1, mean)
  
  index_to_mean <- function(my_index, my_mean) {
    return(my_mean[my_index])
  }
  
  quantile_normalized_data <- apply(data_rank, 2, index_to_mean, data_mean)
  
  # Putting in a uniform interval
  scaled_and_normalized <- apply(quantile_normalized_data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Reordering data
  new_data <- cbind(scaled_and_normalized[order((scaled_and_normalized[, 1])), , drop = FALSE], y)
  
  return(new_data[, -1, drop = FALSE])
}

normalizing_and_scaling_x_test <- function(x) {
  
  # Putting the reference values
  x <- cbind(matrix(1:nrow(x), nrow = nrow(x)), x)
  x <- x[order(x[, 1]), , drop = FALSE]
  
  # Scaling
  data_rank <- apply(x, 2, rank, ties.method = "min")
  
  data_sort <- apply(x, 2, sort)
  data_mean <- rowMeans(x)
  
  index_to_mean <- function(my_index, my_mean) {
    return(my_mean[my_index])
  }
  
  quantile_normalized_data <- apply(data_rank, 2, index_to_mean, data_mean)
  
  # Putting in a uniform interval
  scaled_and_normalized <- apply(quantile_normalized_data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Reordering data
  x_scaled <- scaled_and_normalized[order((scaled_and_normalized[, 1])), , drop = FALSE]
  
  return(x_scaled[, -1, drop = FALSE])
}


# Scaling the data
unit_cube_scale <- function(x) {
  
  # Scaling to -1 to 1 function
  scaling <- function(z) {
    (2 * z - (max(z) + min(z))) / (max(z) - min(z))
  }
  
  # Applying on all covariates
  x_scaled <- apply(x, 2, scaling)
  return(x_scaled)
}


# ==================================#
# Objects to test the stump function
# ==================================#

# Create an stump
stump <- function(x, tau, mu) {
  node <- list()
  node[["node_0"]] <- list(
    node_number = 0, observations_index = 1:nrow(x),
    depth_node = 0,
    node_var = NA,
    node_var_split = NA,
    left = 0,
    right = 0,
    parent_node = NA, # Retrieving the number of parent node
    terminal = 1,
    tau = tau
  )
  node[["node_0"]]$mu <- mu
  return(node)
}

# ==================================#
# Objects to test the grow_tree_verb function
# ==================================#
# tree<-stump(x,tau)
# node_min_size<-5

grow_projection_tree <- function(tree, x, node_min_size, theta = NULL) {
  
  # Controlling the "bad trees"
  bad_trees <- TRUE
  count_bad_trees <- 0
  
  # Try to not get bad_trees
  while (bad_trees) {
    
    # Creating a new tree object to be grown
    new_tree <- tree
    
    # Returning the Tree if there's no node to grow
    if (sum(unlist(lapply(new_tree, function(x) {
      # This condition will see if the terminal nodes have at least more than 2* the min node size
      x$terminal == 1 & length(x$observations_index) > 2 * node_min_size
    }))) == 0) {
      return(tree)
    }
    
    # Selecting the terminal nodes
    nodes_to_grow_names <- names(new_tree[unlist(lapply(new_tree, function(x) {
      (x$terminal == 1) & (length(x$observations_index) > 2 * node_min_size)
    }))])
    
    
    # Selecting one terminal node randomly 
    node_to_grow <- sample(nodes_to_grow_names,size = 1)
    
    
    # Selecting the current node to grow
    current_node <- new_tree[[node_to_grow]]
    
    # Selecting the covariate in case of projection
    node_var_pair <- sample(colnames(x)[apply(x,2,is.numeric)], 2) # selecting the covariate WITH NO ROTATION
    
    # Selecting the node_var that it will be sample
    node_var <- sample(node_var_pair,size = 1)
    
    
    # Creating a node_var list for the projection
    node_var_list <- list(node_var_pair = node_var_pair, node_var = node_var)
    
    # Determining that this node is no longer a terminal one
    new_tree[[node_to_grow]]$terminal <- 0
    
    # Selecting randomly (OR NOT) the angle to be rotated
    if (is.null(theta)) {
      # Selecting a uniform theta
      theta<-runif(n = 1,min = 0 ,max = pi)
      
      # Selecting a grid in theta values
      # theta <- sample(x = seq(0, pi, length.out = 25)[-1], size = 1)
    }
    
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
    rotated_x <- tcrossprod(A(theta), x[, node_var_pair])
    rownames(rotated_x) <- node_var_pair
    
    # Selecting the rotated var split
    node_var_split <- sample(sort(rotated_x[node_var, ]), size = 1)
    
    
    # Creating the left observations
    left_node_split <- current_node$observations_index[which(rotated_x[node_var, current_node$observations_index] < node_var_split)] # Selecting the left node
    
    # Creating the right observations
    right_node_split <- current_node$observations_index[which(rotated_x[node_var, current_node$observations_index] >= node_var_split)] # Selecting the right node
    
    # Updating the depth node
    depth_node <- current_node$depth_node + 1
    
    # Adding to the next node
    i <- max(unlist(lapply(new_tree, function(x) {
      x$node_number
    }))) + 1
    
    # Informing the next node
    new_tree[[paste0("node_", i)]] <- list(
      node_number = i, observations_index = left_node_split,
      depth_node = depth_node,
      node_var = node_var_list,
      node_var_split = node_var_split,
      theta = theta,
      left = 1,
      right = 0,
      parent_node = current_node$node_number, # Retrieving the number of parent node
      terminal = 1,
      mu = current_node$mu,
      tau = current_node$tau
    ) # Adding the node to the new_tree
    
    i <- i + 1 # Adding to the next node
    
    new_tree[[paste0("node_", i)]] <- list(
      node_number = i, observations_index = right_node_split,
      depth_node = depth_node,
      node_var = node_var_list,
      node_var_split = node_var_split,
      theta = theta,
      left = 0,
      right = 1,
      parent_node = current_node$node_number, # Retrieving the number of parent node
      terminal = 1,
      mu = current_node$mu,
      tau = current_node$tau
    ) # Adding the node to the new_tree
    
    
    # Verifying if it is a good or bad tree
    if (any(unlist(lapply(new_tree, function(x) {
      (length(x$observations_index) < node_min_size)
    })))) { # Veryfing if any terminal node is lower than node_min_size
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE # Return the new good tree
    }
    
    # Verfying the limit of counting 2 bad trees
    if (count_bad_trees == 2) {
      return(tree) # Return the original tree
    }
  }
  
  return(new_tree)
}


# ==================================#
# Objects to test the grow__projection_tree_verb function
# ==================================#

grow_tree <- function(tree, x, node_min_size, rotation =  TRUE, theta = NULL){
  
  
  # Controlling the "bad trees"
  bad_trees <- TRUE
  count_bad_trees <- 0
  
  # Try to not get bad_trees
  while (bad_trees) {
    
    # Creating a new tree object to be grown
    new_tree <- tree
    
    # Returning the Tree if there's no node to grow
    if (sum(unlist(lapply(new_tree, function(x) {
      # This condition will see if the terminal nodes have at least more than 2* the min node size
      x$terminal == 1 & length(x$observations_index) > 2 * node_min_size
    }))) == 0) {
      return(tree)
    }
    
    # Selecting the terminal nodes
    nodes_to_grow_names <- names(new_tree[unlist(lapply(new_tree, function(x) {
      (x$terminal == 1) & (length(x$observations_index) > 2 * node_min_size)
    }))])
    
    
    # Selecting one terminal node randomly 
    node_to_grow <- sample(nodes_to_grow_names,size = 1)
    
    
    # Selecting the current node to grow
    current_node <- new_tree[[node_to_grow]]
    
    # Select only the covariates without the rotation
    node_var <- sample(colnames(x), 1) 
    
    
    # Determining that this node is no longer a terminal one
    new_tree[[node_to_grow]]$terminal <- 0
    
    
    # Conditioning if the covariate is continuous or categorical
    if (is.numeric(x[, node_var])) {
      
      node_var_split <- sample(x[current_node$observations_index, node_var], size = 1) # Selecting the splitting point on the variable
      # node_var_split <- sample(c(0.3,0.5,0.7),size = 1)
      # node_var_split <- sample(c(-0.4,0,0.4),size = 1)
      
      left_node_split <- current_node$observations_index[which(x[current_node$observations_index, node_var] < node_var_split)] # Selecting the left node
      right_node_split <- current_node$observations_index[which(x[current_node$observations_index, node_var] >= node_var_split)] # Selecting the right node
      
      # Updating the depth node
      depth_node <- current_node$depth_node + 1
      
      # Adding to the next node
      i <- max(unlist(lapply(new_tree, function(x) {
        x$node_number
      }))) + 1
      
      
      # Informing the next node
      new_tree[[paste0("node_", i)]] <- list(
        node_number = i, observations_index = left_node_split,
        depth_node = depth_node,
        node_var = node_var,
        node_var_split = node_var_split,
        left = 1,
        right = 0,
        parent_node = current_node$node_number, # Retrieving the number of parent node
        terminal = 1,
        mu = current_node$mu,
        tau = current_node$tau
      ) # Adding the node to the new_tree
      
      i <- i + 1 # Adding to the next node
      
      new_tree[[paste0("node_", i)]] <- list(
        node_number = i, observations_index = right_node_split,
        depth_node = depth_node,
        node_var = node_var,
        node_var_split = node_var_split,
        left = 0,
        right = 1,
        parent_node = current_node$node_number, # Retrieving the number of parent node
        terminal = 1,
        mu = current_node$mu,
        tau = current_node$tau
      ) # Adding the node to the new_tree
      
    } else { # Condition to categorical variables
      
      node_var_split <- sample(
        levels(x[current_node$observations_index, node_var]), # Selecting from the observations
        1
      ) # One sample to split
      
      left_node_split <- current_node$observations_index[which(x[current_node$observations_index, node_var] == node_var_split)] # Selecting the left node
      right_node_split <- current_node$observations_index[which(x[current_node$observations_index, node_var] != node_var_split)] # Selecting the right node
      
      # Updating depthnode
      depth_node <- current_node$depth_node + 1
      
      i <- max(unlist(lapply(new_tree, function(x) {
        x$node_number
      }))) + 1 # Adding to the next node
      
      # Informing the next node
      new_tree[[paste0("node_", i)]] <- list(
        node_number = i, observations_index = left_node_split,
        depth_node = depth_node,
        node_var = node_var,
        node_var_split = node_var_split,
        left = 1,
        right = 0,
        parent_node = current_node$node_number,
        terminal = 1,
        mu = current_node$mu,
        tau = current_node$tau
      ) # Adding the node to the new_tree
      
      i <- i + 1 # Adding to the next node
      
      new_tree[[paste0("node_", i)]] <- list(
        node_number = i, observations_index = right_node_split,
        depth_node = depth_node,
        node_var = node_var,
        node_var_split = node_var_split,
        left = 0,
        right = 1,
        parent_node = current_node$node_number,
        terminal = 1,
        mu = current_node$mu,
        tau = current_node$tau
      ) # Adding the node to the new_tree
    }
    
    # Verifying if it is a good or bad tree
    if (any(unlist(lapply(new_tree, function(x) {
      (length(x$observations_index) < node_min_size)
    })))) { # Veryfing if any terminal node is lower than node_min_size
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE # Return the new good tree
    }
    
    # Verfying the limit of counting 2 bad trees
    if (count_bad_trees == 2) {
      return(tree) # Return the original tree
    }
  }
  
  return(new_tree)
  
}

# ==================================#
# Objects to test the prune_tree_verb function
# ==================================#

# Prune a tree verb
prune_tree_verb <- function(tree, x) {
  
  # Returning if the tree is just a single node
  if(length(tree)==1){
    return(tree)
  }
  
  # Selecting parent nodes
  parent_nodes <- unique(unlist(lapply(tree, function(x) {
    x$parent_node
  })))
  parent_nodes <- parent_nodes[!is.na(parent_nodes) & parent_nodes >= 0]
  
  # Choose the pairs nodes
  pairs_terminal_nodes <- t(sapply(parent_nodes, function(x) {
    lapply(tree, function(y) {
      y$parent_node == x
    })
  }))
  
  # NOmes of terminals nodes
  names_pairs_terminal_nodes <- apply(pairs_terminal_nodes, 1, function(x) {
    names(which(unlist(x)))
  })
  
  # Getting BOOLEAN of the pairs of ONLY terminal nodes
  terminal_both_nodes <- unlist(lapply(apply(names_pairs_terminal_nodes, 2, function(x) {
    lapply(tree[x], function(z) {
      z$terminal == 1
    })
  }), function(w) {
    all(unlist(w))
  }))
  
  # Selecting the name of children terminal nodes.
  keep_only_terminal_nodes <- names_pairs_terminal_nodes[, terminal_both_nodes]
  
  # Children pair (if to case of just one to be pruned)
  if (is.matrix(keep_only_terminal_nodes)) {
    child_to_be_pruned <- keep_only_terminal_nodes[, sample(1:ncol(keep_only_terminal_nodes), size = 1)]
  } else {
    child_to_be_pruned <- keep_only_terminal_nodes
  }
  # Getting the parent node
  parent_of_the_pruned_node <- tree[[child_to_be_pruned[1]]]$parent_node
  
  tree[[child_to_be_pruned[1]]] <- NULL # Pruning child node 1
  
  tree[[child_to_be_pruned[2]]] <- NULL # Pruning child node 2
  
  # Transforming the parent in a terminal node
  tree[[paste0("node_", parent_of_the_pruned_node)]]$terminal <- 1
  
  
  return(tree)
}

# tree <- prune_tree_verb(tree = tree)

# ==================================#
# Objects to test the change_tree_verb function
# ==================================#
# tree<-grow_tree(tree = stump(x = x,tau = 1),x = x,node_min_size = 5,rotation = FALSE) %>%
# grow_tree(x=x,node_min_size = 5,rotation = FALSE)


# Change a tree verb
change_tree_verb <- function(tree, x, node_min_size, rotation = FALSE, theta = NULL) {
  
  # Controlling the "bad trees"
  bad_trees <- TRUE
  count_bad_trees <- 0
  
  while (bad_trees) {
    
    # Craeting the dummy for the new tree
    new_tree <- tree
    
    # Randomly select a internal node
    nodes_to_change_names <- names(new_tree[unlist(lapply(new_tree, function(x) {
      x$terminal == 0
    }))])
    
    # Sampling the node to change
    node_to_change <- sample(nodes_to_change_names,size = 1)
    
    # setting the node to be changed as default
    current_node <- new_tree[[node_to_change]]
    
    # Selecting the covariate
    node_var <- sample(colnames(x), 1)
    
    # Analyzing the case where isn't a rotated cov.
    if (is.numeric(x[, node_var])) {
      node_var_split <- sample(sort(x[current_node$observations_index, node_var]), size = 1) # Selecting the splitting point on the variable
      # node_var_split <- sample(c(0.3,0.5,0.7),size = 1)
      # node_var_split <- sample(c(-0.4,0,0.4),size = 1)
      
      
    } else {
      node_var_split <- sample(
        levels(x[current_node$observations_index, node_var]), # Selecting from the observations
        1
      )
    }
    
    
    # Function to get all children nodes from that which was changed.
    get_all_children <- function(new_tree, current_node) {
      aux <- names(which(unlist(lapply(new_tree, function(x) {
        x$parent_node == current_node$node_number
      })))) # Selecting the children nodes.
      
      if (new_tree[[aux[1]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[1]]]))
      }
      if (new_tree[[aux[2]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[2]]]))
      }
      return(aux)
    }
    
    # Get childrens
    children_from_change_node <- get_all_children(new_tree, current_node = current_node)
    
    
    new_tree[[children_from_change_node[1]]]$node_var <- node_var # Changing the node var
    new_tree[[children_from_change_node[2]]]$node_var <- node_var # Changing the node var
    
    new_tree[[children_from_change_node[1]]]$node_var_split <- node_var_split # Changing the node split
    new_tree[[children_from_change_node[2]]]$node_var_split <- node_var_split # Changing the node split
    
    
    
    # Updating all nodes
    for (i in 1:length(children_from_change_node)) {
      
      # Iterating over each one of the tre noees
      current_node_aux <- new_tree[[children_from_change_node[i]]]
      
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
          new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] < current_node_aux$node_var_split)] # Updating the left node
        } else {
          new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] >= current_node_aux$node_var_split)]
        }
        
      } else { # Checking the case where there is no rotated variable
        
        # To continous covariates
        if (is.numeric(x[, current_node_aux$node_var])) {
          
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] >= current_node_aux$node_var_split)] # Updating the right node
          }
          
          
          # To categorical covariates
        } else {
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
          }
        }
      }
    }
    
    
    # Verifying if it is a good or bad tree
    if (any(unlist(lapply(new_tree, function(x) {
      (length(x$observations_index) < node_min_size)
    })))) { # Veryfing if any terminal node is lower than node_min_size
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE # Return the new good tree
    }
    
    # Verfying the limit of counting 2 bad trees
    if (count_bad_trees == 2) {
      return(tree) # Return the original tree
    }
    
  } # Stopping the while
  
  return(new_tree)
}

# ==================================#
# Objects to test the change_tree_verb function
# ==================================#
# tree<-grow_tree(tree = stump(x = x,tau = 1),x = x,node_min_size = 5,rotation = FALSE) %>%
# grow_tree(x=x,node_min_size = 5,rotation = FALSE)


# Change a tree verb
change_projection_tree_verb <- function(tree, x, node_min_size, theta = NULL) {
  
  # Controlling the "bad trees"
  bad_trees <- TRUE
  count_bad_trees <- 0
  
  while (bad_trees) {
    
    # Craeting the dummy for the new tree
    new_tree <- tree
    
    # Randomly select a internal node
    nodes_to_change_names <- names(new_tree[unlist(lapply(new_tree, function(x) {
      x$terminal == 0
    }))])
    
    # Sampling the node to change
    node_to_change <- sample(nodes_to_change_names,size = 1)
    
    # setting the node to be changed as default
    current_node <- new_tree[[node_to_change]]
    
    # Selecting the covariate
    node_var_pair <- sample(colnames(x)[apply(x,2,is.numeric)], 2)
    
    node_var <- sample(node_var_pair,size = 1)
    
    
    # Selecting randomly (OR NOT) the angle to be rotated
    if (is.null(theta)) {
      # Selecting a uniform theta
      theta<-runif(n = 1,min = 0 ,max = pi)
      
      # Selecting a grid in theta values
      # theta <- sample(x = seq(0, pi, length.out = 25)[-1], size = 1)
    }
    
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
    rotated_x <- tcrossprod(A(theta), x[, node_var_pair])
    rownames(rotated_x) <- node_var_pair
    
    # Selecting the rotated var split
    node_var_split <- sample(sort(rotated_x[node_var, ]), size = 1)
    
    # Function to get all children nodes from that which was changed.
    get_all_children <- function(new_tree, current_node) {
      aux <- names(which(unlist(lapply(new_tree, function(x) {
        x$parent_node == current_node$node_number
      })))) # Selecting the children nodes.
      
      if (new_tree[[aux[1]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[1]]]))
      }
      if (new_tree[[aux[2]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[2]]]))
      }
      return(aux)
    }
    
    # Get childrens
    children_from_change_node <- get_all_children(new_tree, current_node = current_node)
    
    # Create a node_var list
    node_var_list <- list(node_var_pair = node_var_pair, node_var = node_var)
    
    new_tree[[children_from_change_node[1]]]$node_var <- node_var_list # Changing the node var
    new_tree[[children_from_change_node[2]]]$node_var <- node_var_list # Changing the node var
    
    new_tree[[children_from_change_node[1]]]$node_var_split <- node_var_split # Changing the node split
    new_tree[[children_from_change_node[2]]]$node_var_split <- node_var_split # Changing the node split
    
    # For a rotated tree, adding the theta parameter
    new_tree[[children_from_change_node[1]]]$theta <- theta
    new_tree[[children_from_change_node[2]]]$theta <- theta
    
    
    # Updating all nodes
    for (i in 1:length(children_from_change_node)) {
      
      # Iterating over each one of the tre noees
      current_node_aux <- new_tree[[children_from_change_node[i]]]
      
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
          new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] < current_node_aux$node_var_split)] # Updating the left node
        } else {
          new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] >= current_node_aux$node_var_split)]
        }
        
      } else { # Checking the case where there is no rotated variable
        
        # To continous covariates
        if (is.numeric(x[, current_node_aux$node_var])) {
          
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] >= current_node_aux$node_var_split)] # Updating the right node
          }
          
          
          # To categorical covariates
        } else {
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[children_from_change_node[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
          }
        }
        
        
      }
    }
    
    
    # Verifying if it is a good or bad tree
    if (any(unlist(lapply(new_tree, function(x) {
      (length(x$observations_index) < node_min_size)
    })))) { # Veryfing if any terminal node is lower than node_min_size
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE # Return the new good tree
    }
    
    
    # Verfying the limit of counting 2 bad trees
    if (count_bad_trees == 2) {
      return(tree) # Return the original tree
    }
    
  } # Counting the tree while
  
  
  return(new_tree)
}


# ==================================#
# Objects to test the swap_tree_verb function
# ==================================#
# tree<-grow_tree(tree = stump(x = x,tau = 1),x = x,node_min_size = 5,rotation = FALSE) %>%
#       grow_tree(x=x,node_min_size = 5,rotation = FALSE) %>%
#       grow_tree(x=x,node_min_size = 5, rotation = FALSE) %>%
#       grow_tree(x=x,node_min_size = 5, rotation = FALSE)


# Testing the tree changes updates
# for(i in 1:10000){
# tree <- change_tree_verb(tree = tree,x = x,node_min_size = 15,rotation = FALSE)
# }


# Swap a tree verb
swap_tree_verb <- function(tree, x, node_min_size) {
  
  # Auxiliary variables to trim and make sure that I getting the true values
  bad_trees <- TRUE
  count_bad_trees <- 0
  
  # Condition to not choose terminal nodes with less than node_min_size
  while (bad_trees) {
    
    # Updating the new tree
    new_tree <- tree
    
    # If less than 3 internal nodes, return the tree
    count_internal <- sum(unlist(lapply(new_tree, function(x) {
      x$terminal == 0
    })))
    
    if (count_internal <= 3) {
      return(tree)
    }
    
    # Selecting internal nodes names
    names_internal_nodes <- names(which(sapply(tree, function(x) {
      x$terminal == 0
    })))
    
    # )) # Randomly select a internal node
    node_to_swap_child_internal_CANDIDATES <- intersect(names(new_tree)[! (names(new_tree) %in% c("node_0","node_1","node_2")) ],
                                                        names_internal_nodes)
    
    # If there is no internal child node
    if(length(node_to_swap_child_internal_CANDIDATES)==0){
      return(tree)
    }
    # Certifying that's a parent and child internal
    
    node_to_swap_child_internal <- sample(node_to_swap_child_internal_CANDIDATES,size = 1)
    
    # Certifying that's a parent and child internal
    # Saving as current node
    current_node <- new_tree[[node_to_swap_child_internal]]
    
    # Selecting the node to swap
    node_to_swap_parent <- names(which(unlist(lapply(new_tree, function(x) {
      x$node_number == current_node$parent_node
    })))) # Selecting the children from swap parent
    
    # Getting all children function
    get_all_children <- function(new_tree, current_node) {
      
      # Selecting the children nodes.
      aux <- names(which(unlist(lapply(new_tree, function(x) {
        x$parent_node == current_node$node_number
      }))))
      
      if (new_tree[[aux[1]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[1]]]))
      }
      if (new_tree[[aux[2]]]$terminal == 0) {
        aux <- c(aux, get_all_children(new_tree, current_node = new_tree[[aux[2]]]))
      }
      return(aux)
    }
    
    # Getting the both next children from the swapped parent node.
    parent_swap_imediate_children <- get_all_children(new_tree, current_node = new_tree[[node_to_swap_parent]])[1:2]
    
    # Getting the both next children from the swapped children node.
    child_swap_imediate_children <- get_all_children(new_tree, current_node = new_tree[[node_to_swap_child_internal]])[1:2]
    
    # Selecting the node_var and node_var_split rule for parent and child internal respectively
    #
    # For the case where a rot. lat/lon is swap
    if ( is.list(new_tree[[parent_swap_imediate_children[1]]]$node_var) ) {
      node_and_split_1 <- list(
        new_tree[[parent_swap_imediate_children[1]]]$node_var, new_tree[[parent_swap_imediate_children[1]]]$node_var_split,
        new_tree[[parent_swap_imediate_children[1]]]$theta
      )
    } else {
      node_and_split_1 <- c(new_tree[[parent_swap_imediate_children[1]]]$node_var, new_tree[[parent_swap_imediate_children[1]]]$node_var_split)
    }
    # For the case where a rot. lat/lon is swaped
    if (is.list(new_tree[[child_swap_imediate_children[1]]]$node_var) ) {
      node_and_split_2 <- list(
        new_tree[[child_swap_imediate_children[1]]]$node_var, new_tree[[child_swap_imediate_children[1]]]$node_var_split,
        new_tree[[child_swap_imediate_children[1]]]$theta
      )
    } else {
      node_and_split_2 <- c(new_tree[[child_swap_imediate_children[1]]]$node_var, new_tree[[child_swap_imediate_children[1]]]$node_var_split)
    }
    
    
    
    # Chaging the first swapped internal child from parent (node one)
    if (is.list(node_and_split_2[[1]])) { # For the case where a rot. lat/lon is swap
      new_tree[[parent_swap_imediate_children[1]]]$node_var <- node_and_split_2[[1]]
      new_tree[[parent_swap_imediate_children[1]]]$node_var_split <- as.numeric(node_and_split_2[[2]])
      new_tree[[parent_swap_imediate_children[1]]]$theta <- as.numeric(node_and_split_2[[3]])
      
    } else {
      new_tree[[parent_swap_imediate_children[1]]]$node_var <- node_and_split_2[1]
      new_tree[[parent_swap_imediate_children[1]]]$node_var_split <- as.numeric(node_and_split_2[2])
    }
    
    
    # Chaging the first swapped internal child from parent (node two)
    if (is.list(node_and_split_2[[1]]) ) { # For the case where a rot. lat/lon is swap
      new_tree[[parent_swap_imediate_children[2]]]$node_var <- node_and_split_2[[1]]
      new_tree[[parent_swap_imediate_children[2]]]$node_var_split <- as.numeric(node_and_split_2[[2]])
      new_tree[[parent_swap_imediate_children[2]]]$theta <- as.numeric(node_and_split_2[[3]])
      
    } else {
      new_tree[[parent_swap_imediate_children[2]]]$node_var <- node_and_split_2[1]
      new_tree[[parent_swap_imediate_children[2]]]$node_var_split <- as.numeric(node_and_split_2[2])
    }
    
    # Chaging the first swapped internal parent from child (node one)
    if (is.list(node_and_split_1[[1]]) ) { # For the case where a rot. lat/lon is swap
      new_tree[[child_swap_imediate_children[1]]]$node_var <- node_and_split_1[[1]]
      new_tree[[child_swap_imediate_children[1]]]$node_var_split <- as.numeric(node_and_split_1[[2]])
      new_tree[[child_swap_imediate_children[1]]]$theta <- as.numeric(node_and_split_1[[3]])
      
    } else {
      new_tree[[child_swap_imediate_children[1]]]$node_var <- node_and_split_1[1]
      new_tree[[child_swap_imediate_children[1]]]$node_var_split <- as.numeric(node_and_split_1[2])
    }
    
    
    # Chaging the first swapped internal parent from child (node two)
    if (is.list(node_and_split_1[[1]]) ) { # For the case where a rot. lat/lon is swap
      new_tree[[child_swap_imediate_children[2]]]$node_var <- node_and_split_1[[1]]
      new_tree[[child_swap_imediate_children[2]]]$node_var_split <- as.numeric(node_and_split_1[[2]])
      new_tree[[child_swap_imediate_children[2]]]$theta <- as.numeric(node_and_split_1[[3]])
      
    } else {
      new_tree[[child_swap_imediate_children[2]]]$node_var <- node_and_split_1[1]
      new_tree[[child_swap_imediate_children[2]]]$node_var_split <- as.numeric(node_and_split_1[2])
    }
    
    # Updating the new_tree observations (Getting all children from that swapped)
    list_nodes <- get_all_children(new_tree = new_tree, current_node = new_tree[[node_to_swap_parent]])
    
    # Updating all nodes
    for (i in 1:length(list_nodes)) {
      current_node_aux <- new_tree[[list_nodes[i]]]
      
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
          new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] < current_node_aux$node_var_split)] # Updating the left node
        } else {
          new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index] >= current_node_aux$node_var_split)]
        }
        
      } else { # Checking the case where there is no rotated variable
        
        # To continous covariates
        if (is.numeric(x[, current_node_aux$node_var])) {
          
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] >= current_node_aux$node_var_split)] # Updating the right node
          }
          
          
          # To categorical covariates
        } else {
          # Updating observations from the left node
          if (current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$observations_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$observations_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
          }
        }
      }
    }
    
    # Verifying if it is a good or bad tree
    if (any(unlist(lapply(new_tree, function(x) {
      (length(x$observations_index) < node_min_size)
    })))) { # Veryfing if any terminal node is lower than node_min_size
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE # Return the new good tree
    }
    
    # Verfying the limit of counting 2 bad trees
    if (count_bad_trees == 2) {
      return(tree) # Return the original tree
    }
  }
  
  return(new_tree)
}



# Sort an action to update the tree
update_tree_verb <- function(tree, x, node_min_size, verb, rotation = TRUE, theta = NULL) {
  
  # Update the tree by a verb
  updated_tree <- switch(verb,
                         grow = grow_tree(tree,
                                          x = x, node_min_size = node_min_size, rotation = rotation,
                                          theta = theta
                         ), # Grow verb
                         grow_projection = grow_projection_tree(tree = tree,
                                                                x = x,
                                                                node_min_size = node_min_size,theta = theta),
                         prune = prune_tree_verb(tree, x = x), # Prune verb
                         change = change_tree_verb(tree, x = x, node_min_size = node_min_size, rotation = rotation, theta = theta), # Change verb
                         change_projection = change_projection_tree_verb(tree = tree,
                                                                         node_min_size = node_min_size,theta = theta,x = x),
                         swap = swap_tree_verb(tree, x = x, node_min_size = node_min_size) # Swap verb
                         
  )
  return(updated_tree)
}

# ==================================#
# Objects to test the tree_complete_conditional function
# ==================================#

# Go to RBART FUNCTION

# Calculating the full conditional probability with a new \nu value

tree_complete_conditional <- function(tree, x, residuals, nu = 1, phi = 1, 
                                      tau_mu,
                                      c_value, number_trees= number_trees,
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
    (sum(z$Omega_plus_I_inv) + c(tau_mu / c_value))#/(tau_mu)
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
  
  # aux <- c()
  # aux2 <- c()
  # for(i in 1:length(residuals_terminal_nodes)){
  #   aux[i] =+ sum(residuals_terminal_nodes[[i]]^2)
  #   aux2[i] =+ sum(residuals_terminal_nodes[[i]])
  # }
  # 
  # sum(aux)
  # sum(aux2)
  
  
  # # Residuals times tau
  # 0.5*n*log(tau) # First term
  # +0.5*sum(log(tau_mu/(tau_mu+nodes_size*tau))) # sencond element
  # -0.5*tau*sum(aux) # third element
  # 0.5*(tau^2)*sum((aux2^2)/(tau_mu+nodes_size*tau))
  
  return(list(log_posterior = log_posterior,
              S = S,
              RTR = RTR,
              R_Omega_I_one =  R_Omega_I_one))
}


# Generate mu_j values
update_mu <- function(tree,
                      x,
                      nu,
                      phi,
                      tau,
                      residuals,
                      c_value,
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
  
  # cat("Mu mean:",unlist(mu_mean),"\n")
  # cat("Mu sd:",unlist(mu_var),"\n")
  
  
  # Adding the mu values calculated
  for (i in 1:length(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- c(mu[[names_terminal_nodes[i]]])
  }
  
  return(tree)
}

# Update tau_j values
update_tau <- function(x, 
                       y,
                       a_tau,
                       d_tau, 
                       predictions) {
  # set.seed(seed)
  
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

# Update tau_j values
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
      distance_matrix =  nodes$distance_matrix,
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
    
    # Acessing trees
    # trees<-trees[[1]]
    
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



# Function to update nu values
update_nu <- function(x, y, current_tree_iter,residuals,
                      seed = NULL,
                      a_nu, d_nu, tau, tau_mu,
                      c_value, phi, nu,number_trees,
                      likelihood_object,p,gp_variables) {
  
  
  # Calculating the values of a and d
  a = a_nu 
  d = d_nu
  
  
  # Proposing a new \nu value
  # nu_proposal <- sample(seq(from = 1e-4,to = 1,length.out = 1000),size = 1) # Sampling new value
  # nu_proposal <- sample(seq(from = 0.1,to = 5,length.out = 1000),size = 1) # Sampling new value
  # nu_proposal <- sample(exp(seq(from = -1,to = 1,length.out = 1000)),size = 1) # Sampling new value
  # nu_proposal <- sample(seq(from = exp(-4),to = exp(4),length.out = 1000),size = 1) # Sampling new value
  # New nu_proposal based over the range of \nu
  nu_proposal <-runif(n = 1,min = 0.1*(0.25*sd(y)^(-2)), max = 20*(0.25*sd(y)^(-2)))
  
  # print(nu_proposal)
  
  # nu_proposal <- 1000
  
  # Getting the old likelihood object
  # likelihood_object<- tree_complete_conditional_from_mu(tree = current_tree_iter,
  #                                                             x = x, tau = tau,
  #                                                             residuals = residuals,
  #                                                             nu = nu,
  #                                                             tau_mu = tau_mu, c_value = c_value,
  #                                                             phi = phi,
  #                                                             number_trees = number_trees,tau_multiplier = tau_multiplier)
  
  # calculating the likelihood from the new step
  tree_from_nu_proposal <- inverse_omega_plus_I(tree = current_tree_iter,
                                                x = x,nu = nu_proposal, tau = tau,
                                                phi = phi,gp_variables = gp_variables)
  
  # likelihood_nu_proposal <- tree_complete_conditional_from_mu(tree = tree_from_nu_proposal,
  #                                                     x = x, tau = tau,
  #                                                     residuals = residuals,
  #                                                     nu = nu_proposal,
  #                                                     tau_mu = tau_mu, c_value = c_value,
  #                                                     phi = phi,
  #                                                     number_trees = number_trees,tau_multiplier = tau_multiplier)
  
  likelihood_nu_proposal <- tree_complete_conditional(tree = tree_from_nu_proposal,
                                                      x = x,
                                                      residuals = residuals,
                                                      nu = nu_proposal,
                                                      tau_mu = tau_mu, c_value = c_value,
                                                      phi = phi,
                                                      number_trees = number_trees,tau_multiplier = tau_multiplier)
  # Calculating the new likelihood
  
  # Old phi likelhood
  l_old_nu <- likelihood_object$log_posterior + 
    sum(dgamma(x = nu,shape = a,rate = d,log = TRUE))
  
  
  
  # Proposal likelihood
  l_proposal_nu <- likelihood_nu_proposal$log_posterior +
    sum(dgamma(x = nu_proposal,shape = a,rate = d,log = TRUE))
  
  
  
  
  # Probability of accept the new proposed tree
  acceptance_nu <- exp(l_proposal_nu - l_old_nu)
  
  # If storage for phi
  
  if (runif(1) < acceptance_nu) { #
    
    # Nu boolean to see if was accepted or not
    nu_boolean <- TRUE
    
    return(list(nu_boolean = nu_boolean,
                likelihood_object = likelihood_nu_proposal,
                tree = tree_from_nu_proposal,
                nu_proposal = nu_proposal)) # Returning the proposal value for nu
  }else{
    nu_boolean <- FALSE
    return( list(nu_boolean = nu_boolean)) # Returning the old value for nu
  } #
}

# ==============#
# rBart-GP FUNCTION
my_rBart_gp <- function(x, y, x_test, # Setting the training data matrix
                        number_trees = 2, # Setting the number of trees
                        node_min_size = round(0.05 * length(y)), # Min node size,
                        mu = 0,
                        alpha = 0.95, # Alpha from prior
                        beta = 2, # Beta from prior
                        tau = 1, # Tau from prior,
                        nu_vector = NULL, # This will be defining the nu the default value
                        phi_vector = rep(0.1 / (sqrt(number_trees)), number_trees), # This parameter is a "scale paramter" to the GP
                        c_value = number_trees,
                        n_iter = 1250, # Number of iterations
                        burn = n_iter/2, # Number of burn
                        thin = 1, # Number of thin
                        rotation = TRUE, # If rotated lon and lat will be used in tree building
                        theta = NULL, # If theta is NULL, then the rotation angle will be randomly selected
                        seed = NULL, # Alpha vector values from the Dirichlet prior
                        scale_boolean = TRUE,
                        a_tau = 1, # Prior from a_v_ratio gamma
                        d_tau = 1, # Prior from d_v_ratio gamma,
                        a_nu = 1, # Prior from a_tau
                        d_nu = 1, # Prior from d_tau
                        predictions=matrix(0, nrow = number_trees, ncol = nrow(x)),
                        predictions_list =  NULL,
                        discrete_phi_boolean = FALSE,
                        nu_multiplier = 1,tau_multiplier = 1, scale_multiplier = 1,
                        x_scale =  TRUE,
                        nu_update = TRUE,
                        phi_update = TRUE,
                        gp_variables = colnames(x),   # Selecting the GP-Variables
                        p = 1, # Shrink main parameter from GP-BART
                        K_bart,
                        error_handling_residuals = FALSE,
                        kappa = 0.5, bart_boolean = TRUE, bart_number_iter = (n_iter/2-1)
                        
) {
  
  
  # Scaling the train data
  data_scaled <- unit_cube_scale(rbind(x, x_test))
  
  # Adjusting the kappa (avoiding the Infinity error)
  if(kappa == 1 ){
    kappa <- kappa - 2*.Machine$double.eps
  }
  
  if(kappa == 0 ){
    kappa <- kappa + 2*.Machine$double.eps
    # print(kappa)
  }
  
  #==
  
  # Getting the maximum and minimum values from a distance matrix
  distance_matrix_x <- symm_distance_matrix(m1 = x)
  distance_range <- range(distance_matrix_x[upper.tri(distance_matrix_x)])
  distance_min <- sqrt(distance_range[1])
  distance_max <- sqrt(distance_range[2])
  
  # Option to scale X 
  if(x_scale){
    x <- data_scaled[1:nrow(x), , drop = FALSE]
    x_test_scaled <- data_scaled[(nrow(x) + 1):nrow(data_scaled), , drop = FALSE]
  } else {
    x_test_scaled <- x_test
  }
  
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
      nu_vector <- rep((4*number_trees*K_bart^2)/((max(y_scale)-min(y_scale))^2),number_trees)
    }
    
    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- (4*(number_trees)*K_bart^2)/((1-kappa)*(max(y_scale)-min(y_scale))^2)
    
  } else {
    
    # Not scaling the y
    y_scale <- y
    
    if( is.null(nu_vector) ){
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4*number_trees*K_bart^2)/((max(y_scale)-min(y_scale))^2),number_trees)
    }
    
    
    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- (4*(number_trees)*K_bart^2)/((1-kappa)*(max(y_scale)-min(y_scale))^2)
    
  }
  
  # Getting the number of observations 
  n <- length(y)
  
  # Adjusting the \nu vector with the kappa values
  nu_vector <- nu_vector / kappa
  
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
    stop("Inset a valid NAMED matrix")
  }
  
  if (node_min_size == 0) {
    stop("Node Minimum Size need to be greater than 0")
  }
  
  if (length(nu_vector) != number_trees) {
    stop("Inset a valid \\nu vector for the number of trees")
  }
  
  if (length(phi_vector) != number_trees) {
    stop("Inset a valid \\phi vector for the number of trees")
  }
  
  # Set the minimum node size to 20
  # if (node_min_size <= round(0.05 * length(y))) {
  #   warning("\n The default value for min_node_size was settled to 15.")
  #   min_node_size <- 15
  # }
  
  # Recomendation about the min_node_size
  if (node_min_size < 15) {
    warning("\n It is recommended that the min_node_size should be of at least 15 observations.")
  }
  
  # Storage containers
  store_size <- (n_iter - burn) / thin
  tree_store <- vector("list", store_size)
  tau_store <- c()
  nu_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  signal_pc_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  
  
  
  y_hat_store <- matrix(NA, ncol = length(y), nrow = store_size)
  y_hat_store_proposal <- matrix(NA, ncol = length(y), nrow = store_size)
  
  # Storing the likelihoods
  log_lik_store <- rep(NA, store_size)
  log_lik_store_fixed_tree <- rep(NA,store_size)
  log_like_nu_store_proposal <- matrix(NA, nrow = (n_iter), ncol = 2)
  
  
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
  
  # if(number_trees==1) current_trees<-current_trees[[1]] #Just to simplify in case of just one tree
  
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
  
  # Getting the verb list
  # verb_store <- data.frame(verb = NA, accepted = NA)
  
  
  # Setting initial values for phi vector
  for (i in 1:n_iter) {
    
    # print(verb_store)
    
    utils::setTxtProgressBar(progress_bar, i)
    
    # Returning to the GP-BART approach
    # if(i >= bart_number_iter){
    #   bart_boolean <- FALSE
    # }
    
    if ((i > burn) & ((i %% thin) == 0)) {
      
      curr <- (i - burn) / thin
      tree_store[[curr]] <- current_trees
      if(scale_boolean){
        
        # Returning to the original
        # tau_store[[curr]] <- (tau)/sqrt((b-a))
        tau_store[[curr]] <- tau
      } else{
        tau_store[[curr]] <- tau
      }
      y_hat_store[curr, ] <- if (number_trees == 1) {
        
        # In case of normalizing
        if(scale_boolean){
          # Keep it the regular scale
          # unnormalize_y(predictions, a = a, b = b)
          predictions 
        }else{
          predictions
        }
      } else {
        # colSums(unnormalize_y(predictions, a = a, b = b))
        colSums(predictions)
      }
      
      current_partial_residuals_list[[curr]] <- 
        if(scale_boolean){
          
          # I need to keep the same because of the \nu values 
          # unnormalize_y(current_partial_residuals_matrix, a = a, b = b)
          current_partial_residuals_matrix
        } else {
          current_partial_residuals_matrix
        } 
      
      # Saving the predictions
      current_predictions_list[[curr]] <- predictions
      
      
      # y_hat_store_proposal[curr,] = predictions_proposal
      phi_store[curr, ] <- phi_vector
      phi_proposal_store[curr, ] <- phi_vector_proposal
      nu_store[curr, ] <- nu_vector
      # print(verb_store)
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
        
        # print(verb_store)
        
        # # # To update the mu values
        current_trees[[j]] <- update_mu_bart(
          tree = current_trees[[j]],
          x = x,
          residuals = current_partial_residuals, 
          tau = tau,
          tau_mu = tau_mu)
        
        
        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j, ] <- update_residuals_bart(
          tree = current_trees[[j]], x = x,
          residuals = current_partial_residuals
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
        likelihood_new <- tree_complete_conditional(
          tree = new_trees[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          nu = nu_vector[j], phi = phi_vector[j], c_value = c_value,
          number_trees = number_trees,
          tau_multiplier = tau_multiplier
        )
        
        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional(
          tree = current_trees[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          nu = nu_vector[j], phi = phi_vector[j], c_value = c_value,
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
          residuals = current_partial_residuals,
          likelihood_object = likelihood_object,
          phi = phi_vector[j], nu = nu_vector[j], c_value = c_value)
        
        
        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j, ] <- update_residuals(
          tree = current_trees[[j]], x = x,
          residuals = current_partial_residuals,
          phi = phi_vector[j], nu = nu_vector[j],tau = tau,
          error_handling_residuals = error_handling_residuals # Boolean to check if handle the residuals or not
        )
        
        # predictions[j, ] <- update_g(
        #   tree = current_trees[[j]], x = x,
        #   residuals = current_partial_residuals,
        #   phi = phi_vector[j], nu = nu_vector[j],p = p
        # )
        
        # Plot the current residuals, and the prediction from it
        # Updating the \nu value
        if(nu_update){
          mh_update_nu<-update_nu(x = x, y = y_scale,
                                  current_tree_iter = current_trees[[j]],
                                  likelihood_object = likelihood_object,
                                  residuals = current_partial_residuals,
                                  phi = phi_vector[j],
                                  nu = nu_vector[j],
                                  a_nu = a_nu,
                                  d_nu = d_nu,
                                  c_value = c_value,
                                  number_trees = number_trees,
                                  gp_variables = gp_variables,
                                  tau = tau,
                                  tau_mu = tau_mu)
          
          # In case of accept the update over \nu update everything
          if(mh_update_nu$nu_boolean){
            
            # Updating the tree and the \nu object from the tree
            current_trees[[j]] <- mh_update_nu$tree
            
            # Updating the likelihood objects
            likelihood_object <- mh_update_nu$likelihood_object
            
            # Updating the nu value
            nu_vector[j] <- mh_update_nu$nu_proposal
            
          } # If doesn't accept, nothing changes.
        }
        
        # To update phi
        if(phi_update){
          mh_update_phi<- update_phi_marginal(current_tree_iter = current_trees[[j]],
                                              residuals = current_partial_residuals,
                                              x = x,nu = nu_vector[j],phi = phi_vector[j],
                                              gp_variables = gp_variables,
                                              c_value = c_value,
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
    nu_store = nu_store,
    y = y_scale,
    X = x,
    scale_boolean = scale_boolean,
    acc_ratio = acc_ratio,
    acc_ratio_nu = acc_ratio_nu,
    acc_ratio_phi = acc_ratio_phi,
    signal_pc_store = signal_pc_store,
    iter = n_iter,
    burn = burn,
    thin = thin,
    store_size = store_size,
    number_trees = number_trees, node_min_size = node_min_size, c_value = c_value,
    log_like_nu_store_proposal = log_like_nu_store_proposal,
    a_min = a_min,
    b_max = b_max,
    a_tau = a_tau,
    d_tau = d_tau,
    a_nu = a_nu,
    d_nu = d_nu,
    y_hat_store_proposal = y_hat_store_proposal, current_partial_residuals_list = current_partial_residuals_list,
    beta = beta, x_test = x_test_scaled,
    current_predictions_list = current_predictions_list, nu_multiplier = nu_multiplier , scale_multiplier = scale_multiplier,
    p = p, tau_mu = tau_mu, kappa = kappa,
    verb_store_list = verb_store_list
  ))
}


# #Do a MH for PHI
update_phi_marginal <- function(x, current_tree_iter,residuals,
                                seed = NULL,
                                tau, 
                                c_value,tau_mu,
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
  
  
  likelihood_phi_proposal <- tree_complete_conditional(tree = tree_from_phi_proposal,
                                                       x = x,
                                                       residuals = residuals,
                                                       nu = nu, tau_mu = tau_mu,
                                                       phi = phi_proposal,c_value = c_value,
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
                       c_value,tau_mu,
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
  
  
  likelihood_phi_proposal <- tree_complete_conditional(tree = tree_from_phi_proposal,
                                                       x = x,
                                                       residuals = residuals,
                                                       nu = nu, tau_mu = tau_mu,
                                                       phi = phi_proposal,c_value = c_value,
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


# #Do a MH for PHI
update_signal_pc <- function(current_iter_tree, # Current j tree from the iteration
                             x = x, 
                             residuals = current_partial_residuals,
                             nu,
                             phi,
                             signal_pc,
                             alpha_signal_pc,
                             a_v_ratio,
                             d_v_ratio,
                             beta_signal_pc,
                             c_value = c_value,
                             v_ratio,number_trees,
                             tau_multiplier
){
  
  # Calculating the values of a and d (FOR \TAU)
  a_tau = tau_multiplier/( (1-signal_pc)*a_v_ratio )
  d_tau = 1/d_v_ratio
  
  # Calculating the values of a and d for \nu
  a_nu = 1*(1-signal_pc)/( signal_pc*a_v_ratio )
  d_nu = 1/d_v_ratio
  
  
  # Increased the range of tree proposal
  signal_pc_proposal <- sample(seq(from = 0.3001, to = 1, by = 0.001), size = 1) # Sampling new value
  
  # Calculating the proposal a and d for tau
  a_proposal_tau = tau_multiplier/( (1-signal_pc_proposal)*a_v_ratio )
  d_proposal_tau = 1/d_v_ratio
  
  
  # Calculating the proposal of a and d for \nu
  a_proposal_nu = 1*(1-signal_pc_proposal)/( signal_pc_proposal*a_v_ratio )
  d_proposal_nu = 1/d_v_ratio
  
  # Selecting the terminal nodes
  terminal_nodes <- current_iter_tree[names(which(sapply(current_iter_tree, function(x) {
    x$terminal == 1
  })))]
  
  # Seleciting the tau values
  tau_j <- sapply(terminal_nodes, function(x) {
    x$tau
  })
  
  
  
  # Calculating the likelihood for the proposal value
  l_proposal_signal_pc <- sum(dgamma(tau_j,shape = a_proposal_tau,rate=d_proposal_tau,log = TRUE))+
    sum(dgamma(nu,shape = a_proposal_nu,rate=d_proposal_nu,log = TRUE))+
    dbeta(x = signal_pc_proposal, shape1 = alpha_signal_pc, shape2 = beta_signal_pc, log = TRUE)
  
  
  # Calculating the likelihood for the old value
  l_old_signal_pc<- sum(dgamma(tau_j,shape = a_tau,rate=d_tau,log = TRUE))+
    sum(dgamma(nu,shape = a_nu,rate=d_nu,log = TRUE))+
    dbeta(x = signal_pc, shape1 = alpha_signal_pc, shape2 = beta_signal_pc, log = TRUE)
  
  
  # Probability of accept the new proposed tree
  acceptance_signal_pc <- exp(l_proposal_signal_pc - l_old_signal_pc)
  
  # If storage for phi
  
  if (runif(1) < acceptance_signal_pc) { #
    return(signal_pc_proposal) # Returning the proposal value for phi
  }else{
    return(signal_pc) # Returning the old value for phi
  } #
  
}


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

tree_count_terminals <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)
  
  for (k in 1:length(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in 1:gpbart_model$number_trees) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- lapply(tree, function(x) {
        x$terminal
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


get_tau_ratio_values <- function(gpbart_model) {
  
  # n_iter
  n_iter <- length(gpbart_model$trees)
  
  tau_ratio_iters <- rep(list(NA), gpbart_model$number_trees)
  
  num_trees <- gpbart_model$number_trees
  
  for (j in 1:n_iter) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for tau_ratio terminals
    tau_ratio_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(tau_ratio_terminals) <- paste0("node_", 0:49)
    
    for (k in 1:num_trees) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)
      
      
      for (nodes in 1:length(all_nodes)) {
        if (tree[[all_nodes[nodes]]]$terminal == 1) {
          # tree[[nodes]]$tau_ratio <- sample(c(0.001,seq(0.005,10,by = 0.005),10,50,100), 1)
          tau_ratio_terminals[k, as.numeric(stringr::str_extract(all_nodes[nodes], pattern = "\\(?[0-9,.]+\\)?")) + 1] <- tree[[all_nodes[nodes]]]$tau_ratio
        }
      }
    }
    
    # Saving all
    tau_ratio_iters[[j]] <- tau_ratio_terminals
  }
  return(tau_ratio_iters)
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

my_predict_rBART <- function(rBart_model, type = c("all", "mean", "median"),
                             pred_bart_only = FALSE) {
  
  # Loading x_test
  x_test <- rBart_model$x_test
  x_train <- rBart_model$X
  
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
    
    # Testing the 100 MH run
    # i=800
    
    # Selecting one tree from BART model
    current_tree <- rBart_model$trees[[i]]
    
    # Getting the predictions from the test observations
    y_pred_aux <- predict_gaussian_from_multiple_trees(
      multiple_trees = current_tree, x_train = rBart_model$X,
      x_new = x_test, partial_residuals = rBart_model$current_partial_residuals_list[[i]],
      phi_vector = rBart_model$phi_store[i, ], nu_vector = rBart_model$nu_store[i,],
      tau = rBart_model$tau_store[i],pred_bart_only = pred_bart_only
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
    # Calculating the diagonal from all matrix
    # cov_hat_matrix[i, ] <- sqrt(diag(Reduce("+", cov_pred_final)))
    
    # cov_hat_matrix[i,]<-sqrt(diag(do.call("mean",model_each_aux$variance_matrix)))
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
                mean = mean(rBart_model$tau_store^(-1/2)),
                median = median(rBart_model$tau_store^(-1/2))
    )
  )
  
  
  
  return(list(out = out, list_matrix_pred = y_list_matrix,
              list_matrix_cov = cov_hat_matrix_list,
              list_matrix_cov_sqrt_sum = cov_hat_matrix_sqrt_sum,
              mcmc_pi_mean =  y_hat_matrix,
              mcmc_pi_sd = y_sd_pi_matrix))
  
}




rmse <- function(obs, pred) {
  return(sqrt(mean((obs - pred)^2)))
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

# Getting a fixed and true TREE
fixed_tree <- function(x_train,x_test){
  
  # Scaling the train data
  data_scaled <- unit_cube_scale(rbind(x, x_test))
  # data_scaled <- (rbind(x, x_test))
  
  x <- data_scaled[1:nrow(x), , drop = FALSE]
  
  x_test_scaled <- data_scaled[(nrow(x) + 1):nrow(data_scaled), , drop = FALSE]
  
  # -0.1575854  0.5167328
  
  
  
  node_1 <- list(node_number = 1,
                 observations_index = which(x_train <= -0.157),
                 depth_node = 1,
                 node_var = "x",
                 node_var_split  = -0.158,
                 left = 1,
                 right = 0,
                 parent_node = NA,
                 terminal = 1,
                 mu = 1,
                 tau = 1)
  
  node_2 <- list(node_number = 1,
                 observations_index = which(x_train <= -0.157 & x_train <= 0.517),
                 depth_node = 1,
                 node_var = "x",
                 node_var_split  = 0.517,
                 left = 1,
                 right = 0,
                 parent_node = 1,
                 terminal = 1,
                 mu = 1,
                 tau = 1)
  
  node_3 <- list(node_number = 1,
                 observations_index = which(x_train <= -0.157 & x_train <= 0.517),
                 depth_node = 1,
                 node_var = "x",
                 node_var_split  = 0.517,
                 left = 0,
                 right = 1,
                 parent_node = 1,
                 terminal = 1,
                 mu = 1,
                 tau = 1)
  
  return(list(node_1 = node_1, node_2 = node_2, node_3 = node_3))
  
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
      distance_matrix =  dist_m,
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
      distance_matrix = nodes$distance_matrix,
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


# Getting the function to normalize the bart
normalize_y <- function(y) {
  a <- min(y)
  b <- max(y)
  y <- (y - a) / (b - a) - 0.5
  return(y)
}

unnormalize_y <- function(z, a, b) {
  z <- (b - a) * (z + 0.5) + a
  return(z)
}


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
      (kernel_function(distance_matrix = nodes$distance_matrix,
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




# Naive sigma_estimation
naive_sigma <- function(x,y){
  
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
  
  return(sigma)
  
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
  zero_tau_prob_squared <- function(x,naive_tau_value, prob , shape){
    
    # Find the zero to the function P(tau < tau_ols ) = 0.1, for a defined   
    return ((pgamma(naive_tau_value,
                    shape = shape,
                    rate = x)-(1-prob))^2 )
    
  }
  
  
  # Getting the root
  min_root <-  try(uniroot(f = zero_tau_prob,interval = c(1e-12,100),
                           naive_tau_value = tau_ols,
                           prob = prob, shape = shape)$root,silent = TRUE)
  
  if( class(min_root) =="try-error"){
    # Verifying the squared version
    min_root <- optim(par = runif(1),fn = zero_tau_prob_squared,
                      method = "L-BFGS-B",lower = c(0), naive_tau_value = tau_ols,
                      prob = prob, shape = shape)$par
    
  }
  
  
  return(min_root)
  
}


# New PI calibration dataset
pi_calibration <- function(gpbart_mod,prob = 0.5, quantile = FALSE){
  
  # Decodomg of I will get it from the quantile or not
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

# This function scale the data to be between -0.5 and 0.5
new_scale <- function(x) {
  
  # Scaling to -1 to 1 function
  scaling <- function(z) {
    (2 * z - (max(z) + min(z))) / (max(z) - min(z))
  }
  
  # Applying on all covariates
  x_scaled <- apply(x, 2, scaling)
  return(x_scaled*0.5)
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


