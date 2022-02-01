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
# Objects to test the grow_project_tree_verb function
# ==================================#

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
# Objects to test the grow__tree_verb function
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
      
      
      # Selecting the splitting point on the variable, unformly on the observed instances
      node_var_split <- sample(x[current_node$observations_index, node_var], size = 1) 
      
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
      
      # Selecting from the observations
      node_var_split <- sample(
        levels(x[current_node$observations_index, node_var]), 
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


# ==================================#
# Objects to test the change_tree_verb function
# ==================================#

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
# Objects to test the change_projction_tree_verb function
# ==================================#

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


# Sort an action to update the tree
update_tree_verb_bart <- function(tree, x, node_min_size, verb) {
  
  # Update the tree by a verb
  updated_tree <- switch(verb,
                         grow = grow_tree(tree,
                                          x = x, node_min_size = node_min_size, rotation = rotation,
                                          theta = theta
                         ), 
                         prune = prune_tree_verb(tree, x = x), # Prune verb
                         change = change_tree_verb(tree, x = x, node_min_size = node_min_size, rotation = rotation, theta = theta), # Change verb
                         swap = swap_tree_verb(tree, x = x, node_min_size = node_min_size) # Swap verb
                         
  )
  return(updated_tree)
}


