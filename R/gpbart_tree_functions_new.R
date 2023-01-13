# Calculating the node likelihood
node_loglikelihood_gpbart <- function(node,
                                x_train,
                                res_vec,
                                tau,
                                tau_mu,
                                nu,
                                phi_vector,
                                gp_variables){

        # Slicing the current res_vec
        res_node <- res_vec[node$obs_train]
        n_obs <- length(res_node)


        # Subsetting x_train
        x_train_node <- x_train[node$obs_train,gp_variables,drop = FALSE]
        sq_dist_matrix <- symm_distance_matrix(m1 = x_train_node,phi_vector = phi_vector)
        cov_function <- kernel_function(squared_distance_matrix_phi = sq_dist_matrix,nu = nu)+diag(1/tau,nrow = nrow(x_train_node)) + 1/tau_mu

        # Getting the likelihood for that node
        log_likelihood <- mvnfast::dmvn(X = res_node,mu = rep(0,nrow(cov_function)),sigma = cov_function,log = TRUE)

        return(log_likelihood)

}

# Creating a function to grow a tree
grow_gpbart <- function(res_vec,
                 tree,
                 x_train,
                 x_test,
                 xcut,
                 tau,
                 tau_mu,
                 alpha,
                 beta,
                 node_min_size,
                 # Passing by nu arguments
                 nu,
                 phi_vector_p,
                 cov_gp,
                 cat_var = NULL){


        # Accepting or not the verb
        accepted <- FALSE

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree)

        # Sampling one terminal node
        g_node_position <- sample(1:length(terminal_nodes),size = 1)
        g_node <- terminal_nodes[[g_node_position]]
        g_node_name <- names(terminal_nodes[g_node_position])
        g_node_position_orig <- which(names(tree) ==g_node_name)

        # Initializing the sample
        split_var_candidates <- colnames(x_train)
        good_tree_index <- 0


        # GETTING CATEGORICAL VARIABLES
        if(is.null(cat_var)){
                cat_var <- names(x_train)[!apply(x_train,2,function(x){length(unique(x))>round(nrow(x_train)*0.4)})]
        }

        while(good_tree_index==0){
                # Selecting a valid split
                split_var <- sample(split_var_candidates,size = 1)


                # For the case of categorical variables
                if(split_var %in% cat_var){

                        # Selecting the split rule for categorical variables
                        split_var_sampled_rule <- sample(unique(x_train[[split_var]]),size = 1)

                        if((sum(x_train[g_node$obs_train,split_var]==split_var_sampled_rule)>=node_min_size) & (sum(x_train[g_node$obs_train,split_var]!=split_var_sampled_rule)>=node_min_size)){
                                good_tree_index <- 1
                        } else {
                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                                }
                        }


                }       else       {

                                # Avoiding invalid max.
                                if((length(x_train[g_node$obs_train,split_var])-node_min_size)<1){
                                        return(list(tree = tree, accepted = accepted))
                                }

                                # Getting the min and maximum observed value within the terminal node
                                min_node_obs <- sort(x_train[g_node$obs_train,split_var])[node_min_size]
                                max_node_obs <- sort(x_train[g_node$obs_train,split_var])[length(x_train[g_node$obs_train,split_var])-node_min_size]


                                # Getting the column from xcut
                                xcut_valid <- xcut[which(xcut[,split_var]>=min_node_obs & xcut[,split_var]<=max_node_obs),split_var]


                                # No valid tree found
                                if(length(xcut_valid) == 0 ){

                                        split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                        if(length(split_var_candidates)==0){
                                                return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                                        }

                                } else {
                                        good_tree_index <- 1

                                        # Getting only unique values of xcut_valid
                                        xcut_valid <- unique(xcut_valid)
                                        split_var_sampled_rule <- sample(xcut_valid,size = 1)
                                }
                }# End for no categorical variables
        }



        # Creating the left and the right nodes
        max_index <- max(get_indexes(tree))
        max_tree_size <- length(tree)

        if(split_var %in% cat_var){
                # Creating the vector of the splitting rules for the categorical
                left_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]==split_var_sampled_rule)]
                right_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]!=split_var_sampled_rule)]

                left_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]==split_var_sampled_rule)]
                right_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]!=split_var_sampled_rule)]

        } else {
                # Creating the vector of new train and test index
                left_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]<=split_var_sampled_rule)]
                right_train_id <- g_node$obs_train[which(x_train[g_node$obs_train,split_var]>split_var_sampled_rule)]

                left_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]<=split_var_sampled_rule)]
                right_test_id <- g_node$obs_test[which(x_test[g_node$obs_test,split_var]>split_var_sampled_rule)]
        }
        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(list(tree = tree, accepted = accepted))
        }

        # Creating the left node
        left_node <- list(index = max_index+1,
                          obs_train = left_train_id,
                          obs_test  = left_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$depth+1,
                          var = split_var,
                          var_split_rule = split_var_sampled_rule,
                          mu = 0)

        right_node <- list(index = max_index+2,
                           obs_train = right_train_id,
                           obs_test  = right_test_id,
                           left = NA,
                           right = NA,
                           parent = g_node$index,
                           terminal = 1,
                           nog = 0,
                           depth = g_node$depth+1,
                           var = split_var,
                           var_split_rule = split_var_sampled_rule,
                           mu = 0)

        # Modifying the new g_node
        new_g_node <-g_node
        new_g_node$left <- max_index+1
        new_g_node$right <- max_index+2
        new_g_node$terminal <- 0
        new_g_node$nog = 1 # It always be a a nog since is given origin to two children


        # Get nog counter ( FOR THE NEW TREE )
        nog_counter <- count_nog(tree = tree[-g_node_position]) + 1

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood_gpbart(res_vec = res_vec,node = left_node,
                                                    x_train = x_train,
                                                    tau = tau,tau_mu = tau_mu,nu = nu,
                                                    phi_vector =  phi_vector_p,gp_variables = cov_gp) +
                node_loglikelihood_gpbart(res_vec = res_vec,node = right_node,
                                          x_train = x_train,
                                          tau = tau,tau_mu = tau_mu,nu = nu,
                                          phi_vector =  phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec,node = g_node,
                                          x_train = x_train,
                                   tau = tau, tau_mu = tau_mu,
                                   nu = nu,phi_vector =  phi_vector_p,gp_variables = cov_gp)

        # Calculate the transition
        transition_loglike <- log(0.3/nog_counter)-log(0.3/length(terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- 2*log(1-alpha*(1+(g_node$depth+1))^(-beta)) + (-beta)*log(1+g_node$depth)+log(alpha) - log(1-alpha*(1+g_node$depth)^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){
                # Maybe use append to make everything easier
                tree[[g_node_name]] <- new_g_node
                # Transform the parent into a nog
                if(!is.na(g_node$parent)){
                        tree[[paste0("node_",g_node$parent)]]$nog <- 0
                }

                new_nodes <- list(left_node,right_node)
                names(new_nodes) <- c(paste0("node_",c(new_nodes[[1]]$index,new_nodes[[2]]$index)))
                tree <- append(tree,new_nodes,after = g_node_position_orig)
                accepted <- TRUE
        }

        return(list(tree = tree, accepted = accepted))

}


# Creating a function to grow a tree
grow_rotation_gpbart <- function(res_vec,
                        tree,
                        x_train,
                        x_test,
                        xcut,
                        tau,
                        tau_mu,
                        alpha,
                        beta,
                        node_min_size,
                        # Passing by nu arguments
                        nu,
                        phi_vector_p,
                        cov_gp,
                        rotation_variables){

        # ACC
        accepted <- FALSE
        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree)

        # Sampling one terminal node
        g_node_position <- sample(1:length(terminal_nodes),size = 1)
        g_node <- terminal_nodes[[g_node_position]]
        g_node_name <- names(terminal_nodes[g_node_position])
        g_node_position_orig <- which(names(tree) ==g_node_name)

        # Initializing the sample
        split_var_candidates <- 1:length(rotation_variables)
        good_tree_index <- 0


        while(good_tree_index==0){

                # Selecting the pair node that will be selected
                split_var_pair <- sample(rotation_variables,2)

                # Selecting a valid split
                split_var <- sample(split_var_pair,size = 1)

                # Selecting an angle to rotate my coordinates
                theta <- stats::runif(n = 1,min = 0,max = pi)

                # Creating the rotated coordinates
                rotated_x <- tcrossprod(A(theta), x_train[,split_var_pair])
                rownames(rotated_x) <- split_var_pair

                # Getting the rotation for the test
                rotated_x_test <- tcrossprod(A(theta), x_test[,split_var_pair])
                rownames(rotated_x_test) <- split_var_pair

                if((length(rotated_x[split_var,g_node$obs_train])-node_min_size)<1){
                        return(list(tree = tree, accepted = accepted))
                }

                # Getting the min and maximum observed value within the terminal node
                min_node_obs <- sort(rotated_x[split_var,g_node$obs_train])[node_min_size]
                max_node_obs <- sort(rotated_x[split_var,g_node$obs_train])[length(rotated_x[split_var,g_node$obs_train])-node_min_size]


                # Getting the x_cut matrix rotated
                xcut_rotated <- tcrossprod(A(theta), xcut[,split_var_pair])
                rownames(xcut_rotated) <- split_var_pair

                # Getting the column from xcut
                xcut_valid <- xcut_rotated[split_var,which(xcut_rotated[split_var,]>=min_node_obs & xcut_rotated[split_var,]<=max_node_obs)]


                # No valid tree found
                if(length(xcut_valid) == 0 ){

                        rotation_variables <-  rotation_variables[!(rotation_variables %in% split_var_pair)]

                        if(length(rotation_variables)==0 || length(rotation_variables)==1){
                                return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                        }

                } else {
                        good_tree_index <- 1
                }
        }
        # Sampling a x_cut_rule
        xcut_valid <- unique(xcut_valid)
        split_var_sampled_rule_rotation <- sample(xcut_valid,size = 1)

        # Creating the left and the right nodes
        max_index <- max(get_indexes(tree))
        max_tree_size <- length(tree)

        # Creating the vector of new train and test index
        left_train_id <- g_node$obs_train[which(rotated_x[split_var,g_node$obs_train]<=split_var_sampled_rule_rotation)]
        right_train_id <- g_node$obs_train[which(rotated_x[split_var,g_node$obs_train]>split_var_sampled_rule_rotation)]

        left_test_id <- g_node$obs_test[which(rotated_x_test[split_var,g_node$obs_test]<=split_var_sampled_rule_rotation)]
        right_test_id <- g_node$obs_test[which(rotated_x_test[split_var,g_node$obs_test]>split_var_sampled_rule_rotation)]

        # No valid tree
        if((length(left_train_id) < node_min_size) || (length(right_train_id)<node_min_size)){
                return(list(tree = tree, accepted = accepted))
        }

        # Creating the left node
        left_node <- list(index = max_index+1,
                          obs_train = left_train_id,
                          obs_test  = left_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$depth+1,
                          var = list(split_var_pair = split_var_pair, split_var = split_var , theta = theta),
                          var_split_rule = split_var_sampled_rule_rotation,
                          mu = 0)

        right_node <- list(index = max_index+2,
                           obs_train = right_train_id,
                           obs_test  = right_test_id,
                           left = NA,
                           right = NA,
                           parent = g_node$index,
                           terminal = 1,
                           nog = 0,
                           depth = g_node$depth+1,
                           var = list(split_var_pair = split_var_pair, split_var = split_var, theta = theta),
                           var_split_rule = split_var_sampled_rule_rotation,
                           mu = 0)

        # Modifying the new g_node
        new_g_node <-g_node
        new_g_node$left <- max_index+1
        new_g_node$right <- max_index+2
        new_g_node$terminal <- 0
        new_g_node$nog = 1 # It always be a a nog since is given origin to two children


        # Get nog counter ( FOR THE NEW TREE )
        nog_counter <- count_nog(tree = tree[-g_node_position]) + 1

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood_gpbart(res_vec = res_vec,node = left_node,
                                                    x_train = x_train,
                                                    tau = tau,tau_mu = tau_mu,nu = nu,
                                                    phi_vector =  phi_vector_p,gp_variables = cov_gp) +
                node_loglikelihood_gpbart(res_vec = res_vec,node = right_node,
                                          x_train = x_train,
                                          tau = tau,tau_mu = tau_mu,nu = nu,
                                          phi_vector =  phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec,node = g_node,
                                          x_train = x_train,
                                          tau = tau, tau_mu = tau_mu,
                                          nu = nu,phi_vector =  phi_vector_p,gp_variables = cov_gp)

        # Calculate the transition
        transition_loglike <- log(0.3/nog_counter)-log(0.3/length(terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- 2*log(1-alpha*(1+(g_node$depth+1))^(-beta)) + (-beta)*log(1+g_node$depth)+log(alpha) - log(1-alpha*(1+g_node$depth)^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){
                # Maybe use append to make everything easier
                tree[[g_node_name]] <- new_g_node
                # Transform the parent into a nog
                if(!is.na(g_node$parent)){
                        tree[[paste0("node_",g_node$parent)]]$nog <- 0
                }

                new_nodes <- list(left_node,right_node)
                names(new_nodes) <- c(paste0("node_",c(new_nodes[[1]]$index,new_nodes[[2]]$index)))
                tree <- append(tree,new_nodes,after = g_node_position_orig)
                accepted <- TRUE
        }

        return(list(tree = tree, accepted = accepted))

}


# Pruning a tree
prune_gpbart <- function(tree,
                  x_train,
                  res_vec,
                  tau,
                  tau_mu,
                  alpha,
                  beta,
                  nu,
                  phi_vector_p,
                  cov_gp){

        #Accept bool
        accepted <- FALSE
        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        t_nodes <- (get_terminals(tree = tree))

        n_terminal_nodes <- length(t_nodes)
        n_nogs <- length(nog_nodes)

        # Returning the  a simple tree
        if(length(nog_nodes)==0){
                return(list(tree = tree, accepted = accepted))
        }

        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                p_node <- nog_nodes[[nog_nodes_index]]
        }


        # Name node to be pruned
        p_node_name <- paste0("node_",p_node$index)

        # Getting the name of non terminals
        names_non_terminals <- names(tree[!(names(tree) %in% names(t_nodes))])
        names_non_terminals <- names_non_terminals[names_non_terminals!=p_node_name] # Removing the current pruned node
        parents_non_t_nodes <- sapply(tree[names_non_terminals],function(x){x$parent})

        left_node_name <- paste0("node_",p_node$left)
        left_node <- tree[[paste0("node_",p_node$left)]]
        right_node_name <- paste0("node_",p_node$right)
        right_node <- tree[[paste0("node_",p_node$right)]]

        # Calculating the loglikelihood of the new tree
        tree_loglikeli <- node_loglikelihood_gpbart(res_vec = res_vec,node = p_node,tau = tau, tau_mu = tau_mu,
                                                    x_train = x_train,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec,node = left_node,tau = tau,tau_mu = tau_mu,
                                          x_train = x_train, nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec,node = right_node, tau = tau,tau_mu = tau_mu,
                                          x_train = x_train, nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp)

        # Calculate the transition
        transition_loglike <- log(0.3/(n_nogs))-log(0.3/length(n_terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- log(1-alpha*(1+p_node$depth)^(-beta)) - ((-beta)*log((1+p_node$depth))+log(alpha)) - 2*log(1-alpha*(1+(p_node$depth+1))^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Deciding weather accept or not
        if(stats::runif(n = 1) < exp(log_acceptance)){
                tree[[p_node_name]]$terminal <- 1
                tree[[p_node_name]]$left <- NA
                tree[[p_node_name]]$right <- NA
                tree[[p_node_name]]$nog <- 0

                # Removing the pruned nodes
                tree[[left_node_name]] <- NULL
                tree[[right_node_name]] <- NULL

                # The parent of the pruned node need to return to be a nog

                # Checking if AFTER pruning this node its parent become a NOG
                if(p_node_name!="node_0"){
                        new_t_nodes_names <- c(names(t_nodes),p_node_name)
                        pruned_node_parent <- tree[[paste0("node_",tree[[p_node_name]]$parent)]]
                        if((paste0("node_",pruned_node_parent$left) %in% new_t_nodes_names) & (paste0("node_",pruned_node_parent$right) %in% new_t_nodes_names)){
                                tree[[paste0("node_",tree[[p_node_name]]$parent)]]$nog <- 1
                        }
                }

                accepted <- TRUE
        }

        return(list(tree = tree, accepted = accepted))

}

# Changing the tree
change_gpbart <- function(res_vec,
                   tree,
                   x_train,
                   x_test,
                   xcut,
                   tau,
                   tau_mu,
                   alpha,
                   beta,
                   node_min_size,
                   nu,
                   phi_vector_p,
                   cov_gp,
                   cat_var = NULL){

        # ACC
        accepted <- FALSE
        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        n_terminal_nodes <- length(get_terminals(tree = tree))

        if(length(tree)==1 || length(tree) == 0){
                return(list(tree = tree, accepted = accepted))
        }


        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                c_node <- nog_nodes[[nog_nodes_index]]
        }

        # Case of just one split
        if( length(tree)==3){
                c_node <- tree[[1]] ## Getting the root node
        }

        if((length(nog_nodes) == 0) & length(tree)!=3){
                return(list(tree = tree, accepted = accepted))
        }

        good_tree_index <- 0

        # Getting the name of the changed node
        split_var_candidates <- colnames(x_train)

        if(is.null(cat_var)){
                cat_var <- names(x_train)[!apply(x_train,2,function(x){length(unique(x))>round(nrow(x_train)*0.4)})]
        }

        while(good_tree_index==0){

                # Selecting a valid split
                split_var <- sample(split_var_candidates,size = 1)

                if(split_var %in% cat_var){

                        # Selecting the split rule for categorical variables
                        split_var_sampled_rule <- sample(unique(x_train[[split_var]]),size = 1)

                        if((sum(x_train[c_node$obs_train,split_var]==split_var_sampled_rule)>=node_min_size) & (sum(x_train[c_node$obs_train,split_var]!=split_var_sampled_rule)>=node_min_size)){
                                good_tree_index <- 1
                        } else {
                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                                }
                        }

                } else {


                        # Case of invalid max
                        if((length(x_train[c_node$obs_train,split_var])-node_min_size)<1){
                                return(list(tree = tree, accepted = accepted))
                        }


                        # Getting the min and maximum observed value within the terminal node
                        min_node_obs <- sort(x_train[c_node$obs_train,split_var])[node_min_size]
                        max_node_obs <- sort(x_train[c_node$obs_train,split_var])[length(x_train[c_node$obs_train,split_var])-node_min_size]

                        # Getting the column from xcut
                        xcut_valid <- xcut[which(xcut[,split_var]>=min_node_obs & xcut[,split_var]<=max_node_obs),split_var]


                        # No valid tree found
                        if(length(xcut_valid) == 0 ){

                                split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                                if(length(split_var_candidates)==0){
                                        return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                                }

                        } else {
                                # Sampling a x_cut_rule
                                # xcut_valid <- unique(xcut_valid)
                                split_var_sampled_rule <- sample(xcut_valid,size = 1)
                                good_tree_index <- 1
                        }

                }



        }


        if(split_var %in% cat_var){

                # Creating the vector of new train and test index
                left_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]==split_var_sampled_rule)]
                right_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]!=split_var_sampled_rule)]

                left_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]==split_var_sampled_rule)]
                right_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]!=split_var_sampled_rule)]

        } else {
                # Creating the vector of new train and test index
                left_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]<=split_var_sampled_rule)]
                right_train_id <- c_node$obs_train[which(x_train[c_node$obs_train,split_var]>split_var_sampled_rule)]

                left_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]<=split_var_sampled_rule)]
                right_test_id <- c_node$obs_test[which(x_test[c_node$obs_test,split_var]>split_var_sampled_rule)]
        }

        # Getting the left and the right
        new_left_name <- paste0("node_",c_node$left)
        new_right_name <- paste0("node_",c_node$right)

        # Creating a new left node and changing it
        old_left_node <- tree[[new_left_name]]
        new_left_node <- tree[[new_left_name]]
        new_left_node$obs_train <- left_train_id
        new_left_node$obs_test <- left_test_id
        new_left_node$var <- split_var
        new_left_node$var_split_rule <- split_var_sampled_rule
        new_left_node$parent <- c_node$index

        # Creating a new right node and changing it
        old_right_node <- tree[[new_right_name]]
        new_right_node <- tree[[new_right_name]]
        new_right_node$obs_train <- right_train_id
        new_right_node$obs_test <- right_test_id
        new_right_node$var <- split_var
        new_right_node$var_split_rule <- split_var_sampled_rule
        new_left_node$parent <- c_node$index

        # No valid tree
        if((length(new_left_node$obs_train) < node_min_size) || (length(new_right_node$obs_train)<node_min_size)){
                return(list(tree = tree, accepted = accepted))
        }

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood_gpbart(res_vec = res_vec,x_train = x_train,node = new_left_node,tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) +
                node_loglikelihood_gpbart(res_vec = res_vec,x_train = x_train, node = new_right_node,tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec, x_train = x_train ,node = old_left_node, tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec, x_train = x_train, node = old_right_node,tau = tau, tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp)

        log_acceptance <- tree_loglikeli

        # Accepting the tree ornot
        # if(stats::runif(n = 1)<=exp(log_acceptance)){
        if(stats::runif(n = 1)<=exp(log_acceptance)){

                # Maybe use append to make everything easier
                tree[[new_left_name]] <- new_left_node
                tree[[new_right_name]] <- new_right_node
                accepted <- TRUE
        }

        return(list(tree = tree, accepted = accepted))

}

# Changing the tree
change_rotation_gpbart <- function(res_vec,
                          tree,
                          x_train,
                          x_test,
                          xcut,
                          tau,
                          tau_mu,
                          alpha,
                          beta,
                          node_min_size,
                          nu,
                          phi_vector_p,
                          # Select the rotation
                          cov_gp,
                          rotation_variables){


        # ACC
        accepted <- FALSE

        # Getting the node
        nog_nodes <- get_nog(tree = tree)
        n_terminal_nodes <- length(get_terminals(tree = tree))

        if(length(tree)==1 || length(tree)==0){
                return(list(tree = tree, accepted = accepted))
        }

        # Sample a node to be pruned
        nog_nodes_index <- sample(1:length(nog_nodes),size = 1)

        if(nog_nodes_index!=0 & length(nog_nodes)!=0){
                c_node <- nog_nodes[[nog_nodes_index]]
        }

        # Case of just one split
        if( length(tree)==3){
                c_node <- tree[[1]] ## Getting the root node
        }


        if((length(nog_nodes) == 0) & length(tree)!=3){
                return(list(tree = tree, accepted = accepted))
        }

        good_tree_index <- 0

        while(good_tree_index==0){

                # Getting the name of the changed node
                split_var_pair <- sample(rotation_variables,size = 2)

                # Selecting a valid split
                split_var <- sample(split_var_pair,size = 1)

                # Selecting an angle to rotate my coordinates
                theta <- stats::runif(n = 1,min = 0,max = pi)


                # Creating the rotated coordinates
                rotated_x <- tcrossprod(A(theta), x_train[,split_var_pair])
                rownames(rotated_x) <- split_var_pair

                # Getting the rotation for the test
                rotated_x_test <- tcrossprod(A(theta), x_test[,split_var_pair])
                rownames(rotated_x_test) <- split_var_pair


                # Getting the min and maximum observed value within the terminal node
                min_node_obs <- sort(rotated_x[split_var,c_node$obs_train])[node_min_size]
                max_node_obs <- sort(rotated_x[split_var,c_node$obs_train])[length(rotated_x[split_var,c_node$obs_train])-node_min_size]


                # Getting the x_cut matrix rotated
                xcut_rotated <- tcrossprod(A(theta), xcut[,split_var_pair])
                rownames(xcut_rotated) <- split_var_pair

                # Getting the column from xcut
                xcut_valid <- xcut_rotated[split_var,which(xcut_rotated[split_var,]>=min_node_obs & xcut_rotated[split_var,]<=max_node_obs)]


                # Case of invalid max
                if((length(rotated_x[split_var,c_node$obs_train])-node_min_size)<1){
                        return(list(tree = tree, accepted = accepted))
                }



                # No valid tree found
                if(length(xcut_valid) == 0 ){

                        rotation_variables <-  rotation_variables[!(rotation_variables %in% split_var_pair)]

                        if(length(rotation_variables)==0 || length(rotation_variables)==1){
                                return(list(tree = tree, accepted = accepted)) # There are no valid candidates for this node
                        }



                } else {
                        good_tree_index <- 1
                }
        }

        # Sampling a x_cut_rule
        xcut_valid <- unique(xcut_valid)
        split_var_sampled_rule_rotation <- sample(xcut_valid,size = 1)

        # Creating the vector of new train and test index
        left_train_id <- c_node$obs_train[which(rotated_x[split_var,c_node$obs_train]<=split_var_sampled_rule_rotation)]
        right_train_id <- c_node$obs_train[which(rotated_x[split_var,c_node$obs_train]>split_var_sampled_rule_rotation)]

        left_test_id <- c_node$obs_test[which(rotated_x_test[split_var,c_node$obs_test]<=split_var_sampled_rule_rotation)]
        right_test_id <- c_node$obs_test[which(rotated_x_test[split_var,c_node$obs_test]>split_var_sampled_rule_rotation)]

        # Getting the left and the right
        new_left_name <- paste0("node_",c_node$left)
        new_right_name <- paste0("node_",c_node$right)

        # Creating a new left node and changing it
        old_left_node <- tree[[new_left_name]]
        new_left_node <- tree[[new_left_name]]
        new_left_node$obs_train <- left_train_id
        new_left_node$obs_test <- left_test_id
        new_left_node$var <- list(split_var_pair = split_var_pair, split_var = split_var, theta = theta )
        new_left_node$var_split_rule <- split_var_sampled_rule_rotation
        new_left_node$parent <- c_node$index

        # Creating a new right node and changing it
        old_right_node <- tree[[new_right_name]]
        new_right_node <- tree[[new_right_name]]
        new_right_node$obs_train <- right_train_id
        new_right_node$obs_test <- right_test_id
        new_right_node$var <- list(split_var_pair = split_var_pair, split_var = split_var, theta = theta )
        new_right_node$var_split_rule <- split_var_sampled_rule_rotation
        new_right_node$parent <- c_node$index

        # No valid tree
        if((length(new_left_node$obs_train) < node_min_size) || (length(new_right_node$obs_train)<node_min_size)){
                return(list(tree = tree, accepted = accepted))
        }

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood_gpbart(res_vec = res_vec,x_train = x_train,node = new_left_node,tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) +
                node_loglikelihood_gpbart(res_vec = res_vec,x_train = x_train, node = new_right_node,tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec, x_train = x_train ,node = old_left_node, tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp) -
                node_loglikelihood_gpbart(res_vec = res_vec, x_train = x_train, node = old_right_node,tau = tau, tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp)

        log_acceptance <- tree_loglikeli

        # Accepting the tree ornot
        if(stats::runif(n = 1)<=exp(log_acceptance)){

                # Maybe use append to make everything easier
                tree[[new_left_name]] <- new_left_node
                tree[[new_right_name]] <- new_right_node
                accepted <- TRUE

        }

        return(list(tree = tree, accepted = accepted))


}
