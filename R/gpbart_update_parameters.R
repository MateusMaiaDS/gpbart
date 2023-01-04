update_phi_gpbart <- function(tree,
                              x_train,
                              res_vec,
                              nu,
                              phi_vector_p,
                              tau,
                              tau_mu,
                              cov_gp,
                              prior_phi,
                              proposal_phi){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        for(i in 1:length(phi_vector_p)){
                # phi_proposal <- sample(x = 1/(2*pi*up_crossings),size = 1)
                # phi_proposal <- stats::runif(n = 1,min = (3/4)*phi_vector_p[i],max = (4/3)*phi_vector_p[i])

                # Setting a proposal given by the list element "proposal phi"
                if(proposal_phi[["proposal_mode"]]=="discrete_grid"){
                        if(is.null(proposal_phi[["grid"]])){
                                phi_proposal <- sample(c(0.1,seq(0,10,by=0.5)[-1],seq(10,20,by = 1), 25,30,50,75,100,125),size = 1)
                        } else {
                                phi_proposal <- sample(proposal_phi[["grid"]],size = 1)
                        }
                } else if(proposal_phi[["proposal_mode"]]=="sliding_window"){
                        phi_proposal <- stats::runif(n = 1,min = (3/4)*phi_vector_p[i],max = (4/3)*phi_vector_p[i])
                } else {
                        stop("Insert a valid proposal_mode")
                }
                # phi_proposal <- sample(exp(seq(log(0.05),log(100),length.out = 50)),size = 1)
                # phi_proposal <- stats::runif(0,100,n = 1)

                new_phi_vector_p <- phi_vector_p
                new_phi_vector_p[i] <- phi_proposal

                old_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train[,cov_gp, drop = FALSE],
                                                                                          tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p,gp_variables = cov_gp)}))
                new_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train[,cov_gp, drop = FALSE],
                                                                                                     tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = new_phi_vector_p,gp_variables = cov_gp)}))


                # Transition loglikelihood
                # transition_log <- stats::dunif(x = phi_vector_p[i],min = (3/4)*phi_proposal,max = (4/3)*phi_proposal,log = TRUE) - stats::dunif(x = phi_proposal,min = (3/4)*phi_vector_p[i],max = (4/3)*phi_vector_p[i],log = TRUE)

                if(proposal_phi[["proposal_mode"]] == "sliding_window"){
                        transition_log <- stats::dunif(x = phi_vector_p[i],min = (3/4)*phi_proposal,max = (4/3)*phi_proposal,log = TRUE) - stats::dunif(x = phi_proposal,min = (3/4)*phi_vector_p[i],max = (4/3)*phi_vector_p[i],log = TRUE)
                }

                # Calculating the prior log value
                if(is.null(prior_phi[["type"]])){
                         prior_log <- log(0.7*stats::dgamma(x = phi_proposal,shape = 5000,rate = 100)+0.3*stats::dgamma(x = phi_proposal,shape = 3,rate = 2.5))-log(0.7*stats::dgamma(x = phi_vector_p[i],shape = 5000,rate = 100)+0.3*stats::dgamma(x = phi_vector_p[i],shape = 3,rate = 2.5))
                } else if(prior_phi[["type"]]=="gamma_mixture"){
                        if(any(is.null(prior_phi$prob_1),is.null(prior_phi$prob_2),is.null(prior_phi$shape_1),is.null(prior_phi$shape_2),is.null(prior_phi$rate_1),is.null(prior_phi$rate_2)) ){
                                stop("Insert valid prior parameters")
                        } else {
                                prior_log <- log(prior_phi$prob_1*stats::dgamma(x = phi_proposal,shape = prior_phi$shape_1,rate = prior_phi$rate_1)+prior_phi$prob_2*stats::dgamma(x = phi_proposal,shape = prior_phi$shape_2,rate = prior_phi$rate_2))-log(prior_phi$prob_1*stats::dgamma(x = phi_vector_p[i],shape = prior_phi$shape_1,rate = prior_phi$rate_1)+prior_phi$prob_2*stats::dgamma(x = phi_vector_p[i],shape = prior_phi$shape_2,rate = prior_phi$rate_2))
                        }
                } else {
                        stop("Insert a valid list of prior_phi")
                }


                # Calculating acceptance
                if(proposal_phi[["proposal_mode"]]=="sliding_window"){
                        acceptance <- exp(new_log_like-old_log_like + prior_log + transition_log ) #+ stats::dgamma(x = phi_proposal,shape = 5,rate = 1,log = TRUE) - stats::dgamma(x = phi_vector_p[i],shape = 5,rate = 1,log = TRUE) )
                } else {
                        acceptance <- exp(new_log_like-old_log_like + prior_log ) #+ stats::dgamma(x = phi_proposal,shape = 5,rate = 1,log = TRUE) - stats::dgamma(x = phi_vector_p[i],shape = 5,rate = 1,log = TRUE) )

                }
                acceptance <- exp(new_log_like-old_log_like)

                if(stats::runif(n = 1,min = 0,max = 1)<=acceptance){
                        phi_vector_p <- new_phi_vector_p
                }
        }

        return(phi_vector_p)
}


update_nu_gpbart <- function( tree,
                              x_train,
                              res_vec,
                              nu,
                              phi_vector_p,
                              tau,
                              tau_mu){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        nu_proposal <- stats::runif(n = 1,min = 1,max = 100)

        old_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train,
                                                                                             tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p)}))
        new_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec, x_train = x_train,
                                                                                             tau = tau,tau_mu = tau_mu,nu = nu_proposal,phi_vector = phi_vector_p)}))

        # Calculating acceptance
        acceptance <- exp(new_log_like-old_log_like)

        if(stats::runif(n = 1,min = 0,max = 1)<=acceptance){
                return(nu_proposal)
        } else {
                return(nu)
        }

}



# Update mu node
update_mu_node <- function(node,
                           x_train,
                           x_test,
                           res_vec, # The complete res_vector
                           tau,
                           tau_mu,
                           nu,
                           phi_vector_p,
                           gp_variables){

        # Calculating the x_train from node
        x_train_node <- x_train[node$obs_train,gp_variables,drop = FALSE]
        res_node <- res_vec[node$obs_train]
        # Calculating the v factor from equation 11
        distance_sq_matrix <- symm_distance_matrix(m1 = x_train_node,phi_vector = phi_vector_p)
        omega_plus_tau_diag <- (nu^(-1))*exp(-distance_sq_matrix)+diag(1/tau,nrow = nrow(distance_sq_matrix))
        inv_omega_plus_tau <- chol2inv(chol(omega_plus_tau_diag))
        # Calculating S
        S <- sum(inv_omega_plus_tau)+tau_mu

        # Calculating the mean
        mu_mean <- (S^(-1))*crossprod(res_node,inv_omega_plus_tau)
        mu_var <- S^(-1)

        node$mu <- stats::rnorm(n = 1,mean = mu_mean,sd = sqrt(mu_var))

        return(node)
}


update_mu_gpbart <- function(tree,
                             x_train,
                             res_vec,
                             nu,
                             phi_vector_p,
                             tau,
                             tau_mu,
                             cov_gp){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        new_t_nodes_tree <-lapply(t_nodes, function(node){ update_mu_node(node = node,x_train = x_train,res_vec = res_vec,
                                                                          tau = tau,tau_mu = tau_mu,nu = nu,phi_vector_p = phi_vector_p,
                                                                          gp_variables = cov_gp)})

        # Updating all nodes
        tree[names(t_nodes)] <- new_t_nodes_tree

        return(tree)

}

# Update mu node
update_g_node <- function(node,
                           x_train,
                           x_test,
                           res_vec, # The complete res_vector
                           tau,
                           tau_mu,
                           nu,
                           phi_vector_p,
                           test_only = FALSE,
                           gp_variables){

        # Calculating the x_train from node
        x_train_node <- x_train[node$obs_train,gp_variables,drop = FALSE]
        x_test_node <- x_test[node$obs_test,gp_variables,drop = FALSE]

        # Creating the vetors
        g_sample <- numeric(nrow(x_train_node))
        g_sample_test <- numeric(nrow(x_test_node))

        res_node <- res_vec[node$obs_train]
        # Calculating the v factor from equation 11
        distance_sq_matrix <- symm_distance_matrix(m1 = x_train_node,phi_vector = phi_vector_p)
        omega <- (nu^(-1))*exp(-distance_sq_matrix)
        omega_plus_tau_diag <- omega+ diag(1/tau,nrow = nrow(distance_sq_matrix))
        inv_omega_plus_tau <- chol2inv(chol(omega_plus_tau_diag))



        # Calculating S
        S <- sum(inv_omega_plus_tau)+tau_mu

        if(!test_only){
                # Calculating the mean
                g_mean <- node$mu + crossprod(omega,crossprod(inv_omega_plus_tau,(res_node-node$mu)))
                g_var <- omega - crossprod(omega,crossprod(inv_omega_plus_tau,omega))
        }

        # Diagnostics
        # plot(res_node,g_mean)

        # Calculating test quantities
        distance_sq_matrix_test_star <- distance_matrix(m1 = x_train_node,m2 = x_test_node,phi_vector = phi_vector_p)
        distance_sq_matrix_test_star_star <- symm_distance_matrix(m1 = x_test_node,phi_vector = phi_vector_p)

        omega_star <- (nu^(-1))*exp(-distance_sq_matrix_test_star)
        omega_star_star <- (nu^(-1))*exp(-distance_sq_matrix_test_star_star)


        g_test_mean <- node$mu+crossprod(omega_star,crossprod(inv_omega_plus_tau,res_node-node$mu))
        g_test_var <- omega_star_star - crossprod(omega_star,crossprod(inv_omega_plus_tau,omega_star))
        # ====

        # Returning g sample
        if(!test_only){
                g_sample <- mvnfast::rmvn(n = 1,mu = g_mean,sigma = g_var+diag(1e-8,nrow = nrow(g_var)))
        }

        g_sample_test <- mvnfast::rmvn(n = 1,mu = g_test_mean,sigma = g_test_var+diag(1e-6,nrow = nrow(g_test_var)))

        return(list(train_sample = g_sample, test_sample = g_sample_test))
}


update_g_gpbart <- function(tree,
                            x_train,
                            x_test,
                            res_vec,
                            tau,
                            tau_mu,
                            nu,
                            phi_vector_p,
                            cov_gp){

        # new g
        g_sample_train <- numeric(nrow(x_train))
        g_sample_test <- numeric(nrow(x_test))
        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        # Iterate over terminal nodes
        for(i in 1:length(t_nodes)){
               g_aux <- update_g_node(node = t_nodes[[i]],x_train = x_train,x_test = x_test,
                                      res_vec = res_vec,tau = tau,
                                      tau_mu = tau_mu,nu = nu,phi_vector_p = phi_vector_p,gp_variables = cov_gp)

               g_sample_train[t_nodes[[i]]$obs_train] <- g_aux$train_sample
               g_sample_test[t_nodes[[i]]$obs_test] <- g_aux$test_sample
        }


        if(any(is.na(g_sample_train)) || any(is.na(g_sample_test))){
                stop("error somehwere in the prediction")
        }
        return(list(g_sample = g_sample_train, g_sample_test = g_sample_test))
}

# Tree permutation matrix
s_permutation <- function(tree,
                          x_train){
        t_nodes <- get_terminals(tree)
        mu_vec <- numeric(nrow(x_train))
        permutation_matrix <- matrix(0, nrow = nrow(x_train) , ncol = nrow(x_train))

        # Get train obs index
        all_index <- do.call(c,lapply(t_nodes,function(x){x$obs_train}))
        all_mu <- do.call(c,lapply(t_nodes, function(x){rep(x$mu,length(x$obs_train))}))
        for(i in 1:nrow(permutation_matrix)){
                permutation_matrix[i,all_index[i]] <- 1
                mu_vec[all_index[i]] <- all_mu[i]
        }

        return(list(permutation_matrix = permutation_matrix,
                    mu_vec = mu_vec))

}


# Verification of the thing that I want to do
# node_one_index <- tree$node_1$obs_train
# permutation_test <- s_permutation(tree = current_trees[[1]],x_train = x_train)
# cov_node_one <- symm_distance_matrix(m1 = x_train[node_one_index,,drop = FALSE],phi_vector = rep(1,5))
# cov_node_one_perm <- permutation_test%*%symm_distance_matrix(m1 = x_train,phi_vector = rep(1,5))%*%t(permutation_test)
# cov_node_one==cov_node_one_perm[1:24,1:24]

# NEED TO DO A TEST TO SEE IF MY ORIGINAL DISTANCE MATRIX WILL BE RETURNED BY THE PRODUCT OF THE PERMUTATION MATRIX


# Updating a single nu
update_single_nu <- function(current_trees,
                             phi_matrix,
                             x_train,
                             current_nu,
                             y_train,
                             tau,
                             K_bart = 2){


        # Getting n_tree
        n_tree <- length(current_trees)

        # Getting the current_nu
        proposal_nu <- stats::runif(n = 1,min = 0.5*current_nu,max = 1.5*current_nu)

        # Creating the list of terminal nodes
        list_terminal_nodes_cov <- vector("list",length = length(current_trees))
        list_mus <- vector("list",length = length(current_trees))
        list_distance_matrix <- vector("list",length = length(current_trees))
        # list_terminal_nodes_cov <- list()
        list_mus <- list()

        for(t in 1:length(current_trees)){

                # Permutation and nu elements
                perm_aux <- s_permutation(tree = current_trees[[t]],x_train = x_train)
                list_terminal_nodes_cov[[t]] <- perm_aux$permutation_matrix
                list_mus[[t]] <- perm_aux$mu_vec
                list_distance_matrix[[t]] <- symm_distance_matrix(m1 = x_train,phi_vector = phi_matrix[t,])
        }

        # Creating the big distance block matrix
        big_block_permutation <- do.call(cbind,list_terminal_nodes_cov)
        big_block_distance <- Matrix::bdiag(list_distance_matrix)
        big_long_mean_matrix <- do.call(c,list_mus)

        # Calcualting the dmn.mean
        y_post_mean <- big_block_permutation%*%big_long_mean_matrix
        y_post_cov <- big_block_permutation%*%big_block_distance%*%t(big_block_permutation)
        kernel_post_cov <- kernel_function(as.matrix(y_post_cov),nu = 1)
        old_nu_log <- mvnfast::dmvn(X = y_train,mu = as.matrix(y_post_mean),sigma = (current_nu^(-1))*kernel_post_cov+diag(tau^-1,nrow = nrow(x_train)),log = TRUE)
        new_nu_log <- mvnfast::dmvn(X = y_train,mu = as.matrix(y_post_mean),sigma = (proposal_nu^(-1))*kernel_post_cov+diag(tau^-1,nrow = nrow(x_train)),log = TRUE)

        # Addin gthe prior term
        prior_old <- stats::dgamma(x = current_nu,shape = 4*(K_bart^2)*n_tree*0.1,rate = 0.1,log = TRUE)
        prior_new <- stats::dgamma(x = proposal_nu,shape = 4*(K_bart^2)*n_tree*0.1,rate = 0.1,log = TRUE)

        acceptance <- new_nu_log-old_nu_log+prior_new-prior_old

        if(stats::runif(n = 1) < exp(acceptance)){
                return(proposal_nu)
        } else {
                return(current_nu)
        }
}

