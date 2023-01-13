## GP-Bart
#' @useDynLib gpbart
#' @importFrom Rcpp sourceCpp
#'
# Building the GP-BART function
#' @export
#'
gp_bart <- function(x_train,
                 y_train,
                 x_test,
                 n_tree = 1,
                 n_mcmc = 2000,
                 n_burn = 500,
                 node_min_size = 1,
                 tau = 1,
                 alpha = 0.5, beta = 2,
                 df = 3, sigquant = 0.9,
                 numcut = 100,
                 scale_boolean = TRUE,
                 K_bart = 2,
                 bart_boolean = FALSE,
                 keeptrees = FALSE,
                 bart_warmup = 250,
                 x_scale = TRUE,
                 update_nu = FALSE,
                 store_verb = TRUE,
                 gp_variables_,
                 rotation_variables_ = NULL,
                 cat_var_ = NULL,
                 prior_phi_ = list(type = NULL,prob_1 = NULL, prob_2 = NULL, shape_1 = NULL, shape_2 = NULL,
                                   rate_1 = NULL, rate_2 = NULL),
                 proposal_phi_ = list(proposal_mode = "discrete_grid", grid = NULL)){


        # Verifying if x_train and x_test are matrices
        if(!is.matrix(x_train) || !is.matrix(x_test)){
                x_train <- as.matrix(x_train)
                x_test <- as.matrix(x_test)
        }

        if(is.null(colnames(x_train)) || is.null(colnames(x_test))) {
                stop("Insert valid NAMED matrix for x_train or x_test")
        }

        # Saving a_min and b_max
        a_min <- NULL
        b_max <- NULL

        a_min <- min(y_train)
        b_max <- max(y_train)


        # Scaling the x
        if(x_scale){

                # Getting the x values
                x_min <- apply(rbind(x_train),2,min)
                x_max <- apply(rbind(x_train),2,max)

                # Getting the training original
                x_train_original <- x_train
                x_test_original <- x_test

                # # Creating the xscale_train/test
                xscale_train <- x_train
                xscale_test <- x_test

                # Normalize all the columns
                for(i in 1:ncol(x_train)){
                        xscale_train[,i] <- normalize_covariates_bart(y = x_train[,i],a = x_min[i],b = x_max[i])
                        xscale_test[,i] <- normalize_covariates_bart(y = x_test[,i],a = x_min[i],b = x_max[i])
                }


                # The result of scaling
                x_train <- as.matrix(xscale_train)
                x_test <- as.matrix(xscale_test)

        } else{
                x_train_original <- x_train
                x_test_original <- x_test
        }

        # Cut matrix
        xcut <- matrix(NA,ncol = ncol(x_train),nrow = numcut)

        # Getting possible x values
        for(j in 1:ncol(x_train)){
                xs <- stats::quantile(x_train[ , j], type=7,
                                      probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]

                xcut[,j] <-xs
        }# Error of the matrix

        # Adding a small pertubation over the xcut var to not match the current observed values
        xcut <- xcut + 1e-8

        if(is.null(colnames(x_train)) || is.null(colnames(x_test)) ) {
                stop("Insert a valid NAMED matrix")
        }

        # Giving the names for xcut
        colnames(xcut) <- colnames(x_train)

        if(!is.vector(y_train)) {
                stop("Insert a valid y vector")
        }

        # Scale values
        if(scale_boolean) {
                # Normalizing y
                y_scale <- normalize_bart(y = y_train)

                # Calculating \tau_{\mu} based on the scale of y
                tau_mu <- (8 * n_tree * (K_bart^2))

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

                cat("a_tau is " , round(a_tau,digits = 3)," \n")
                cat("d_tau is " , round(d_tau,digits = 3)," \n")

        } else {

                # Not scaling the y
                y_scale <- y_train

                # Calculating \tau_{\mu} based on the scale of y
                # Need to change this value in case of non-scaling
                tau_mu <- (8 * n_tree * (K_bart^2))/((b_max-a_min)^2)
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

                cat("a_tau is " , round(a_tau,digits = 3)," \n")
                cat("d_tau is " , round(d_tau,digits = 3)," \n")

        }

        # Defining other quantities
        n_train <- nrow(x_train)
        n_test <- nrow(x_test)
        n_post <- n_mcmc-n_burn

        # Getting the y_hat for train and test
        y_train_hat_post <- matrix(0, ncol = n_train,nrow = n_post)
        y_test_hat_post <- matrix(0, ncol = n_test,nrow = n_post)
        curr <- 0
        y_train_hat_trees <- matrix(rep(y_scale/n_tree,n_tree), nrow = n_tree, ncol = n_train)
        y_test_hat_trees <- matrix(0, nrow = n_tree, ncol = n_test)


        # Storing the current partial residuals
        current_partial_residuals_matrix <- matrix(NA,nrow = n_tree, ncol = n_train,byrow = TRUE)

        current_partial_residuals_list <- vector("list",length = n_post)
        all_tree_prediction <- vector("list", length = n_post)

        # Initialising values for phi_vec, and nu
        phi_vec_matrix <- matrix(1, nrow = n_tree,ncol = ncol(x_train[,gp_variables_, drop = FALSE]))
        phi_post <- list(n_post)
        nu <- 4*(K_bart^2)*n_tree
        nu_post <- numeric(n_post)

        tau_post <- numeric(n_post)
        post_trees <- vector("list",length = n_tree)

        # Getting initial trees
        current_trees <- vector("list",length = n_tree)

        for(i in 1:n_tree){
                current_trees[[i]] <- new_tree(x_train = x_train,x_test = x_test)
        }

        # Creating a manual progress bart
        progress_bart_limits <- round(seq(1,n_mcmc,length.out=100))

        # Checking if the correct number of bart_warmup iterations were selected
        if(n_burn < bart_warmup){
                warning(" BART iterations are being used to predict the model")
        }


        # Tree verb ratio acceptance
        df_verb <- as.data.frame(matrix(nrow = 0,ncol = 2,dimnames = list(NULL,c("verb","accepted_verb"))))

        for(i in 1:n_mcmc){

                # Small progress bar
                if(i %in% progress_bart_limits){
                        cat("|")
                }

                for(t in 1:n_tree){

                        # Adding BART warmup iterations
                        if(i > bart_warmup) {
                                bart_boolean <- FALSE
                        }

                        # USING BART BOOLEAN OR NOT
                        if(bart_boolean){
                                partial_residuals <- y_scale - colSums(y_train_hat_trees[-t,,drop = FALSE])

                                # Selecting one verb
                                verb <- sample(c("grow","grow_rotation","prune","change","change_rotation"), prob = c(0.15,0.15,0.3,0.2,0.2),size = 1)

                                if(length(current_trees[[t]])==1){
                                        verb <- sample(c("grow","grow_rotation"),size = 1)
                                }

                                # For the case where the name is null
                                if( (verb == "grow_rotation") & is.null(rotation_variables_) ){
                                        verb <- "grow"
                                }

                                if( (verb == "change_rotation") & is.null(rotation_variables_) ){
                                        verb <- "change"
                                }



                                # Selecting a new tree
                                current_trees[[t]]  <- if(verb=="grow"){
                                        grow(res_vec = partial_residuals,tree = current_trees[[t]],
                                             x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                             tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,cat_var = cat_var_)
                                } else if(verb == "grow_rotation") {
                                        grow_rotation(res_vec = partial_residuals,tree = current_trees[[t]],
                                                      x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                      tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                      rotation_variables = rotation_variables_)
                                } else if(verb=="prune"){
                                        prune(tree = current_trees[[t]],res_vec = partial_residuals,
                                              tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta)
                                } else if(verb=="change"){
                                        change(res_vec = partial_residuals,tree = current_trees[[t]],
                                               x_train = x_train,x_test = x_test,xcut = xcut,
                                               tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                               cat_var = cat_var_)
                                } else if(verb == "change_rotation") {
                                        change_rotation(res_vec = partial_residuals,tree = current_trees[[t]],
                                                        x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                        tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                        rotation_variables = rotation_variables_)
                                } else {
                                        stop("No valid verb for BART approach")
                                }

                                # Updating mu
                                current_trees[[t]] <- update_mu(tree = current_trees[[t]],
                                                                partial_residuals = partial_residuals,
                                                                tau = tau,tau_mu = tau_mu)

                                # Prediction aux
                                pred_obj <- getPrediction(tree = current_trees[[t]],x_train = x_train, x_test = x_test)

                                y_train_hat_trees[t,] <- pred_obj$train_pred
                                y_test_hat_trees[t,] <- pred_obj$test_pred

                                # =================================
                        } else {# CALCULATING THE GP-ITERATIONS
                                # =================================

                                # Calculating partial residuals
                                partial_residuals <- y_scale - colSums(y_train_hat_trees[-t,,drop = FALSE])

                                # Storing current partial
                                current_partial_residuals_matrix[t,] <- partial_residuals

                                verb <- sample(x = c("grow","grow_rotate","prune","change","change_rotate"),size = 1,prob = c(0.15,0.15,0.3,0.2,0.2))

                                if(length(current_trees[[t]])==1){
                                        verb <- sample(c("grow","grow_rotate"),size = 1)
                                }


                                # Restricting the case where there are no rotation
                                if(verb == "grow_rotate" & is.null(rotation_variables_)){
                                        verb <- "grow"
                                }

                                if(verb == "change_rotate" & is.null(rotation_variables_)){
                                        verb <- "change"
                                }

                                # Selecting one verb movement
                                if(verb == "grow"){
                                        current_tree_aux <- grow_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                          x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                          tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                          nu = nu,phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_,
                                                                          cat_var = cat_var_)
                                } else if( verb == "grow_rotate"){
                                        current_tree_aux <- grow_rotation_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                                   x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                                   tau_mu = tau_mu, alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                                   nu = nu, phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_,
                                                                                   rotation_variables = rotation_variables_)

                                }else if( verb == "prune"){
                                        current_tree_aux <- prune_gpbart(res_vec = partial_residuals,
                                                                           x_train = x_train,
                                                                           tree = current_trees[[t]],
                                                                           tau = tau, tau_mu = tau_mu, alpha = alpha, beta = beta,
                                                                           nu = nu, phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_)
                                } else if( verb == "change"){
                                        current_tree_aux <- change_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                          x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                          tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                          nu = nu,phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_,
                                                                          cat_var = cat_var_)
                                } else if (verb == "change_rotate"){
                                        current_tree_aux <- change_rotation_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                                   x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                                   tau_mu = tau_mu, alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                                   nu = nu, phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_,
                                                                                   rotation_variables = rotation_variables_)
                                } else {
                                        stop("Error no valid-verb")
                                }


                                # Getting the tree
                                current_trees[[t]] <- current_tree_aux$tree

                                #  Choosing if store the verbs or not
                                if(store_verb){
                                        df_verb <- rbind(df_verb,data.frame( verb = verb,
                                                                             accepted_bool = current_tree_aux$accepted))
                                }
                                # Changing for the current tree
                                # Updating the phi
                                phi_vec_matrix[t,] <- update_phi_gpbart(tree = current_trees[[t]],x_train = x_train,res_vec = partial_residuals,
                                                          phi_vector_p = phi_vec_matrix[t,],nu = nu,tau = tau,tau_mu = tau_mu,cov_gp = gp_variables_,
                                                          proposal_phi = proposal_phi_,prior_phi = prior_phi_)

                                # Update the mu values
                                # current_trees[[t]] <- update_mu_gpbart(tree = current_trees[[t]],x_train = x_train,res_vec = partial_residuals,nu = nu,
                                #                                   phi_vector_p = phi_vec_matrix[t,],tau = tau,tau_mu = tau_mu,cov_gp = gp_variables_)

                                # This one is the most complicated, I need to update the predictions based on the tree structure
                                sample_g_aux <- update_g_gpbart(tree = current_trees[[t]],
                                                                x_train = x_train,
                                                                x_test = x_test,
                                                                res_vec = partial_residuals,
                                                                tau =  tau,
                                                                tau_mu = tau_mu,
                                                                nu = nu,
                                                                phi_vector_p = phi_vec_matrix[t,],cov_gp = gp_variables_)

                                y_train_hat_trees[t,] <- sample_g_aux$g_sample
                                y_test_hat_trees[t,] <- sample_g_aux$g_sample_test
                        }
                }

                # Storing tau and getting new tau
                tau <- update_tau(y = y_scale,y_hat = colSums(y_train_hat_trees),a_tau = a_tau,d_tau = d_tau)


                # Adding the option of updating nu or not
                if(update_nu & !(bart_boolean)){
                        nu <- update_single_nu(current_trees = current_trees,y_train = colSums(y_train_hat_trees),
                                               current_nu = nu,phi_matrix = phi_vec_matrix,x_train = x_train,tau = tau)
                }

                tau_post[i] <- tau

                # Storing the posterior elements
                if(i>n_burn){
                        curr = curr + 1
                        current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix
                        all_tree_prediction[[curr]] <- y_train_hat_trees
                        y_train_hat_post[curr,] <- colSums(y_train_hat_trees)
                        y_test_hat_post[curr,] <- colSums(y_test_hat_trees)
                        nu_post[curr] <- nu
                        phi_post[[curr]] <- phi_vec_matrix

                        # Saving the trees
                        if(keeptrees) {
                                post_trees[[curr]] <- current_trees
                        }

                        # Saving the posterior of the hyperparameters
                        phi_post[[curr]] <- phi_vec_matrix

                }


        }


        # Adjusting tau and y_hat for the scale factor
        if(scale_boolean){
                tau_post <- tau_post/((b_max-a_min)^2)
                y_train_hat_post <- unnormalize_bart(z = y_train_hat_post,a = a_min,b = b_max)
                y_test_hat_post <- unnormalize_bart(z = y_test_hat_post,a = a_min,b = b_max)
        }

        # Diagnostic
        # plot(y_train,colMeans(y_train_hat_post))
        # crossprod(y_train-colMeans(y_train_hat_post))

        # Returning the posterior objets
        if(keeptrees){
                post_obj <- list(tau_post = tau_post,
                     y_hat_post = y_train_hat_post,
                     y_test_hat_post = y_test_hat_post,
                     last_trees = current_trees,
                     data = list(x_train = x_train,
                                 y_train = y_train,
                                 x_test = x_test),
                     prior = list(tau_mu = tau_mu,
                                  a_tau = a_tau,
                                  d_tau = d_tau,
                                  scale_boolean = scale_boolean,
                                  K_bart = 2,
                                  bart_boolean = bart_boolean,
                                  bart_warmup = bart_warmup,
                                  x_scale = x_scale,
                                  update_nu = update_nu,
                                  gp_variables_ = gp_variables_,
                                  rotation_variables_ = rotation_variables_,
                                  n_mcmc = n_mcmc,
                                  n_burn = n_burn),
                     posterior = list(phi_post = phi_post,
                                      nu_post = nu_post,
                                      partial_residuals = current_partial_residuals_list,
                                      all_tree_prediction = all_tree_prediction,
                                      trees = post_trees,
                                      post_verbs = df_verb))
        } else {
                post_obj <- list(tau_post = tau_post,
                     y_hat_post = y_train_hat_post,
                     y_test_hat_post = y_test_hat_post,
                     last_trees = current_trees,
                     data = list(x_train = x_train,
                                 y_train = y_train,
                                 x_test = x_test),
                     prior = list(tau_mu = tau_mu,
                                  a_tau = a_tau,
                                  d_tau = d_tau,
                                  scale_boolean = scale_boolean,
                                  K_bart = 2,
                                  bart_boolean = bart_boolean,
                                  bart_warmup = bart_warmup,
                                  x_scale = x_scale,
                                  update_nu = update_nu,
                                  gp_variables_ = gp_variables_,
                                  rotation_variables_ = rotation_variables_,
                                  n_mcmc = n_mcmc,
                                  n_burn = n_burn),
                     posterior = list(phi_post = phi_post,
                                      nu_post = nu_post,
                                      partial_residuals = current_partial_residuals_list,
                                      all_tree_prediction = all_tree_prediction,
                                      post_verbs = df_verb))

        }

        return(post_obj)

}

# GP-BART predict ---  a function to predict given new observations
#' @export
#'
gp_bart_predict <- function(gpbart_mod_example,
                            x_new){

        if(!is.matrix(x_new) || !all(colnames(x_new)==colnames(gpbart_mod_example$data$x_train))){
                stop ("Insert a valid new X matrix")
        }

        # Storing the new trees
        new_trees <- gpbart_mod_example$posterior$trees

        y_hat_new_test <- matrix(NA,nrow = length(new_trees),ncol = nrow(x_new))
        y_hat_trees_new_test <- matrix(NA,nrow = length(new_trees[[1]]),ncol = nrow(x_new))

        # All trees prediction
        all_tree_hat_list <- vector("list",length = length(new_trees))

        # Returning tau to the scaled version
        gpbart_mod_example$tau_post <- gpbart_mod_example$tau_post*(max(gpbart_mod_example$data$y_train)-min(gpbart_mod_example$data$y_train))^2

        # Creating the new test set divisions
        for(i in 1:length(new_trees)){

                for(t in 1:length(gpbart_mod_example$posterior$trees[[1]])){

                        tree_aux <- new_trees[[i]][[t]]

                        for(leaf in 1:length(tree_aux)){

                                # Setting the new root
                                if(is.na(tree_aux[[leaf]]$var_split_rule)){
                                        tree_aux[[leaf]]$obs_test <- 1:nrow(x_new)
                                }

                                # Get the current node
                                current_node_obs <- tree_aux[[leaf]]$obs_test

                                # Updating the left node
                                if(!is.na(tree_aux[[leaf]]$left)){


                                        left_node <- tree_aux[[paste0("node_",tree_aux[[leaf]]$left)]]

                                        if(is.list(left_node$var)){

                                                # Rotate the the new observations
                                                rotated_x_new <- tcrossprod(A(left_node$var$theta),x_new[,left_node$var$split_var_pair])
                                                rownames(rotated_x_new) <- left_node$var$split_var_pair
                                                tree_aux[[paste0("node_",tree_aux[[leaf]]$left)]]$obs_test <- current_node_obs[rotated_x_new[left_node$var$split_var,current_node_obs]<left_node$var_split_rule]

                                        } else {

                                                tree_aux[[paste0("node_",tree_aux[[leaf]]$left)]]$obs_test <- current_node_obs[x_new[current_node_obs, left_node$var]<left_node$var_split_rule]
                                        }
                                }

                                # Updating the left node
                                if(!is.na(tree_aux[[leaf]]$right)){
                                        right_node <- tree_aux[[paste0("node_",tree_aux[[leaf]]$right)]]

                                        if(is.list(right_node$var)){

                                                rotated_x_new <- tcrossprod(A(right_node$var$theta),x_new[,right_node$var$split_var_pair])
                                                rownames(rotated_x_new) <- right_node$var$split_var_pair
                                                tree_aux[[paste0("node_",tree_aux[[leaf]]$right)]]$obs_test <- current_node_obs[rotated_x_new[right_node$var$split_var,current_node_obs]>=right_node$var_split_rule]

                                        } else {
                                                tree_aux[[paste0("node_",tree_aux[[leaf]]$right)]]$obs_test <- current_node_obs[x_new[current_node_obs, right_node$var]>=right_node$var_split_rule]
                                        }
                                }


                                # Predicting for the current node
                                if(tree_aux[[leaf]]$terminal==1){
                                       y_hat_trees_new_test[t,tree_aux[[leaf]]$obs_test] <-  update_g_node(node = tree_aux[[leaf]],
                                                      x_train = gpbart_mod_example$data$x_train,
                                                      x_test = x_new,
                                                      res_vec = gpbart_mod_example$posterior$partial_residuals[[i]][t,],
                                                      tau = gpbart_mod_example$tau_post[i],tau_mu = gpbart_mod_example$prior$tau_mu,
                                                      nu = gpbart_mod_example$posterior$nu_post[i],
                                                      phi_vector_p = gpbart_mod_example$posterior$phi_post[[i]][t,],test_only = TRUE)$test_sample
                                }

                        }
                } # Ending the iterations over trees
                all_tree_hat_list[[i]] <- if(gpbart_mod_example$prior$scale_boolean){
                        y_hat_trees_new_test
                } else {
                        y_hat_trees_new_test
                }
                y_hat_new_test[i,] <- colSums(y_hat_trees_new_test)
        }

        # Getting the GPBART mod
        if(gpbart_mod_example$prior$scale_boolean){
                return(list( all_trees = all_tree_hat_list,
                        y_hat_new_post = unnormalize_bart(z = y_hat_new_test,a = min(gpbart_mod_example$data$y_train),b = max(gpbart_mod_example$data$y_train))))
        } else {
                return(list( all_trees_post = all_tree_hat_list,
                             y_hat_new_post = y_hat_new_test))
        }

}

# x_test_m <- as.matrix(x_test)
# x_test_m
# gp_pred_test <- gp_bart_predict(gpbart_mod_example = gpbart_mod_example,x_new = x_test_m)
# plot(y_test,gp_pred_test$y_hat_new_post %>% colMeans())
#
# #Calculating a bart model
# bart_mod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)
# plot(bart_mod$yhat.test.mean,y_test)
#
# crossprod((y_test-bart_mod$yhat.test.mean))
# crossprod((y_test-gp_pred_test$y_hat_new_post %>% colMeans()))
