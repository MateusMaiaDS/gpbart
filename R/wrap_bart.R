## GP-Bart
#' @useDynLib gpbart
#' @importFrom Rcpp sourceCpp

# A function to retrive the number which are the factor columns
base_dummyVars <- function(df) {
        num_cols <- sapply(df, is.numeric)
        factor_cols <- sapply(df, is.factor)

        if((sum(num_cols)==0) & (sum(factor_cols)==0)){
             stop("Invalid data.frame")
        }
        return(list(continuousVars = names(df)[num_cols], facVars = names(df)[factor_cols]))
}

# Translation phi into the grid
translate_phi_matrix <- function(phi_mat) {
        phi_grid = c(0.1, 0.5, 1.0, 2.5, 3.0, 4.0, 5.0, 10.0, 25.0, 50.0)
        phi_mat <- apply(phi_mat, c(1, 2), function(x) phi_grid[x+1])
        return(phi_mat)
}

#' GP-BART: Gaussian Processes Bayesian Additive Regression Trees
#'
#'  GP-BART is an extension to to the Bayesian Additive Regression Trees (BART)
#'  
#' @param x_train Set of explanatory variables of the training data. It must be a `data.frame()`
#' @param y Response variable for the training data
#' @param x_test Set of explanatory variables of the test data. It must be a `data.frame()`
#' @param n_tree Number of Trees used in the GP-BART model.
#' @param node_min_size Node minimum size of observations within a terminal node
#' @param n_mcmc The number of MCMC iterations
#' @param n_burn The number of MCMC iterations to be trated as burn in
#' @param alpha Base parameter for the tree prior
#' @param beta Power parameter for the tree prior
#' @param df Degrees of freedom for the residual precision prior
#' @param sigquant The quantile of the residual precision prior;
#' @param kappa The number of prior standard deviations away from the from the range of the response.
#' @param tau Initial value for the residual precision
#' @param scale_bool A Boolean to choose if will be scaled or not.
#' @param nu Value for the GP precision. The default value is \eqn{\nu = 4\kappa^{2}T}
#' @param rand_tau_init A Boolean to let the initial value be initialised or not. The default is \code{TRUE}.
#' @param verbose Verbosity flag for printing progress. The default is \code{TRUE}
#'
#' @export
gpbart <- function(x_train,
                  y,
                  x_test,
                  n_tree = 20,
                  node_min_size = 2,
                  n_mcmc = 3500,
                  n_burn = 1500,
                  alpha = 0.95,
                  beta = 2,
                  df = 3,
                  sigquant = 0.9,
                  kappa = 2,
                  tau = 100,
                  scale_bool = TRUE,
                  nu = 1,
                  rand_tau_init = TRUE,
                  verbose = TRUE
                  ) {


     # Test parameters settled as false
     stump <- FALSE
     sample_phi <- TRUE
     sample <- TRUE
     no_rotation_bool = FALSE
     only_rotation_bool = TRUE
     # ===========
     
     # Verifying if x_train and x_test are matrices
     if(!is.data.frame(x_train) || !is.data.frame(x_test)){
          stop("Insert valid data.frame for both data and xnew.")
     }

     if(no_rotation_bool & only_rotation_bool) {
             stop("Insert valid rotation boolean setting")
     }

     # Getting the valid
     dummy_x <- base_dummyVars(x_train)

     # Create a data.frame aux

     # Create a list
     if(length(dummy_x$facVars)!=0){
             for(i in 1:length(dummy_x$facVars)){
                     # See if the levels of the test and train matches
                     if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
                        stop("Levels differs")
                        levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
                     }
                     df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
                     formula_aux <- stats::aggregate(y~x,df_aux,mean)
                     formula_aux$y <- rank(formula_aux$y)
                     x_train[[dummy_x$facVars[i]]] <- factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y))
                     levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[i]]])
                     # Doing the same for the test set
                     x_train[[dummy_x$facVars[i]]] <- as.numeric(x_train[[dummy_x$facVars[[i]]]])-1
                     x_test[[dummy_x$facVars[i]]] <- as.numeric(x_test[[dummy_x$facVars[[i]]]])-1
             }
     }
     x_train_scale <- as.matrix(x_train)
     x_test_scale <- as.matrix(x_test)

     # Scaling x
     x_min <- apply(as.matrix(x_train_scale),2,min)
     x_max <- apply(as.matrix(x_train_scale),2,max)

     # Storing the original
     x_train_original <- x_train
     x_test_original <- x_test


     # Normalising all the columns
     for(i in 1:ncol(x_train)){
             x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
             x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
     }

     # Selecting the variables to be only used into the GP's
     if(length(dummy_x$facVars)!=0){
        x_train_scale_gp <- x_train_scale[,dummy_x$continuousVars, drop = FALSE]
        x_test_scale_gp <- x_test_scale[,dummy_x$continuousVars, drop = FALSE]
     } else {
        x_train_scale_gp <- x_train_scale
        x_test_scale_gp <- x_test_scale
     }

     # Checking if
     if(!(length(dummy_x$facVars)==0) & (identical(x_train_scale,x_train_scale_gp))){
             stop("No numerical variables to be used here. Enter a new dataset.")
     }

     # Getting some y values
     if(!is.vector(y)){
               stop("The y should be a numeric vector.")
     }

     # Scaling the y
     min_y <- min(y)
     max_y <- max(y)

     # Getting the min and max for each column
     min_x <- apply(x_train_scale,2,min)
     max_x <- apply(x_train_scale, 2, max)

     # Scaling "y"
     if(scale_bool){
        y_scale <- normalize_bart(y = y,a = min_y,b = max_y)
        nu <- tau_mu <- (8*n_tree*(kappa^2))

     } else {
        y_scale <- y
        nu <- tau_mu <- (8*n_tree*(kappa^2))/((max_y-min_y)^2)
     }

     # Getting the naive sigma value
     nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

     # Calculating tau hyperparam
     a_tau <- df/2

     # Calculating lambda
     qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
     lambda <- (nsigma*nsigma*qchi)/df
     d_tau <- (lambda*df)/2


     # Call the bart function
     # tau_init <- tau
     tau_init <- nsigma^(-2)
     mu_init <- mean(c(unlist(y_scale)))

     if(rand_tau_init) {
             # tau_init <- stats::runif(n = 1,min = 0.1*tau_init,max = 4*tau_init)
             tau_init <- stats::rgamma(n = 1,shape =  a_tau,rate = d_tau)
     }

     # Creating the vector that stores all trees
     all_tree_post <- vector("list",length = round(n_mcmc-n_burn))

     # Generating the BART obj
     bart_obj <- cppgpbart(x_train_scale,
                          x_train_scale_gp,
                          y_scale,
                          x_test_scale,
                          x_test_scale_gp,
                          n_tree,
                          node_min_size,
                          n_mcmc,
                          n_burn,
                          tau_init,
                          mu_init,
                          tau_mu,
                          alpha,
                          beta,
                          a_tau,d_tau,
                          nu,
                          stump,
                          sample,
                          sample_phi,
                          verbose,
                          no_rotation_bool,
                          only_rotation_bool)


     if(scale_bool){
             # Tidying up the posterior elements
             y_train_post <- unnormalize_bart(z = bart_obj[[1]],a = min_y,b = max_y)
             y_test_post <- unnormalize_bart(z = bart_obj[[2]],a = min_y,b = max_y)
             tau_post <- bart_obj[[3]]/((max_y-min_y)^2)
             # y_train_sd_post <- unnormalize_bart(z = sqrt(bart_obj[[8]]),a = min_y,b = max_y)
             # y_test_sd_post <- unnormalize_bart(z = sqrt(bart_obj[[9]]),a = min_y,b = max_y)
             y_train_sd_post <- sqrt(bart_obj[[8]])*((max_y-min_y))
             y_test_sd_post <- sqrt(bart_obj[[9]])*((max_y-min_y))
             tau_post_all <- bart_obj[[10]]/((max_y-min_y)^2)
             tau_init <- tau_init/((max_y-min_y)^2)


             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  unnormalize_bart(z = bart_obj[[4]][,,i],a = min_y,b = max_y)
             }
     } else {
             y_train_post <- bart_obj[[1]]
             y_test_post <- bart_obj[[2]]
             y_train_sd_post <- sqrt(bart_obj[[8]])
             y_test_sd_post <- sqrt(bart_obj[[9]])
             tau_post <- bart_obj[[3]]
             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  bart_obj[[4]][,,i]
             }
             tau_post_all <- bart_obj[[10]]

     }

     # Return the list with all objects and parameters
     return(list(y_hat = y_train_post,
                 y_hat_test = y_test_post,
                 tau_post = tau_post,
                 tau_post_all = tau_post_all,
                 all_tree_post = all_tree_post,
                 all_phi_post = translate_phi_matrix(bart_obj[[5]]),
                 prior = list(n_tree = n_tree,
                              alpha = alpha,
                              beta = beta,
                              tau_mu = tau_mu,
                              a_tau = a_tau,
                              d_tau = d_tau,
                              tau_init = tau_init),
                 mcmc = list(n_mcmc = n_mcmc,
                             n_burn = n_burn),
                 data = list(x_train = x_train,
                             y = y,
                             x_test = x_test)))
}

