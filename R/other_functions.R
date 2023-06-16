# Normalize BART function (Same way ONLY THE COVARIATE NOW)
normalize_covariates_bart <- function(y, a = NULL, b = NULL) {

     # Defining the a and b
     if( is.null(a) & is.null(b)){
          a <- min(y)
          b <- max(y)
     }
     # This will normalize y between -0.5 and 0.5
     y  <- (y - a)/(b - a)
     return(y)
}


# Normalize BART function (Same way ONLY THE COVARIATE NOW)
normalize_bart <- function(y, a = NULL, b = NULL) {

     # Defining the a and b
     if( is.null(a) & is.null(b)){
          a <- min(y)
          b <- max(y)
     }
     # This will normalize y between -0.5 and 0.5
     y  <- (y - a)/(b - a) - 0.5
     return(y)
}

# Getting back to the original scale
unnormalize_bart <- function(z, a, b) {
     # Just getting back to the regular BART
     y <- (b - a) * (z + 0.5) + a
     return(y)
}


# Naive sigma_estimation
naive_sigma <- function(x,y){

     # Getting the valus from n and p
     n <- length(y)

     # Getting the value from p
     p <- ifelse(is.null(ncol(x)), 1, ncol(x))

     # Adjusting the df
     df <- data.frame(x,y)
     colnames(df)<- c(colnames(x),"y")

     # Naive lm_mod
     lm_mod <- stats::lm(formula = y ~ ., data =  df)

     # Getting sigma
     sigma <- summary(lm_mod)$sigma
     return(sigma)

}



# Function to create a vector of variables that being categorical will
#have the same code
recode_vars <- function(x_train, dummy_obj){

        vars <- numeric()
        j <- 0
        i <- 0
        c <- 1
        while(!is.na(colnames(x_train)[c])){
                if(colnames(x_train)[c] %in% dummy_obj$facVars){
                        curr_levels <- dummy_obj$lvls[[colnames(x_train)[c]]]
                        for(k in 1:length(curr_levels)){
                             i = i+1
                             vars[i] <- j
                        }
                } else {

                     i = i+1
                     vars[i] <- j
                }
                j = j+1
                c = c+1
        }

        return(vars)
}

#' RMSE
#' 
#' Calculating the Root Mean Squared Error 
#'  
#' @param x A numerical vector of observed variables
#' @param y A numerical vector of predicted variables
#'
#' @export
rmse <- function(x,y){
     return(sqrt(mean((y-x)^2)))
}

# Calculating CRPS from (https://arxiv.org/pdf/1709.04743.pdf)
crps <- function(y,means,sds){

     # scaling the observed y
     z <- (y-means)/sds

     crps_vector <- sds*(z*(2*stats::pnorm(q = z,mean = 0,sd = 1)-1) + 2*stats::dnorm(x = z,mean = 0,sd = 1) - 1/(sqrt(pi)) )

     return(list(CRPS = mean(crps_vector), crps = crps_vector))
}


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



