# # Update nu
update_nu_bgp <- function(x_train,
                          y_train,
                          phi_vec,
                          nu,
                          tau){

        # Calculating the distance matrix
        sq_dist_matrix_x <- symm_distance_matrix(m1 = x_train,phi_vector = phi_vec)

        # Setting a proposal for \nu
        proposal_nu <- stats::runif(n = 1,min = 0,max = 100)

        # cat("Proposal nu equal to: ", round(proposal_nu,digits = 3)," ")

        # Getting the covariance matrix
        K_y <- kernel_function(squared_distance_matrix_phi = sq_dist_matrix_x,nu = nu)


        # Setting the cov_old
        cov_old <- K_y + diag(1/tau,nrow = nrow(K_y))

        # Setting the cov_new
        cov_new <- (proposal_nu^-1)*K_y + diag(1/tau, nrow = nrow(K_y))

        # Getting the loglikelihood for the old and for the new nu proposal
        log_like_old <- mvnfast::dmvn(X = y_train,mu = rep(0,nrow(K_y)),sigma = cov_old,log = TRUE)
        log_like_new <- mvnfast::dmvn(X = y_train,mu = rep(0,nrow(K_y)),sigma = cov_new,log = TRUE)

        # Setting a uniform prior for nu
        log_acceptance <- log_like_new - log_like_old

        if(stats::runif(n = 1,min = 0,max = 1) < exp(log_acceptance)){
                return(proposal_nu)
        } else {
                return(nu)
        }
}

# Updating phi ( phi has dimension 1xp)
update_phi_bgp <- function(x_train,
                       y_train,
                       phi_vec,
                       nu,
                       tau){


  # Creating a function to update phi
  for(i in 1:length(phi_vec)){
    phi_proposal <- stats::runif(n = 1.,min = 0,max = 10)
    new_phi_vec <- phi_vec
    new_phi_vec[i] <- phi_proposal
    new_dist_m1_x <- symm_distance_matrix(m1 = x_train,phi_vector = new_phi_vec)
    old_dist_m1_x <- symm_distance_matrix(m1 = x_train,phi_vector = phi_vec)


    new_cov <- kernel_function(squared_distance_matrix_phi = new_dist_m1_x,nu = nu)+diag(1/tau, nrow = nrow(new_dist_m1_x))
    old_cov <- kernel_function(squared_distance_matrix_phi = old_dist_m1_x,nu = nu)+diag(1/tau, nrow = nrow(new_dist_m1_x))


    old_loglike <- mvnfast::dmvn(X = y_train,mu = rep(0,nrow(old_cov)),sigma = old_cov,log = TRUE)
    new_loglike <- mvnfast::dmvn(X = y_train,mu = rep(0,nrow(new_cov)), sigma = new_cov,log = TRUE)

    log_acceptance <- new_loglike - old_loglike


    # Accepting the new vector
    if(stats::runif(n = 1,min = 0,max = 1) < log_acceptance) {
      phi_vec <- new_phi_vec
    }

  }
  return(phi_vec)

}

# Sampling the y_hat
sample_y_hat_bgp <- function(y_train,
                         x_train,
                         phi_vec,
                         nu,
                         tau){

  # Creating the distance matrix
  sq_dist_matrix <- symm_distance_matrix(m1 = x_train,phi_vector = phi_vec)
  K_y <- kernel_function(squared_distance_matrix_phi =  sq_dist_matrix,nu = nu)

  # Creating the posterior mean
  post_mean <- crossprod(K_y,solve(K_y+diag(1/tau,nrow = nrow(K_y)),y_train))
  post_var <- K_y - crossprod(K_y,solve(K_y+diag(1/tau,nrow = nrow(K_y)),K_y))

  y_sample <- mvnfast::rmvn(n = 1,mu = post_mean,sigma = post_var+diag(1e-12,nrow = nrow(post_var)))

  return(y_sample)

}
