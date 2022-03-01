library(Rcpp)
sourceCpp(file = "gpbart/distmatrix_rcpp/src/dist_matrix.cpp")

# GP-function main
gp_main <- function(x_train, y_train, x_star, tau, phi, nu, distance_matrix_train) {
  
  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train,m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star) 
  
  
  # Getting the distance matrix from x_train and x_star (SPATIAL VALUES)
  # distance_matrix_K_star <- distance_matrix(m1 = x_train[,c("lat","lon")],m2 = x_star[,c("lat","lon")])
  # distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star[,c("lat","lon")])
  
  # Calculating the K elements from the covariance structure
  K_y <- (kernel_function(distance_matrix = distance_matrix_train,
                          nu = nu,
                          phi = phi) + diag(x = tau^-1,nrow = nrow(x_train)))
  K_star <- (kernel_function(distance_matrix = distance_matrix_K_star,
                             nu = nu, phi = phi))
  
  # ===
  # Commented because is not necessary to calculate it here
  # === 
  # K_star_star <- (kernel_function(distance_matrix = distance_matrix_K_star_star,
  #                                 nu = nu, phi = phi)) #+ diag(x = tau^-1,nrow = nrow(x_star)) )
  
  # Calculating \alpha
  L <- chol(K_y)
  alpha <- backsolve(L,backsolve(L,y_train, transpose = TRUE))
  mu_star <- crossprod(K_star, alpha)
  
  # Here the abs is because the smallest values that are coming from here are due to numerical approximations.
  # v <- backsolve(L,K_star,transpose = TRUE)
  # cov_star <-K_star_star - crossprod(v)
  cov_star <- c()
  
  # ===============#
  return(list(mu_pred = mu_star, cov = cov_star))
}


# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(distance_matrix, nu, phi) {
  
  # Calculating the square matrix
  kernel_matrix <- exp(-distance_matrix / (2 * phi^2)) / nu
  
  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0, nrow = dim(distance_matrix)[1],
                            ncol = dim(distance_matrix)[2])
  }
  # Getting the kernel matrix
  return(kernel_matrix)
}


# Testing the matrix
# matrix_test <- matrix( c(1,4,1,3,1,3),
#                        nrow=3)
# symm_distance_matrix(matrix_test) -> aux_dist
# aux_dist_default <- (as.matrix(dist(matrix_test)))^2
# aux_dist_double <- distance_matrix(matrix_test,matrix_test)
# aux_dist
# aux_dist_default
# aux_dist_double


# Test over kernel_function
# dist_x_test <- matrix(rnorm(5))
# dist_matrix_test <- as.matrix(dist(dist_x_test))
# omega <- kernel_function(distance_matrix = dist_matrix_test, nu = 1, phi = 0.0001)
# 
# # Testing the standard BART settings
# mu <- rep(0.2,5)
# tau <- 0.1
# residuals <- rnorm(5)
# diagonal_d <- diag(100000,5)
# 
# 
# # What is approximately equals to mu
# mean_r_star <- mu + crossprod(omega,solve((omega+diagonal_d),(residuals-mu))) # Return the \mu vec
# var_r_star <- tau*omega - crossprod(omega,solve((omega+diagonal_d),omega)) # Return a diagonal matrix with \tau inverse values
