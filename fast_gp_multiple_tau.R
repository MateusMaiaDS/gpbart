library(Rcpp)
sourceCpp(file = "dist_matrix.cpp")

# GP-function main
gp_main <- function(x_train, y_train, x_star, tau, phi, nu, distance_matrix_train) {
  
  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train,m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star) 
  
  # Calculating the K elements from the covariance structure
  K_y <- (tau^-1)*(kernel_function(distance_matrix = distance_matrix_train,
                                   nu = nu,
                                   phi = phi) + diag(x = 1,nrow = nrow(x_train)))
  K_star <- (tau^-1)*(kernel_function(distance_matrix = distance_matrix_K_star,
                                      nu = nu, phi = phi))
  K_star_star <- (tau^-1)*(kernel_function(distance_matrix = distance_matrix_K_star_star,
                                           nu = nu, phi = phi))
  
  # Calculating \alpha
  L <- chol(K_y)
  alpha <- backsolve(L,backsolve(L,y_train, transpose = TRUE))
  mu_star <- crossprod(K_star, alpha)
  
  # Here the abs is because the smallest values that are coming from here are due to numerical approximations.
  v <- backsolve(L,K_star,transpose = TRUE)
  cov_star <-K_star_star - crossprod(v)
  
  # ===============#
  return(list(mu_pred = mu_star, cov = cov_star))
}


# Function to create the the function K that will be used
# in a gaussian process (Andrew Version)
kernel_function <- function(distance_matrix, nu, phi) {
  
  # Calculating the square matrix
  kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * distance_matrix^2)
  
  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0,nrow = dim(distance_matrix)[1],
                            ncol = dim(distance_matrix)[2])
  }
  # Getting the kernel matrix
  return(kernel_matrix)
}

