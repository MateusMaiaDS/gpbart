# GP-function main
gp_main <- function(x_train, y_train, x_star, tau, phi, nu, distance_matrix_train, get_cov_star = FALSE) {
  
  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train, m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star)
  
  # Getting the distance matrix from x_train and x_star (SPATIAL VALUES)
  # distance_matrix_K_star <- distance_matrix(m1 = x_train[,c("lat","lon")],m2 = x_star[,c("lat","lon")])
  # distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star[,c("lat","lon")])
  
  # Calculating the K elements from the covariance structure
  n_train <- nrow(x_train)
  K_y <- kernel_function(squared_distance_matrix = distance_matrix_train,
                         nu = nu,
                         phi = phi) + diag(x = 1e-18, nrow = n_train)
  K_diag <- is_diag_matrix(K_y)
  K_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star,
                            nu = nu, phi = phi)
  
  # Calculating \alpha
  if(K_diag) {
    L <- diag(K_y)
    alpha <- y_train/L
  } else {
    L <- chol(K_y)
    alpha <- backsolve(L, backsolve(L, y_train, transpose = TRUE, k = n_train), k = n_train)
  }
  mu_star <- crossprod(K_star, alpha)
  
  # Here the abs is because the smallest values that are coming from here are due to numerical approximations.
  if(isTRUE(get_cov_star)) {
    K_star_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star_star,
                                   nu = nu, phi = phi) + diag(x = 1/tau, nrow = nrow(x_star))
    v <- if(K_diag) K_star/L else backsolve(L, K_star, transpose = TRUE, k = n_train)
    cov_star <- K_star_star - crossprod(v)
    results <- list(mu_pred = mu_star, cov_pred = cov_star)
  } else {
    results <- list(mu_pred = mu_star)
  }
  
  # ===============#
  return(results)
}

# GP-function main
gp_main_sample <- function(x_train, y_train, x_star, tau, phi, nu, distance_matrix_train, get_cov_star = TRUE) {
  
  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train, m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star)
  
  # Calculating the K elements from the covariance structure
  n_train <- nrow(x_train)
  K_y <- kernel_function(squared_distance_matrix = distance_matrix_train,
                         nu = nu,
                         phi = phi) + diag(x = 1/(tau), nrow = n_train)
  K_diag <- is_diag_matrix(K_y)
  K_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star,
                            nu = nu, phi = phi)
  
  # Calculating \alpha
  if(K_diag) {
    L <- diag(K_y)
    alpha <- y_train/L
  } else {
    L <- chol(K_y)
    alpha <- backsolve(L, backsolve(L, y_train, transpose = TRUE, k = n_train), k = n_train)
  }
  mu_star <- crossprod(K_star, alpha)
  
  # Here the abs is because the smallest values that are coming from here are due to numerical approximations.
  if(isTRUE(get_cov_star)) {
    K_star_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star_star,
                                   nu = nu, phi = phi) + diag(x = 1/tau, nrow = nrow(x_star))
    v <- if(K_diag) K_star/L else backsolve(L, K_star, transpose = TRUE, k = n_train)
    cov_star <- K_star_star - crossprod(v)
    
    residuals_sample <- mapply(FUN=rMVN_var, mu_star, cov_star, SIMPLIFY = FALSE)
    
    results <- list(mu_pred = residuals_sample, cov_pred = cov_star)
  } else {
    results <- list(mu_pred = mu_star)
  }
  
  # ===============#
  return(results)
}

# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix, nu, phi) {
  
  # Calculating the square matrix
  kernel_matrix <- exp(-squared_distance_matrix / (2 * phi^2)) / nu
  
  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix)[1],
                            ncol = dim(squared_distance_matrix)[2])
  }
  # Getting the kernel matrix
  return(kernel_matrix)
}
