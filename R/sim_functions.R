# Simulating functions
#' @export
two_d_sim_two_node_gp_sum_each_tree <- function(n, # Number of observations
                                                mu1 =  c(-10,0,10), # Mean of the first terminal node
                                                mu2  = c(5,20,-15) , # Mean of the second terminal nodde
                                                tau1 = c(10,20,100), # First \tau variance value
                                                tau2 =  c(5,10,10), # Tau variance value
                                                nu = 0.1, # Getting the \nu parameter
                                                phi = 0.1, # Getting the \phi parameter
                                                seed = NULL # Setting the seed
){
  #
  # Setting the seed 
  set.seed(seed)
  
  # Defining the kernel function structure
  omega_function <- function(x, x_star = NULL, nu, phi) {
    
    # Calculating the square matrix
    if (is.null(x_star)) {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x))^2)  +
        diag(1e-8, nrow = nrow(x))
    } else {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x, x_star) )^2)
    }
    
    # Getting the kernel matrix
    return(kernel_matrix)
  }
  
  
  # Generating the x axis
  x <- expand.grid(lat = seq(-1,1,length.out = round(sqrt(n))),
                   lon = seq(-1,1,length.out = round(sqrt(n))))
  
  # Creating the true response
  y <- numeric(nrow(x))
  
  # Getting observation from different tree rules
  tree_one_n1 <- length(which(x[,1]< x[,2]))
  tree_one_n2 <- length(which(x[,1]>= x[,2]))
  
  tree_two_n1 <- length(which(x[,1] < -x[,2]))
  tree_two_n2 <- length(which(x[,1] >= - x[,2]))
  
  tree_three_n1 <- length(which(x[,1] > 0))
  tree_three_n2 <- length(which(x[,1] <= 0))
  
  
  # ==== Sampling for the first treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_one_node_one <-  ((tau1[1])^-1)*omega_function(x = x[x[,1] < x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) + diag((tau1[1]^-1),
                                                                            nrow = tree_one_n1) 
  # Defining the mean
  mean_tree_one_node_one <- rep(mu1[1],tree_one_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < x[,2] ] <-  y [ x[,1] < x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_one,
                                                              Sigma = var_tree_one_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_one_node_two <-  ((tau2[1])^-1)*omega_function(x = x[x[,1] >= x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) + diag((tau2[1]^-1),
                                                                            nrow = tree_one_n2) 
  # Defining the mean
  mean_tree_one_node_two <- rep(mu2[1],tree_one_n2)
  
  y [ x[,1] >= x[,2] ] <-  y [ x[,1] >= x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_two,
                                                                Sigma = var_tree_one_node_two)
  
  
  # ==== Sampling for the second treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_two_node_one <-  ((tau1[2])^-1)*omega_function(x = x[x[,1] < - x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) + diag((tau1[2]^-1),
                                                                            nrow = tree_two_n1) 
  # Defining the mean
  mean_tree_two_node_one <- rep(mu1[2],tree_two_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < -x[,2] ] <-  y [ x[,1] < -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_one,
                                                                Sigma = var_tree_two_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_two_node_two <-  ((tau2[2])^-1)*omega_function(x = x[x[,1] >= -x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) + diag((tau2[2]^-1),
                                                                            nrow = tree_two_n2) 
  # Defining the mean
  mean_tree_two_node_two <- rep(mu2[2],tree_two_n2)
  
  y [ x[,1] >= -x[,2] ] <-  y [ x[,1] >= -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_two,
                                                                  Sigma = var_tree_two_node_two)
  
  
  # ==== Sampling for the third treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_three_node_one <-  ((tau1[3])^-1)*omega_function(x = x[x[,1] < 0, ,drop =FALSE],
                                                            nu = nu, 
                                                            phi = phi) + diag((tau1[3]^-1),
                                                                              nrow = tree_three_n1) 
  # Defining the mean
  mean_tree_three_node_one <- rep(mu1[3],tree_three_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < 0 ] <-  y [ x[,1] < 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_one,
                                                      Sigma = var_tree_three_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_three_node_two <-  ((tau2[3])^-1)*omega_function(x = x[x[,1] >= 0, ,drop =FALSE],
                                                            nu = nu, 
                                                            phi = phi) + diag((tau2[3]^-1),
                                                                              nrow = tree_three_n2) 
  # Defining the mean
  mean_tree_three_node_two <- rep(mu2[3],tree_three_n2)
  
  y [ x[,1] >= 0 ] <-  y [ x[,1] >= 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_two,
                                                        Sigma = var_tree_three_node_two)
  
  
  
  return(as.matrix(data.frame(x, y = y)))
  
}
