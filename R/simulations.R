#' @export
sim_friedman <- function(n,p = 10,seed,sd){

     # Setting a seed if necessary
     set.seed(seed)
     
     # Getting the X matrix
     x <- matrix(stats::runif(p* n), ncol = p)
     
     y <- 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
     
     # Adding the noise
     y <- y + stats::rnorm(n = n,sd = sd)
     
     sim_data <- cbind(x,y)
     colnames(sim_data) <- c(paste0("x.",1:p),"y")
     
     # Returning the simulated data 
     return(sim_data)
}


