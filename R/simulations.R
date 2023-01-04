#' @export
sim_friedman <- function(n,seed,sd){

     #Setting a seed
     set.seed(seed)
     return(as.matrix(as.data.frame(mlbench::mlbench.friedman1(n = n,sd = sd))))
}


#' @export
sim_friedman_andrew <- function(n,seed,sd){
        
        # Setting a seed if necessary
        set.seed(seed)
        
        # Getting the X matrix
        x <- matrix(stats::runif(2* n), ncol = 2)
        
        y <- 20*(x[,1]-0.5)^2
        
        # Adding the noise
        y <- y + stats::rnorm(n = n,sd = sd)
        
        sim_data <- cbind(x,y)
        colnames(sim_data) <- c(paste0("x.",1:2),"y")
        
        # Returning the simulated data 
        return(sim_data)
}