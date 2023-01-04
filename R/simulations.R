#' @export
sim_friedman <- function(n,seed,sd){

     #Setting a seed
     set.seed(seed)
     return(as.matrix(as.data.frame(mlbench::mlbench.friedman1(n = n,sd = sd))))
}


