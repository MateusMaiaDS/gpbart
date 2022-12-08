# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix_phi, nu) {

        # Calculating the square matrix
        kernel_matrix <- (exp(-squared_distance_matrix_phi)) / nu

        # Case nu = 0
        if(nu == 0 || nu > 1e13){
                kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix_phi)[1],
                                        ncol = dim(squared_distance_matrix_phi)[2])
        }
        # Getting the kernel matrix
        return(kernel_matrix)
}


# Calculating RMSE
#' @export
rmse <- function(obs, pred) {
        return(sqrt(stats::mean((obs - pred)^2)))
}

# Calculating CRPS from (https://arxiv.org/pdf/1709.04743.pdf)
#' @export
crps <- function(y,means,sds){

        # scaling the observed y
        z <- (y-means)/sds

        crps_vector <- sds*(z*(2*stats::pnorm(q = z,mean = 0,sd = 1)-1) + 2*stats::dnorm(x = z,mean = 0,sd = 1) - 1/(sqrt(pi)) )

        return(list(CRPS = stats::mean(crps_vector), crps = crps_vector))
}


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

# Getting a half cauchy
dhalfcauchy <- function(x,mu,sigma, log = FALSE) {

        if(!log){
        # Vectorised version of dhalfcauchy
                ifelse(x>mu,(2/(pi*sigma))*(1/(1+((x-mu)^2)/(sigma^2))),0)
        } else {# Vectorised version of dhalfcauchy
                log(ifelse(x>mu,(2/(pi*sigma))*(1/(1+((x-mu)^2)/(sigma^2))),0))
        }

}

# Define the function
progress_bar <- function(i, n) {
        # Calculate the progress value
        progress <- i / n
        
        # Calculate the number of hashes to display
        hashes <- round(progress * 40)
        
        # Create the progress bar string
        bar <-(c(cat(rep("#", hashes)), cat(rep("-", 40 - hashes))))
        
        # Add the progress value to the progress bar
        progress_string <- paste(round(progress * 100), "%")
        bar <- paste0(bar, "  ", progress_string)
        
        # Return the progress bar string
        return(bar)
}
