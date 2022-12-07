# Exporting the kfold function
#' @export
k_fold <- function(data,dependent_variable="y",k_partitions,seed=NULL,as_data_frame = FALSE){

     # Setting the seed.
     set.seed(seed)

     # Creating the object k_fold
     k_fold_validation <- list()

     # Creating the test_ratio
     partitions_index <- kfold(x = data,k = k_partitions)

     # Without use any package
     # if(nrow(data)%%k_partitions!=0){
     # partitions_index <- sample(rep(1:k_partitions,floor(nrow(data)/k_partitions)+1),
     #                            size = nrow(data),replace = FALSE)
     # }else{
     #   partitions_index <- sample(rep(1:k_partitions,k_partitions,
     #                              size = nrow(data),replace = FALSE))
     # }
     for(k in 1:k_partitions){
          # Training values
          cat("Running repetition number ",k,"\n")

          # Training
          if(as_data_frame){
               x_train = data[partitions_index!=k,colnames(data)!="y",drop=FALSE] #%>% as.matrix
               y_train = data[partitions_index!=k,colnames(data)=="y",drop=FALSE] #%>% as.matrix

               # Test
               x_test = data[partitions_index==k,colnames(data)!="y",drop=FALSE] #%>% as.matrix
               y_test = data[partitions_index==k,colnames(data)=="y",drop=FALSE] #%>% as.matrix
          } else {
               x_train = as.matrix(data[partitions_index!=k,colnames(data)!="y",drop=FALSE]) #%>% as.matrix
               y_train = as.matrix(data[partitions_index!=k,colnames(data)=="y",drop=FALSE]) #%>% as.matrix

               # Test
               x_test = as.matrix(data[partitions_index==k,colnames(data)!="y",drop=FALSE]) #%>% as.matrix
               y_test = as.matrix(data[partitions_index==k,colnames(data)=="y",drop=FALSE]) #%>% as.matrix

          }

          # Saving the data split
          k_fold_validation[[k]]<-list(x_train=x_train,y_train=y_train,x_test=x_test,y_test=y_test)
     }

     return(k_fold_validation)

}

# PROBLEM TO LOAD
# Author: Robert Hijmans
# January 2010
# License GPL3


# FROOMM PACKAGE DISMOO

kfold <- function(x, k=5, by=NULL) {

     singlefold <- function(obs, k) {
          if (k==1) {
               return(rep(1, obs))
          } else {
               i <- obs / k
               if (i < 1) {
                    stop('insufficient records:', obs, ', with k=', k)
               }
               i <- round(c(0, i * 1:(k-1), obs))
               times = i[-1] - i[-length(i)]

               group <- c()
               for (j in 1:(length(times))) {
                    group <- c( group, rep(j, times=times[j]) )
               }

               r <- order(stats::runif(obs))
               return(group[r])
          }
     }

     if (is.vector(x)) {
          if (length(x) == 1) {
               if (x > 1) {
                    x <- 1:x
               }
          }
          obs <- length(x)
     }  else {
          obs <- nrow(x)
     }
     if (is.null(by)) {
          return(singlefold(obs, k))
     }

     by = as.vector(as.matrix(by))
     if (length(by) != obs) {
          stop('by should be a vector with the same number of records as x')
     }
     un <- unique(by)
     group <- vector(length=obs)
     for ( u in un ) {
          i = which(by==u)
          kk = min(length(i), k)
          if (kk < k) warning('lowered k for by group: ', u  ,'  because the number of observations was  ',  length(i))
          group[i] <- singlefold(length(i), kk)
     }
     return(group)
}
