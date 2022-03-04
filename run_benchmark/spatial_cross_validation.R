# This file generate spatial cross-validation, i.e: use the number "k" of partitions as the number of clusters in k-means partition
spatial_k_fold<-function(data,dependent_variable="y",k_partitions,seed=NULL){
  
  # Create the list of the k_fold_validation
  k_fold_validation<-list()
  
  # Defining a seed
  if(!is.null(seed)){
    set.seed(seed)
  }
  #
  # Selecting spatial coordinates
  spatial_coordinates<-data[,c("lat","lon"),drop=FALSE]
  
  # Calling the clusters related with the number of partitions
  cluster_model<-kmeans(spatial_coordinates,centers = k_partitions)
  
  # Saving partitions index
  partitions_index<-cluster_model$cluster
  
  for(k in 1:k_partitions){
    # Training values
    cat("Running repetition number ",k,"\n")
    
    # Training 
    x_train = data[partitions_index!=k,colnames(data)!="y",drop=FALSE] %>% as.matrix
    y_train = data[partitions_index!=k,colnames(data)=="y",drop=FALSE] %>% as.matrix
    
    # Test
    x_test = data[partitions_index==k,colnames(data)!="y",drop=FALSE] %>% as.matrix
    y_test = data[partitions_index==k,colnames(data)=="y",drop=FALSE] %>% as.matrix
    
    
    # Saving the data split
    k_fold_validation[[k]]<-list(x_train=x_train,y_train=y_train,x_test=x_test,y_test=y_test)
  }
  
  return(k_fold_validation)
}


# Creating the cross-validation scenario
cross_validation <- function(data, training_ratio = 0.7,
                             seed = NULL){
  # Setting the seed
  set.seed(seed)
  
  # Creating the test ratio
  test_ratio <- 1 - sum(training_ratio)
  
  # Verifying if the sum of ratios return to 1
  if(sum(training_ratio,test_ratio) != 1) {stop("The sum of the ratios must be equal to 1.")}
  
  # Creating the train, validation and test
  set.seed(seed)
  
  # Sampling the training samples
  training_index <- sample(x = 1:nrow(data),
                           size = round(nrow(data)*training_ratio))
  
  # Setting the training sample
  training_sample <-  data[training_index,]
  
  # Gathering the validation index
  training_index <- sample(x = (1:nrow(data)),
                           size = round(nrow(data)*training_ratio))
  
  # Getting the test sample
  test_sample <- data[-c(training_index),]
  
  # Return a list with the training, validation, and test sample
  return(list(train_sample = training_sample,
              test_sample = test_sample))
}


k_fold <- function(data,dependent_variable="y",k_partitions,seed=NULL,as_data_frame = FALSE){
  
  
  # Setting the seed.
  set.seed(42)
  
  # Creating the object k_fold
  k_fold_validation <- list()
  
  # Creating the test_ratio 
  partitions_index <- dismo::kfold(x = data,k = k_partitions)
  
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
      x_train = data[partitions_index!=k,colnames(data)!="y",drop=FALSE] %>% as.matrix
      y_train = data[partitions_index!=k,colnames(data)=="y",drop=FALSE] %>% as.matrix
      
      # Test
      x_test = data[partitions_index==k,colnames(data)!="y",drop=FALSE] %>% as.matrix
      y_test = data[partitions_index==k,colnames(data)=="y",drop=FALSE] %>% as.matrix
    }
    
    # Saving the data split
    k_fold_validation[[k]]<-list(x_train=x_train,y_train=y_train,x_test=x_test,y_test=y_test)
  }
  
  return(k_fold_validation)
  
}
