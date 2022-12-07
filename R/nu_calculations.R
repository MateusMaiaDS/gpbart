# # Getting the x index
# # Getting the current tree
# tree <- current_trees[[1]]
# t_nodes <- get_terminals(tree)
# n_t_nodes <- length(t_nodes)
#
# obs_train_vec <- t(replicate(numeric(nrow(train_data)),n = length(t_nodes)))
# n_tree <- 1
# big_matrix <- matrix(0,nrow = nrow(x_train),ncol = nrow(x_train)*n_tree)
#
# # Getting the values
# for(i in 1:nrow(x_train)){
#         tree_vector <- numeric(nrow(x_train))
#         tree_vector[i] <- 1
#         big_vector<- rep(tree_vector,n_tree)
#         big_matrix[i,] <- big_vector
# }
#
# # Or doing the a fast way with diagonal
# big_m <- diag(nrow = nrow(x_train))
# y <- rnorm(n = nrow(x_train))
#
# y_calculation <- big_m%*%y
#
# # Getting more more trees
# n_tree <- 10
#
# for(j in 1)
#
# for(j in 1:length(t_nodes)){
#   for(i in 1:length(t_nodes[[j]]$obs_train)){
#       obs_train_vec[j,t_nodes[[j]]$obs_train[i]] <- 1
#   }
# }
# colSums(obs_train_vec)
