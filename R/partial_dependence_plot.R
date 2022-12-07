# Create a function to return a partial dependence plot
# partial_dependence_plot <- function(gp_bart_obj,
#                                     x_train,
#                                     y_train,
#                                     var) {
#
#         # Creating partial dependence plots
#         N <- N_new <- gp_bart_obj$y_hat_post %>% ncol()
#         partial_dependence_x3 <- numeric()
#         partial_dependence_each_tree <- matrix(NA,nrow = nrow(x_train), ncol = 1)
#         K <- gp_bart_obj$last_trees %>% length()
#
#         nu <- gp_bart_obj$posterior$nu_post %>% apply(2,mean)
#         tau <- gp_bart_obj$tau_post %>% median()
#         partial_dependence_each_tree <- matrix(NA,nrow = nrow(x_train), ncol = K)
#
#         all_tree_phi <-list()
#         for(j in 1:nrow(gp_bart_obj$posterior$phi_post[[1]])){
#                 all_phi_mcmc_aux <- matrix(NA,ncol = ncol(x_train))
#                 for(i in 1:length(gp_bart_obj$posterior$phi_post)){
#                         all_phi_mcmc_aux <- rbind(all_phi_mcmc_aux,gp_bart_obj$posterior$phi_post[[i]][j,])
#                 }
#
#                 all_tree_phi[[j]] <- all_phi_mcmc_aux[-1,]
#
#         }
#         gp_bart_obj$posterior$phi_post[[1]]
#         phi <- do.call(rbind,gp_bart_obj$posterior$phi_post)
#
#         h <- normalize_bart(y = y_train) %>% matrix(ncol = 1)
#         h <- Reduce("+",gp_bart_obj$posterior$partial_residuals)/2500
#
#         for(i in 1:nrow(x_train)){
#                 x_test <- x_train
#                 x_test[,var] <- x_train[,var][i]
#
#                 # x_new <- seq(-pi, pi, length = N_new)
#                 JAGS_Sigma_new <- array(NA, dim = c(K, N, N_new))
#                 JAGS_Sigma_star <- JAGS_cov_h_new <- array(NA, dim = c(K, N_new, N_new))
#                 JAGS_Sigma <- array(NA, dim = c(K, N, N))
#                 JAGS_h_new <- matrix(NA, ncol = K, nrow = N_new)
#                 JAGS_cov_h <- array(NA, dim = c(K, N, N))
#
#
#                 for (k in 1:K) {
#                         JAGS_Sigma_new[k, , ] <- c((nu[k]^(-1))) * exp(-distance_matrix(m1 = x_train,m2 = x_test,phi_vector = all_tree_phi[[k]] %>% colMeans()))
#                         JAGS_Sigma_star[k, , ] <- c((nu[k]^(-1))) * exp(-symm_distance_matrix(m1 = x_test,phi_vector = all_tree_phi[[k]] %>% colMeans())) #+ diag((1/tau), N_new)
#                         JAGS_Sigma[k, , ] <- (nu[k]^(-1)) * exp(-symm_distance_matrix(m1 = x_train,phi_vector = all_tree_phi[[k]] %>% colMeans())) + diag((1/(tau*(max(y_train)-min(y_train))^2)), nrow = N)
#
#                         JAGS_h_new[, k] <- t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], h[k, ] )
#
#                         JAGS_cov_h_new[k, , ] <- JAGS_Sigma_star[k, , ] - t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], JAGS_Sigma_new[k, , ])
#                 }
#
#
#                 partial_dependence_x3[i] <- mean( unnormalize_bart(rowSums(JAGS_h_new),a = min(y_train),b = max(y_train)))
#                 partial_dependence_each_tree[i,] <- unnormalize_bart(JAGS_h_new %>% as.matrix %>% apply(2,mean),a = min(y_train),b = max(y_train))
#
#                 print(i)
#         }
#
#         plot(x_train[,var],partial_dependence_x3)
#
#         # Plotting each tree dataset.
#         plot_each_tree <- partial_dependence_each_tree %>% cbind( x_train[,"x.3"]) %>% as.data.frame()
#         colnames(plot_each_tree) <- c(paste0("tree_",1:K),"x")
#         plot_each_tree <- plot_each_tree %>% pivot_longer(cols = starts_with("tree"),names_to = "tree") %>% rename(pred = value)
#
#         # Plot output
#         # Plotting the data
#         ggplot()+
#                 geom_line(data = plot_each_tree, mapping = aes(x = x,
#                                                                y = pred, col = tree),show.legend = FALSE,
#                           alpha = 0.5)+
#                 geom_line(data = plot_each_tree, mapping = aes(x = x,
#                                                                y = pred, col = tree),show.legend = FALSE,
#                           alpha = 0.5)+
#
#                 # geom_ribbon(data = data.frame(x = x_new,
#                 #                               ci_up = JAGS_pred_mean + qnorm(0.75)*JAGS_pred_sd,
#                 #                               ci_low = JAGS_pred_mean - qnorm(0.75)*JAGS_pred_sd),
#                 #             mapping = aes(x = x,
#                 #                           ymin = ci_up,
#                 #                           ymax = ci_low), alpha = 0.2, fill = "red")+
#                 # geom_line(data = data.frame(x = x_new,
#                 #                             y = JAGS_pred_mean),
#                 #           mapping = aes( x = x,
#                 #                          y = y), col = "red")+
#         geom_point(data = data.frame(x = x_train[,"x.3"],
#                                      y = partial_dependence_x3),
#                    mapping = aes(x = x, y = y), )+
#                 # geom_line(data = data.frame(x = x_new,
#                 #                             y = sqrt(true_nu^(-1)*(true_tau^-1))*sin(x_new)),
#                 #           mapping = aes(x = x, y = y), col = "blue")+
#                 # ylim(c(-2,2))+
#                 ggtitle(paste("SUM-GP"))+
#                 theme_classic()
#
#
# }


# Getting new test value
new_partial_dependence_plot <- function(gp_bart_mod,
                                        var){

        # Saving the mod
        x_train <- gp_bart_mod$data$x_train
        partial_dependence_plot <- numeric(nrow(x_train))
        partial_dependence_each_tree <- matrix(NA,nrow = nrow(x_train), ncol = length(gp_bart_mod$last_trees))

        for(i in 1:nrow(x_train)){
                x_test <- x_train
                x_test[,var] <- x_test[,var][i]

                new_pred <- gp_bart_predict(gpbart_mod_example = gp_bart_mod,x_new = x_test)

                partial_dependence_plot[i] <- stats::mean(colMeans(new_pred$y_hat_new_post))
                # partial_dependence_each_tree[i,] <- unnormalize_bart(rowMeans(Reduce("+",new_pred$all_trees)/length(new_pred$all_trees)),
                #                                                      a = min(gp_bart_mod$data$y_train),
                #                                                      b = max(gp_bart_mod$data$y_train))
        }

        return(list(mod = gp_bart_mod,
                    partial_dependence_plot = partial_dependence_plot,
                    var = var))

}


# new_partial_d_plot <- new_partial_dependence_plot(gp_bart_mod = gpbart_mod_example,var = "x.3")

# plot(new_partial_d_plot$mod$data$x_train[,"x.3"],new_partial_d_plot$partial_dependence_plot)
