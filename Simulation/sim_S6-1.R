library(rTensor)
library(MatTransMix)
source("libMatTransFull.R")
source("functions.R")
# data(iris)
# data <- as.matrix(iris[,-5])
# n <- nrow(data)
# p <- 2
# T <- 2
# X1 <- array(NA, dim = c(p, T, n))
# 
# for(i in 1:n){
#   X1[1,1,i] <- data[i,1]
#   X1[1,2,i] <- data[i,2]
#   X1[2,1,i] <- data[i,3]
#   X1[2,2,i] <- data[i,4]
# }
# 
# 
# M <- MatTrans.EM(X1, initial = MatTrans.init(X1, K = 2, n.start = 1), trans = "None", silent = TRUE, size.control = 10)
# M$best.result[[1]]$id


set.seed(12345)
sim_lam_varying_two_comp <- function(x, M_list_list, n, r){
  one_sim <- function(M_list, n, r){
    K <- length(M_list)+1
    true_label <- 2*rbinom(n,1,0.5)-1
    X_array <- generate_data_X(M_list[[1]],true_label, n,d)
    ##======= Low-rank Lloyd with spectral initialization =======##
    rU <- K*r; rV <- K*r;
    X_tsr <- as.tensor(X_array);
    s_spec_init <- lr_spec_init(X_tsr,rU,rV,K)
    lrlloyds_res <- lrlloyds(X_tsr, s_spec_init, rep(r,K), TRUE)
    mattrans_res <- MatTrans.EM(X_array, initial = MatTrans.init(X_array, K = 2, n.start = 1), trans = "None", silent = TRUE, size.control = 10)
    return(c(effective_hamming_error(s_spec_init,true_label),
             effective_hamming_error(lrlloyds_res$shat,true_label+2),
             effective_hamming_error(mattrans_res$best.result[[1]]$id,true_label+2)))
  }
  error_mat <- matrix(NA, nrow = length(M_list_list), ncol = 3)
  for (i in 1:length(M_list_list)){
    error_mat[i,] <- one_sim(M_list_list[[i]], n, r)
  }
  return(error_mat)
}
##======= low-dimensional setting =======##
d = 50; n=200; r = 2; K = 2; sigma = 1;  num <- 20;
lambda_list <- seq(0.5,3, length.out=num)
M_list_list <- vector(mode='list', length=num)
for (i in 1:num) {
  M_list_list[[i]] <- generate_M_list_lam(d,r,lambda_list[i], kappa=1.3)
}
error_mat <- lapply(1:1, sim_lam_varying_two_comp, M_list_list, n, r)
# save(error_mat, file = "error_mat_lam_two_component_compare.RData")
