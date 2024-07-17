library(rTensor)
source("functions.R")
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
    return(c(effective_hamming_error(s_spec_init,true_label),
             effective_hamming_error(lrlloyds_res$shat,true_label+2)))
  }
  error_mat <- matrix(NA, nrow = length(M_list_list), ncol = 2)
  for (i in 1:length(M_list_list)){
    error_mat[i,] <- one_sim(M_list_list[[i]], n, r)
  }
  return(error_mat)
}
##======= high-dimensional setting =======##
d = 100; n=50; r = 2; K = 2; sigma = 1;  num <- 20;
lambda_list <- seq(1.5,3.8, length.out=num)
M_list_list <- vector(mode='list', length=num)
for (i in 1:num) {
  M_list_list[[i]] <- generate_M_list_lam(d,r,lambda_list[i], kappa=1.3)
}
error_mat <- lapply(1:100, sim_lam_varying_two_comp, M_list_list, n, r)
save(error_mat, file = "error_mat_lam_two_component_high.RData")
