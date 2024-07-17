library(rTensor)
source("functions.R")
set.seed(12345)
sim_n_varying <- function(x, M_list, n_list, r){
  one_sim <- function(M_list, n, r){
    K <- length(M_list)+1
    true_label <- 2*rbinom(n,1,0.5)-1
    X_array <- generate_data_X(M_list[[1]],true_label, n,d)
    ##======= Low-rank Lloyd with spectral initialization =======##
    rU <- K*r; rV <- K*r;
    X_tsr <- as.tensor(X_array)
    s_spec_init <- lr_spec_init(X_tsr,rU,rV,K)
    lrlloyds_res <- lrlloyds(X_tsr, s_spec_init, rep(r,K), FALSE)
    return(c(est_err, effective_hamming_error(s_spec_init,true_label),
             effective_hamming_error(lrlloyds_res$shat,true_label+2)))
  }
  error_mat <- matrix(NA, nrow = length(n_list), ncol = 5)
  for (i in 1:length(n_list)){
    error_mat[i,] <- one_sim(M_list, n_list[i], r)
  }
  return(error_mat)
}
##======= Example where clustering is more challenging =======##
d = 5; n=1000; r = 2; K = 2; sigma = 1; lambda <- 0.5;
M_list <- generate_M_list_lam(d,r,lambda, kappa=1.3)
error_mat <- lapply(1:100, sim_n_varying, M_list,
                    n_list = round(seq(1000, 10000, length.out = 20)), r)
save(error_mat, file = "error_mat_n.RData")
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(10,5,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(10,5,100)), c(1,2), var))

