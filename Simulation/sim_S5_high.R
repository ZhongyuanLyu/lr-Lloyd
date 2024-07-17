library(rTensor)
source("functions.R")
set.seed(12345)
sim_r_varying <- function(x, M_list, n, r, rc_list){
  one_sim <- function(M_list, n, r, rc){
    K <- length(M_list)+1
    true_label <- 2*rbinom(n,1,0.5)-1
    X_array <- generate_data_X(M_list[[1]],true_label, n,d)
    ##======= Low-rank Lloyd with spectral initialization =======##
    rU <- K*min(r,rc); rV <- K*min(r,rc);
    X_tsr <- as.tensor(X_array)
    s_spec_init <- lr_spec_init(X_tsr,rU,rV,K)
    lrlloyds_res <- lrlloyds(X_tsr, s_spec_init, rep(rc,K), FALSE)
    return(c(effective_hamming_error(s_spec_init,true_label),
             effective_hamming_error(lrlloyds_res$shat,true_label+2)))
  }
  error_mat <- matrix(NA, nrow = length(rc_list), ncol = 2)
  for (i in 1:length(rc_list)){
    error_mat[i,] <- one_sim(M_list, n, r, rc_list[i])
  }
  return(error_mat)
}
##======= Robustness w.r.t overspecified r_k =======##
d = 100; n=50; r = 3; K = 2; sigma = 1; lambda <- 2.5;
M_list <- generate_M_list_lam(d,r,lambda, kappa=1.3)
error_mat <- lapply(1:1, sim_r_varying, M_list, n, r,
                    rc_list = seq(r,10))
save(error_mat, file = "error_mat_r_high.RData")
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(8,2,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(8,2,100)), c(1,2), var))

