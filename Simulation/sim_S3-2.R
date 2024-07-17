library(rTensor)
source("functions.R")
set.seed(12345)
sim_lam_varying <- function(x, M_list_list, n, r){
  one_sim <- function(M_list, n, r){
    K <- length(M_list)
    true_label <- 2*rbinom(n,1,0.5)-1
    X_tsr <- generate_data(M_list,true_label)
    # est_err <- spec_aggregation(M_list[[1]],X_array, n,d,r)
    ##======= Low-rank Lloyd with spectral initialization =======##
    rU <- K*r; rV <- K*r;
    s_spec_init <- lr_spec_init(X_tsr,rU,rV,K)
    lrlloyds_res <- lrlloyds(X_tsr, s_spec_init, rep(r,K), TRUE)
    est_err_K <- rep(NA,K)
    for (k in 1:K){
      tmp <- apply(X_tsr@data[,,lrlloyds_res$shat==k],c(1,2),mean)
      svd_tmp <- svd(tmp)
      M_lr <- svd_tmp$u[,1:r]%*%diag(svd_tmp$d[1:r])%*%t(svd_tmp$v[,1:r])
      est_err_K[k] <- min(sum((M_lr-M_list[[1]])^2)/sum(M_list[[1]]^2),
                          sum((M_lr-M_list[[2]])^2)/sum(M_list[[2]]^2))
    }
    return(c(max(est_err_K), effective_hamming_error(s_spec_init,true_label),
             effective_hamming_error(lrlloyds_res$shat,true_label+2)))
  }
  error_mat <- matrix(NA, nrow = length(M_list_list), ncol = 3)
  for (i in 1:length(M_list_list)){
    error_mat[i,] <- one_sim(M_list_list[[i]], n, r)
  }
  return(error_mat)
}
##======= Example where clustering is more challenging =======##
d = 5; n=1000; r = 2; K = 2; sigma = 1;
num <- 20;
lambda_mat <- matrix(NA, num, 2)
lambda_mat[,1] <- 0.1
lambda_mat[,2] <- seq(0.5,10, length.out=num)
M_list_list <- vector(mode='list', length=num)
for (i in 1:num) {
  M_list_list[[i]] <- generate_M_list_lam(d,r,lambda_mat[i,], kappa=1.3)
}
error_mat <- lapply(1:100, sim_lam_varying, M_list_list, n, r)
save(error_mat, file = "error_mat_lam.RData")


