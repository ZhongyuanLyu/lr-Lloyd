library(rTensor)
source("functions.R")
sim_n_varying <- function(x, M_list, n_list, r, sigma){
  one_sim <- function(M_list, n, r, sigma){
    K <- length(M_list)
    true_label <- sample(c(1:K), n, replace = TRUE)
    data_tsr <- generate_data(M_list, true_label, sigma)
    # print(svd(k_unfold(data_tsr, m=1)@data)$d)
    # rU <- K*r; rV <- K*r
    rU <- r; rV <- r
    ## Vectorized Lloyd with  spectral initialization
    vec_res <- vec_lloyd(data_tsr,K)
    ## Vectorized Lloyd with HOSVD
    m3_tsr <- k_unfold(data_tsr, m=3)@data
    print(min(svd(m3_tsr)$d))
    hosvd_res <- kmeans(m3_tsr,  K, nstart = 30)
    s_hosvd <- hosvd_res$cluster
    center_hosvd <- hosvd_res$centers
    s_hosvd_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, center_hosvd, iter.max = 3*log(n),
                           algorithm = "Lloyd", nstart = 30)$cluster
    ## Low-rank Lloyd with spectral initialization
    s_spec <- lr_spec_init(data_tsr,rU,rV,K)
    lrlloyds_spec_res <- lrlloyds(data_tsr, s_spec, rep(r,K), FALSE)
    ## Low-rank Lloyd with HOSVD
    lrlloyds_hosvd_res <- lrlloyds(data_tsr, s_hosvd, rep(r,K), FALSE)
    return(c(effective_hamming_error(vec_res$lloyd,true_label),
             effective_hamming_error(lrlloyds_spec_res$shat,true_label),
             # effective_hamming_error(s_hosvd,true_label),
             effective_hamming_error(s_hosvd_vec_lloyds,true_label),
             effective_hamming_error(lrlloyds_hosvd_res$shat,true_label)))
  }
  error_mat <- matrix(NA, nrow = length(n_list), ncol = 4)
  for (i in 1:length(n_list)){
    error_mat[i,] <- one_sim(M_list, n_list[i], r, sigma)
  }
  return(error_mat)
}
set.seed(123456)
# d = 100; r = 3; K = 2; sigma = 1; lambda <- 10; Delta <- 40
## Setting K
# gen_M <- generate_M_list(d,K,r,lambda)
# gen_M$Delta
# error_mat <- lapply(1:100, sim_n_varying, gen_M$M_list,
#                     n_list = c(100,200), r, sigma)
## Setting (I)
d = 50; r = 3; K = 2; sigma = 1; lambda <- 10; Delta <- 10
M_list <- generate_two_M(d,r,lambda, Delta)
norm(M_list[[1]]-M_list[[2]],"F")
error_mat <- lapply(1:1, sim_n_varying, M_list,
                    n_list = c(100,200), r, sigma)
## Setting (II)
# gen_M <- generate_two_M_o(d,r,lambda)
# gen_M$Delta
# error_mat <- lapply(1:1, sim_n_varying, gen_M$M_list,
#                     n_list = c(100, 200), r, sigma)
## Setting (relaxed)
# gen_M <- generate_M_list_relaxed(d,K,r,lambda)
# gen_M$Delta
# error_mat <- lapply(1:100, sim_n_varying, gen_M$M_list,
#                     n_list = c(100, 200), r, sigma)


# apply(array(unlist(error_mat), dim = c(2,4,100)), c(1,2), mean)
# sqrt(apply(array(unlist(error_mat), dim = c(2,4,100)), c(1,2), var))

