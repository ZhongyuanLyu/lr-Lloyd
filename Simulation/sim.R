library(rTensor)
source("functions.R")
set.seed(12345)
d = 100; n=100; r = 3; K = 2; sigma = 1; lambda <- 2.5
gen_M <- generate_M_list(d,K,r,lambda)
true_label <- sample(c(1:K), n, replace = TRUE)
data_tsr <- generate_data(gen_M$M_list, true_label, sigma)
rU <- K*r;rV <- K*r;
## Vectorized Lloyd
s_vec_lloyds <- vec_lloyd(data_tsr,K)
## Low-rank Lloyd with spectral Initialization
s_spec_init <- lr_spec_init(data_tsr,rU,rV,K)
lrlloyds_res <- lrlloyds(data_tsr, s_spec_init, rep(r,K), TRUE)
effective_hamming_error(s_spec_init,true_label)
effective_hamming_error(s_vec_lloyds,true_label)
effective_hamming_error(lrlloyds_res$shat,true_label)
