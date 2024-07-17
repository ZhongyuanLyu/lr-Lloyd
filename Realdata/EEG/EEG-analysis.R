library(rTensor)
library(sparcl)
library(tensorsparse)
source("functions.R")
source("DTC.R")
# load("data_tsr_S1nonzero.Rdata")
load("data_tsr_S2match_nonzero.Rdata")
# load("data_tsr_S2nonmatch_nonzero.Rdata")
# load("data_tsr_all_nonzero.Rdata")
scale_tsr <- aperm(apply(data_tsr,c(1,2), scale),c(2,3,1))
data_tsr <- as.tensor(scale_tsr)
screeplot(princomp(t(k_unfold(data_tsr, m=1)@data)), npcs = 5, 
          cex.lab=1.5, cex.axis=1, cex.main=5, cex.sub=2,
          type = "lines", main = "Scree plot of EEG: Mode 1")

screeplot(princomp(t(k_unfold(data_tsr, m=2)@data)), npcs = 5, 
          cex.lab=1.5, cex.axis=1, cex.main=5, cex.sub=2,
          type = "lines", main = "Scree plot of EEG: Mode 2")

n <- data_tsr@modes[3]; rU <- 3;rV <- 3; K <- 2;
## Spectral Initialization
HOSVD_init <- hosvd(data_tsr, ranks = c(rU, rV, K))
U_int <- HOSVD_init$U[[1]]
V_int <- HOSVD_init$U[[2]]
Ghat <- ttl(data_tsr, list(U_int%*%t(U_int),V_int%*%t(V_int)), c(1,2))
M3Ghat <- k_unfold(Ghat, m=3)@data
s_spec_init <- kmeans(M3Ghat, K, nstart = 30)$cluster
effective_hamming_error(s_spec_init,true_label+1)
s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, nstart = 30)$cluster
# s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, iter.max = 60,
                       # algorithm = "Lloyd", nstart = 50)$cluster
effective_hamming_error(s_vec_lloyds,true_label+1)
## lr-Lloyd's algorithm
lrlloyds_res <- lrlloyds(data_tsr, s_spec_init, c(2,1), TRUE)
effective_hamming_error(lrlloyds_res$shat,true_label+1)

## Determine the rank
# Mean_mat_list <- vector(mode = "list", length = K)
# for (k in 1:K) {
#   Mean_mat_list[[k]] <- apply(data_tsr@data[,,s_spec_init==k], c(1,2), mean)
# }
# screeplot(princomp(Mean_mat_list[[1]]), npcs = 10, 
#           type = "lines", main = "scree plot")
# screeplot(princomp(Mean_mat_list[[2]]), npcs = 10, 
#           type = "lines", main = "scree plot")




## SKM
km.perm <- KMeansSparseCluster.permute(M3Ghat,K=2, wbounds=seq(2,100,len=15), nperms=5)
effective_hamming_error(KMeansSparseCluster(M3Ghat, K, wbounds = km.perm$bestw)[[1]]$Cs,true_label+1)

## TBM
TBM_res <- tbmClustering(data_tsr@data,3,3,2)
effective_hamming_error(TBM_res$Es,true_label+1)

## DTC
# tune_res <- mytune(data_tsr@data, rank_list = 3, sparse_list = seq(0.3, 0.7, 0.2), fuse_list = seq(0.3, 0.7, 0.2), tune_method = "sequential")
tune_res <- mytune(data_tsr@data, rank_list = 3, tune_method = "sequential")
temp <- stdtruncate(data_tsr@data, sparse_para = tune_res$sparse.opt, smooth_para=tune_res$fuse.opt ,Rank = tune_res$Rank.opt, is_first2modes_symmetric = FALSE)
DTC_res <- mydtc(temp, K)
effective_hamming_error(DTC_res$membership,true_label+1)
tune_res$sparse.opt
tune_res$fuse.opt
## Plot True mean mat
# Mean_true_list <- vector(mode = "list", length = K)
# r_list=c(2,2)
# for (k in 1:K) {
#   mean_mat <- apply(data_tsr@data[,,(true_label+1)==k], c(1,2), mean)
#   mean_svd <- svd(mean_mat, nu = r_list[k], nv = r_list[k])
#   if (r_list[k]==1){
#     Mean_true_list[[k]] <- mean_svd$d[1]* mean_svd$u %*% t(mean_svd$v)
#   } else{
#     Mean_true_list[[k]] <- mean_svd$u %*% diag(mean_svd$d[1:r_list[k]]) %*% t(mean_svd$v)
#   }
# }

# library(rgl)
# persp3D(x = c(1:256), y=seq(1/64, 1, length.out = 64), z = Mean_true_list[[1]], col="burlywood1",
#         shade = .4, theta = 130, phi = 20, along = "xy",border = "black", lwd=0.2,
#         zlim = c(-1.5,1.5),  axes=TRUE, ticktype="detailed",bty = "b2",
#         xlab = "Channels", ylab = "Time (sec)", zlab= "Standardized Voltage")
# 
# persp3D(x = c(1:256), y=seq(1/64, 1, length.out = 64), z = Mean_true_list[[2]], col="burlywood1",
#         shade = .4, theta = 130, phi = 20, along = "xy",border = "black", lwd=0.2,
#         zlim = c(-1.5,1.5),  axes=TRUE, ticktype="detailed",bty = "b2",
#         xlab = "Channels", ylab = "Time (sec)", zlab= "Standardized Voltage")

