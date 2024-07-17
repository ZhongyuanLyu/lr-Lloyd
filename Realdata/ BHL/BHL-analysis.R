rm(list=ls())
library(rTensor)
library(rgl)
# library(GEOquery)
source("functions.R")
load("data_tsr.Rdata")
scale_tsr <- aperm(apply(data_tsr,c(1,2), scale, FALSE),c(2,3,1))
data_tsr <- as.tensor(scale_tsr)
plot(scale_tsr[,1,1], "p")

persp3d(x = c(1:1124), y=seq(1,4), z = scale_tsr[,,1], col="burlywood1",
        shade = .4, theta = 130, phi = 20, along = "xy",border = "black", lwd=0.2,
        zlim = c(0,2),  axes=TRUE, ticktype="detailed",bty = "b2",
        xlab = "Channels", ylab = "Time (sec)", zlab= "Standardized Voltage")




screeplot(princomp(k_unfold(data_tsr, m=1)@data), npcs = 4, 
          cex.lab=1.5, cex.axis=1, cex.main=5, cex.sub=2,
          type = "lines", main = "Scree plot of BHL: Mode 1")

screeplot(princomp(t(k_unfold(data_tsr, m=2)@data)), npcs = 4, 
          cex.lab=1.5, cex.axis=1, cex.main=5, cex.sub=2,
          type = "lines", main = "Scree plot of BHL: Mode 2")

n <- data_tsr@modes[3]; K <- 3;
## Spectral Initialization
HOSVD_init <- hosvd(data_tsr, ranks = c(1, 1, K))
U_int <- HOSVD_init$U[[1]]
V_int <- HOSVD_init$U[[2]]
Ghat <- ttl(data_tsr, list(U_int%*%t(U_int),V_int%*%t(V_int)), c(1,2))
M3Ghat <- k_unfold(Ghat, m=3)@data
s_init <- kmeans(M3Ghat, K, nstart = 50)$cluster
effective_hamming_error(s_init,true_labels)
s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, nstart = 50)$cluster
# s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, 2, iter.max = 30,
#                        algorithm = "Lloyd", nstart = 50)$cluster
effective_hamming_error(s_vec_lloyds,true_labels)

## Low-rank Lloyd's Algorithm
## Determine the rank
Mean_mat_list <- vector(mode = "list", length = K)
for (k in 1:K) {
  Mean_mat_list[[k]] <- apply(data_tsr@data[,,s_init==k], c(1,2), mean)
}
screeplot(princomp(Mean_mat_list[[1]]), npcs = 4, 
          type = "lines", main = "scree plot")
screeplot(princomp(Mean_mat_list[[2]]), npcs = 4, 
          type = "lines", main = "scree plot")
## Iterations
shat <- lrlloyds(data_tsr, s_init, c(1,1,1))
effective_hamming_error(shat,true_labels)
shat_vec <- lloyds(data_tsr, s_init)
effective_hamming_error(shat_vec,true_labels)




