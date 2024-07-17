library(rTensor)
source("functions.R")
load("adj_tsr_directed_weight.Rdata")
load("commodites.Rdata")
load("countryname.Rdata")
commodities <- temp2
remove(temp2)
continent_list <- c("Asia","North America","Europe","Asia",
                    "Europe","Europe","Europe","Asia",
                    "Asia","Asia","Europe","North America",
                    "Asia","Europe","North America","Europe",
                    "Europe","Europe","Europe","Australia",
                    "Europe","Asia","Asia","South America",
                    "Europe","Europe","Europe","Asia",
                    "Europe","Asia","Europe","Europe",
                    "Europe","Europe","Asia","Europe",
                    "Africa","Europe","Europe","Europe",
                    "Asia","Asia","Europe",
                    "South America","South America","Asia","Africa","Asia")
data_tsr <- as.tensor(log(A_weight+1))
screeplot(princomp(t(k_unfold(data_tsr, m=1)@data)), npcs = 10,
          type = "lines", main = "scree plot")
screeplot(princomp(t(k_unfold(data_tsr, m=2)@data)), npcs = 10,
          type = "lines", main = "scree plot")
screeplot(princomp(t(k_unfold(data_tsr, m=3)@data)), npcs = 10,
          type = "lines", main = "scree plot")

n <- data_tsr@modes[3]; rU <- 3;rV <- 3; K <- 2;
## Spectral Initialization
HOSVD_init <- hosvd(data_tsr, ranks = c(rU, rV, K))
U_int <- HOSVD_init$U[[1]]
V_int <- HOSVD_init$U[[2]]
Ghat <- ttl(data_tsr, list(U_int%*%t(U_int),V_int%*%t(V_int)), c(1,2))
M3Ghat <- k_unfold(Ghat, m=3)@data
s_spec_init <- kmeans(M3Ghat, K, nstart = 30)$cluster
# s_init <- kmeans(M3Ghat, K, iter.max = 20, algorithm = "Lloyd", nstart = 10)$cluster
# true_label <- rep(c(0,1,2), each = 45)

s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, nstart = 30)$cluster
# s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, iter.max = 60,
# algorithm = "Lloyd", nstart = 50)$cluster
effective_hamming_error(s_vec_lloyds,s_spec_init)


## Determine the rank
Mean_mat_list <- vector(mode = "list", length = K)
for (k in 1:K) {
  Mean_mat_list[[k]] <- apply(data_tsr@data[,,s_spec_init==k], c(1,2), mean)
}
screeplot(princomp(Mean_mat_list[[1]]), npcs = 10,
          type = "lines", main = "scree plot")
screeplot(princomp(Mean_mat_list[[2]]), npcs = 10,
          type = "lines", main = "scree plot")
# screeplot(princomp(Mean_mat_list[[3]]), npcs = 10,
#           type = "lines", main = "scree plot")
## spectral initialization
lrlloyds_res <- lrlloyds(data_tsr, s_spec_init, c(2,3,2), FALSE)
effective_hamming_error(lrlloyds_res$shat,s_spec_init)
effective_hamming_error(lrlloyds_res$shat,s_vec_lloyds)


class_1 <- commodities[which(lrlloyds_res$shat==1),]
class_2 <- commodities[which(lrlloyds_res$shat==2),]
# class_3 <- commodities[which(lrlloyds_res$shat==3),]

# sum(A_weight[,,c(6:14)])/sum(A_weight[,,c(6:15)])
# sum(A_weight[,,c(16:18,23:24)])/sum(A_weight[,,c(16:24)])
# sum(A_weight[,,c(26)])/sum(A_weight[,,c(25:27)])
# sum(A_weight[,,c(31,36:37)])/sum(A_weight[,,c(28:38)])
# sum(A_weight[,,c(41,43)])/sum(A_weight[,,c(41:43)])
# sum(A_weight[,,c(45:47)])/sum(A_weight[,,c(44:49)])
# sum(A_weight[,,c(50:55,57:58,60)])/sum(A_weight[,,c(50:63)])
# sum(A_weight[,,c(65:67)])/sum(A_weight[,,c(64:67)])
# sum(A_weight[,,which(commodities$Commodity.Code %in% c(75,78:81))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(72:83))])
# sum(A_weight[,,which(commodities$Commodity.Code %in% c(86,89))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(86:89))])
# sum(A_weight[,,which(commodities$Commodity.Code %in% c(91:93,97))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(90:99))])
# 
# 1-sum(A_weight[,,c(6:14)])/sum(A_weight[,,c(6:15)])
# 1-sum(A_weight[,,c(16:18,23:24)])/sum(A_weight[,,c(16:24)])
# 1-sum(A_weight[,,c(26)])/sum(A_weight[,,c(25:27)])
# 1-sum(A_weight[,,c(31,36:37)])/sum(A_weight[,,c(28:38)])
# 1-sum(A_weight[,,c(41,43)])/sum(A_weight[,,c(41:43)])
# 1-sum(A_weight[,,c(45:47)])/sum(A_weight[,,c(44:49)])
# 1-sum(A_weight[,,c(50:55,57:58,60)])/sum(A_weight[,,c(50:63)])
# 1-sum(A_weight[,,c(65:67)])/sum(A_weight[,,c(64:67)])
# 1-sum(A_weight[,,which(commodities$Commodity.Code %in% c(75,78:81))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(72:83))])
# 1-sum(A_weight[,,which(commodities$Commodity.Code %in% c(86,89))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(86:89))])
# 1-sum(A_weight[,,which(commodities$Commodity.Code %in% c(91:93,97))])/sum(A_weight[,,which(commodities$Commodity.Code %in% c(90:99))])


