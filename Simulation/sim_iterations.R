library(rTensor)
library(ggplot2)
library(latex2exp)
source("functions.R")
change_label <- function(x){
  return(ifelse(x==1,2,1))
}
set.seed(12345)
one_line <- function(n = 200, d = 50, r = 2, K = 2, sigma = 1, lambda = 2.5){
  gen_M <- generate_M_list(d,K,r,lambda)
  true_label <- sample(c(1:K), n, replace = TRUE)
  pertub_label <- true_label
  ind <- sample(c(1:n),0.45*2/K*n)
  pertub_label[ind] <- change_label(true_label[ind])
  res_mat <- matrix(NA, nrow = ceiling(3*log(n))+1, ncol = 30)
  for (i in 1:30){
    data_tsr <- generate_data(gen_M$M_list, true_label, sigma)
    res_mat[,i] <- lrlloyds_iter_err(data_tsr, pertub_label, rep(r,K), true_label)
  }
  return(list("res_mat"=res_mat,"Delta"=gen_M$Delta))
}
res_K2 <- one_line()
# res_K3 <- one_line(K=3)
# res_K4 <- one_line(K=4)
# res_K5 <- one_line(K=5)
# res_K <- log(c(rowMeans(res_K2$res_mat),rowMeans(res_K3$res_mat),
#                rowMeans(res_K4$res_mat),rowMeans(res_K5$res_mat)))
# res_K <- cbind(rep(1:(ceiling(3*log(200))+1),4),res_K)
# res_K <- cbind(res_K, rep(c("K=2","K=3","K=4","K=5"), each = (ceiling(3*log(200))+1)))
# res_K_df <- as.data.frame(res_K)
# colnames(res_K_df) <- c("Iterations","Error","Ksize")
# res_K_df$Ksize <- as.factor(res_K_df$Ksize)
# res_K_df$Iterations <- as.integer(res_K_df$Iterations)
# res_K_df$Error <- as.numeric(res_K_df$Error)
# ggplot(res_K_df,aes(x=Iterations,y=Error,group=Ksize,color=Ksize)) +
#   geom_line(aes(x=Iterations,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=Iterations,y=Error, shape=Ksize),size=5) +
#   ylab("log(clustering error)") +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size=40),
#         legend.position = "right",
#         # legend.position = c(0.08, 0.2),
#         legend.background = element_rect(fill='transparent'),
#         axis.text.x = element_text(size=40),
#         axis.text.y = element_text(size=40),
#         axis.title.y = element_text(size=40, face="bold"),
#         axis.title.x = element_text(size=40, face="bold"))

# one_line <- function(n = 200, d = 50, r = 2, K = 3){
#   gen_P <- generate_P_list(d = 50, r = 2, K = 3)
#   true_label <- sample(c(1:K), n, replace = TRUE)
#   pertub_label <- true_label
#   ind <- sample(c(1:n),0.45*2/K*n)
#   pertub_label[ind] <- change_label(true_label[ind])
#   res_mat <- matrix(NA, nrow = ceiling(3*log(n))+1, ncol = 30)
#   for (i in 1:1){
#     data_tsr <- generate_network_data(gen_P$P_list, true_label)
#     res_mat[,i] <- lrlloyds_iter_err(data_tsr, pertub_label, rep(r,K), true_label)
#   }
#   return(list("res_mat"=res_mat,"Delta"=gen_M$Delta))
# }
# 
# temp <- one_line()

# res_lam1 <- one_line(lambda = 1.9)
# res_lam2 <- one_line(lambda = 2.1)
# res_lam3 <- one_line(lambda = 2.3)
# res_lam4 <- one_line(lambda = 2.5)

# res_lam <- log(c(rowMeans(res_lam1$res_mat),rowMeans(res_lam2$res_mat),
#                rowMeans(res_lam3$res_mat),rowMeans(res_lam4$res_mat)))
# res_lam <- cbind(rep(1:(ceiling(3*log(200))+1),4),res_lam)
# res_lam <- cbind(res_lam, rep(c("∆=4.23","∆=4.66","∆=5.11","∆=5.45"), each = (ceiling(3*log(200))+1)))
# res_lam_df <- as.data.frame(res_lam)
# colnames(res_lam_df) <- c("Iterations","Error","Lambda")
# res_lam_df$Lambda <- as.factor(res_lam_df$Lambda)
# res_lam_df$Iterations <- as.integer(res_lam_df$Iterations)
# res_lam_df$Error <- as.numeric(res_lam_df$Error)
# ggplot(res_lam_df,aes(x=Iterations,y=Error,group=Lambda,color=Lambda)) +
#   geom_line(aes(x=Iterations,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=Iterations,y=Error, shape=Lambda),size=5) +
#   ylab("log(clustering error)") +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size=40),
#         legend.position = "right",
#         # legend.position = c(0.08, 0.2),
#         legend.background = element_rect(fill='transparent'),
#         axis.text.x = element_text(size=40),
#         axis.text.y = element_text(size=40),
#         axis.title.y = element_text(size=40, face="bold"),
#         axis.title.x = element_text(size=40, face="bold"))

