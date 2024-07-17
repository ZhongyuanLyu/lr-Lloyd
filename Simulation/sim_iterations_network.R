library(rTensor)
library(ggplot2)
library(latex2exp)
source("functions.R")
change_label <- function(x){
  return(ifelse(x==1,2,1))
}
set.seed(12345)
one_line <- function(n = 200, d = 50, r = 2, K = 3, maxp=0.05, gap_ratio = 0.5){
  gen_P <- generate_P_list(d, r, K, maxp, gap_ratio)
  true_label <- sample(c(1:K), n, replace = TRUE)
  pertub_label <- true_label
  ind <- sample(c(1:n),0.45*2/K*n)
  pertub_label[ind] <- change_label(true_label[ind])
  res_mat <- matrix(NA, nrow = ceiling(3*log(n))+1, ncol = 30)
  for (i in 1:30){
    data_tsr <- generate_network_data(gen_P$P_list, true_label)
    res_mat[,i] <- lrlloyds_iter_err(data_tsr, pertub_label, rep(r,K), true_label)
  }
  return(list("res_mat"=res_mat,"Delta"=gen_P$Delta))
}


res_p1 <- one_line(maxp=0.05)
res_p2 <- one_line(maxp=0.08)
res_p3 <- one_line(maxp=0.1)
res_p4 <- one_line(maxp=0.15)

res_p <- log(c(rowMeans(res_p1$res_mat),rowMeans(res_p2$res_mat),
               rowMeans(res_p3$res_mat),rowMeans(res_p4$res_mat)))
res_p <- cbind(rep(1:(ceiling(3*log(200))+1),4),res_p)
res_p <- cbind(res_p, rep(c("∆=0.72","∆=1.19","∆=1.46","∆=2.15"), each = (ceiling(3*log(200))+1)))
res_p_df <- as.data.frame(res_p)
colnames(res_p_df) <- c("Iterations","Error","Delta")
res_p_df$Delta <- as.factor(res_p_df$Delta)
res_p_df$Iterations <- as.integer(res_p_df$Iterations)
res_p_df$Error <- as.numeric(res_p_df$Error)
ggplot(res_p_df,aes(x=Iterations,y=Error,group=Delta,color=Delta)) +
  geom_line(aes(x=Iterations,y=Error), linetype="twodash",size=1.5) +
  geom_point(aes(x=Iterations,y=Error, shape=Delta),size=5) +
  ylab("log(clustering error)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=40),
        legend.position = "right",
        # legend.position = c(0.08, 0.2),
        legend.background = element_rect(fill='transparent'),
        axis.text.x = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.y = element_text(size=40, face="bold"),
        axis.title.x = element_text(size=40, face="bold"))

