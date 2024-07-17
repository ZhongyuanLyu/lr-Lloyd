library(rTensor)
library(MatTransMix)
source("libMatTransFull.R")
source("functions.R")

##======= low-dimensional setting =======##
d = 5; n=200; r = 2; K = 2; sigma = 1;  num <- 20;
time_mat <- matrix(NA, nrow = 100, ncol = 10)
rU <- K*r; rV <- K*r;
true_label <- 2*rbinom(n,1,0.5)-1
M <- generate_M_list_lam(d,r,1.5, kappa=1.3)
X_array <- generate_data_X(M[[1]],true_label, n,d)
X_tsr <- as.tensor(X_array);
s_spec_init <- lr_spec_init(X_tsr,rU,rV,K)
s_mat_init <- MatTrans.init(X_array, K = 2, n.start = 1)
mattrans_res <- MatTrans.EM(X_array, initial = s_mat_init, trans = "None", silent = TRUE, size.control = 10)
for (i in 1:100) {
  time_mat[i,1:5] <- system.time(lrlloyds(X_tsr, s_spec_init, rep(r,K), TRUE))
  time_mat[i,6:10] <- system.time(MatTrans.EM(X_array, initial = s_mat_init, trans = "None", silent = TRUE, size.control = 30))
}

# save(time_mat, file = "time_mat_compare.RData")

time_df <- matrix(NA, nrow(time_mat)*2, 2)
time_df[1:nrow(time_mat),1] <- time_mat[,3]
time_df[(nrow(time_mat)+1):(2*nrow(time_mat)),1] <- time_mat[,8]
time_df <- as.data.frame(time_df)
time_df[,2] <- c(rep("lr-Lloyds",nrow(time_mat)), rep("MatTransMix",nrow(time_mat)))
colnames(time_df) <- c("Time","Method")

get_box_stats <- function(y, upper_limit = max(time_df$Time) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "s\n",
      "Median =", round(median(y), 2), "s\n"
    )
  ))
}
ggplot(time_df, aes(x=Method, y=Time, fill=Method)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 2, size = 15) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=40),
        legend.position = "right",
        legend.background = element_rect(fill='transparent'),
        axis.text.x = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.y = element_text(size=40, face="bold"),
        axis.title.x = element_text(size=40, face="bold"))
  
