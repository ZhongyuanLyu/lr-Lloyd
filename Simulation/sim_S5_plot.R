library(rTensor)
library(ggplot2)
library(latex2exp)
source("functions.R")

load(file = "error_mat_r.RData")
numrow <- nrow(error_mat[[1]])
mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), var))
num <- nrow(mean_err_mat)
err_df <- matrix(NA, num*4, 3)
err_df[1:num,2] <- mean_err_mat[,2]
load(file = "error_mat_r_SNR_med.RData")
mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
load(file = "error_mat_r_SNR_medhigh.RData")
mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
err_df[(2*num+1):(num*3),2] <- mean_err_mat[,2]
load(file = "error_mat_r_SNR_low.RData")
mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
err_df[(3*num+1):(num*4),2] <- mean_err_mat[,2]
err_df[,1] <-rep(round(seq(3,10),digits = 1),times = 4)
err_df[,3] <- c(rep(c("SNR-0.8"),num),rep("SNR-0.7",num),rep("SNR-0.65",num), rep("SNR-0.6",num))
err_df <- as.data.frame(err_df)


colnames(err_df) <- c("rc","Error","Type")
err_df$Type <- as.factor(err_df$Type)
err_df$Error <- as.numeric(err_df$Error)
err_df$rc <- as.numeric(err_df$rc)
ggplot(err_df,aes(x=rc,y=Error,group=Type,color=Type)) +
  geom_line(aes(x=rc,y=Error), linetype="twodash",size=1.5) +
  geom_point(aes(x=rc,y=Error, shape=Type),size=5) +
  xlab(TeX("$\\r_c$"))+
  ylab("Error") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=40),
        legend.position = "right",
        legend.background = element_rect(fill='transparent'),
        axis.text.x = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.y = element_text(size=40, face="bold"),
        axis.title.x = element_text(size=40, face="bold"))


# numrow <- nrow(error_mat[[1]])
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), var))
# num <- nrow(mean_err_mat)
# err_df <- matrix(NA, num*2, 5)
# err_df[1:num,2] <- mean_err_mat[,1]
# err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
# err_df[,1] <-rep(round(seq(3,10),digits = 1),times = 2)
# err_df[,3] <- c(rep("TS-Init",num),rep("lr-Lloyds",num))
# err_df[,4] <- c(mean_err_mat[,1]-2*sd_err_mat[,1],mean_err_mat[,2]-2*sd_err_mat[,2])
# err_df[,5] <- c(mean_err_mat[,1]+2*sd_err_mat[,1],mean_err_mat[,2]+2*sd_err_mat[,2])
# err_df <- as.data.frame(err_df)
# colnames(err_df) <- c("rc","Error","Type","Lower","Upper")
# err_df$Type <- as.factor(err_df$Type)
# err_df$Error <- as.numeric(err_df$Error)
# err_df$rc <- as.numeric(err_df$rc)
# ggplot(err_df,aes(x=rc,y=Error,group=Type,color=Type)) +
#   geom_line(aes(x=rc,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=rc,y=Error, shape=Type),size=5) +
#   xlab(TeX("$\\r_c$"))+
#   ylab("Error") +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size=40),
#         legend.position = "right",
#         legend.background = element_rect(fill='transparent'),
#         axis.text.x = element_text(size=40),
#         axis.text.y = element_text(size=40),
#         axis.title.y = element_text(size=40, face="bold"),
#         axis.title.x = element_text(size=40, face="bold"))
