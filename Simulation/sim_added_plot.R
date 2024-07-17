library(rTensor)
library(ggplot2)
library(latex2exp)
source("functions.R")
cal_err <- function(mat){
  mat_new <- matrix(NA, nrow = nrow(mat), ncol = 2)
  mat_new[,1]<- mat[,1]/mat[,2]
  mat_new[,2]<- mat[,5]
  return(mat_new)
}

# load(file = "error_mat_n.RData")
# error_matnew <- lapply(error_mat, cal_err)
# mean_err_mat <- apply(array(unlist(error_matnew), dim = c(20,2,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_matnew), dim = c(20,2,100)), c(1,2), var))
# mean_err_mat <- mean_err_mat[2:20,]
# sd_err_mat <- sd_err_mat[2:20,]
# num <- nrow(mean_err_mat)
# err_df <- matrix(NA, num*2, 5)
# err_df[1:num,2] <- mean_err_mat[,1]
# err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
# err_df[,1] <-rep(round(seq(1000,10000, length.out = 20))[-1],times = 2)
# err_df[,3] <- c(rep("Estimation",num),rep("Clustering",num))
# err_df[,4] <- c(mean_err_mat[,1]-2*sd_err_mat[,1],mean_err_mat[,2]-2*sd_err_mat[,2])
# err_df[,5] <- c(mean_err_mat[,1]+2*sd_err_mat[,1],mean_err_mat[,2]+2*sd_err_mat[,2])
# err_df <- as.data.frame(err_df)
# colnames(err_df) <- c("n","Error","Type","Lower","Upper")
# err_df$Type <- as.factor(err_df$Type)
# err_df$n <- as.integer(err_df$n)
# err_df$Error <- as.numeric(err_df$Error)
# ggplot(err_df,aes(x=n,y=Error,group=Type,color=Type)) +
#   geom_line(aes(x=n,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=n,y=Error, shape=Type),size=5) +
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
# 



# load(file = "error_mat_lam.RData")
# numrow <- nrow(error_mat[[1]])
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,3,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,3,100)), c(1,2), var))
# mean_err_mat <- mean_err_mat[2:numrow,]
# sd_err_mat <- sd_err_mat[2:numrow,]
# mean_err_mat <- mean_err_mat[,-2]
# sd_err_mat <- sd_err_mat[,-2]
# num <- nrow(mean_err_mat)
# err_df <- matrix(NA, num*2, 5)
# err_df[1:num,2] <- mean_err_mat[,1]
# err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
# err_df[,1] <-rep((seq(0.5,10, length.out=20)/0.13)[-1],times = 2)
# err_df[,3] <- c(rep("Estimation",num),rep("Clustering",num))
# err_df[,4] <- c(mean_err_mat[,1]-2*sd_err_mat[,1],mean_err_mat[,2]-2*sd_err_mat[,2])
# err_df[,5] <- c(mean_err_mat[,1]+2*sd_err_mat[,1],mean_err_mat[,2]+2*sd_err_mat[,2])
# err_df <- as.data.frame(err_df)
# colnames(err_df) <- c("Lambda","Error","Type","Lower","Upper")
# err_df$Type <- as.factor(err_df$Type)
# err_df$Lambda <- as.integer(err_df$Lambda)
# err_df$Error <- as.numeric(err_df$Error)
# ggplot(err_df,aes(x=Lambda,y=Error,group=Type,color=Type)) +
#   geom_line(aes(x=Lambda,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=Lambda,y=Error, shape=Type),size=5) +
#   xlab(TeX("$\\sigma_{r_1}(M_1)/sigma_{1}(M_2)$"))+
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



# load(file = "error_mat_lam_two_component.RData")
# numrow <- nrow(error_mat[[1]])
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), var))
# mean_err_mat <- mean_err_mat[1:11,]
# sd_err_mat <- sd_err_mat[1:11,]
# num <- nrow(mean_err_mat)
# err_df <- matrix(NA, num*2, 5)
# err_df[1:num,2] <- mean_err_mat[,1]
# err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
# err_df[,1] <-rep(round(seq(0.5,3, length.out=numrow)[1:11],digits = 1),times = 2)
# err_df[,3] <- c(rep("TS-Init",num),rep("lr-Lloyds",num))
# err_df[,4] <- c(mean_err_mat[,1]-2*sd_err_mat[,1],mean_err_mat[,2]-2*sd_err_mat[,2])
# err_df[,5] <- c(mean_err_mat[,1]+2*sd_err_mat[,1],mean_err_mat[,2]+2*sd_err_mat[,2])
# err_df <- as.data.frame(err_df)
# colnames(err_df) <- c("Lambda","Error","Type","Lower","Upper")
# err_df$Type <- as.factor(err_df$Type)
# err_df$Error <- as.numeric(err_df$Error)
# ggplot(err_df,aes(x=Lambda,y=Error,group=Type,color=Type)) +
#   geom_line(aes(x=Lambda,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=Lambda,y=Error, shape=Type),size=5) +
#   xlab(TeX("$\\lambda$"))+
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



# load(file = "error_mat_lam_two_component_high.RData")
# numrow <- nrow(error_mat[[1]])
# mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), mean)
# sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,2,100)), c(1,2), var))
# mean_err_mat <- mean_err_mat[9:20,]
# sd_err_mat <- sd_err_mat[9:20,]
# num <- nrow(mean_err_mat)
# err_df <- matrix(NA, num*2, 5)
# err_df[1:num,2] <- mean_err_mat[,1]
# err_df[(num+1):(num*2),2] <- mean_err_mat[,2]
# err_df[,1] <-rep(round(seq(1.5,3.8, length.out=numrow)[9:20],digits = 1),times = 2)
# err_df[,3] <- c(rep("TS-Init",num),rep("lr-Lloyds",num))
# err_df[,4] <- c(mean_err_mat[,1]-2*sd_err_mat[,1],mean_err_mat[,2]-2*sd_err_mat[,2])
# err_df[,5] <- c(mean_err_mat[,1]+2*sd_err_mat[,1],mean_err_mat[,2]+2*sd_err_mat[,2])
# err_df <- as.data.frame(err_df)
# colnames(err_df) <- c("Lambda","Error","Type","Lower","Upper")
# err_df$Type <- as.factor(err_df$Type)
# err_df$Error <- as.numeric(err_df$Error)
# ggplot(err_df,aes(x=Lambda,y=Error,group=Type,color=Type)) +
#   geom_line(aes(x=Lambda,y=Error), linetype="twodash",size=1.5) +
#   geom_point(aes(x=Lambda,y=Error, shape=Type),size=5) +
#   xlab(TeX("$\\lambda$"))+
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
# 
# seq(1.5,3.8, length.out=numrow)[9:20]



load(file = "error_mat_lam_two_component_compare.RData")
numrow <- nrow(error_mat[[1]])
mean_err_mat <- apply(array(unlist(error_mat), dim = c(numrow,3,100)), c(1,2), mean)
sd_err_mat <- sqrt(apply(array(unlist(error_mat), dim = c(numrow,3,100)), c(1,2), var))
mean_err_mat <- mean_err_mat[1:20,]
sd_err_mat <- sd_err_mat[1:20,]
num <- nrow(mean_err_mat)
err_df <- matrix(NA, num*2, 5)
err_df <- as.data.frame(err_df)
err_df[1:num,2] <- mean_err_mat[,2]
err_df[(num+1):(num*2),2] <- mean_err_mat[,3]
# err_df[(num*2+1):(num*3),2] <- mean_err_mat[,3]
err_df[,1] <-rep(round(seq(0.5,3, length.out=numrow)[1:20],digits = 1),times = 2)
err_df[,3] <- c(rep("lr-Lloyds",num), rep("MatTransMix",num))

colnames(err_df) <- c("Lambda","Error","Type","Lower","Upper")
err_df$Type <- as.factor(err_df$Type)
# err_df$Error <- as.numeric(err_df$Error)
ggplot(err_df,aes(x=Lambda,y=Error,group=Type,color=Type)) +
  geom_line(aes(x=Lambda,y=Error), linetype="twodash",size=1.5) +
  geom_point(aes(x=Lambda,y=Error, shape=Type),size=5) +
  xlab(TeX("$\\lambda$"))+
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


      