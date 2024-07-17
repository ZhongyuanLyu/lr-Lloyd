rm(list=ls())
library(rTensor)
library(rgl)
library(GEOquery)
source("functions.R")
full_data <- getGEO(filename='GDS1083.soft')
data_df <- Table(full_data)
sample_info <- full_data@header[["sample_id"]]
sapply(strsplit(sample_info[1], split=",")[[1]], grepl, sample_info[10])
sample_IDs<- vector(mode = "list", length = 27)
for (i in 1:9){
  for (j in 1:3){
    sample_IDs[[3*(i-1)+j]] <- names(which(sapply(strsplit(sample_info[i], split=",")[[1]], grepl, sample_info[9+j])==TRUE))
  }
}
true_labels <- rep(c(1,2,3),9)
data_tsr <- array(NA, dim = c(1124,4,27))
for (i in 1:27){
  for (j in 1:4){
    data_tsr[,j,i] <- data_df[,sample_IDs[[i]][j]]
  }
}
# save(data_tsr,true_labels, file = "data_tsr.Rdata")




