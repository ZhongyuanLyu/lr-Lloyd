library(tidyr)
extract_sample <- function(file_directory, keyword = "S1 obj"){
  file_list <- list.files(file_directory)
  size <- length(file_list)
  check_name <- vector(mode = "character", length = size)
  for (i in 1:size){
    check_name[i] <- readLines(paste(file_directory,file_list[i], sep = "/"), n=4)[4]
  }
  file_S1_list <-file_list[grepl(keyword, check_name)]
  n0 <- length(file_S1_list); skip <- 0
  average_df <- read.table(paste(file_directory,file_S1_list[1], sep = "/"))
  average_df$V4 <- 0
  for (i in 1:length(file_S1_list)){
    temp <- read.table(paste(file_directory,file_S1_list[i], sep = "/"))
    if (mean(temp$V4==0)>0.95){
      n0 <- n0 - 1; skip <- skip + 1; next
    } else{
      average_df$V4 <- average_df$V4 + temp$V4
    }
  }
  # for (i in 1:length(file_S1_list)){
  #   temp <- read.table(paste(file_directory,file_S1_list[i], sep = "/"))
  #   average_df$V4 <- average_df$V4 + temp$V4
  # }
  average_df$V4 <-average_df$V4/n0
  avg_mat<- spread(average_df[,-1], key = V2, value = V4)
  avg_mat <- as.matrix(avg_mat[,-1])
  cat("We skip:",skip,"samples \n")
  return(avg_mat)
}

file_directory_list <- paste("eeg_full",list.files("eeg_full"), sep ="/")
n <- length(file_directory_list)
data_tsr<-array(0, dim=c(256, 64, n))
for (i in 1:n){
  data_tsr[,,i] <- extract_sample(file_directory_list[i], "trial ")
}
true_label <- grepl("a0000", file_directory_list)
# save(data_tsr,true_label, file = "data_tsr_S2match_nonzero.Rdata")
# save(data_tsr,true_label, file = "data_tsr_all_nonzero.Rdata")

# file_directory_list <- paste("eeg_full",list.files("eeg_full"), sep ="/")
# n <- length(file_directory_list)
# file_list <- list.files(file_directory_list[103])
# size <- length(file_list)
# check_name <- vector(mode = "character", length = size)
# for (i in 1:size){
#   check_name[i] <- readLines(paste(file_directory_list[103],file_list[i], sep = "/"), n=4)[4]
# }
