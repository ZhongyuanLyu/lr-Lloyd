rm(list=ls())
library(igraph)
library(rTensor)
source("functions.R")
included=c(1,2,3,4,5,6,7,8,9)
for(i in included){
  oname = paste("HVR_", i, sep="")
  netname= paste("Net_HVR_", i, sep="")
  adjname= paste("A_HVR_", i, sep="")
  assign(oname, read.csv(paste(oname, ".txt", sep=""),stringsAsFactors=F))
  temp=get(oname)
  G=graph_from_data_frame(temp,directed=F)
  assign(netname,G )
  assign(adjname, as_adjacency_matrix(G,type="both",sparse=F))
} # load 9 edge list

net_nodes=row.names(A_HVR_5)
for(i in included)
{
  adjname= paste("A_HVR_", i, sep="")
  net_nodes=intersect(net_nodes,row.names(get(adjname)))
}
n=length(net_nodes)

for(i in included){
  adjname= paste("A_HVR_", i, sep="")
  inds=match(net_nodes, row.names(get(adjname)))
  A_temp=get(adjname)[inds,]
  assign(adjname,A_temp[,inds])
  print(get(adjname)[1:4,1:4])
  netname= paste("Net_HVR_care_", i, sep="")
  assign(netname, graph_from_adjacency_matrix(get(adjname),mode="undirected"))
}

for(i in included){
  netname= paste("Net_HVR_", i, sep="")
  pdf(file=paste0(netname,".pdf"))
  plot(get(netname), vertex.size=4, vertex.color="red", vertex.label=NA)
  dev.off()
}
D=rep(0,length(included))
for (i in included){
  adjname= paste("A_HVR_", i, sep="")
  D[i]=sum(get(adjname))
}
included=c(1,2,3,4,5,6,7,8,9)
t=length(included)
indices <- c(n,n,t)
variable_list=paste0('A_HVR_',included)
mylist<-lapply(variable_list,get)
arr<-array(as.numeric(unlist(mylist)), dim=c(n, n,t))
arr_six <- arr[,,1:6]
# data_tsr <- as.tensor(arr_six)
data_tsr <- as.tensor(arr)
princomp(t(k_unfold(data_tsr, m=1)@data))
screeplot(princomp(t(k_unfold(data_tsr, m=1)@data)), npcs = 25,
          type = "lines", main = "scree plot")
screeplot(princomp(t(k_unfold(data_tsr, m=2)@data)), npcs = 25,
          type = "lines", main = "scree plot")


n <- data_tsr@modes[3]; rU <- 15;rV <- 15; K <- 6;

## Spectral Initialization
HOSVD_init <- hosvd(data_tsr, ranks = c(rU, rV, K))
U_int <- HOSVD_init$U[[1]]
V_int <- HOSVD_init$U[[2]]
Ghat <- ttl(data_tsr, list(U_int%*%t(U_int),V_int%*%t(V_int)), c(1,2))
M3Ghat <- k_unfold(Ghat, m=3)@data
screeplot(princomp(t(M3Ghat)), npcs = 25, type = "lines", main = "scree plot")
princomp(t(M3Ghat))
s_spec_init <- kmeans(M3Ghat, K, nstart = 50)$cluster
s_spec_init
s_vec_lloyds <- kmeans(k_unfold(data_tsr, m=3)@data, K, nstart = 50)$cluster
s_vec_lloyds
effective_hamming_error(s_spec_init,s_vec_lloyds)
## Low-rank Lloyd's Algorithm
## Determine the rank
Mean_mat_list <- vector(mode = "list", length = K)
for (k in 1:K) {
  Mean_mat_list[[k]] <- apply(data_tsr@data[,,s_spec_init==k], c(1,2), mean)
}
screeplot(princomp(Mean_mat_list[[1]]), npcs = 20, 
          type = "lines", main = "scree plot")
screeplot(princomp(Mean_mat_list[[2]]), npcs = 20, 
          type = "lines", main = "scree plot")
## Iterations
## spectral initialization
rc <- 6
lrlloyds_res <- lrlloyds(data_tsr, s_spec_init, c(rc,rc,rc,rc,rc,rc), TRUE)
lrlloyds_res$shat
effective_hamming_error(lrlloyds_res$shat,s_vec_lloyds)
# shat_vec <- lloyds(data_tsr, s_init)
# effective_hamming_error(shat_vec,true_label+1)


###Functions
Initialzation<-function(tnsr, ranks=NULL)
{
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  temp_mat <-matrix(0,ncol=tnsr@modes[1], nrow=tnsr@modes[2])
  for(i in 1:tnsr@modes[3])
  {
    temp_mat=temp_mat+tnsr@data[,,i]
  }
  U_list[[1]] <- eigen(temp_mat,symmetric = T)$vector[,c(1:ranks[1])]
  U_list[[2]] <- U_list[[1]]
  outer_production=NULL
  for(i in 1:dim(U_list[[1]])[1])
  {
    row_matrix=NULL
    for (j in 1:dim(U_list[[1]])[2])
    {
      temp=U_list[[1]][i,j]*U_list[[2]]
      row_matrix=cbind(row_matrix,temp)
    }
    outer_production=rbind(outer_production,row_matrix)
  }
  # outer_production <- outer(U_list[[1]],U_list[[2]])
  temp_mat <- rs_unfold(tnsr, m = 3)@data %*% outer_production
  U_list[[3]] <- svd(temp_mat, nu = ranks[3])$u
  return(U_list)
}

PowerIteration<- function(tnsr, ranks=NULL, U_0_list, max_iter = 25, tol = 1e-05)
{
  stopifnot(is(tnsr, "Tensor"))
  if (is.null(ranks)) 
    stop("ranks must be specified")
  if (sum(ranks > tnsr@modes) != 0) 
    stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks <= 0) != 0) 
    stop("ranks must be positive")
  num_modes <- tnsr@num_modes
  U_list <- U_0_list
  tnsr_norm <- fnorm(tnsr@data)
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z, U_list) {
    est <- ttl(Z, U_list, ms = 1:num_modes)
    curr_resid <- fnorm(tnsr@data - est@data)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1) 
      return(FALSE)
    if (abs(curr_resid - fnorm_resid[curr_iter - 1])/tnsr_norm < 
        tol) 
      return(TRUE)
    else {
      return(FALSE)
    }
  }
  pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  while ((curr_iter < max_iter) && (!converged)) {
    setTxtProgressBar(pb, curr_iter)
    modes <- tnsr@modes
    modes_seq <- 1:num_modes
    for (m in modes_seq) {
      X <- ttl(tnsr, lapply(U_list[-m], t), ms = modes_seq[-m])
      U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
    }
    Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)
    if (CHECK_CONV(Z, U_list)) {
      converged <- TRUE
      setTxtProgressBar(pb, max_iter)
    }
    else {
      curr_iter <- curr_iter + 1
    }
  }
  close(pb)
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1)/tnsr_norm)) * 
    100
  est <- ttl(Z, U_list, ms = 1:num_modes)
  invisible(list(Z = Z, U = U_list, conv = converged, est = est, 
                 norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid, 
                                                                 1), all_resids = fnorm_resid))
}

U_0=Initialzation(data_tsr,ranks=c(15,15,3))
decomp=PowerIteration(data_tsr,ranks=c(15,15,3),U_0,max_iter = 10000,tol=1e-05)

network_embedding_Our=decomp[["U"]][[3]]
s_twist <- kmeans(network_embedding_Our, K, nstart = 50)$cluster
s_twist
