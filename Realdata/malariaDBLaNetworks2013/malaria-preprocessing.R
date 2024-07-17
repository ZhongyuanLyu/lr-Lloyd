rm(list=ls())
library(igraph)
library(rTensor)

included=c(1,2,3,4,5,6,7,8,9)
for(i in included)
{
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

for(i in included)
{
  adjname= paste("A_HVR_", i, sep="")
  inds=match(net_nodes, row.names(get(adjname)))
  A_temp=get(adjname)[inds,]
  assign(adjname,A_temp[,inds])
  print(get(adjname)[1:4,1:4])
  netname= paste("Net_HVR_care_", i, sep="")
  assign(netname, graph_from_adjacency_matrix(get(adjname),mode="undirected"))
}


for(i in included)
{
  netname= paste("Net_HVR_", i, sep="")
  pdf(file=paste0(netname,".pdf"))
  plot(get(netname), vertex.size=4, vertex.color="red", vertex.label=NA)
  dev.off()
}



D=rep(0,length(included))
for (i in included)
{
  adjname= paste("A_HVR_", i, sep="")
  D[i]=sum(get(adjname))
}

included=c(1,2,3,4,5,6,7,8,9)
t=length(included)
indices <- c(n,n,t)
variable_list=paste0('A_HVR_',included)
mylist<-lapply(variable_list,get)
arr<-array(as.numeric(unlist(mylist)), dim=c(n, n,t))
arrT <- as.tensor(arr)



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
  tnsr_norm <- fnorm(tnsr)
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z, U_list) {
    est <- ttl(Z, U_list, ms = 1:num_modes)
    curr_resid <- fnorm(tnsr - est)
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

U_0=Initialzation(arrT,ranks=c(15,15,3))
decomp=PowerIteration(arrT,U_0,ranks=c(15,15,3),max_iter = 10000,tol=1e-05)
#nodes_embedding_Our=decomp[["U"]][[1]]
network_embedding_Our=decomp[["U"]][[3]]

plot(network_embedding_Our[,c(2,3)], cex=2, xlab="second eigen vector", ylab="third eigen vector")
text(network_embedding_Our[,2],network_embedding_Our[,3], as.character(c(1:9)), cex=2.5, pos=c(4,2,1,4,3,rep(1,3),3), col=c(rep("red",6),"black","green","blue"))

plot(network_embedding_Our[,c(2,3)]/network_embedding_Our[,1], cex=2, xlab="second eigen vector", ylab="third eigen vector")
text(network_embedding_Our[,2]/network_embedding_Our[,1],network_embedding_Our[,3]/network_embedding_Our[,1], as.character(c(1:9)), cex=2.5, pos=c(4,2,1,4,3,rep(1,3),3), col=c(rep("red",6),"black","green","blue"))


plot(network_embedding_Our[1:6,c(2,3)], cex=3,xlab="second eigen vector", ylab="third eigen vector")
text(network_embedding_Our[1:6,2],network_embedding_Our[1:6,3], as.character(c(1:6)), cex=3, pos=c(rep(3,5),1), col="red")

included=c(1,2,3,4,5,6) #1,,6,85,,8,6,7,9)
t=length(included)
indices <- c(n,n,t)
variable_list=paste0('A_HVR_',included)
mylist<-lapply(variable_list,get)
arr<-array(as.numeric(unlist(mylist)), dim=c(n, n,t))
arrT <- as.tensor(arr)
U_0=Initialzation(arrT,ranks=c(4,4,1))
decomp=PowerIteration(arrT,U_0,ranks=c(4,4,1),max_iter = 10000,tol=1e-05)
membership=decomp[["U"]][[1]]
change_point=decomp[["U"]][[3]]
plot(membership[,c(2,3)], main="Node embedding of t_(1)")
plot(membership[,c(3,4)], main="Node embedding of t_(1)")
plot(membership[,c(4,2)], main="Node embedding of t_(1)")
metadata_CysPoLv=read.table("metadata_CysPoLv.txt")
metadata_UPS=read.table("metadata_UPS.txt")
true_label=metadata_CysPoLv[which(row.names(metadata_CysPoLv) %in% net_nodes),]
true_label_UPS=metadata_UPS[which(row.names(metadata_UPS) %in% net_nodes),]

LP_7=eigen(A_HVR_6)
plot(LP_7$vectors[,c(2,3)], main="Node embedding of t_(1)")
plot(LP_7$vectors[,c(3,4)], main="Node embedding of t_(1)")
plot(LP_7$vectors[,c(4,2)], main="Node embedding of t_(1)")

### Visualization
library(plotly)
library(EMMIXskew)
fit_LP=kmeans(LP_7$vectors[,c(2,3,4)],4,iter.max = 10000)
fit_Tesor=kmeans(membership[,c(2,3,4)],4,iter.max = 10000)
plot_ly( x = membership[,2], y =membership[,3], z =membership[,4],size=2
         ,color=as.character(fit_Tesor$cluster))
plot_ly( x = LP_7$vectors[,2], y =LP_7$vectors[,3], z =LP_7$vectors[,4],size=2
         ,color=as.character(5-fit_LP$cluster))

library(dbscan)
dbscan_tensor=dbscan(membership[,c(2,3,4)],eps=.05,minPts = 5)
dbscan_LP6=dbscan(LP_7$vectors[,c(2,3,4)],eps=.04,minPts = 5)
plot_ly( x = membership[,2], y =membership[,3], z =membership[,4],size=2
         ,color=as.character(dbscan_tensor$cluster))
plot_ly( x = LP_7$vectors[,2], y =LP_7$vectors[,3], z =LP_7$vectors[,4],size=2
         ,color=as.character(dbscan_LP6$cluster))

node_membership=cbind(net_nodes,fit_Tesor$cluster)
rask_et_al_DBLa <- read.table("~/Google Drive/tensor for dynamic networks/real data/malariaDBLaNetworks2013/rask_et_al_DBLa.fasta", quote="\"", comment.char="")
n=dim(rask_et_al_DBLa)[1]
gene_inform=cbind(as.character(rask_et_al_DBLa[seq(from=1, to=n, by=2),1]),
                   as.character(rask_et_al_DBLa[seq(from=2, to=n, by=2),1]))

gene_membership=cbind(gene_inform[as.integer(node_membership[,1]),],node_membership)
gene_membership=gene_membership[,-3]
gene_membership=gene_membership[,c(1,3,2)]
colnames(gene_membership)=c("gene","group","sequence")
#write.csv(gene_membership,file="gene_membership.csv")


#Heatmap
A1=A_HVR_1
rownames(A1)=gene_membership[,1]
colnames(A1)=fit_Tesor$cluster
for(i in 2:6)
{
  adjname=paste0("A_HVR_",i)
  A1=A1+get(adjname)
}
ind_sort1=sort(fit_Tesor$cluster,index.return=T)$ix

A_temp=A1[ind_sort1,]
A1_sorted=A_temp[,ind_sort1]

palf <- colorRampPalette(c("white", "black")) 
heatmap(A1_sorted[,212:1], Rowv = NA, Colv = NA, col = palf(100),scale="none", margins=c(10,10))
