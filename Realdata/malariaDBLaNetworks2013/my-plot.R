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
  pdf(file=paste0(netname,"blue.pdf"))
  plot(get(netname), vertex.size=3, vertex.color="SkyBlue2", vertex.label=NA)
  dev.off()
}


