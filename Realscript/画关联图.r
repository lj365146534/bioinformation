library(igraph)
cells<-c(0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,1,1,0,3,0,3,3,3,0,0,0,0,0,0,0,0,3,0,3,1,1,1,0,0,0,0,0,0,1,1,0,3,0,0,0,0,1,0,0,0,0,0,1,1,3,1,0,0,3,0,0,0,0,0,0,0,0,0,3,1,0,3,0,0,3,1,0,3,0,0,1,1,3,1,0,0,0,0,0,3,0,3,1,1,0,0,0,0,1,3,3,0,0,3,1,3,0,0,0,0,0,0,0,0,1,3,3,0,0,3,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0)
cells=matrix(cells,14,14,byrow=T) #创建邻接矩阵
myCoord<-c(1, 2, 7.5, 5, 3, 6, 6, 8, 8,11, 8, 10, 11, 13,2, 1, 2, 4, 5.5, 1, 6, 4, 9,8, 14, 5.5, 2.2, 4)
myCoord<-matrix(myCoord,14,2,byrow=F) #创建顶点坐标
cnames=paste("e",1:14,sep="")
g=graph.adjacency(cells,mode="undirected",weighted=T)
plot(g,vertex.color="green",layout=myCoord,vertex.shape="square",
     vertex.label=cnames,vertex.label.font=2,vertex.label.dist=-1,
     vertex.label.degree=-pi/2,vertex.label.color="black",
     vertex.frame.color="gray",
  edge.width=E(g)$weight,edge.color="gray")
  