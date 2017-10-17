x <-read.table("aa.txt")
y <- as.matrix(x)
for ( i in 1:8){
	for (j in 1:8){
		if(i >j){
			y[i,j] <- y[j,i]
		}
	}
}

for ( i in 1:8){
	for (j in 1:8){
		if(y[i,j] < 95){
			y[i,j] <- 95
		}
	}
}
p <- ggplot(melt(y),aes(Var2,Var1))+geom_tile(aes(fill = value),colour= "white")+ scale_fill_gradient(low ="red",high="green",limits=c(95,100))+
	scale_x_discrete(expand = c(0, 0),labels=colnames(y),breaks=1:8)+scale_y_discrete(expand = c(0, 0),labels=rownames(y),breaks=1:8)+labs(x="", y="")
pdf(file = "C:/Users/trainrun/Desktop/yang/MP.pdf",width = 8,height=8)
print(p)
dev.off();
