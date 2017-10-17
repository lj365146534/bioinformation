library(maps)
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(phangorn)
library(sysfonts)
library(showtext)
font.add("Calibri","Calibri.ttf")
font.add("Arial","Arial.ttf")

color <- c("black","orange","black,green","#FF00FF","#871F78")
month <- c("","_MILD","654654654654","79797","")
rootNode <- "A/CHICKEN/SHANDONG/01/2010"
size <- 6
m <- function(tr,node,a=0){
	i <- 0
	while (i<a){
		node <- tr$edge[which(tr$edge[,2] == node),1]
		i <- i+1
	}
	return(node)
}
seq <- read.dna(file= "C:/Users/Ja/Desktop/pdd/pb2_a.fas",format="fasta")
seq.name <- row.names(seq)
seq.colorframe <- data.frame(name=seq.name, color=length(color))

con <- function(y,x=seq.colorframe,w=1){

	for(i in y){
		a <- grep(i,x$name)
		x$color[a] <- w
		w <- w+1
	}
	return(x)
}

dd <- con(month)
dd$color <- as.character(dd$color)
fitseq <- phyDat(seq)
dm <- dist.ml(fitseq)
tree <- NJ(dm)
nodeNumber <- which(tree$tip == rootNode)
nodeNumber <- m(tree,nodeNumber,size)
nodeNumber
#tree <-reroot(tree,nodeNumber,0)
#tree <- reroot(tree,which(tree$tip==rootNode),0)
# NJtrees <- bootstrap.phyDat(fitseq, FUN=function(x)NJ(dist.logDet(x)), bs=100)
# tre <- plotBS(tree, NJtrees, phylogram)
# fit <- pml(tre, fitseq, k=4)
fit <- pml(tree, fitseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")
phangorn <- phyPML(fit, type="ml")
plot(1:4)
###一下参数为画大概100-150个序列进化树，氨基酸为蓝色，_MILD和其他为黑，_SEVERE为红
showtext.begin()
(ggtree(phangorn) + geom_text(aes(x=branch, label=AA_subs, vjust=.5),size=2.0,family="Calibri",color='blue')
	) 	%<+% dd + geom_tiplab(size=2.0,family="Arial",vjust=.5,aes(color=color))+ scale_color_manual(values=color)
showtext.end()


