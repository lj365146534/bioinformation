###H7N9进化树标记位点##
library(ggplot2)
library(ggtree)
library(phangorn)

###ha##
#treefile <- system.file("extdata", "ha.nwk", package="ggtree")
#tipseqfile <- system.file("extdata", "ha.fas", package="ggtree")
tree <- read.tree("C:/Users/Ja/Desktop/pdd/pb2_RAxML_bestTree.result")
tipseq <- read.phyDat("C:/Users/Ja/Desktop/pdd/pb2_a.fas",format="fasta")

#root枝的名字
#rootNode <- "A/SHANGHAI/05/2013"
rootNode <- "A/CHICKEN/ZHEJIANG/HJ/2007"
#算出root枝的node
#nodeNumber <- which(tree$tip == rootNode)
#tree <-reroot(tree,nodeNumber,0)
tree <- reroot(tree,which(tree$tip==rootNode),0)


fit <- pml(tree, tipseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")
p <- phyPML(fit, type="ml")


pdf("pb2.pdf", width= 15, height = 150)
ggtree(p) + geom_text(aes(x=branch, label=gsub("/","",AA_subs), vjust=-.5), color='blue', size=1.5)+ geom_tiplab(size=2) + scale_x_continuous(limits=c(0,0.08))
dev.off()
