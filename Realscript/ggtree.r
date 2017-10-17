
library(ggplot2)
library(ggtree)
library(phangorn)


#读取文件方法1#注意进化树不能有bootstrap值
	treefile <- "D:\\工作\\H9N2\\20160809meeting\\RAxML_bestTree.H9.nwk"
	tipseqfile <- "D:\\工作\\H9N2\\20160809meeting\\H9N2-HA.fas"
#读取文件方法2
	treefile <- system.file("extdata", "ha.nwk", package="ggtree")
#读取nwk进化树
tre <- read.tree(treefile)


tipseqfile <- system.file("extdata", "ha.fas", package="ggtree")

#读取fas文件
tipseq <- read.phyDat(tipseqfile,format="fasta")


fit <- pml(tre, tipseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")
p <- phyPML(fit, type="ml")
ggtree(phangorn) + geom_text(aes(x=branch, label=AA_subs, vjust=-.5))

ggtree(phangorn) + geom_text(aes(x=branch, label=gsub("/","",AA_subs), vjust=-.5), color='blue', size=1.5)+ geom_tiplab(grep("Jul",label),color='green', size=2)



###基因型,进化树与genotype.txt中的毒株名称必须一致,不多不少。
library(ggplot2)
library(ggtree)
treefile <- system.file("extdata", "gH5/2344.nwk", package="ggtree")
tre <- read.tree(treefile)

genotype_file <- system.file("extdata/gH5/genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)

library(phytools)
#显示所有毒株名称
tre$tip.label
#以其中一个为根
tree <- reroot(tre,120)
#绘制进化树
p <- ggtree(tree) + geom_tiplab(size=1)
p <- ggtree(tree)
gheatmap(p, genotype, offset = 0, width=0.5)


###H7N9进化树标记位点##
library(ggplot2)
library(ggtree)
library(phangorn)

###ha##
treefile <- system.file("extdata", "ha.nwk", package="ggtree")
tipseqfile <- system.file("extdata", "ha.fas", package="ggtree")
tre <- read.tree(treefile)
tipseq <- read.phyDat(tipseqfile,format="fasta")

###na##
treefile <- system.file("extdata", "na.nwk", package="ggtree")
tre <- read.tree(treefile)
tipseqfile <- system.file("extdata", "na.fas", package="ggtree")
tipseq <- read.phyDat(tipseqfile,format="fasta")

tre$tip.label
tree <- reroot(tre,32)

fit <- pml(tree, tipseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")
p <- phyPML(fit, type="ml")

pdf("test.pdf", width= 10, height = 100)
ggtree(p) + geom_text(aes(x=branch, label=gsub("/","",AA_subs), vjust=-.5), color='blue', size=1.5)+ geom_tiplab(size=2)
dev.off()
