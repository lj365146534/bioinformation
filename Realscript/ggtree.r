
library(ggplot2)
library(ggtree)
library(phangorn)


#��ȡ�ļ�����1#ע�������������bootstrapֵ
	treefile <- "D:\\����\\H9N2\\20160809meeting\\RAxML_bestTree.H9.nwk"
	tipseqfile <- "D:\\����\\H9N2\\20160809meeting\\H9N2-HA.fas"
#��ȡ�ļ�����2
	treefile <- system.file("extdata", "ha.nwk", package="ggtree")
#��ȡnwk������
tre <- read.tree(treefile)


tipseqfile <- system.file("extdata", "ha.fas", package="ggtree")

#��ȡfas�ļ�
tipseq <- read.phyDat(tipseqfile,format="fasta")


fit <- pml(tre, tipseq, k=4)
fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
                 optInv=T, optGamma=T, optEdge=TRUE,
                 optRooted=FALSE, model = "GTR")
p <- phyPML(fit, type="ml")
ggtree(phangorn) + geom_text(aes(x=branch, label=AA_subs, vjust=-.5))

ggtree(phangorn) + geom_text(aes(x=branch, label=gsub("/","",AA_subs), vjust=-.5), color='blue', size=1.5)+ geom_tiplab(grep("Jul",label),color='green', size=2)



###������,��������genotype.txt�еĶ������Ʊ���һ��,���಻�١�
library(ggplot2)
library(ggtree)
treefile <- system.file("extdata", "gH5/2344.nwk", package="ggtree")
tre <- read.tree(treefile)

genotype_file <- system.file("extdata/gH5/genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)

library(phytools)
#��ʾ���ж�������
tre$tip.label
#������һ��Ϊ��
tree <- reroot(tre,120)
#���ƽ�����
p <- ggtree(tree) + geom_tiplab(size=1)
p <- ggtree(tree)
gheatmap(p, genotype, offset = 0, width=0.5)


###H7N9���������λ��##
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
