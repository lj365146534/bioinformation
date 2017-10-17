library(ggplot2)
library(reshape2)
a <- read.table("all.txt",header = T)
#a1 <- melt(a, varnames=c(SitesNumber))
a2 <- ggplot(a, aes(variable, SitesNumber)) + geom_tile()
a3 <- a2 + facet_grid(. ~ Segment, scale = "free_x", space = "free")
a4 <- a3 + aes(fill=factor(value)) + scale_fill_brewer(palette = "RdPu")
a5 <- a4 + theme(axis.text.x =element_blank())
###去掉X轴内容
pdf("test.pdf", width= 100, height = 12)
a5
dev.off()

##首先图中不显示横纵坐标签即plot(...,xaxt="n",yaxt="n")，再使用axis(1,c(1,2,3),labels=c("a","b","c")
## brewer.pal(7,"RdPu")
##[1] "#FEEBE2" "#FCC5C0" "#FA9FB5" "#F768A1" "#DD3497" "#AE017E" "#7A0177"
##> display.brewer.pal(7,"RdPu")           显示色彩颜色


#colour = c("white","#FCC5C0","#F768A1","#AE017E","#7A0177")