##########################################################################
# 画两组间GSEA分析的NES得分的棒棒糖图
library(ggplot2)
df = read.table("./GSEA_GOBP_barplot.txt", header=T, sep="\t")
p = ggplot(df, aes(x=NES, y=reorder(Description, NES), color=pvalue))+
geom_point(aes(size=abs(NES)))+
geom_segment(aes(x=0, xend=NES, y=Description, yend=Description), color="black", linewidth=1)+
scale_color_gradient(low="hotpink1", high="steelblue1")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))
ggsave("GSEA_GOBP_barplot.pdf", p, width=6, height=6)
