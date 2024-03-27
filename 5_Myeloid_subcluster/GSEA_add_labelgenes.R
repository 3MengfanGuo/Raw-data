# 在GSEA中标注候选基因
# R包GseaVis能轻松在GSEA标准图上标注候选基因，同时也支持点阵图样式
# 部分所需R包安装与载入：
devtools::install_github("junjunlab/GseaVis")
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GseaVis)

# 1.先使用clusterProfiler包完成GSEA富集
data = read.table("D:/data.txt", header=T, row.names=1, sep="\t")
head(data)
#            baseMean log2FoldChange     lfcSE       stat       pvalue         padj
# MARCH1   258.031664    -0.60749636 0.2377749 -2.5549225 1.062114e-02 2.951119e-02
# MARCH10    5.590111    -0.97812598 1.4549549 -0.6722724 5.014103e-01 6.326400e-01
# MARCH2   443.304767    -0.09801133 0.2223848 -0.4407285 6.594096e-01 7.674181e-01
# MARCH3   116.590427     2.28195037 0.3795243  6.0126595 1.825041e-09 5.024491e-08
# MARCH4     4.369855     0.85402077 1.8312338  0.4663636 6.409553e-01 7.522414e-01
# MARCH5  1728.756453     0.49617999 0.1397650  3.5501017 3.850823e-04 1.850531e-03

# symbol转entrez ID：
symbol <- rownames(data)
entrez <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(entrez)
#    SYMBOL ENTREZID
# 11   A1BG        1
# 12   A1CF    29974
# 13    A2M        2
# 14 A4GALT    53947
# 15  A4GNT    51146
# 16   AAAS     8086

# 准备genelist文件(entrez+log2FC)：
genelist <- data$avg_logFC
names(genelist) <- rownames(data)
#过滤掉ID转换中缺失部分基因：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]

KEGG_ges <- gseKEGG(
geneList=genelist,
organism = "hsa",
keyType = "kegg",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
pAdjustMethod="BH",
verbose=FALSE,
eps=0
)

##############################################################################################
df = read.table("./RSA_FC_sorted.txt", header=T, row.names=1, sep="\t")
genelist = df$FC
names(genelist) = rownames(df)

# GSEA_KEGG富集分析：
R.utils::setOption("clusterProfiler.download.method", "auto") ##如果富集时报错就加上这句代码
KEGG_ges <- gseGO(
geneList=genelist,
ont="BP",
OrgDb=org.Hs.eg.db,
keyType="SYMBOL",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
pAdjustMethod="BH",
verbose=FALSE,
eps=0
)
#################################################################################
#将entrez重新转换为symbol：
KEGG_ges2 <- setReadable(KEGG_ges, OrgDb=org.Hs.eg.db, keyType="ENTREZID")

core <- result$core_enrichment[15]
core_genes <- str_split(core ,'/')[[1]]
gseaNb(
    object = KEGG_ges,
    geneSetID = result$ID[15],
    addPval = T,
    pvalX = 0.95,
    pvalY = 0.8,
    newGsea = T,
    addPoint = F,
    newCurveCol = c("dodgerblue1","#b9b4ad", "darkorange1"),
    newHtCol = c("darkorange1", "#b9b4ad", "dodgerblue1"),
    addGene = core_genes, #选择添加我们定义的基因集
    kegg = F,
    geneCol = 'black'
)

##########################################################################
library(ggplot2)
df = read.table("./GO_GSEA_barplot.txt", header=T, sep="\t")
ggplot(df, aes(x=NES, y=reorder(Description, NES), color=pvalue))+
geom_point(aes(size=abs(NES)))+
geom_segment(aes(x=0, xend=NES, y=Description, yend=Description), color="black", linewidth=1)+
scale_color_gradient(low="darkorange1", high="dodgerblue1")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))




