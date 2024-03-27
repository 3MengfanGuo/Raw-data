# 在GSEA中标注候选基因
# R包GseaVis能轻松在GSEA标准图上标注候选基因，同时也支持点阵图样式
# 部分所需R包安装与载入：
devtools::install_github("junjunlab/GseaVis")
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GseaVis)

# 1.计算两组间的基因的FoldChange
mydata = readRDS("../1_preprocessing/mydata_cluster.rds")
mydata = subset(mydata, cell_type=="Cytotoxic NK/T cells")
data = FindMarkers(mydata, group.by="Type", ident.1="metastasis", ident.2="primary tumor", logfc.threshold=0, min.pct=0)
write.table(data, "./Cytotoxic_NKT_metastasis_primary_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")

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
genelist <- data$log2FoldChange
names(genelist) <- rownames(data)
#过滤掉ID转换中缺失部分基因：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]
genelist = sort(genelist, decreasing=T)

##############################################################################################
# GSEA_GO富集分析：
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
# 画带有基因标签的GSEA分析图
KEGG_ges2 <- setReadable(KEGG_ges, OrgDb=org.Hs.eg.db, keyType="SYMBOL")
core <- result$core_enrichment[78]
core_genes <- str_split(core ,'/')[[1]]
p = gseaNb(
    object = KEGG_ges,
    geneSetID = result$ID[78],
    addPval = T,
    pvalX = 0.95,
    pvalY = 0.8,
    newGsea = T,
    addPoint = F,
    newCurveCol = c("steelblue1","lightgray", "hotpink1"),
    newHtCol = c("hotpink1", "lightgray", "steelblue1"),
    addGene = core_genes, #选择添加我们定义的基因集
    kegg = F,
    geneCol = 'black'
)
ggsave("GSEA_Tcell_Activation.pdf", p, width=6, height=6)




