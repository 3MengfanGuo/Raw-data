# https://cloud.tencent.com/developer/article/1692240
# https://zhuanlan.zhihu.com/p/358986392
library(Seurat)
library(doParallel)
library(SCopeLoomR)
mydata = readRDS("./ADGRL4_cells.rds")
exprMat <- as.matrix(mydata@assays$RNA@counts)
cellInfo = mydata@meta.data[, c(14, 3, 2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')

### Initialize settings
library(SCENIC)
# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
# 设置分析环境
mydbDIR <- "./database"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", dbDir=mydbDIR, dbs=mydbs, nCores=4, datasetTitle="Epithelial cells") 
saveRDS(scenicOptions, "int/scenicOptions.rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
## 计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
## TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts=20)
## 需要注意nParts参数，它的作用是把表达矩阵分成n份分开计算，目的是防止数据量大时内存不够。

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top10perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")        # choose settings

####################################################################################################################
# Seurat可视化SCENIC结果
# 把SCENIC结果中最重要的regulonAUC矩阵导入Seurat，这样就与之前的分析联系起来

##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
mydata <- AddMetaData(mydata, AUCmatrix) #这样就把每个细胞的regulon的AUC值转到meta.data里了
##################################################################
## 展示每个细胞的regulon的AUC值。可以用UMAP图、小提琴图、也可以计算细胞类型的均值画热图。

## output文件夹内，Step2_regulonTargetsInfo.tsv，就是高度可信的TF-target关系对
## 请特别注意这个文件，后续分析找到了有价值的regulon，需要回到这个文件找对应的转录因子和靶基因 。
## SCENIC中regulon的名称有两种，一种是TF名称+extended+靶基因数目，另一种是TF名称+靶基因数目。第一种是转录因子与所有靶基因组成的基因调控网络，第二种是转录因子与高可信靶基因（即highConfAnnot=TRUE的基因）组成的基因调控网络。
## 所以，挑选出没有extended的regulon

## regulon在每个细胞中AUC值，以regulon为行细胞为列的矩阵。int/3.4_regulonAUC.Rds
## 计算每个regulon的AUC阈值，细胞中regulonAUC值>阈值，代表此regulon在细胞中处于激活状态，否则代表沉默状态。int/3.5_AUCellThresholds.Rds

df = read.table("./TF_AUC_heatmap.txt", header=T, row.names=1, sep="\t")
pheatmap(df, scale="row", angle_col="45", color=colorRampPalette(c("#0da9ce", "white", "#e74a32"))(500), clustering_method="median", border_color="white")




