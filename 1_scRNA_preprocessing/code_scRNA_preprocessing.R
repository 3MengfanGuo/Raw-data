library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
###############################################################################
counts <- Read10X(data.dir="./0_rawdata/", gene.column=1)
scRNA = CreateSeuratObject(counts, min.cells=3, min.features=200)
scRNA = subset(scRNA, Type %in% c("Ca", "Li", "Ov"))
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<6000&percent.mito<15) # 线粒体比例<15%, 基因个数在200~6000之间
p = VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
ggsave("nFeature_RNA_violin.pdf", p, width=6, height=6)
colors = c4a("carto.pastel", 8)
#########################################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
scRNA = SCTransform(scRNA, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### 划分所有细胞的亚群时，resolution=0.1
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 用小提琴图来检测marker基因在细胞亚群之间的分布
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0: Regulatory T cells: ICOS, CTLA4, CD28, TNFRSF4, CCR7
# 1: Cytotoxic NK/T cells 1: CD8B, GZMK, 
# 2: Plasma B cells: JCHAIN, IGHA1, MZB1
# 3: Myeloid cells: LYZ, CST3, C1QC, S100A8
# 4: Cytotoxic NK/T cells 2: GNLY, TRDC, XCL1, FGFBP2, KLRD1
# 5: B cells: MS4A1, BANK1, VPREB3
# 6: Myofibroblasts: TAGLN, COL1A1, BGN, DCN
# 7: Endothelial cells: VWF, PLVAP, CLDN5
# 8: Mast cells: TPSAB1, CPA3
cell_label = c(
"Regulatory T cells", "Cytotoxic NK/T cells", "Plasma B cells", "Myeloid cells", "Cytotoxic NK/T cells",
"B cells", "Myofibroblasts", "Endothelial cells", "Mast cells")
#################################################################################
## 给细胞的标签命名
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
saveRDS(mydata, "./mydata_cluster.rds")
colors = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977", "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A", "#39737C", "#86B4A9", "#82853B", "#CCC94D")
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
genes = c("IL7R", "ICOS", "CTLA4", "CD28", "NKG7", "GNLY", "CD8B", "GZMK", "TRDC", "IGHA1", "MZB1", "LYZ", "CST3", "C1QC", "S100A8", "MS4A1", "BANK1", "TAGLN", "COL1A1", "VWF", "PLVAP", "TPSAB1", "CPA3")
p = DotPlot(mydata, features=genes, cols=c("snow", "chartreuse4"))+coord_flip()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

# marker基因的UMAP映射图
set = c("CTLA4", "GNLY", "MZB1", "CST3", "MS4A1", "COL1A1", "VWF", "TPSAB1")
FeaturePlot(mydata, features=set, cols=c("snow", "red"), ncol=4)

set = c("IL7R", "GNLY", "MZB1", "LYZ", "MS4A1", "COL1A1", "VWF", "TPSAB1")
VlnPlot(mydata, features=set, pt.size=0, cols=colors, ncol=4)+NoLegend()+theme(axis.title.x=element_blank())


#####################################################################  配对柱状图
Type_label = c("Control", "24h Treated", "72h Treated")
bar$Type = factor(bar$Type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("skyblue", "orange"))+theme_classic()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())


########################## 画配对秩和检验的箱线图
library(ggplot2)
library(ggpubr)

df = read.table("./celltype_number_percent.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -percent, sum), y=percent, fill=Type))+
scale_fill_manual(values=c("orange1", "greenyellow"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="t.test")+
theme(axis.text.x=element_text(angle=15, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.x=element_blank())


	



