library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(glmGamPoi)

NKT = SCTransform(NKT, method="glmGamPoi", vars.to.regress="percent.mito", verbose=FALSE)
NKT = RunPCA(NKT, verbose=FALSE)
NKT = RunHarmony(NKT, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(NKT)
NKT <- FindNeighbors(NKT, dims=1:20, reduction="harmony")
NKT <- RunUMAP(NKT, dims=1:20, reduction="harmony")

colors = c("#FFCC00", "#009999", "#CC3366")
colors = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A")
######################################  画图展示批次效应去除的结果
DimPlot(NKT, reduction="umap", group.by="patient_ID", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)
#################################################################################
library(clustree)
obj = FindClusters(mydata, resolution = seq(0.1, 0.3,by=0.05))
clustree(obj)
#############################  NK/T细胞分类的resolution=0.1
mydata <- FindClusters(NKT, resolution=0.1)
UMAPPlot(mydata, pt.size=1.5, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
#################################################################################
# 发现一共有13个亚群
FeaturePlot(mydata, features=c("Ly6d"), cols=c("lightgray", "red"))+NoLegend()
VlnPlot(mydata, features=c("AGER"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())

# 0: CD8+ T cells: TRAC, CD8B, CAMK4, CD27, IL7R
# 1: NK cells: KLRF1, NCR3
# 2: Regulatory T cells: KIR2DL4, GEM, KRT86, TNFRSF18, ENTPD1

cell_label = c("CD8+ T cells", "NK cells", "Regulatory T cells")
#################################################################################
## 给细胞的标签命名
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
UMAPPlot(mydata, pt.size=1.5, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)

set = c("TRAC", "CD8B", "CAMK4", "CD27", "KLRF1", "NCR3", "KIR2DL4", "GEM", "KRT86", "TNFRSF18", "ENTPD1")
DotPlot(mydata, features=set)+coord_flip()+theme_bw()+scale_color_gradient2(low="steelblue1", mid="white", high="hotpink1")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
############################################################################ 画柱状图
bar = mydata@meta.data %>% group_by(Type, cell_type) %>% count()
bar$cell_type = factor(bar$cell_type, levels=cell_label)
bar$Type = factor(bar$Type, levels=c("Normal", "Tumor", "Liver_metastasis"))
ggplot(data=bar, aes(x=Type, y=n, fill=cell_type))+ 
geom_bar(stat="identity", position=position_fill())+
scale_fill_manual(values=colors)+theme_classic()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=30, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")

#####################################################################  配对柱状图
Tissue_label = c("primary tumor", "metastasis")
bar2$Tissue = factor(bar2$Tissue, levels=Tissue_label)
ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Tissue))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("dodgerblue1", "darkorange1"))+theme_classic()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())

############################################################################ MMP12做差异表达分析
MMP12 = subset(mydata, cell_type=="MMP12+ NKT")
Idents(MMP12) = MMP12@meta.data$Type
DEG <- FindAllMarkers(MMP12, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(DEG, "MMP12_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
saveRDS(mydata, "./NKT_sub.rds")






