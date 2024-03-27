library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(glmGamPoi)

Myeloid = SCTransform(Myeloid, method="glmGamPoi", vars.to.regress="percent.mito", verbose=FALSE)
Myeloid = RunPCA(Myeloid, verbose=FALSE)
Myeloid = RunHarmony(Myeloid, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(Myeloid)
Myeloid <- FindNeighbors(Myeloid, dims=1:20, reduction="harmony")
Myeloid <- RunUMAP(Myeloid, dims=1:20, reduction="harmony")

colors = c("#66CCCC", "#FF99CC", "#CCFF66")

######################################  画图展示批次效应去除的结果
DimPlot(Myeloid, reduction="umap", group.by="patient_ID", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)
#################################################################################
library(clustree)
obj = FindClusters(mydata, resolution = seq(0.1, 0.3,by=0.05))
clustree(obj)
#############################  髓系细胞分类的resolution=0.1
mydata <- FindClusters(Myeloid, resolution=0.1)
UMAPPlot(mydata, pt.size=1.5, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
#################################################################################
# 发现一共有13个亚群
FeaturePlot(mydata, features=c("Ly6d"), cols=c("lightgray", "red"))+NoLegend()
VlnPlot(mydata, features=c("AGER"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())

# 0: Monocytes: CD300E, LILRA5, S100A12, CDC42EP2, NLRP3, EHD1, IRAK3, RXRA, HES4
# 1: Macrophages: SELENOP, SLC40A1, GPNMB, PLTP
# 2: Dendritic cells: FCER1A, CD1C

cell_label = c("Monocytes", "Macrophages", "Dendritic cells")
#################################################################################
## 给细胞的标签命名
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
UMAPPlot(mydata, pt.size=1.5, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)

set = c("CD300E", "LILRA5", "S100A12", "NLRP3", "EHD1", "IRAK3", "RXRA", "HES4", "SELENOP", "SLC40A1", "GPNMB", "PLTP", "FCER1A", "CD1C")
DotPlot(mydata, features=set)+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
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
MMP12 = subset(mydata, cell_type=="MMP12+ Myeloid")
Idents(MMP12) = MMP12@meta.data$Type
DEG <- FindAllMarkers(MMP12, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(DEG, "MMP12_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
saveRDS(mydata, "./Myeloid_sub.rds")






