#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
#求求您别盗版了

library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(ggpubr)
library(reticulate)
devtools::install_github("jinworks/CellChat")
#数据标准化----
af=readRDS("osf.rds")
afgs="Ferroptosis.txt"

# QC指标使用小提琴图可视化,ncol为图片排列列数
pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)+scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))
dev.off()

# 指标之间的相关性
plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#组图
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

##标准化,使用LogNormalize方法
af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)

## 鉴定高变基因
# 高变基因：在一些细胞中表达高，另一些细胞中表达低的基因
# 变异指标： mean-variance relationship
# 返回2000个高变基因，用于下游如PCA降维分析。
af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因
top10 <- head(VariableFeatures(af), 10)
top10

# 展示高变基因
plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "01.topgene.pdf",width =7,height = 6)
plot2                   #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

## 归一化
# 归一化处理：每一个基因在所有细胞中的均值变为0，方差标为1，对于降维来说是必需步骤
# 归一化后的值保存在：af[["RNA"]]@scale.data
af <- ScaleData(af)

# 可以选择全部基因归一化
all.genes <- rownames(af)
af <- ScaleData(af, features = all.genes)

##########0.3.part3 降维(绘制原始分布)##########################
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
af <- Seurat::RunPCA(af, features = VariableFeatures(object = af))
af <- Seurat::RunTSNE(af,dims = 1:20)
pdf(file = "02.rawtsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right、left、bottom(意思同英文)
dev.off()
pdf(file = "02.rawpca.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "pca",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
colaa=distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf",width =8,height =5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) #选择展示的主成分，这边是PC2与PC1
dev.off()
#（选做 harmony 去批次与降维）
af <- RunHarmony(af, group.by.vars = "Type")

# 矫正后结果可视化
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "harmony",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af <- Seurat::RunTSNE(af,dims = 1:20,reduction ='harmony')
pdf(file = "03.tsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
pdf(file = "03.tsne.split.pdf",width =12,height = 7.5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) #选择展示的主成分，这边是PC2与PC1
dev.off()

#####作图----
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
# 前两个PC特征基因可视化
VizDimLoadings(af, dims = 1:2, reduction = "pca")
#热图可视化前15个PC
pdf(file = "04.pc_heatmap.pdf",width =7.5,height = 9)
DimHeatmap(af, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()
##确定使用PC个数
# each PC essentially representing a ‘metafeature’
af <- JackStraw(af, num.replicate = 100)
af <- ScoreJackStraw(af, dims = 1:20)
pdf(file = "04.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(af, dims = 1:20)
dev.off()
pdf(file = "04.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30,reduction = "harmony")
dev.off()

#选择PC----
afPC=7
##对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
#af <- FindNeighbors(af, dims = 1:15)
#设置不同的分辨率，观察分群效果，dim为PCA选择的主成分数
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

p2_tree=clustree(af@meta.data, prefix = "RNA_snn_res.")
pdf(file = "04.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()
# 聚类并最优化----
# resolution参数：值越大，细胞分群数越多，根据前面进行选择
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 
#选择分辨率进行降维
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
af <- FindClusters(af, resolution = 0.8)#key point

# 查看聚类数ID
head(Idents(af), 5)

# 查看每个类别多少个细胞
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# 鉴定各个细胞集群的标志基因only.pos：只保留上调差异表达的基因
af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "05.cluster_markers.csv")

## 将细胞在低维空间可视化UMAP/tSNE
af <- RunUMAP(af, dims = 1:afPC, reduction = "harmony")
af <- RunTSNE(af, dims = 1:afPC, reduction = "harmony")

# 可视化UMAP/tSNE
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()



# AddModuleScore计算基因集评分----
afgenes=read.table("ppi.hub.txt",header = F,sep = "\t")[,1]
afgss=as.character(read.table(afgs,header = F,sep = "\t")[,1])
af <- AddModuleScore(
  object = af,
  features =list(afgss[afgss %in% rownames(af)]) ,
  ctrl = 100, #默认值是100
  name = gsub(".txt","",afgs)
)
colnames(af@meta.data)[length(colnames(af@meta.data))] <- gsub(".txt","",afgs) 
colsa = distinctColorPalette(100) #随机颜色

###显著性分析----
#group diff
source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)
comparisons <- list(names(table(af$Type)))
afvp(af, gene_signature = gene_sig, file_name = "06-group.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.format",group="Type",widplot=6,heiplot=5.5,ak=0.9) #
#报错添加点
comparisons <- list(names(table(af$Type)))
afvp(af, gene_signature = intersect(afgenes,rownames(af)), file_name = "06-group.hub_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=15,heiplot=40,ak=0.9) 
#cell diff
comparisons <- list()
comp=combn(names(table(af$seurat_clusters)),2)
names(table(af$seurat_clusters))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af, gene_signature = gene_sig, file_name = "06-cluster.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="seurat_clusters",widplot=20,heiplot=8,ak=0.01,split = "Type")
#cell diff
pdf(file = "06-cluster.hub_VlnPlot.pdf",width =10,height = 6)
VlnPlot(af, features = afgenes,group.by = "seurat_clusters", stack=TRUE,cols = colsa, slot = "data")+ NoLegend()   
dev.off()


#0.5.part5 singleR细胞注释----
refdata=celldex::HumanPrimaryCellAtlasData()    #网络为公用可能加载不好，用手机开热点吧

refdata
head(colnames(refdata))
head(rownames(refdata))

# 查看共有多少种细胞类型
unique(refdata@colData@listData[["label.main"]])
# 使用的数据为标化后的数据
testdata <- GetAssayData(af, layer ="data")
dim(testdata)
testdata[1:30,1:4]
clusters <- af@meta.data$seurat_clusters
table(clusters)
cellpred <- SingleR(test = testdata,  
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

str(cellpred,max.level = 3)
metadata <- cellpred@metadata
head(metadata)

celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = F)
celltype
write.csv(celltype, "07.singleR.celltype_anno_SingleR.csv")
# 打分热图上面的注释结果需要校正
pdf(file = "07-singleR.pdf",width =7.5,height = 6.5)
p = plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")
p
dev.off()
#########sigleR注释后结果可视化
newLabels=cellpred$labels
names(newLabels)=levels(af)
af=RenameIdents(af, newLabels)
# 可视化UMAP/tSNE
pdf(file = "07-scRNA.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "07-scRNA.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af$cellType=Idents(af)

##人工注释----
#实用网址：
#https://www.thermofisher.cn/cn/zh/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work.html
#http://xteam.xbio.top/CellMarker/
#https://www.jianshu.com/p/15dddefc7038
#toppgene.cchmc.org
####################################################################
genes <- list("Cancer stem cells" = c("PROM1","CD34","CD90"),
"Monocyte"=c("CD14"),
"M1 macrophages" = c("CD16", "FCGR3B","FCGR1A"),
"M2 macrophages" = c("MSR1","CD163","MRC1","CSF1R"),
"CD8+ T cells" = c("GZMK", "CD8A","CD8B"),
"CD4+ memory cells" = c("IL7R", "CD27","CCR7"),
"B cells" = c("CD79A", "CD37", "MS4A1"),
"Regulatory T cells" = c("LAG3", "ITGA2","FOXP3","HELIOS","NRP1"),
"NK cells" = c("CD160","NKG7", "GNLY", "CD247", "CCL3", "GZMB"),
"Fibroblasts" = c("FGF7", "MME", "COL3A1"),
"Endothelial cells" = c("PECAM1", "VWF"),
"Neurons" = c("ENO2"),
"Epithelial_cells"=c("cD24","CDH1","CLDN4"),
"myeloid" = c("CD163","AIF1"),
"plasma cell" = c("JCHAIN"),
"Mucosal-associated invariant T cells" = c( "KLRB1","SLC4A10"),
"myeloid Dendritic" = c("CD1C"),
"Plasmacytoid Dendritic Cells" = c("IL3RA"),
"Keratinocytes (KC)" = c("KRT14", "KRT1", "DMKN"),
"Melanocytes (MLNC)" = c("DCT", "TYRP1", "PMEL"),
"Eccrine gland cells (ECG)" = c("PIP", "DCD", "MUCL1"),
"Smooth muscle cells (SMC)" = c("ACTA2", "TAGLN", "MYL9"),
"Nerve cells (Nerve)" = c("MPZ", "PLP1", "S100B"),
"Mast cells (Mast)" = c("CPA3", "TPSAB1", "CTSG")
)
#)
################
af <- FindClusters(af, resolution = 0.8)
save.image(file="tem.Rdata")
#for Epithelial and Keratinocytes
genes1 <- list("Epithelial_cells"=c("CDH1","CLDN4"),
               "Eccrine gland cells (ECG)" = c("PIP", "MUCL1"),
              "Plasmacytoid Dendritic Cells" = c("IL3RA"),
              "Keratinocytes" = c("KRT14", "KRT1", "DMKN"),
              "Melanocytes" = c("DCT", "PMEL")
)
pdf(file = "08.ann_cluster_marker_main.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = genes1,dot.scale = 15,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

BC <- list("Basal Cells" = c("KRT5", "KRT14", "TP63",
                             "CD44",  "SOX2", "EGFR")
)
pdf(file = "08.ann_cluster_marker_BC.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = BC,dot.scale = 15,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

# for others
genes2 <- list("Cancer stem cells" = c("PROM1","CD34"),
              "M2 macrophages" = c("MRC1","CSF1R"),
              "CD8+ T cells" = c("GZMK", "CD8A"),
              "CD4+ memory cells" = c("IL7R", "CD27"),
              "B cells" = c("CD79A", "CD37", "MS4A1"),
              "Regulatory T cells" = c("LAG3"),
              "NK cells" = c("NKG7", "GNLY", "GZMB"),
              "Fibroblasts" = c("FGF7", "COL3A1"),
              "Endothelial cells" = c("PECAM1","VWF"),
              "myeloid Dendritic" = c("CD1C"),
              "Plasmacytoid Dendritic Cells" = c("IL3RA")
)
pdf(file = "08.ann_cluster_marker_others.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = genes2,dot.scale = 15,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

M2 <-list("M2" = c(
  "IDO",          # 分泌
  "IL10",         # 分泌
  "TGFB",         # 分泌
  "CSF1R",        # 表面 (CD115)
  "CD204",        # 表面
  "CD163",        # 表面
  "MRC1",         # 表面 (CD206)
  "CD209",        # 表面 (DC-SIGN)
  "FCER1",        # 表面
  "VSIG4",        # 表面
  "IRF4",         # 胞内/转录因子
  "STAT6"         # 胞内/转录因子
))
pdf(file = "08.ann_cluster_marker_M2.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = M2,dot.scale = 15,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

T_Cell <-list("T Cells" = c("CD3D", "CD3E"),
              "Tc" = c( "CD8A", "CD8B")
              )
pdf(file = "08.ann_cluster_marker_Tcell.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = T_Cell,dot.scale = 15,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()


#人工注释
table(af@active.ident)
ann.ids <- c("Basal_cell",           #cluster0
             "Epithelial_cells",     #cluster1
             "Epithelial_cells",     #cluster2
             "Keratinocytes",        #cluster3
             "Epithelial_cells",   #cluster4
             "Keratinocytes",         #cluster5
             "Keratinocytes",         #cluster6
             "Epithelial_cells",   #cluster7
             "Epithelial_cells",         #cluster8
             "Epithelial_cells",   #cluster9
             "Epithelial_cells",   #cluster10
             "Keratinocytes",         #cluster11
             "Basal_cell",            #cluster12
             "Killer_T_cell",         #cluster13
             "Basal_cell",            #cluster14
             "Basal_cell",            #cluster15
             "Fibroblasts",           #cluster16
             "M2c_macrophages",       #cluster17
             "Melanocytes"            #cluster18
)


afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens
af$cellType=Idents(af)
#########人工注释后结果可视化
# 可视化UMAP/tSNE
pdf(file = "08-ann.scRNA.UMAP.pdf",width =9,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "08-ann.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


# 查看当前SingleR输出有哪些字段
colnames(cellpred)

# 每行是一个cluster的注释诊断信息
annotation_info <- data.frame(
  ClusterID = rownames(cellpred),
  SingleR_label = as.character(cellpred$labels),
  Delta_next = cellpred$delta.next,
  Pruned_label = as.character(cellpred$pruned.labels),
  Flagged_uncertain = is.na(cellpred$pruned.labels)
)

write.csv(
  annotation_info,
  "SingleR_annotation_diagnostics.csv",
  row.names = FALSE
)

# 各cluster对所有参考细胞类型的匹配分数
write.csv(
  as.matrix(cellpred$scores),
  "SingleR_score_matrix.csv"
)


meta <- af[[]]

# 行：供体；列：最终细胞类型
cell_counts <- table(
  Donor = meta$orig.ident,
  CellType = meta$cellType,
  useNA = "ifany"
)

# 每位供体内部计算比例
cell_proportions <- prop.table(cell_counts, margin = 1)

write.csv(
  as.data.frame.matrix(cell_counts),
  "Cell_counts_per_donor.csv"
)

write.csv(
  100 * as.data.frame.matrix(cell_proportions),
  "Cell_percentages_per_donor.csv"
)






# only.pos：只保留上调差异表达的基因----
af.markers <- FindAllMarkers(af, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "08.cell_markers.csv")
# get top 10 genes
top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# plot (需要采取抽样)
pdf(file = "09-cell_marker.hetmap.pdf",width =15,height = 10)
DoHeatmap(af,features = top5af.markers$gene,
          group.colors = colsa) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()
# 抽样plot
af_sub=subset(af, downsample = 50)
library(pheatmap)
colanno=af_sub@meta.data 
colanno$barcode=rownames(colanno)
colanno=colanno%>%arrange(cellType)
rownames(colanno)=colanno$barcode
colanno$barcode=NULL


colanno$cellType=factor(colanno$cellType,levels = unique(colanno$cellType))
colanno1=colanno[,23]
colanno1=as.data.frame(colanno1)
rownames(colanno1)=rownames(colanno)
colnames(colanno1)="celltype"

rowanno=top5af.markers
colnames(top5af.markers)
rowanno=rowanno%>%arrange(cluster)

mat4 <- MinMax(scale(as.matrix(af[["RNA"]]@data)), -2, 2)[rowanno$gene,rownames(colanno1)]
pdf(file = "09-cell_marker.hetmap_50samples.pdf",width =15,height = 10)
pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) 
         
)
dev.off()
#基因分布图
colaa=distinctColorPalette(100)
afgenes=read.table("ppi.hub.txt",header = F,sep = "\t")[,1]
pdf(file = "09-cell_FeaturePlot.pdf",width =12,height = 10)
FeaturePlot(af, features = afgenes, cols = c("grey", "red"),min.cutoff = 0.1, max.cutoff = 1,ncol=4,pt.size = 0.5, slot = "counts")    #min.cutoff与max.cutoff修改截断以更好可视化结果，通过颜色强调基因的分布
dev.off()
pdf(file = "09-cell_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = colaa, slot = "counts")+ NoLegend()   
dev.off()
pdf(file = "09-cell_Group_VlnPlot.pdf",width =20,height = 5)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = c("red3","blue3"), slot = "counts", split.by ="Type" )
dev.off()
colsa = distinctColorPalette(100)
pdf(file = "09-cell_GS_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = gsub(".txt","",afgs),group.by = "cellType",cols=colsa, split.by ="Type")
dev.off()
#####################显著性分析----
#group diff
source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)
#cell diff
comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af, gene_signature = gene_sig, file_name = "09-cell.GS.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=10,heiplot=10,ak=0.9)
#cell diff
comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af, gene_signature = intersect(afgenes,rownames(af)), file_name = "09-cell.hub.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=20,heiplot=100,ak=5)

save.image(file="OSF_mannual_ann.Rdata")

######特定基因的细胞亚群比例-----
#差异基因有bulk里的TF，AKR1C2,SAT1,
#甲基化里的SLC40A1, TF, ACSL6, ALOX15, ATG5, STEAP3, AKR1C2, SLC1A5；
#PPi里的GPX4，HMOX1，ACSL4，FTH1,NCOA4, SLC40A1,TP53

geneselect="TF"
cellselect="Epithelial_cells"
af$geneType=ifelse(as.numeric(as.matrix(af@assays$RNA@scale.data)[geneselect,])>median(sort(as.numeric(as.matrix(af[,which(af$cellType %in% c(cellselect))]@assays$RNA@scale.data)[geneselect,]))),paste0("High ",geneselect," ",cellselect),paste0("Low ",geneselect," ",cellselect))
table(
  af$geneType[af$cellType %in% "Epithelial_cells"]
)
epi_cells <- colnames(af)[af$cellType %in% "Epithelial_cells"]
tf_epi <- as.numeric(af@assays$RNA@scale.data["TF", epi_cells])
cutoff <- median(tf_epi)

c(
  Median = cutoff,
  Below = sum(tf_epi < cutoff),
  Equal = sum(tf_epi == cutoff),
  Above = sum(tf_epi > cutoff)
)
Cellratio <- prop.table(table( af[,which(af$cellType %in% c(cellselect))]$geneType,af[,which(af$cellType %in% c(cellselect))]$Type), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]="Celltype"
colourCount = length(unique(Cellratio$Celltype))
colaa=distinctColorPalette(100)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Type',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # 图例："left" 左, "right" 右,  "bottom" 下, "top" 上
  scale_fill_manual(values=colaa)
ggsave("10-cell_geneselect_ration_TF_EC.pdf",width = 6,height = 3.5)

afvp(af,gene_signature = geneselect, file_name = "06-cell.geneselect._VlnPlot_TF", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=12,heiplot=8,ak=0.9,split = "cellType")

####
afc=af[,which(af$cellType %in% c(cellselect))]
Idents(afc)=afc$geneType
# 可视化UMAP/tSNE
pdf(file = "10-sg.scRNA.UMAP.pdf",width =9,height = 5.5)
DimPlot(afc, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "10-sg.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(afc, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

###循环-----
# 定义 geneselect 和 cellselect 列表
geneselect_list <- c("TF", "AKR1C2", "SAT1", "SLC40A1", "ALOX15", "ATG5", "STEAP3", "SLC1A5", "GPX4", 
                     "HMOX1", "ACSL4", "FTH1", "NCOA4", "TP53")
cellselect_list <- c("Basal_cell", "Epithelial_cells", "Keratinocytes", "Eccrine_gland_cells", 
                     "Killer_T_cell", "Fibroblasts", "M2c_macrophages", "Melanocytes")

# 加载必要的库
library(ggplot2)
library(Seurat)
library(randomcoloR) # 用于生成 distinctColorPalette

# 遍历基因和细胞类型
for (gene in geneselect_list) {
  for (cell in cellselect_list) {
    # 创建以基因为名的文件夹（如果不存在）
    dir.create(gene, showWarnings = FALSE)
    
    # 创建 geneType 分组
    af$geneType <- ifelse(
      as.numeric(as.matrix(af@assays$RNA@scale.data)[gene, ]) > 
        median(sort(as.numeric(as.matrix(af[, which(af$cellType %in% c(cell))]@assays$RNA@scale.data)[gene, ]))),
      paste0("High ", gene, " ", cell),
      paste0("Low ", gene, " ", cell)
    )
    
    # 计算细胞群比例
    Cellratio <- prop.table(
      table(
        af[, which(af$cellType %in% c(cell))]$geneType,
        af[, which(af$cellType %in% c(cell))]$Type
      ), 
      margin = 2
    )
    Cellratio <- as.data.frame(Cellratio)
    colnames(Cellratio)[1] <- "Celltype"
    
    # 生成颜色
    colourCount <- length(unique(Cellratio$Celltype))
    colaa <- distinctColorPalette(100)
    
    # 绘制比例图
    p <- ggplot(Cellratio) + 
      geom_bar(
        aes(x = Var2, y = Freq, fill = Celltype),
        stat = "identity",
        width = 0.7,
        size = 0.5,
        colour = '#222222'
      ) + 
      theme_classic() +
      labs(x = 'Type', y = 'Ratio') +
      coord_flip() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right"
      ) +
      scale_fill_manual(values = colaa)
    
    # 保存比例图
    ggsave(
      filename = paste0(gene, "/cell_geneselect_ratio_", gene, "_", cell, ".pdf"),
      plot = p,
      width = 6,
      height = 3.5
    )
    
    # 绘制小提琴图
    afvp(
      af,
      gene_signature = gene,
      file_name = paste0(gene, "/VlnPlot_", gene, "_", cell),
      test_sign = comparisons,
      pta = 0.1,
      cols = colsa,
      label = "p.signif",
      group = "Type",
      widplot = 12,
      heiplot = 8,
      ak = 0.9,
      split = "cellType"
    )
  }
}


#使用irGSEA进行基因集分析：https://github.com/chuiqin/irGSEA----
library(irGSEA)
library(GSVA)
library(GSEABase)
#irGSEA调用msigdbr中的基因集进行后续分析，注意哈，更改category参数即可指定基因集进行分析，支持的基因集有如下:
#"【category"】	 【description】
#"H"	        hallmark gene sets
#"C1"	        positional gene sets
#"C2"	        curated gene sets
#"C3"	        motif gene sets
#"C4"	        computational gene sets
#"C5"	        GO gene sets
#"C6"	        oncogenic signatures
#"C7"	        immunologic signatures
af.final <- irGSEA.score(object = af[,which(af$cellType %in% c(cellselect))], 
                         assay = "RNA", 
                         slot = "data", 
                         seeds = 123, 
                         ncores = 1,
                         min.cells = 3, 
                         min.feature = 0,
                         custom = F, 
                         msigdb = T, 
                         species = "Homo sapiens", 
                         category = "C7",  
                         geneid = "symbol",
                         method = c("AUCell", "UCell", "singscore", 
                                    "ssgsea"),
                         kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = af.final, 
                               group.by = "geneType",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, heatmap.width =30)
pdf(file = "11-immunologic_gene.heatmap.pdf",width =10,height = 8)
irGSEA.heatmap.plot
dev.off()

#分析循环----
# 定义类别及对应描述
categories <- list(
  H = "hallmark gene sets",
  C1 = "positional gene sets",
  C2 = "curated gene sets",
  C3 = "motif gene sets",
  C4 = "computational gene sets",
  C5 = "GO gene sets",
  C6 = "oncogenic signatures",
  C7 = "immunologic signatures"
)

# 遍历每个类别
for (cat in names(categories)) {
  # 打印当前处理的类别信息
  cat(sprintf("Processing category: %s (%s)\n", cat, categories[[cat]]))
  
  # 计算 irGSEA 分数
  af.final <- irGSEA.score(
    object = af[, which(af$cellType %in% c(cellselect))], 
    assay = "RNA", 
    slot = "data", 
    seeds = 123, 
    ncores = 1, 
    min.cells = 3, 
    min.feature = 0, 
    custom = FALSE, 
    msigdb = TRUE, 
    species = "Homo sapiens", 
    category = cat, 
    geneid = "symbol",
    method = c("AUCell", "UCell", "singscore", "ssgsea"),
    kcdf = 'Gaussian'
  )
  
  # 整合 irGSEA 结果
  result.dge <- irGSEA.integrate(
    object = af.final, 
    group.by = "geneType",
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  
  # 生成热图
  irGSEA.heatmap.plot <- irGSEA.heatmap(
    object = result.dge, 
    method = "RRA", 
    top = 50
  )
  
  # 保存热图到 PDF
  output_file <- sprintf("%s-%s.heatmap.pdf", cat, gsub(" ", "_", categories[[cat]]))
  pdf(file = output_file, width = 10, height = 8)
  print(irGSEA.heatmap.plot)  # 需要显式调用 print() 才能保存 ggplot 对象
  dev.off()
}



#######高低细胞组的差异分析----
af.markers <- FindAllMarkers(afc, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "11.geneGroup_Diff.csv")
###enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(af.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
af.markers=merge(af.markers,ids,by.x='gene',by.y='SYMBOL')
View(af.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(af.markers$ENTREZID, af.markers$cluster) 
## KEGG,12,15
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_KEGG.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_KEGG.pdf",width =8,height = 8)
p +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+ theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     #ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_GO.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_GO.pdf",width =8,height = 8)
p+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60)) + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()

# get top 10 genes
top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

comparisons <- list()
comp=combn(names(table(afc$geneType)),2)
names(table(afc$geneType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af=afc,gene_signature = top5af.markers$gene, file_name = "11-geneGroup_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="geneType",widplot=16,heiplot=49,ak=0.9)

save.image(file = "Before_time_analysis.Rdata")
#######高低细胞组的拟时序分析----
logFCfilter=1           
adjPvalFilter=0.05
af.markers=af.markers[(abs(as.numeric(as.vector(af.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(af.markers$p_val_adj))<adjPvalFilter),]
monocle.matrix=as.matrix(afc@assays$RNA@counts, 'sparseMatrix')
afmetadata=afc@meta.data
monocle.sample=afmetadata[,8,drop=F]
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn$gene_kk_name=monocle.geneAnn$gene_short_name
monocle.markers=af.markers

#将seurat对象转化为monocle输入格式
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- new_cell_data_set(expression_data = as.matrix(monocle.matrix), cell_metadata = pd, gene_metadata = fd    
)

names(pData(cds))[names(pData(cds))=="Size_Factor"]="Cluster"
ssss=as.data.frame(Idents(afc))
pData(cds)[,"geneType"]=paste0(ssss$`Idents(afc)`)
pData(cds)$Type=afc$Type
source("order_cells.R")
source("beam.R")
#devtools::install_version("igraph", version = "2.0.3")
library(igraph)
library(SingleCellExperiment)
library(ggsci)
#开始细胞轨迹分析
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#monocle选择高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree',auto_param_selection = F)
cds <- orderCells(cds)
save(af, cds, disp.genes, file="TF_EC_time.Rdata")
plot(cds$State, cds$Pseudotime)
plot_complex_cell_trajectory(cds, color_by = "State")
#树枝的细胞轨迹图
pdf(file = "12.cds_geneType.pdf",width =5,height = 5)
m1=plot_cell_trajectory(cds,color_by = "geneType", cell_size = 0.5)#+facet_wrap(~cell_type2,ncol=5)  #适当利用分面
m1
dev.off()
#时间的细胞轨迹图
pdf(file = "12.cds_time.pdf",width =5,height = 5)
m2=plot_cell_trajectory(cds,color_by = "Pseudotime", cell_size = 0.5) #+facet_wrap(~cell_type2,ncol=5)  #适当利用分面
m2
dev.off()
#细胞名称的细胞轨迹图
pdf(file = "12.cds_state.pdf",width =6.5,height = 7)
m3=plot_cell_trajectory(cds,color_by = "State", cell_size = 0.5) #+facet_wrap(~celltype,ncol=3) #适当利用分面
m3
dev.off()

pdf(file = "12.cds_GS.pdf",width =8,height = 7)
pData(cds)[,gene_sig] = afc@meta.data[,gene_sig]
plot_cell_trajectory(cds, color_by = gene_sig)  + scale_color_gsea()
dev.off()

pdf(file = "12.cds_geneselect.pdf",width =8,height = 7)
pData(cds)[,geneselect] = afc@assays$RNA@scale.data[geneselect,]
plot_cell_trajectory(cds, color_by = geneselect)  + scale_color_gsea()
dev.off()

#这里用的是disp.genes，https://www.jianshu.com/p/9995cd707002
#选择10分支
BEAM_res <- BEAM(cds[disp.genes,], branch_point = 10, cores = 2) 
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2)  #也可以对所有基因基因进行排序
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "12.BEAM_res.csv", row.names = F)
pdf(file = "12.BEAM_cluster.pdf",width =14,height = 60)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-6)),],
                            branch_point = 10, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)#展示行名
dev.off()

a=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-6)),],
                              branch_point = 10, #绘制的是哪个分支
                              num_clusters = 4, #分成几个cluster，根据需要调整
                              cores = 2,
                              use_gene_short_name = T,
                              show_rownames = T, return_heatmap = T)

#保存基因及其对应的聚类类型----
clusters <- a$annotation_row
clusters <- data.frame(clusters)
write.csv(clusters,file = "12.BEAM_cluster.csv",quote = F,row.names = T,col.names = T)

#细胞通讯的比较分析----
aflist=SplitObject(af,split.by = "Type") 
#创建CellChat 对象：
#用户可以从数据矩阵、Seurat 或SingleCellExperiment对象创建新的 CellChat 对象。如果输入是 Seurat 或SingleCellExperiment对象，则默认情况下将使用对象中的meta data，用户必须提供该数据来定义细胞分组。例如，group.by=Seurat 对象中默认的细胞标识（我们在此处使用前面放置的细胞类型注释“celltype”）。
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#创建cellchat对象
############################################################################################
cellchat.Control = createCellChat(object = aflist$Control, group.by = "cellType")
levels(cellchat.Control@idents) #展示以下现在的细胞分组
groupSize <- as.numeric(table(cellchat.Control@idents)) #细胞亚群各组数量
#设置配体受体交互数据库 
CellChatDB = CellChatDB.human #如果是老鼠的话使用内置“CellChatDB.mouse”数据
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
###选择使用的细胞通讯集合
#使用所有集合进行细胞通讯分析
CellChatDB.use = CellChatDB
cellchat.Control@DB = CellChatDB.use
#预处理表达数据以进行细胞间通讯分析 
cellchat.Control <- subsetData(cellchat.Control)
# devtools::install_github('immunogenomics/presto')
cellchat.Control <- identifyOverExpressedGenes(cellchat.Control)
cellchat.Control <- identifyOverExpressedInteractions(cellchat.Control)
#细胞通信网络的推断 
cellchat.Control <- computeCommunProb(object=cellchat.Control,raw.use = TRUE)
#在信号通路级别推断细胞-细胞通信 
cellchat.Control <- computeCommunProbPathway(cellchat.Control)
#计算整合的细胞通信网络 
cellchat.Control <- aggregateNet(cellchat.Control)
############################################################################################
cellchat.OSF = createCellChat(object = aflist$OSF, group.by = "cellType")
levels(cellchat.OSF@idents) #展示以下现在的细胞分组
groupSize <- as.numeric(table(cellchat.OSF@idents)) #细胞亚群各组数量
#设置配体受体交互数据库 
CellChatDB = CellChatDB.human #如果是老鼠的话使用内置“CellChatDB.mouse”数据
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
###选择使用的细胞通讯集合
#使用所有集合进行细胞通讯分析
CellChatDB.use = CellChatDB
cellchat.OSF@DB = CellChatDB.use
#预处理表达数据以进行细胞间通讯分析 
cellchat.OSF <- subsetData(cellchat.OSF)
cellchat.OSF <- identifyOverExpressedGenes(cellchat.OSF)
cellchat.OSF <- identifyOverExpressedInteractions(cellchat.OSF)
#细胞通信网络的推断 
cellchat.OSF <- computeCommunProb(object=cellchat.OSF,raw.use = TRUE)
#在信号通路级别推断细胞-细胞通信 
cellchat.OSF <- computeCommunProbPathway(cellchat.OSF)
#计算整合的细胞通信网络 
cellchat.OSF <- aggregateNet(cellchat.OSF)
#############合并两个分组的对象
object.list <- list(Control = cellchat.Control, OSF = cellchat.OSF)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf(file = "13.cellchat_compareInteractions.pdf",width =9,height = 4.5)
gg1 + gg2

dev.off()
# ===== 在这里插入你上面整段细胞数统计代码 =====
object.list <- list(
  Control = cellchat.Control,
  OSF = cellchat.OSF
)

# 各组每种细胞类型的实际细胞数
cell_counts <- do.call(
  rbind,
  lapply(names(object.list), function(group_name) {
    
    tab <- table(object.list[[group_name]]@idents)
    
    data.frame(
      Group = group_name,
      CellType = names(tab),
      N_cells = as.integer(tab),
      row.names = NULL
    )
  })
)

# 汇总总细胞数
group_summary <- aggregate(
  N_cells ~ Group,
  data = cell_counts,
  FUN = sum
)

# 已确认：每组各1位供体
group_summary$N_donors <- 1L
group_summary$Downsampling <- "Not performed"

print(group_summary)
print(cell_counts)

write.csv(
  group_summary,
  "CellChat_group_summary.csv",
  row.names = FALSE
)

write.csv(
  cell_counts,
  "CellChat_celltype_counts.csv",
  row.names = FALSE
)

# 保存实际运行版本，Methods按这里填写
packageVersion("CellChat")
########绘制细胞间通讯强度，包括number、weight-----
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(cellchat.Control@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in Control")
netVisual_circle(cellchat.Control@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights in Control")
netVisual_circle(cellchat.OSF@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in OSF")
netVisual_circle(cellchat.OSF@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights in OSF")
###############################
######分组的细胞通讯差异(OSF vs control)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(2, 1), title = "OSF VS Control")
################################绘制细胞通讯强度的热图
a1=netVisual_heatmap(cellchat.Control, title.name = "Control")
a2=netVisual_heatmap(cellchat.OSF, title.name = "OSF")
a1+a2



#细胞通讯的功能差异
par(mfrow = c(1,1), xpd=TRUE)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
reticulate::install_miniconda()
reticulate::py_install("umap-learn")
install.packages("uwot")
library(uwot)
cellchat <- netEmbedding(cellchat,umap.method = "uwot", type = "functional")
rankSimilarity(cellchat, type = "functional", comparison2=c(2, 1), title = "OSF VS Control")
###############################细胞通讯的差异
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2,1))
##############################细胞通讯的受配体对差异，指定source
levels(cellchat@idents$joint) #查看索引，单核细胞第三个
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5),  comparison = c(2, 1), angle.x = 45)

save.image(file="Final_results.Rdata")
