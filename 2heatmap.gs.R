#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
setwd("C:/Users/csu206/OneDrive/OSF_single cell/Code/2.1.diff_gs")
#引用包
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
#lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(randomcoloR)
library(ComplexHeatmap)
library(GSEABase)
library(GSVA)
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
#           "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
#           "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr', 'ComplexHeatmap')
#  for(i in 1:length(depens)){
#   depen<-depens[i]
#    if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
#  }
#devtools::install_github("IOBR/IOBR")
library(IOBR)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggExtra)

#if (!require('R.utils')) install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )

#参数设置
GSE="GSE20170"      #表达矩阵文件名称，不用后缀
gsname="Ferroptosis Signal" #基因集文件名
C="Healthy"                 #正常控制组名称
P="OSF"                      #疾病实验组的名称
is.rna="T"               #是否为转录组测序数据,TRUE
Ccol= "#0073C2FF"        #热图注释条正常组颜色
Pcol="#EFC000FF"        #热图注释条疾病组颜色
lowcol="#0073C2FF"    #热图格子低表达量颜色
midcol="white"     #热图格子中表达量颜色
highcol="#EFC000FF"      #热图格子高表达量颜色
geneNum=30               #展示p值最显著的前30个
immune.pvaluefi=0.05    #免疫浸润样本的p值过滤,不过滤改成1
cor.method="pearson"      #或spearman
xcol="red3"   #散点图x轴分布颜色
ycol="green2"     #散点图y轴分布颜色 

#rt=read.table(paste0(GSE,".txt"),sep = "\t",header=T,check.names=F)
GSE20170 <- read.delim("GSE20170.txt")
rt=as.matrix(GSE20170)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt) #去重

###1.数据准备----
#分组
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
colnames(rt) <- trimws(colnames(rt))    # 去除列名中的前后空格
rownames(sample) <- trimws(rownames(sample))  # 去除行名中的前后空格
rt=rt[,rownames(sample)]
afcon=sum(sample[,1]==C)
mat=normalizeBetweenArrays(rt)
#判断原始数据是否去了log
max(rt)
if(max(rt, na.rm = TRUE)>50) 
  rt=log2(rt+1)     #rt最大值大于50则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

##2.差异分析----
#如果需要使用未经矫正的数据则将下方数据去除
data=rt1
#data=rt
#保存一下log的表达矩阵
outexp=rbind(ID=colnames(rt1),rt1)
write.table(outexp,"logexp.txt",quote = F,sep = "\t",col.names = F)
conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

if (is.rna %in% c("T","TURE","t","ture") ) {
  #limma-voom  
  Type=c(rep("con",conNum),rep("treat",treatNum))
  design <- model.matrix(~0+factor(Type))
  colnames(design) <- c("con","treat")
  v <- voom(rt,design,normalize="quantile")
  fit <- lmFit(v,design)
  cont.matrix<-makeContrasts(treat-con,levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
}else{
  #limma-trend
  Type=c(rep("con",conNum),rep("treat",treatNum))
  design <- model.matrix(~0+factor(Type))
  colnames(design) <- c("con","treat")
  fit <- lmFit(rt,design)
  cont.matrix<-makeContrasts(treat-con,levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
}
Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
aaaa=Diff
save.image(file= "OSF_bulk.Rdata")


#转录组数据中测到基因集中基因的数目----
afGene=read.table(paste0(gsname,".txt"),sep = "\t",header = F)[,1]
afGene=intersect(rownames(Diff),afGene)
print(length(afGene))
#保存结果
Diff=Diff[afGene,]
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
Diff$abs.logFC=abs(Diff$logFC)
write.csv(Diff,"result.csv",quote = F)
#提取
if (length(rownames(Diff)) > geneNum) {
  Diff=Diff[1:geneNum,]
}
afExp=rt1[rownames(Diff),]
#分组标签
Type=c(rep(C,conNum),rep(P,treatNum))
names(Type)=colnames(rt)
top=as.data.frame(Type)
left=as.data.frame(rownames(afExp))
colnames(left)="geneName"
rownames(left)=left[,1]
right=as.data.frame(Diff$logFC)
colnames(right)="logFC"
right$adj.P=as.vector(Diff$adj.P.Val) #
rownames(right)=rownames(Diff)


############# GSEA分析----
GSinput=gsname
geneset=read.table(paste0(GSinput,".txt"),sep="\t",header=F,check.names=F)
geneset$name=rep(GSinput,nrow(geneset))
geneset=geneset[,c(2,1)]
colnames(geneset)=c("term","gene")
rownames(geneset)=geneset$gene
Ensembl_ID <- bitr(geneset$gene, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
kegmt=data.frame(Ensembl_ID ,geneset[match(Ensembl_ID$SYMBOL,geneset$gene),])
kegmt=kegmt[,c(3,2)]
colnames(kegmt)=c("term","gene")
deg=aaaa
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
anyDuplicated(DEG$symbol) 
DEG$symbol <- gsub("\\s+", "", DEG$symbol)  # 去除基因名中的空格
DEG$symbol <- toupper(DEG$symbol)  # 统一基因名为大写字母
DEG <- DEG[!duplicated(DEG$symbol), ]
data_all_sort <- DEG %>% 
  arrange(desc(logFC))
anyDuplicated(data_all_sort$symbol) 
geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
#geneList = geneList[!is.na(geneList)]
names(geneList) <- data_all_sort$symbol #给上面提取的foldchange加上对应上ENTREZID
table(names(geneList))
head(geneList)
anyDuplicated(geneList) 


#开始GSEA富集分析----
kk2<-GSEA(geneList,TERM2GENE = geneset,pvalueCutoff = 0.5,pAdjustMethod = "fdr") #GSEA分析
#保存GSEA结果
GSEAOUT=as.data.frame(kk2@result)

write.table(GSEAOUT,file="GSEAOUT.xls",sep="\t",quote=F,col.names=T,row.names = F)

#保存
library(enrichplot)
pdf(file="GSEA.pdf",height=6,width=6)
gseaplot2(kk2,geneSetID = rownames(kk2@result),
          title = gsname,  #设置title
          color="red", #线条颜色
          base_size = 12, #基础字体的大小
          subplots = 1:3, #展示上3部分
          pvalue_table = T)# 显示p值
dev.off()
dev.new()

####热图汇总----
afgmt=c(GSinput,"NA",geneset$gene)
afgmt=as.matrix(afgmt)
write.table(t(afgmt),file="afGMT.gmt",sep="\t",quote=F,col.names=F,row.names = F)
geneSet=getGmt("afGMT.gmt",geneIdType=SymbolIdentifier())
mat=rt1[rowMeans(rt1)>0,]

#BiocManager::install("GSVAdata")
library(GSVAdata)
ssgseaParams <- ssgseaParam(mat, geneSet)
ssgseaScore=gsva(ssgseaParams)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut, file="ssgseaOut.xls", sep="\t", quote=F, col.names=F)
top$ssGSEA=t(normalize(ssgseaScore))
colnames(top)[2]="ssGSEA score"
#可视化
kk=Heatmap(scale(afExp),
        name = "scale(exp)",
        row_names_side ="left",
        cluster_columns = F, #关闭列聚类
        cluster_rows = F, #关闭行聚类
        row_names_gp = gpar(fontsize = 11,    #列名字体大小
                               fontface = "bold", #列名字体
                               fill = "white", #列名字体背景色
                               col = "black", # 列名字体颜色
                               border = "grey"# 列名字体边框色
        ),
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), #颜色范围及类型，超过此范围的值的颜色都将使用阈值所对应的颜色,也支持颜色向量：c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837')
        border_gp = gpar(col = "grey", lwd = 3), #设置热图的边框与粗度
        rect_gp = gpar(col = "white", lwd = 1), #设置热图中格子的边框与粗度
        column_split = top$Type,#生成24个模拟分类数据对热图列进行分割
        column_gap = unit(2, "mm"), #设置行列分割宽度
        column_title_gp = gpar(fontsize = 20, # 行标题大小
                            fontface = "bold" , # 行标题字体
                            fill = c("blue3","red3"),  # 行标题背景色,根据分割可分别定义颜色
                            col = "white", # 行标题颜色
                            border = "grey"# 行标题边框色
        ),
        top_annotation = HeatmapAnnotation('ssGSEA score' = top[,2],     #添加顶部注释，原理同上
                                           'Expression of all genes' = anno_boxplot(rt1, height = unit(1, "cm"), gp = gpar(fill = 1:24)),annotation_name_gp =gpar(fontsize = 10,fontface = "bold") ), #添加箱式图、修改其高度为1并为每一个样本的箱式图添加颜色,修改图例字体大小
        right_annotation = rowAnnotation(                        #添加多种类型的行注释，支持分类与数值
                                         
                                         '-log(adj.P)' = anno_simple(-log(right$adj.P),  #添加顶部注释，使用ann_simple函数组合，将mat第一行定义为颜色深浅
                                                              pch = ifelse(right$adj.P>0.05,"ns",ifelse(right$adj.P>0.01,"*",ifelse(right$adj.P>0.001,"**","***"))),  #设置格子内的字体，因为有24个故模拟填充24个数值
                                                              pt_size = unit(2, "mm"),
                                                              pt_gp = gpar(fontsize = 10,fontface = "bold",col = as.character(ifelse(right$adj.P>0.05,"grey","red")))
                                                              ),
                                         logFC = anno_barplot(right$logFC,width=unit(3, "cm"), gp = gpar(fill = colorRampPalette(c("blue", "red"))(length(right$logFC))))),
        )
#展示一下看看
kk
#保存结果
pdf(file="complexheatmap.pdf",  width=8, height=8)
print(kk)
dev.off()

#################cibersort免疫浸润
mat=mat[rowMeans(mat)>0,]
afcibersort<-deconvo_tme(mat, method = "cibersort",perm = 1000)
write.table(afcibersort,"CIBERSORT.txt",sep="\t",quote = F,row.names = F)

###细胞比例
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<immune.pvaluefi,]
data_im=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data_im)=gsub("_CIBERSORT"," ",colnames(data_im))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data_im)=gsub("_"," ",colnames(data_im))
cluster=read.table("sample.txt", header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data_im), row.names(cluster))
data_im=cbind(data_im[sameSample,,drop=F], cluster[sameSample,,drop=F])
data_im=data_im[order(data_im$Type),]
gaps=c(1, as.vector(cumsum(table(data_im$Type))))
xlabels=levels(factor(data_im$Type))
data_im$Type=factor(data_im$Type, levels=c(cluster[1,1],cluster[nrow(cluster),1]))
data_im$GSM=rownames(data_im)
data_im=melt(data_im,id.vars=c("Type","GSM"))
colnames(data_im)=c("Type","GSM","Celltype", "Freq")
Cellratio=data_im
colourCount = length(unique(Cellratio$Celltype))
#定义细胞浸润的颜色
colaa= c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
           "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
           "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
           "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
           "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
           "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f",unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                                                                                               pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                                                                                               pal_locuszoom("default")(7),pal_igv("default")(51),
                                                                                               pal_uchicago("default")(9),pal_startrek("uniform")(7),
                                                                                               pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                                                                                               pal_simpsons("springfield")(16),pal_gsea("default")(12))))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell cycle phase',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # 图例："left" 左, "right" 右,  "bottom" 下, "top" 上
  scale_fill_manual(values=colaa[sample(1:length(colaa), size = 25)])+
  theme_bw()+
  xlab(NULL)+
  theme(axis.text.x  = element_blank())+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))+
  facet_grid(. ~ Type,scales="free")+
  theme(strip.text.x = element_text(size = 30,colour = "black"))+       #分面字体颜色
  theme(strip.background.x = element_rect(fill = c("white"), colour = "black")) #分面颜色
ggsave("immune.ration.pdf",width = 15,height = 8)        #输出图片

###细胞比列差异箱式----
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<immune.pvaluefi,]
data_im=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data_im)=gsub("_CIBERSORT"," ",colnames(data_im))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data_im)=gsub("_"," ",colnames(data_im))
cluster=read.table("sample.txt", header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data_im), row.names(cluster))
data_im=cbind(data_im[sameSample,,drop=F], cluster[sameSample,,drop=F])
data_im=data_im[order(data_im$Type),]
gaps=c(1, as.vector(cumsum(table(data_im$Type))))
xlabels=levels(factor(data_im$Type))
data_im=melt(data_im,id.vars=c("Type"))
colnames(data_im)=c("Type", "Immune", "Expression")
group=levels(factor(data_im$Type))
data_im$Type=factor(data_im$Type, c(cluster[1,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data_im, x="Immune", y="Expression",  fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme_bw()+
  rotate_x_text(50)

pdf(file="immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()

#####免疫相关性----
dir.create("immune.cor")
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<immune.pvaluefi,]
data_im=t(as.matrix(immune[,1:(ncol(immune)-3)]))
GS.ssgsea=as.numeric(top[,2])[which(as.character(top[,1])!=C)]
data_im=data_im[,rownames(top)[which(as.character(top[,1])!=C)]]

outTab=data.frame()
    for(j in row.names(data_im)){
      if(sd(data_im[j,])>0.001){
        x=GS.ssgsea
        y=as.numeric(data_im[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>0) ){
          outTab=rbind(outTab,cbind(genekkk=j,clusterkkk=gsname,cor,pvalue,Regulation="postive"))
          df1=as.data.frame(cbind(x,y))
          df1$gsname=df1[,1]
          p1=ggplot(df1, mapping = aes(x, y)) + 
            xlab(gsname) + ylab(j)+
            geom_point(shape=19,color="cyan",size=3) + geom_smooth(method="lm",formula = y ~ x,color='red',se = T) + theme_bw()+
            stat_cor(method = cor.method, aes(x =x, y =y))+ scale_colour_gradient(low = lowcol, high = highcol)
          p2=ggMarginal(p1, type="density", xparams=list(fill = xcol), yparams=list(fill = ycol))
          
          pdf(file=paste0("immune.cor/",j,".pdf"), width=4.5, height=4)
          print(p2)
          dev.off()
        }
        if((cor< 0) ){
          outTab=rbind(outTab,cbind(genekkk=j,clusterkkk=gsname,cor,pvalue,Regulation="negative"))
          df1=as.data.frame(cbind(x,y))
          df1$gsname=df1[,1]
          p1=ggplot(df1, mapping = aes(x, y)) + 
            xlab(gsname) + ylab(j)+
            geom_point(shape=19,color="cyan",size=3) + geom_smooth(method="lm",formula = y ~ x,color='red',se = T) + theme_bw()+
            stat_cor(method = cor.method, aes(x =x, y =y))+ scale_colour_gradient(low = lowcol, high = highcol)
          p2=ggMarginal(p1, type="density", xparams=list(fill =  xcol), yparams=list(fill = ycol))
          pdf(file=paste0("immune.cor/",j,".pdf"), width=4.5, height=4)
          print(p2)
          dev.off()
        }
      }
    }

write.table(file=paste0("immune.cor/","cortable.xls"),outTab,sep="\t",quote=F,row.names=F)
save.image(file="OSFvsFerroptosis.Rdata")
