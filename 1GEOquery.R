#R包的安装，未安装请去除下面的#号键进行安装
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biobase")
#BiocManager::install("GEOquery")
#BiocManager::install("biomaRt")

#引用包
library(Biobase)
library(GEOquery)

#输入GSE系列号和GPL----
GSE="GSE64216"
GPL="GPL10558"
#输入gene所在列名称,一般为Gene symbol
genename="Symbol"
genename1="ILMN_Gene"

#提取下载GEO数据
gset <- getGEO(GSE, GSEMatrix =T, getGPL = T, AnnotGPL = F) 
if (length(gset) > 1) idx <- grep(GPL, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#查看gset里的信息
str(gset)

#提取表达矩阵
afexp<-data.frame(exprs(gset))
#如果gset里没有gene symbol就需要用其它方法进行注释，但是可以先把探针矩阵整理出来（需要整理探针矩阵可以去除#运行）
#保存探针矩阵
#annmatrix=rbind(ID=colnames(afexp),afexp)
#write.table(annmatrix,file="annmatrix.txt",sep="\t",quote=F,col.names = F)

#简单看看平台的基因名称在哪，顺便提取一下吧
head(gset@featureData@data)
nname=which(colnames(gset@featureData@data)==genename)
afexp$ID=as.character(gset@featureData@data[,nname])

#保存探针信息
ann=cbind(as.character(gset@featureData@data[,nname]),rownames(afexp))
write.table(ann,file="ann.xls",sep="\t",quote=F,col.names = F,row.names = F)

#删除没有gene symbol的探针组
afexp<-afexp[afexp$ID!="",]
afexp=na.omit(afexp)

#对重复基因取平均值并保存好整理的表达矩阵
uniafexp<-aggregate(.~ID,afexp,mean)
uniafexp_adj <- uniafexp[, -c(5, 9)]
uniafexp_adj <- uniafexp_adj[-1,]
write.table(uniafexp_adj,file=paste0(GSE,".txt"),sep="\t",quote=F,col.names = T,row.names = F)

#保存样本临床处理信息
cli=pData(gset)
cliaf=rbind(ID=colnames(cli),cli)
write.table(cliaf,file="clinical.xls",sep="\t",quote=F,col.names = F)


#读取原始数据
## 文件所在目录
library(limma)
preview1 <- read_tsv("GSE64216_non-normalized_data.txt/GSE64216_non-normalized_data.txt")
colnames(preview1)

####1----
GSE64216_non_norm <- read.ilmn(files = "GSE64216_non-normalized_data.txt/GSE64216_non-normalized_data.txt",
                               expr = "AVG_Signal",#这里要根据数据实际的样子来改，有时候又变成了sample全
                               other.columns = "Detection",
                               probeid = "PROBE_ID")

head(GSE64216_non_norm$E)
save(GSE64216_non_norm, file = "GSE64216_non_norm.Rdata")

# 提取基因符号
gene_symbols <- GSE64216_non_norm$genes$SYMBOL

# 提取检测值矩阵
detection_values <- GSE64216_non_norm$E

# 将基因符号与检测值合并为一个数据框
result <- data.frame(Symbol = gene_symbols, Detection = detection_values)

# 查看合并结果
head(result)

#删除没有gene symbol的探针组
afexp<-result[result$Symbol!="",]
afexp=na.omit(afexp)

#对重复基因取平均值并保存好整理的表达矩阵
uniafexp<-aggregate(.~Symbol,afexp,mean)
uniafexp_adj <- uniafexp[, -c(5, 9)]
colnames(uniafexp_adj) <- c("Symbol", "GSM1566477", "GSM1566478", "GSM1566479", "GSM1566481", "GSM1566482", "GSM156683")
write.table(uniafexp_adj,file=paste0(GSE,".txt"),sep="\t",quote=F,col.names = T,row.names = F)



