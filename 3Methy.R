# (0) Preliminaries
################################################################################
# 引用包
# Missing required package GOstats or its dependency
# GhostScript
library(RnBeads)
library(wateRmelon)
library(Gviz)
library(GOstats)
#BiocManager::install("RnBeads.hg38")
#BiocManager::install("Gviz")
Sys.setenv(GS_CMD="\"C:/Program Files/gs/gs10.01.1/bin/gswin64c.exe\"") #文件路径
Sys.setenv(R_ZIPCMD = "\"C:/Program Files/7-Zip/7z.exe\"")
######################### 定义文件路径
# setwd(".")
dataDir <- file.path(getwd(), "data")
resultDir <- file.path(getwd(), "results")
# dataset and file locations
datasetDir <- file.path(dataDir)
idatDir <- file.path(datasetDir, "idat")
sampleSheet <- file.path(datasetDir, "sample_annotation.csv")
reportDir <- file.path(resultDir,"report")
################################################################################
# (1) Set analysis options
################################################################################
rnb.options(
	filtering.sex.chromosomes.removal = TRUE,  #移除性染色体
	identifiers.column                = "Sample_ID" #样本id
)
# optionally disable some parts of the analysis to reduce runtime,??rnb.options
rnb.options(
	exploratory.correlation.qc        = FALSE,
	exploratory.intersample           = FALSE,
	# exploratory.region.profiles       = c("genes"),
	exploratory.region.profiles       = character(0),
	exploratory.clustering            = "top",
	exploratory.clustering.top.sites  = 100,
	# region.types                      = c("promoters", "genes", "tiling"),
	region.types                      = c("promoters", "cpgislands", "genes", "tiling"),   
	differential.report.sites         = TRUE,
	differential.comparison.columns   = c("Type"),
	differential.enrichment.go = T,
	exploratory.gene.symbols = "TF"   #基于Gviz包绘制基因组结构
)
################################################################################
# (2) Run the analysis
################################################################################
#rnb.options(assembly = "hg38")
#BiocManager::install("RnBeads.hg19")
rnb.options(assembly = "hg19")
rnb.run.analysis(
	dir.reports=reportDir,
	sample.sheet=sampleSheet,
	data.dir=idatDir,
	data.type="infinium.idat.dir"
)

dataSource <- c(idatDir, sampleSheet)
#仅执行分析但不输出
rnbs <- rnb.execute.import(data.source = dataSource,
                           data.type="infinium.idat.dir")
#直接从geo导入数据
#rnbs.geo <- rnb.execute.import(data.source="GSE156669", data.type="infinium.GEO")

nsites(rnbs)
samples(rnbs)

#提取并检查CpG水平的甲基化值。
mm <- meth(rnbs)
dim(mm)
boxplot(mm, col="steelblue", las=2)
hist(mm[,5], col="steelblue", breaks=50)
#计算预定义基因组区域的平均甲基化水平。
summarized.regions(rnbs)
mm.p <- meth(rnbs, type="promoters")
dim(mm.p)
boxplot(mm.p, col="steelblue", las=2)
#作为β值的替代，可以获得M值
summary(mval(rnbs))

save.image(file="OSF_methy.Rdata")
################################################################################
# Link to finished analysis
################################################################################
# see the results at:
# http://rnbeads.mpi-inf.mpg.de/reports/tutorial/epigenomics2016/results/report_Ziller2011_vanilla/index.html
