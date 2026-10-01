setwd("C:/Users/csu206/OneDrive/OSF_single_cell/Code/9.immune.cor")
#引用R包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

#参数设置
immfile="CIBERSORT.txt"         #免疫浸润结果   
exp="logexp.txt"       #表达矩阵
C="Health"                     #正常控制组命名
geness= "GSgene.txt"           #基因列表
highcol="red"            #热图高相关性颜色
midcol= "white"                #热图相关性为0颜色
lowcol= "blue"              #热图负相关性颜色
pfi=1                       #样本的p值过滤，根据cibersort结果，最好为0.05，这里是1即对样本不进行过滤

immune=read.table(immfile, header=T, sep="\t", check.names=F, row.names=1)
#删除正常样本
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
immune=immune[rownames(sample)[-(1:sum(sample[,1]==C))],]

colnames(immune)=gsub("_CIBERSORT","",colnames(immune))
immune=immune[immune[,"P-value"]<pfi,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
expdata=read.table(exp, header=T, sep="\t", check.names=F, row.names=1)
af=read.table(file =geness,sep = "\t",header = F,check.names = F)
expdata=expdata[intersect(af[,1],rownames(expdata)),]
sameSample=intersect(rownames(data),colnames(expdata))
data=data[sameSample,]
expdata=expdata[,sameSample]
expdata=t(expdata)
data1=data
#计算cor
outTab=data.frame()
for(immune in colnames(data)){
	for(gene in colnames(expdata)){
		x=as.numeric(data[,immune])
		y=as.numeric(expdata[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value #多重检验需要改
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}
write.table(outTab,quote = F,sep = "\t","cor.xls",row.names = F)
outTab$cor=as.numeric(outTab$cor)
#去除NA值
outTab=outTab[which(is.na(outTab$cor)=="FALSE"),]

pdf(file="cor.pdf", width=15, height=8)
ggplot(outTab, aes(Gene, Immune))  + 
  theme_bw() + 
  theme(panel.grid.major = element_blank()) + 
  theme(legend.key=element_blank())  +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) +   
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  
        axis.text.y = element_text(size = 10, face = "bold"))  + 
  geom_tile(aes(fill=cor), colour = "white", size = 0.5)+
  labs(fill =paste0("Correlation: ","***  p<0.001","  ", "**  p<0.01","  ", " *  p<0.05","       ")) + 
  scale_fill_gradient2(low = lowcol,mid = midcol, high = highcol)+ 
  geom_text(aes(label=text),col ="grey",size = 5) + theme(legend.position="top")
dev.off()



# 附加验证----
# =========================
# MCP-counter + xCell
# =========================

library(IOBR)
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)

#--------------------------
# 1. 文件
#--------------------------
expfile <- "logexp.txt"
samplefile <- "sample.txt"

C <- "Health"   # 正常组名称

#--------------------------
# 2. 读取表达矩阵
#--------------------------
expdata <- read.table(
  expfile,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

expdata <- as.matrix(expdata)
storage.mode(expdata) <- "numeric"

# 查看
dim(expdata)
expdata[1:5, 1:5]

# 去除重复基因
expdata <- avereps(expdata)

# 确认没有重复基因
sum(duplicated(rownames(expdata)))

#--------------------------
# 3. 读取分组
#--------------------------
sample <- read.table(
  samplefile,
  sep = "\t",
  header = FALSE,
  check.names = FALSE,
  row.names = 1
)

colnames(sample) <- "Group"

table(sample$Group)

# 保留表达矩阵和分组共同样本
sameSample <- intersect(colnames(expdata), rownames(sample))

expdata <- expdata[, sameSample, drop = FALSE]
sample <- sample[sameSample, , drop = FALSE]

# 保证顺序完全一致
expdata <- expdata[, rownames(sample), drop = FALSE]

stopifnot(all(colnames(expdata) == rownames(sample)))

# =========================
# 4. MCP-counter
# =========================

mcp <- deconvo_tme(
  eset = expdata,
  method = "mcpcounter"
)

head(mcp)
dim(mcp)
colnames(mcp)

write.table(
  mcp,
  file = "MCPcounter_result.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#==========================
# MCP-counter group comparison
#==========================

# 看第一列样本ID叫什么
colnames(mcp)

# IOBR一般第一列是ID
mcp2 <- mcp

# 如果第一列不是 Sample，统一改名
colnames(mcp2)[1] <- "Sample"

mcp2$Group <- sample[mcp2$Sample, "Group"]

# 转成长数据
mcp_long <- mcp2 %>%
  pivot_longer(
    cols = -c(Sample, Group),
    names_to = "Cell",
    values_to = "Score"
  )

# Wilcoxon
mcp_stat <- mcp_long %>%
  group_by(Cell) %>%
  summarise(
    Health_median = median(Score[Group == C], na.rm = TRUE),
    OSF_median = median(Score[Group != C], na.rm = TRUE),
    pvalue = wilcox.test(
      Score[Group == C],
      Score[Group != C],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH"),
    Direction = ifelse(
      OSF_median > Health_median,
      "Higher in OSF",
      "Lower in OSF"
    )
  ) %>%
  arrange(pvalue)

write.table(
  mcp_stat,
  "MCPcounter_OSF_vs_Health.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mcp_stat


sig_mcp <- mcp_stat %>%
  filter(pvalue < 0.05) %>%
  pull(Cell)

# sig_xcell <- xcell_stat %>%
#   filter(pvalue < 0.05) %>%
#   pull(Cell)

tiff("MCPcounter_boxplot.tiff", width = 12, height = 7,units = "in",
     res = 600,
     compression = "lzw")

ggplot(
  mcp_long,
  aes(x = Group, y = Score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~Cell, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  ylab("MCP-counter score")


dev.off()



# ============================
# MCP-counter × FRG correlation
# ============================

library(dplyr)

# 1. MCP-counter矩阵：sample × immune cell
mcp_cor <- as.data.frame(mcp2)

rownames(mcp_cor) <- mcp_cor$Sample

# 去掉 Sample 和 Group 两列
mcp_cor <- mcp_cor[, !colnames(mcp_cor) %in% c("Sample", "Group"), drop = FALSE]

# 检查
dim(mcp_cor)
head(rownames(mcp_cor))


# 2. 读取你的64个FRGs
af <- read.table(
  "GSgene.txt",
  sep = "\t",
  header = FALSE,
  check.names = FALSE
)

# 找表达矩阵中存在的FRGs
genes <- intersect(as.character(af[,1]), rownames(expdata))

length(genes)
genes
# expdata原来是 gene × sample
# 转为 sample × gene
gene_exp <- t(
  expdata[genes, , drop = FALSE]
)

gene_exp <- as.data.frame(gene_exp)

dim(gene_exp)
head(rownames(gene_exp))
sameSample2 <- intersect(
  rownames(mcp_cor),
  rownames(gene_exp)
)

length(sameSample2)
sameSample2
mcp_cor <- mcp_cor[sameSample2, , drop = FALSE]
gene_exp <- gene_exp[sameSample2, , drop = FALSE]

# 强制相同顺序
gene_exp <- gene_exp[rownames(mcp_cor), , drop = FALSE]

dim(mcp_cor)
dim(gene_exp)

all(rownames(mcp_cor) == rownames(gene_exp))


outTab_mcp <- data.frame()

for (immune_cell in colnames(mcp_cor)) {
  
  for (gene in colnames(gene_exp)) {
    
    x <- as.numeric(mcp_cor[, immune_cell])
    y <- as.numeric(gene_exp[, gene])
    
    # 删除成对缺失值
    keep <- complete.cases(x, y)
    x2 <- x[keep]
    y2 <- y[keep]
    
    # 避免全是同一个数导致cor.test报错
    if (
      length(x2) >= 3 &&
      length(unique(x2)) > 1 &&
      length(unique(y2)) > 1
    ) {
      
      corT <- cor.test(
        x2,
        y2,
        method = "spearman",
        exact = FALSE
      )
      
      outTab_mcp <- rbind(
        outTab_mcp,
        data.frame(
          Gene = gene,
          Immune = immune_cell,
          cor = as.numeric(corT$estimate),
          pvalue = corT$p.value
        )
      )
    }
  }
}
outTab_mcp$FDR <- p.adjust(
  outTab_mcp$pvalue,
  method = "BH"
)

outTab_mcp$text <- ifelse(
  outTab_mcp$pvalue < 0.001, "***",
  ifelse(
    outTab_mcp$pvalue < 0.01, "**",
    ifelse(
      outTab_mcp$pvalue < 0.05, "*", ""
    )
  )
)

write.table(
  outTab_mcp,
  file = "MCPcounter_FRG_correlation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#补充
library(ggpubr)

p <- ggplot(
  mcp_long,
  aes(x = Group, y = Score, fill = Group)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    shape = 21,
    fill = "white"
  ) +
  
  stat_compare_means(
    comparisons = list(c("Health", "OSF")),
    method = "wilcox.test",
    label = "p.format",
    size = 3.5
  ) +
  
  facet_wrap(
    ~Cell,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "Health" = "#4DBBD5",
      "OSF" = "#E64B35"
    )
  ) +
  
  labs(
    x = NULL,
    y = "MCP-counter score"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      size = 9,
      face = "bold"
    ),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    legend.position = "none"
  )
tiff(
  "MCPcounter_boxplot_color_Pvalue.tiff",
  width = 12,
  height = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

print(p)

dev.off()
save.image(file="MCP.Rdata")
