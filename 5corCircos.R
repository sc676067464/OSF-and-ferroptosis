#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
setwd("C:/Users/csu206/OneDrive/OSF_single_cell/Code/8.cor.marker")
#引用包
options(stringsAsFactors=F)
library(corrplot)
library(circlize)
library(limma)
library(PerformanceAnalytics)
#输入文件
expFile="logexp.txt"     #表达矩阵
hub="GSgene.txt"           #核心基因
C="Health"                     #正常控制组命名
higcol="#008837"                   #高正相关性颜色
midcol="#f7f7f7"                #相关性等于0颜色
lowcol="#7b3294"                #高负相关性颜色
showtype="upper"      #一共有三种"full", "lower", "upper",分别为展示全图、下半张图，上半张图
tlcex=0.45          #对称轴字体大小
numbercex=0.5       #相关性系数字体大小
method="pearson"   #相关性检验方法有三种，"pearson", "kendall", "spearman" 

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

#删除正常样本
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
rt=rt[,rownames(sample)[-(1:sum(sample[,1]==C))]]

#提取目的基因
VENN=intersect(rownames(rt),read.table(hub,sep = "\t",header = F,check.names = F)[,1])
rt1=rt[c(VENN),]

rt <- t(rt1) 
# rownames(rt)=rownames(rt1)
# rt=t(rt)  #转置
M=cor(rt,method = method)     #相关型矩阵
res <- cor.mtest(rt)

# ==========================================================
# 审稿补充：表达值、样本数、相关性、BH校正、重复值核查
# 依赖对象：
# rt   = 样本 × 基因的分析矩阵
# data = 原代码中avereps()之前的数值表达矩阵
# ==========================================================

outdir <- "Correlation_verification"
dir.create(outdir, showWarnings = FALSE)

options(digits = 17)

save_csv <- function(x, filename) {
  write.csv(
    x,
    file = file.path(outdir, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

# 1. 导出实际参与分析的样本级表达数据
expression_used <- data.frame(
  SampleID = rownames(rt),
  Group = as.character(sample[rownames(rt), 1]),
  rt,
  check.names = FALSE
)

save_csv(expression_used, "01_expression_used.csv")

# 2. 所有不重复基因对：Pearson r、P值和95%CI
pair_index <- combn(seq_len(ncol(rt)), 2)

all_results <- do.call(rbind, lapply(
  seq_len(ncol(pair_index)),
  function(k) {
    i <- pair_index[1, k]
    j <- pair_index[2, k]
    
    x <- rt[, i]
    y <- rt[, j]
    
    test <- cor.test(
      x, y,
      method = "pearson",
      alternative = "two.sided",
      conf.level = 0.95
    )
    
    ci <- if (is.null(test$conf.int)) {
      c(NA_real_, NA_real_)
    } else {
      as.numeric(test$conf.int)
    }
    
    data.frame(
      Gene1 = colnames(rt)[i],
      Gene2 = colnames(rt)[j],
      N = length(x),
      Pearson_r = unname(test$estimate),
      CI_low = ci[1],
      CI_high = ci[2],
      P_value = test$p.value
    )
  }
))

all_results$P_adj_BH <- p.adjust(
  all_results$P_value,
  method = "BH"
)

# 保留高精度字符列，方便区分1与四舍五入后的1.00
all_results$r_high_precision <- sprintf(
  "%.17g", all_results$Pearson_r
)

save_csv(all_results, "02_all_pair_correlations_BH.csv")

# 3. 提取审稿人要求核查的两个基因对
target_pairs <- list(
  c("AKR1C2", "SAT1"),
  c("TXNRD1", "PCBP1")
)

get_pair_result <- function(g1, g2) {
  all_results[
    (all_results$Gene1 == g1 & all_results$Gene2 == g2) |
      (all_results$Gene1 == g2 & all_results$Gene2 == g1),
    , drop = FALSE
  ]
}

needed_genes <- unique(unlist(target_pairs))

if (!all(needed_genes %in% colnames(rt))) {
  stop(paste(
    "分析矩阵中缺少：",
    paste(setdiff(needed_genes, colnames(rt)), collapse = ", "),
    "。请核对GSgene.txt及基因名，不能直接跳过。"
  ))
}

target_results <- do.call(rbind, lapply(
  target_pairs,
  function(g) get_pair_result(g[1], g[2])
))

save_csv(target_results, "03_target_pair_results.csv")
print(target_results, digits = 12)

# 4. 导出这四个基因在avereps之前的输入行
# 注意：这是logexp.txt中的输入值，不等于平台原始数据
source_rows <- which(rownames(data) %in% needed_genes)

source_values <- data.frame(
  Input_data_row = source_rows,
  Gene = rownames(data)[source_rows],
  data[source_rows, rownames(rt), drop = FALSE],
  check.names = FALSE
)

save_csv(source_values, "04_target_rows_before_averaging.csv")

# 每个基因在输入文件中对应多少行
gene_counts <- table(rownames(data))

gene_audit <- data.frame(
  Gene = colnames(rt),
  Input_row_count = as.integer(gene_counts[colnames(rt)])
)

save_csv(gene_audit, "05_gene_row_counts.csv")

# 5. 检查所选样本的全表达谱是否完全重复
# 用所有输入基因检查，避免只凭4个基因判断样本重复
full_profiles <- t(data[, rownames(rt), drop = FALSE])

duplicate_profile <- duplicated(full_profiles, MARGIN = 1) |
  duplicated(full_profiles, MARGIN = 1, fromLast = TRUE)

sample_audit <- data.frame(
  SampleID = rownames(rt),
  Group = as.character(sample[rownames(rt), 1]),
  Identical_full_profile = duplicate_profile
)

save_csv(sample_audit, "06_sample_profile_audit.csv")

# 6. 检查两基因表达值是否相同、是否有重复坐标、
#    是否接近完全线性关系
pair_audit <- do.call(rbind, lapply(target_pairs, function(g) {
  x <- rt[, g[1]]
  y <- rt[, g[2]]
  
  fit <- lm(y ~ x)
  
  data.frame(
    Gene1 = g[1],
    Gene2 = g[2],
    N = length(x),
    Unique_values_gene1 = length(unique(x)),
    Unique_values_gene2 = length(unique(y)),
    Repeated_XY_rows = sum(duplicated(data.frame(x, y))),
    Exactly_identical_vectors = all(x == y),
    Max_absolute_difference = max(abs(x - y)),
    Regression_intercept = unname(coef(fit)[1]),
    Regression_slope = unname(coef(fit)[2]),
    Max_absolute_residual = max(abs(residuals(fit))),
    Residual_SD_ratio = max(abs(residuals(fit))) / sd(y)
  )
}))

save_csv(pair_audit, "07_target_pair_value_audit.csv")

# 7. 留一法敏感性检查：每次去掉1个样本重新计算r
# 只有去掉后还剩至少3个样本才计算
loo_results <- do.call(rbind, lapply(target_pairs, function(g) {
  do.call(rbind, lapply(seq_len(nrow(rt)), function(i) {
    x <- rt[-i, g[1]]
    y <- rt[-i, g[2]]
    
    valid <- length(x) >= 3 && sd(x) > 0 && sd(y) > 0
    
    data.frame(
      Gene1 = g[1],
      Gene2 = g[2],
      Omitted_sample = rownames(rt)[i],
      N_remaining = length(x),
      Pearson_r = if (valid) cor(x, y) else NA_real_
    )
  }))
}))

save_csv(loo_results, "08_leave_one_out_correlations.csv")

# 8. 创建BH校正P值矩阵，供第六张图使用
P_BH <- matrix(
  NA_real_, ncol(rt), ncol(rt),
  dimnames = list(colnames(rt), colnames(rt))
)

diag(P_BH) <- 0

for (k in seq_len(nrow(all_results))) {
  g1 <- all_results$Gene1[k]
  g2 <- all_results$Gene2[k]
  
  P_BH[g1, g2] <- all_results$P_adj_BH[k]
  P_BH[g2, g1] <- all_results$P_adj_BH[k]
}

save_csv(
  data.frame(
    N_samples = nrow(rt),
    N_genes = ncol(rt),
    N_unique_pair_tests = nrow(all_results),
    Correlation_method = "Pearson",
    Adjustment = "Benjamini-Hochberg"
  ),
  "09_analysis_summary.csv"
)

capture.output(
  sessionInfo(),
  file = file.path(outdir, "sessionInfo.txt")
)

cat("\n核查输出目录：", normalizePath(outdir), "\n")




# ==========================================================
# 散点图：每个点代表一个样本，标出样本ID
# 同时输出PDF和300 dpi TIFF
# ==========================================================

draw_target_scatter <- function(g1, g2) {
  x <- rt[, g1]
  y <- rt[, g2]
  stat <- get_pair_result(g1, g2)
  
  # 给样本标签预留空间
  xpad <- diff(range(x)) * 0.18
  ypad <- diff(range(y)) * 0.18
  
  par(mar = c(5, 5, 5, 2))
  
  plot(
    x, y,
    pch = 21,
    bg = "#008837",
    col = "white",
    cex = 1.7,
    xlim = range(x) + c(-xpad, xpad),
    ylim = range(y) + c(-ypad, ypad),
    xlab = paste0(g1, " expression (input scale)"),
    ylab = paste0(g2, " expression (input scale)"),
    main = paste(g1, "vs", g2),
    cex.lab = 1.1,
    bty = "l"
  )
  
  abline(lm(y ~ x), col = "#7b3294", lwd = 2)
  
  text(
    x, y,
    labels = rownames(rt),
    pos = rep(c(3, 1), length.out = length(x)),
    cex = 0.65
  )
  
  mtext(
    paste0(
      "n = ", stat$N,
      "; Pearson r = ", sprintf("%.8f", stat$Pearson_r),
      "\nP = ", format.pval(stat$P_value, digits = 3, eps = 1e-16),
      "; BH-adjusted P = ",
      format.pval(stat$P_adj_BH, digits = 3, eps = 1e-16)
    ),
    side = 3, line = 0.4, cex = 0.82
  )
}

for (g in target_pairs) {
  file_base <- paste(g, collapse = "_")
  
  pdf(
    file.path(outdir, paste0(file_base, "_scatter.pdf")),
    width = 7, height = 6
  )
  draw_target_scatter(g[1], g[2])
  dev.off()
  
  tiff(
    file.path(outdir, paste0(file_base, "_scatter.tiff")),
    width = 7, height = 6,
    units = "in", res = 300,
    compression = "lzw"
  )
  draw_target_scatter(g[1], g[2])
  dev.off()
}



#绘制相关性图形
tiff(file="1.corpot1.tiff",width=14,height=14,  units = "in",
     res = 300,
     compression = "lzw")
corrplot(M,
         method = "circle",
         order = "hclust", #聚类
         type=showtype,
         col=colorRampPalette(c(lowcol, midcol, higcol))(50)
)
dev.off()

#第二个图
tiff(file="2.corpot2.tiff",width=20,height=20, units = "in",
     res = 300,
     compression = "lzw")
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7, #相关系数
         addCoef.col = "black",
         diag = TRUE,
         type=showtype,
         tl.col="black",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50))
dev.off()

#第三个图
tiff(file="3.corpot3.tiff",width=20,height=20, units = "in",
     res = 300,
     compression = "lzw")
corrplot(M,
         type="upper",
         tl.pos="tp",
         order="AOE",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         tl.cex = tlcex)
corrplot(M,add=TRUE, 
         type="lower",
         method="number",
         order="AOE",
         diag=FALSE,
         tl.pos="n",
         cl.pos="n",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         number.cex=numbercex
         )
dev.off()

#第4个图
tiff(file="4.corpot4.tiff",width=20,height=20, units = "in",
     res = 300,
     compression = "lzw")
corrplot(
  M, 
  order = 'AOE',method = 'pie',
  type = 'lower', 
  tl.pos = 'd',
  col=colorRampPalette(c(lowcol, midcol, higcol))(50),
  tl.cex = tlcex)
corrplot(
  M, 
  add = TRUE, 
  type = 'upper',  
  method = 'number',
  order = 'AOE', 
  diag = FALSE,  
  tl.pos = 'n', 
  cl.pos = 'n',
  col=colorRampPalette(c(lowcol, midcol, higcol))(50),
  number.cex=numbercex)
dev.off()


#第五个图
tiff(file="5.corpot5.tiff",width=25,height=25, units = "in",
     res = 300,
     compression = "lzw")
chart.Correlation(rt,method = method)
dev.off()

#第6个图
tiff(file="6.corpot6.tiff", width=15, height=15, units = "in",
     res = 300,
     compression = "lzw")
corrplot(M,
         order="original",
         method = "circle",
         type=showtype,
         tl.cex=0.8, pch=T,
         p.mat = res$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         tl.col="black")
dev.off()

#第七个图
#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
M[M==1]=0  #删掉相关性=1(即自己与自己做相关性)
c1 = ifelse(c(M)>=0,rgb(1,0,0,abs(M)),rgb(0,1,0,abs(M)))
col1 = matrix(c1,nc=ncol(rt))
tiff(file="7.corpot7.tiff",width=15,height=15, units = "in",
     res = 300,
     compression = "lzw")
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(M)-1)), start.degree = 180)
chordDiagram(M, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))       #????ͼ??
dev.off()
circos.clear()

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
