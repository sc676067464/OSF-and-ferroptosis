# 单基因 ROC：AUC及95%CI、固定训练阈值下的敏感性/特异性及95%CI
# 完整替代原LASSO脚本；无需glmnet，无交叉验证，无拼接Bootstrap样本。
# 表达矩阵：第一列基因名，第一行样本名；每列一个独立受试者。
# 分组表：无表头，恰好两列：样本ID、分组；不要求对照组排在前面。
# GSgene.txt：无表头，第一列候选基因。所有候选基因分别评价，不自动筛选。
# 本脚本未在用户实际数据上运行。一个基因也可以运行。

setwd("C:/Users/csu206/OneDrive/OSF_single_cell/Code/6.lasso")
train.exp <- "logexp.txt"
train.group <- "sample.txt"
gene.file <- "GSgene.txt"
C <- "Health"
D <- "OSF"
have.test <- TRUE
test.exp <- file.path("test", "GSE274203.txt")
test.group <- file.path("test", "sample.txt")
C.test <- "Health"
D.test <- "OSF"
train.name <- "Training"
test.name <- "GSE274203"
outdir <- "Single_gene_ROC_CI"

# TRUE：将训练集表达阈值直接用于测试集，计算测试集敏感性/特异性。
# 前提是两个数据集表达量单位、预处理和尺度可比；仅都取log2并不保证可比。
# 若不满足，改为FALSE：仍计算测试集AUC，但测试集阈值相关指标输出NA。
# 不要按测试集标签重新选方向或阈值，不要各自Z标准化后直接套用阈值。
threshold.transferable <- TRUE

if (!requireNamespace("pROC", quietly = TRUE))
  stop('请先运行 install.packages("pROC")')
use.test <- toupper(as.character(have.test)) %in% c("TRUE", "T")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(z, filename) {
  write.table(z, file.path(outdir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
read_expression <- function(path) {
  a <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, comment.char = "", quote = "")
  if (!nrow(a) || !ncol(a)) stop("表达矩阵为空：", path)
  if (anyDuplicated(rownames(a)) || anyDuplicated(colnames(a)))
    stop("表达矩阵有重复基因或样本名：", path)
  if (!all(vapply(a, is.numeric, logical(1)))) stop("含非数值表达列：", path)
  z <- t(as.matrix(a))
  storage.mode(z) <- "double"
  if (any(!is.finite(z))) stop("表达矩阵存在NA/Inf：", path)
  z
}
read_labels <- function(path, ids, control, disease) {
  a <- read.table(path, header = FALSE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE, comment.char = "", quote = "")
  if (ncol(a) != 2L) stop("分组表必须无表头且恰好两列：", path)
  a[] <- lapply(a, function(v) trimws(as.character(v)))
  if (anyNA(a) || any(a[[1]] == "") || anyDuplicated(a[[1]]))
    stop("分组表有缺失或重复ID：", path)
  if (control == disease || !all(a[[2]] %in% c(control, disease)))
    stop("分组名称不匹配：", path, "；实际分组：", paste(unique(a[[2]]), collapse = ", "))
  if (!setequal(a[[1]], ids)) stop("分组表和表达矩阵样本ID不完全一致：", path)
  labels <- a[[2]][match(ids, a[[1]])]
  y <- as.integer(labels == disease)
  if (min(table(factor(y, levels = 0:1))) < 2L)
    stop("本脚本的AUC区间计算要求每组至少2个独立样本：", path)
  list(y = y, labels = labels)
}

x <- read_expression(train.exp)
tr <- read_labels(train.group, rownames(x), C, D)
y <- tr$y
xt <- NULL
te <- NULL
if (use.test) {
  xt <- read_expression(test.exp)
  te <- read_labels(test.group, rownames(xt), C.test, D.test)
  if (length(intersect(rownames(x), rownames(xt))))
    stop("训练和测试存在相同样本ID，请核实是否重复使用样本。")
}
g <- read.table(gene.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                comment.char = "", quote = "")
genes <- unique(trimws(as.character(g[[1]])))
genes <- genes[!is.na(genes) & nzchar(genes)]
if (!length(genes)) stop("候选基因列表为空。")
audit <- data.frame(Gene = genes, In_training = genes %in% colnames(x),
                    In_test = if (use.test) genes %in% colnames(xt) else NA)
audit$Training_variable <- vapply(genes, function(gene) {
  if (!gene %in% colnames(x)) return(FALSE)
  length(unique(x[, gene])) > 1L
}, logical(1))
audit$Analysed <- audit$In_training & audit$Training_variable
write_tsv(audit, "gene_audit.tsv")
genes <- audit$Gene[audit$Analysed]
if (!length(genes)) stop("没有可分析的训练集候选基因，请检查gene_audit.tsv。")
# 训练集不变基因列入audit但跳过；测试集缺失基因不影响该基因的训练集分析。

roc_fixed <- function(y, score) {
  pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
}
choose_cutoff <- function(r) {
  a <- as.data.frame(pROC::coords(r, x = "all",
    ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE))
  a <- a[is.finite(a$threshold), , drop = FALSE]
  if (!nrow(a)) stop("无法确定有限阈值。")
  j <- a$sensitivity + a$specificity - 1
  best <- a[abs(j - max(j)) < 1e-12, , drop = FALSE]
  # 并列时选择较低的定向评分阈值，即较高敏感性；不看测试集。
  min(best$threshold)
}
fmt_ci <- function(est, low, high) {
  if (anyNA(c(est, low, high))) return("NA")
  sprintf("%.3f (%.3f-%.3f)", est, low, high)
}
evaluate_gene <- function(v, y, ids, labels, gene, dataset, role, sign, cutoff,
                          apply.threshold = TRUE) {
  # score只是定向表达量，不是患病概率。
  score <- sign * v
  r <- roc_fixed(y, score)
  auc <- as.numeric(pROC::auc(r))
  ci.messages <- character()
  ci <- tryCatch(withCallingHandlers(
    as.numeric(pROC::ci.auc(r, method = "delong", conf.level = 0.95)),
    warning = function(w) {
      ci.messages <<- c(ci.messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      ci.messages <<- c(ci.messages, conditionMessage(e))
      rep(NA_real_, 3L)
    })
  degenerate <- all(is.finite(ci)) && abs(ci[3] - ci[1]) < 1e-12
  ci.note <- if (degenerate) {
    "Degenerate DeLong CI; unreliable uncertainty estimate, not certainty"
  } else if (anyNA(ci)) {
    "AUC CI unavailable"
  } else "Approximate DeLong CI; very small sample, exploratory"
  pred <- rep(NA_integer_, length(y))
  tp <- fn <- tn <- fp <- sensitivity <- specificity <- NA_real_
  sens.ci <- spec.ci <- c(NA_real_, NA_real_)
  if (apply.threshold) {
    pred <- as.integer(score >= cutoff)
    tp <- sum(pred == 1 & y == 1); fn <- sum(pred == 0 & y == 1)
    tn <- sum(pred == 0 & y == 0); fp <- sum(pred == 1 & y == 0)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    sens.ci <- as.numeric(stats::binom.test(tp, tp + fn, conf.level = 0.95)$conf.int)
    spec.ci <- as.numeric(stats::binom.test(tn, tn + fp, conf.level = 0.95)$conf.int)
  }
  metrics <- data.frame(
    Gene = gene, Dataset = dataset, Evaluation = role,
    N = length(y), Cases = sum(y == 1), Controls = sum(y == 0),
    Direction_from_training = if (sign == 1) "Higher_expression_is_OSF" else "Lower_expression_is_OSF",
    AUC = auc, AUC_CI_low = ci[1], AUC_CI_high = ci[3],
    AUC_95CI = paste0(fmt_ci(auc, ci[1], ci[3]), if (degenerate) " *" else ""),
    AUC_CI_note = ci.note, AUC_CI_warnings = paste(unique(ci.messages), collapse = "; "),
    Training_expression_threshold = sign * cutoff,
    Disease_rule = if (sign == 1) "Expression >= threshold" else "Expression <= threshold",
    Threshold_applied = apply.threshold,
    Sensitivity = sensitivity, Sensitivity_CI_low = sens.ci[1], Sensitivity_CI_high = sens.ci[2],
    Sensitivity_95CI = fmt_ci(sensitivity, sens.ci[1], sens.ci[2]),
    Specificity = specificity, Specificity_CI_low = spec.ci[1], Specificity_CI_high = spec.ci[2],
    Specificity_95CI = fmt_ci(specificity, spec.ci[1], spec.ci[2]),
    TP = tp, FN = fn, TN = tn, FP = fp,
    Control_mean = mean(v[y == 0]), Case_mean = mean(v[y == 1]),
    Mean_difference = mean(v[y == 1]) - mean(v[y == 0])
  )
  samples <- data.frame(Gene = gene, Dataset = dataset, ID = ids, Group = labels,
    Outcome = y, Expression = v, Oriented_expression = score, Predicted_case = pred)
  list(roc = r, metrics = metrics, samples = samples)
}
plot_result <- function(z, v, y, gene, dataset) {
  grp <- factor(y, levels = 0:1, labels = c("Control", "OSF"))
  stripchart(v ~ grp, vertical = TRUE, method = "stack", pch = 19,
    col = c("#377EB8", "#E64B35"), xlab = "", ylab = "Expression (input scale)",
    main = paste(gene, dataset))
  means <- tapply(v, grp, mean)
  segments(c(0.8, 1.8), means, c(1.2, 2.2), means, lwd = 2)
  pROC::plot.roc(z$roc, col = "#377EB8", lwd = 2, legacy.axes = TRUE,
    main = paste(gene, "ROC"))
  legend("bottomright", bty = "n", cex = 0.72, legend = c(
    paste("AUC (95% CI):", z$metrics$AUC_95CI),
    paste("Sensitivity:", z$metrics$Sensitivity_95CI),
    paste("Specificity:", z$metrics$Specificity_95CI)))
  if (grepl("Degenerate", z$metrics$AUC_CI_note))
    mtext("* Degenerate CI: uncertainty not reliably estimated", side = 1, line = 3, cex = 0.65)
}

run_analysis <- function() {
  rows <- samples <- thresholds <- list()
  pdf(file.path(outdir, "Single_gene_ROC.pdf"), width = 11, height = 5.5)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  for (gene in genes) {
    v <- x[, gene]
    # 预定规则：训练集病例均值较大时sign=1，否则sign=-1；均值相等取1。
    # 若已有预先指定的生物学方向，应在查看结果前修改规则。
    sign <- if (mean(v[y == 1]) >= mean(v[y == 0])) 1 else -1
    cutoff <- choose_cutoff(roc_fixed(y, sign * v))
    a <- evaluate_gene(v, y, rownames(x), tr$labels, gene, train.name,
      "Apparent training; direction and threshold selected here", sign, cutoff)
    rows[[length(rows) + 1L]] <- a$metrics
    samples[[length(samples) + 1L]] <- a$samples
    plot_result(a, v, y, gene, train.name)
    thresholds[[length(thresholds) + 1L]] <- data.frame(Gene = gene, Score_sign = sign,
      Oriented_score_threshold = cutoff, Expression_threshold = sign * cutoff,
      Disease_rule = a$metrics$Disease_rule,
      Selection = "Training Youden among finite thresholds; lowest score cutoff on ties")
    if (use.test && gene %in% colnames(xt)) {
      b <- evaluate_gene(xt[, gene], te$y, rownames(xt), te$labels, gene, test.name,
        "External; training direction and threshold fixed", sign, cutoff, threshold.transferable)
      rows[[length(rows) + 1L]] <- b$metrics
      samples[[length(samples) + 1L]] <- b$samples
      plot_result(b, xt[, gene], te$y, gene, test.name)
    }
  }
  list(metrics = do.call(rbind, rows), samples = do.call(rbind, samples),
       thresholds = do.call(rbind, thresholds))
}
results <- run_analysis()
write_tsv(results$metrics, "Single_gene_performance.tsv")
write_tsv(results$samples, "Sample_predictions.tsv")
write_tsv(results$thresholds, "Training_thresholds.tsv")
write_tsv(results$metrics[, c("Gene", "Dataset", "Cases", "Controls", "AUC_95CI",
  "Training_expression_threshold", "Disease_rule", "Sensitivity_95CI", "Specificity_95CI",
  "TP", "FN", "TN", "FP", "AUC_CI_note")], "Results_for_table.tsv")
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
writeLines(c(
  "单基因探索性分析；每组2例不足以可靠估计诊断性能。",
  "AUC的95%CI使用DeLong近似法。退化CI仍保留原始数值，汇总列加*并附警示。",
  "AUC=1且CI=1-1不代表完美诊断；Bootstrap CI同样可能退化。",
  "敏感性/特异性区间使用Clopper-Pearson精确二项法。数值0-1，不是百分数。",
  "例如2/2=100%，精确95%CI约15.8%-100%；0/2的区间约0%-84.2%。",
  "训练方向按病例和对照均值确定；均值相等取高表达为病例。",
  "训练阈值按有限候选阈值中的最大Youden指数确定；并列取较高敏感性的阈值。",
  "训练集性能属于表观性能。区间未计入方向选择、阈值选择和候选基因筛选的不确定性。",
  "测试集沿用训练方向；AUC<0.5也原样报告，不根据测试结果翻转方向。",
  "测试敏感性/特异性使用固定训练表达阈值，前提是两数据集尺度兼容。",
  "两平台均为log2表达并不保证阈值可迁移；不可比时将threshold.transferable改为FALSE。",
  "此时测试AUC仍输出，测试敏感性/特异性为NA，不用测试标签重新优化阈值。",
  "输入表达量不自动log2或标准化。所有输入样本必须独立，候选基因来源需如实记录。",
  "不输出LASSO系数或联合评分；表达量不是患病概率，不计算校准/Brier。",
  "gene_audit.tsv记录基因缺失和训练集中无变异情况。",
  "缺失测试基因只跳过该基因外部评价。测试集恒定表达可得AUC=0.5且CI退化。",
  "无需按AUC选出最好的基因；任何数据驱动筛选均应另行独立验证。",
  "参考：pROC https://cran.r-project.org/web/packages/pROC/refman/pROC.html",
  "参考：binom.test https://stat.ethz.ch/R-manual/R-devel/library/stats/html/binom.test.html"
), file.path(outdir, "analysis_notes.txt"), useBytes = TRUE)
print(results$metrics[, c("Gene", "Dataset", "AUC_95CI", "Sensitivity_95CI", "Specificity_95CI")])
message("完成。输出目录：", normalizePath(outdir))
