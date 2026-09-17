################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（21Q1 → 26Q1）

### ================================================================
### 主题：拷贝数变异（CNV）下的易感基因筛选流程（In silico CRISPR Screen）
### 灵感来源：https://www.nature.com/articles/s41586-022-04638-9
###
### 分析思路：
###   CCNE1（cyclin E1）扩增的肿瘤细胞系 → 对哪些基因更加依赖？
###   经典发现：CCNE1 扩增 → 依赖 PKMYT1（已进入临床试验）
###
### 分析流程：
###   第一部分：用已有批量结果画火山图（展示 CCNE1 扩增后的已知依赖性）
###   第二部分：正式计算（新版数据从头跑）
###     - GeneEffect 数据（ACH-编号）
###     - CNV 数据（MC-编号 → 映射为 ACH-编号）
###     - CCNE1 扩增组 vs 野生型组
###     - 批量 Wilcoxon + t-test
###   第三部分：火山图复现 + median vs mean 对比
###
### 新版数据格式变化：
###   旧版 Achilles_gene_effect.csv → 新版 CRISPRGeneEffect.csv（格式一致，ACH-编号）
###   旧版 CCLE_gene_cn.csv → 新版 OmicsCNGeneMC_WES.csv
###     关键区别：新版用 ModelConditionID（MC-编号），需通过 ModelCondition.csv 映射为 ACH-编号
### ================================================================

############################################################
## =================== 第一部分：已有结果火山图 ===================
## 使用 07 脚本预处理时保存的 batchData RDS 画图
## （原脚本使用 Nature 论文的 source data，此处改为用自有计算结果）
############################################################
rm(list = ls())

## 读取已有的批量计算结果（如不存在则跳过此部分，直接看第二部分）
batchDataFile <- "TM00/DepMap_TM00/output/CCNE1_amp_Nature_fig1c_batchData.rds"
if(file.exists(batchDataFile)) {
  data <- readRDS(batchDataFile)

  library(ggplot2)       ## 绘图引擎
  library(ggrepel)       ## 标签智能避让

  ## 要标注的关键基因（来自 Nature 论文的发现）
  genes_to_show = c("CCNE1", "PKMYT1", "CDK2", "ANAPC15", "FBXW7", "UBE2S", "UBE2C", "WEE1")

  ## 火山图：dep.FC vs -log10(Wilcoxon p-value)
  ## dep.FC < 0 = 扩增组对该基因更依赖（negative = more dependent）
  p0 <- ggplot(data, aes(dep.FC, -log10(Wilcox_p_value))) +
    geom_point(shape = 21, colour = "black", fill = "lightyellow") +
    geom_point(data = subset(data, GeneSymbol %in% genes_to_show),
               shape = 21, colour = "black", fill = "skyblue3", size = 3) +
    geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), aes(label = GeneSymbol)) +
    theme_bw() +
    xlab(expression(paste(Delta, "FC (CCNE1 amplification vs WT)"))) +
    ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
    labs(title = "DepMap dependencies", subtitle = "CCNE1-amplified tumour cell lines") +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 12))

  print(p0)
  ggsave("TM00/DepMap_TM00/output/12_CCNE1_source_volcano.png", p0, width = 8, height = 6, dpi = 300)
  cat("第一部分：已有结果火山图已保存\n\n")
} else {
  cat("第一部分：跳过（", batchDataFile, "不存在，请先运行第二部分）\n\n")
}


############################################################
## =================== 第二部分：正式计算（新版数据） ===================
## 从头计算：CCNE1 扩增组 vs 野生型组 → 全基因组 Wilcoxon/t-test
############################################################
rm(list = ls())

## ---- 1. 读取基因效应数据（GeneEffect，ACH-编号）----
## CRISPRGeneEffect.csv：行=细胞系(ACH-编号)，列=基因
## 值 = 基因效应得分（越负 = 越依赖）
depData <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)

## 列名去掉 (EntrezID) 后缀，如 "A1BG (1)" → "A1BG"
colnames(depData) <- gsub("\\s+\\(\\d+\\)", "", colnames(depData))

## 第一列是无名行索引（ACH-编号），设为行名
rownames(depData) <- depData[, 1]
depData <- depData[, -1]

## ---- 2. 读取拷贝数数据（CNV，MC-编号）----
## OmicsCNGeneMC_WES.csv：
##   行=细胞系(ModelConditionID, MC-编号)，列=基因
##   值 = log2 拷贝数比值（相对于二倍体）
##   注意：新版用 MC-编号，需要映射为 ACH-编号
cnData_raw <- data.table::fread("data/OmicsCNGeneMC_WES.csv", data.table = F)

## 过滤默认条目（和突变矩阵处理方式一致）
cnData_raw <- cnData_raw[cnData_raw$IsDefaultEntryForMC == "Yes", ]

## 通过 ModelCondition.csv 将 MC-编号映射为 ACH-编号
mcMap <- data.table::fread("data/ModelCondition.csv", data.table = F,
                           select = c("ModelConditionID", "ModelID"))
cnData_raw$ModelID <- mcMap$ModelID[match(cnData_raw$ModelConditionID,
                                           mcMap$ModelConditionID)]

## 去除映射失败的行（NA）
cnData_raw <- cnData_raw[!is.na(cnData_raw$ModelID), ]

## 去除重复的 ModelID（一个 ACH 可能有多个默认 MC，取第一个）
cnData_raw <- cnData_raw[!duplicated(cnData_raw$ModelID), ]

## 设行名为 ACH-编号
rownames(cnData_raw) <- cnData_raw$ModelID

## 去掉元数据列，只保留基因列
cnData <- cnData_raw[, !(colnames(cnData_raw) %in% c("ModelConditionID", "IsDefaultEntryForMC", "ModelID"))]

## 列名去掉 (EntrezID) 后缀
colnames(cnData) <- gsub("\\s+\\(\\d+\\)", "", colnames(cnData))

cat("GeneEffect 细胞系数:", nrow(depData), "\n")
cat("CNV 细胞系数（映射后）:", nrow(cnData), "\n")

## ---- 3. 取交集细胞系 ----
## 只保留同时有 GeneEffect 和 CNV 数据的细胞系
index <- intersect(rownames(depData), rownames(cnData))
depData <- depData[index, ]
cnData <- cnData[index,]

cat("交集细胞系数:", length(index), "\n\n")

## ---- 4. 按 CCNE1 拷贝数分组 ----
## 阈值 1.58（log2 ratio）= 约 3 倍扩增，来自 Nature 论文
cnGene <- "CCNE1"
cnGeneData <- data.frame(DepMap_ID = rownames(cnData), cnGene = cnData[, cnGene])

## 扩增组：CCNE1 log2 CN > 1.58
ampCells <- cnGeneData$DepMap_ID[cnGeneData$cnGene > 1.58]
## 野生型组：CCNE1 log2 CN <= 1.58
wtCells <- cnGeneData$DepMap_ID[cnGeneData$cnGene <= 1.58]

cat("CCNE1 扩增组细胞系数:", length(ampCells), "\n")
cat("CCNE1 野生型组细胞系数:", length(wtCells), "\n\n")

## 提取两组的基因效应数据
ampCellsDepData <- depData[ampCells, ]
wtCellsDepData <- depData[wtCells, ]

## ---- 5. 单基因验证（CDKN1A）----
## 测试阳性对照：CDKN1A 是已知的 CCNE1 扩增敏感基因
cat("=== 阳性对照验证：CDKN1A ===\n")
cat("扩增组 median:", median(ampCellsDepData[, "CDKN1A"], na.rm = TRUE), "\n")
cat("野生组 median:", median(wtCellsDepData[, "CDKN1A"], na.rm = TRUE), "\n")
cat("dep.FC (median):", median(ampCellsDepData[, "CDKN1A"], na.rm = TRUE) -
                          median(wtCellsDepData[, "CDKN1A"], na.rm = TRUE), "\n")

## Wilcoxon 秩和检验（非参数）
test_wilcox <- wilcox.test(ampCellsDepData[, "CDKN1A"],
                           wtCellsDepData[, "CDKN1A"],
                           alternative = "less")
cat("Wilcoxon p-value:", test_wilcox$p.value, "\n")

## t 检验（参数检验）
test_t <- t.test(ampCellsDepData[, "CDKN1A"],
                 wtCellsDepData[, "CDKN1A"],
                 alternative = "less")
cat("t-test p-value:", test_t$p.value, "\n\n")

## ---- 6. 批量遍历所有基因 ----
## 对每个基因计算：扩增组 vs 野生组的 median/mean 差异 + 统计检验
results <- data.frame()
cat("共", ncol(depData), "个基因需要遍历\n")

for (i in 1:ncol(depData)) {
  if(i %% 2000 == 0) cat("进度:", i, "/", ncol(depData), "\n")

  ## 提取该基因在两组中的效应得分
  ampD <- ampCellsDepData[, i]
  wtD  <- wtCellsDepData[, i]

  ## 如果某组有效值太少（<3），跳过统计检验（避免报错）
  if(sum(!is.na(ampD)) < 3 || sum(!is.na(wtD)) < 3) {
    results[i, 1] <- colnames(depData)[i]
    next
  }

  ## 1. GeneSymbol
  results[i, 1] <- colnames(depData)[i]

  ## 2-4. median 差异
  results[i, 2] <- median(ampD, na.rm = TRUE)   ## 扩增组 median
  results[i, 3] <- median(wtD, na.rm = TRUE)    ## 野生组 median
  results[i, 4] <- median(ampD, na.rm = TRUE) - median(wtD, na.rm = TRUE)  ## dep.FC

  ## 5-7. mean 差异
  results[i, 5] <- mean(ampD, na.rm = TRUE)     ## 扩增组 mean
  results[i, 6] <- mean(wtD, na.rm = TRUE)      ## 野生组 mean
  results[i, 7] <- mean(ampD, na.rm = TRUE) - mean(wtD, na.rm = TRUE)      ## dep.FC.mean

  ## 8. Wilcoxon 秩和检验（单侧：扩增组更依赖 = 得分更低）
  results[i, 8] <- tryCatch(
    wilcox.test(ampD, wtD, alternative = "less")$p.value,
    error = function(e) NA
  )

  ## 9. t 检验（单侧）
  results[i, 9] <- tryCatch(
    t.test(ampD, wtD, alternative = "less")$p.value,
    error = function(e) NA
  )
}

## 设置列名
colnames(results) <- c("GeneSymbol",
                       "amp_median_dep",
                       "wt_median_dep",
                       "dep.FC",
                       "amp_mean_dep",
                       "wt_mean_dep",
                       "dep.FC.mean",
                       "Wilcox_p_value",
                       "Ttest_p_value")

## 增加 FDR 校正的 p 值
results$Wilcox_p_value_fdr <- p.adjust(results$Wilcox_p_value, method = "fdr")
results$Ttest_p_value_fdr  <- p.adjust(results$Ttest_p_value, method = "fdr")

## 保存结果
saveRDS(results, file = "TM00/DepMap_TM00/output/CCNE1_amp_Nature_fig1c_batchData.rds")
cat("\n批量计算完成，结果已保存\n")


############################################################
## =================== 第三部分：火山图复现 + 对比 ===================
## 用第二部分的结果画火山图，比较 median vs mean
############################################################
rm(list = ls())

## 读取批量计算结果
data <- readRDS("TM00/DepMap_TM00/output/CCNE1_amp_Nature_fig1c_batchData.rds")

library(ggplot2)
library(ggrepel)
library(patchwork)    ## 图形拼接

## 要标注的关键基因
genes_to_show <- c("CCNE1", "PKMYT1", "CDK2", "ANAPC15", "FBXW7", "UBE2S", "UBE2C", "WEE1")

## ---- 火山图1：median FC ----
## dep.FC < 0 表示扩增组更依赖该基因（是潜在治疗靶点）
p1 <- ggplot(data, aes(dep.FC, -log10(Wilcox_p_value))) +
  geom_point(shape = 21, colour = "black", fill = "lightyellow") +
  geom_point(data = subset(data, GeneSymbol %in% genes_to_show),
             shape = 21, colour = "black", fill = "skyblue3", size = 3) +
  geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), aes(label = GeneSymbol)) +
  theme_bw() +
  xlab(expression(paste(Delta, "FC (CCNE1 amplification vs WT)"))) +
  ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
  labs(title = "DepMap dependencies (median)",
       subtitle = "CCNE1-amplified tumour cell lines") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

print(p1)
ggsave("TM00/DepMap_TM00/output/12_CCNE1_volcano_median.png", p1, width = 8, height = 6, dpi = 300)

## ---- 火山图2：mean FC ----
## 对比 mean 和 median 的差异
p2 <- ggplot(data, aes(dep.FC.mean, -log10(Wilcox_p_value))) +
  geom_point(shape = 21, colour = "black", fill = "lightyellow") +
  geom_point(data = subset(data, GeneSymbol %in% genes_to_show),
             shape = 21, colour = "black", fill = "skyblue3", size = 3) +
  geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), aes(label = GeneSymbol)) +
  theme_bw() +
  xlab(expression(paste(Delta, "FC (CCNE1 amplification vs WT)"))) +
  ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
  labs(title = "DepMap dependencies (mean)",
       subtitle = "CCNE1-amplified tumour cell lines") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

print(p2)
ggsave("TM00/DepMap_TM00/output/12_CCNE1_volcano_mean.png", p2, width = 8, height = 6, dpi = 300)

## ---- 并排对比 median vs mean ----
p_combo <- p1 | p2
print(p_combo)
ggsave("TM00/DepMap_TM00/output/12_CCNE1_volcano_combo.png", p_combo, width = 14, height = 6, dpi = 300)
cat("第三部分：火山图（median + mean + 对比图）已保存\n")

cat("\n=== 脚本全部完成 ===\n")
