################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（24Q4 → 26Q1，Sanger 独立文件 → 统一 DepMap 门户）

### ================================================================
### 主题：Sanger (KY 文库) 与 DepMap Broad (Avana 文库) 两套独立
###       CRISPR 筛选数据的一致性比较
###
### 背景：
###   旧版 Sanger 数据存放在独立的 Sanger_CRISPR_Chronos/ 目录。
###   新版 26Q1 中，Sanger 数据已整合进 DepMap 统一门户：
###     - ScreenGeneEffect.csv 包含所有文库（Avana/KY/Brunello/...）
###     - CRISPRGeneEffect.csv 是 Broad Avana 合并后结果（每模型一行）
###   因此需要从 ScreenGeneEffect.csv 中提取 KY (Sanger) 文库的数据，
###   再与 CRISPRGeneEffect.csv (Broad) 比较。
###
### 分析内容：
###   1. 数据加载与探索（共有细胞系、共有基因）
###   2. 每个细胞系：DepMap vs Sanger 基因效应相关性（Pearson + Spearman）
###   3. 每个基因：跨细胞系的相关性（全部基因 + 高 SD 基因子集）
###   4. 单基因散点图（以 RPS29 为例）
###   5. Top-100 必需基因重合度
### ================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

output_dir <- "TM00/DepMap_TM00/output"

#############################################################
## 第一部分：数据加载
#############################################################

## ---- 1. 读取 DepMap Broad (Avana) 基因效应数据 ----
## CRISPRGeneEffect.csv: 行 = ACH-编号（每模型一行），列 = 基因
## 列名格式 "A1BG (1)"，需要去掉 (EntrezID) 后缀
cat("正在读取 DepMap Broad 数据...\n")
depData <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
colnames(depData) <- gsub("\\s+\\(\\d+\\)", "", colnames(depData))  ## 去掉 " (123)" 后缀
rownames(depData) <- depData[, 1]     ## 第一列 V1 = ACH-编号
depData <- depData[, -1]
depData_depmap <- as.matrix(depData)
rm(depData)
cat("  DepMap Broad:", nrow(depData_depmap), "模型 ×", ncol(depData_depmap), "基因\n")

## ---- 2. 读取 Sanger (KY 文库) 基因效应数据 ----
## ScreenGeneEffect.csv 包含所有文库的 per-screen 数据
## 需要用 ScreenSequenceMap.csv 做 ScreenID → ModelID 映射
cat("正在读取 Sanger KY 数据（从 ScreenGeneEffect 提取）...\n")
screenData <- data.table::fread("data/ScreenGeneEffect.csv", data.table = F)
colnames(screenData) <- gsub("\\s+\\(\\d+\\)", "", colnames(screenData))
rownames(screenData) <- screenData[, 1]
screenData <- screenData[, -1]

## ScreenSequenceMap: ScreenID → ModelID, Library
sm <- data.table::fread("data/ScreenSequenceMap.csv", data.table = F,
                        select = c("ScreenID", "ModelID", "Library"))
sm <- unique(sm)

## 提取 KY (Sanger) 文库的 ScreenID
ky_screens <- sm$ScreenID[sm$Library == "KY"]
ky_screen_to_model <- setNames(sm$ModelID[sm$Library == "KY"], sm$ScreenID[sm$Library == "KY"])

## 从 ScreenGeneEffect 中提取 KY 行
ky_in_data <- rownames(screenData)[rownames(screenData) %in% ky_screens]
cat("  KY screens in ScreenGeneEffect:", length(ky_in_data), "\n")

depData_sanger <- as.matrix(screenData[ky_in_data, ])
rm(screenData)

## 将行名从 ScreenID 映射为 ModelID（ACH-编号）
## 若同一模型有多个 KY screen，取平均值
rownames(depData_sanger) <- ky_screen_to_model[rownames(depData_sanger)]

## 合并同一模型的多个 screen（取行均值）
dup_models <- names(which(table(rownames(depData_sanger)) > 1))
if (length(dup_models) > 0) {
  for (m in dup_models) {
    idx <- which(rownames(depData_sanger) == m)
    depData_sanger[idx[1], ] <- colMeans(depData_sanger[idx, , drop = F], na.rm = T)
  }
  depData_sanger <- depData_sanger[!duplicated(rownames(depData_sanger)), ]
}
cat("  Sanger KY:", nrow(depData_sanger), "模型 ×", ncol(depData_sanger), "基因\n")

## ---- 3. 读取细胞系信息 ----
cellinfor <- readRDS(file = "TM00/DepMap_TM00/output/cellinfor.rds")

#############################################################
## 第二部分：数据探索 —— 共有细胞系和基因
#############################################################

## 共有细胞系：同时被 Sanger KY 和 Broad Avana 筛过的模型
common_cell <- intersect(rownames(depData_sanger), rownames(depData_depmap))
common_cell <- intersect(cellinfor$ModelID, common_cell)
cat("\n共有细胞系:", length(common_cell), "\n")

## 共有基因
common_gene <- intersect(colnames(depData_sanger), colnames(depData_depmap))
cat("共有基因:", length(common_gene), "\n")

## 取交集子集（行 = 共有细胞系，列 = 共有基因）
depData_sanger <- depData_sanger[common_cell, common_gene]
depData_depmap <- depData_depmap[common_cell, common_gene]
cat("对齐后矩阵:", nrow(depData_depmap), "×", ncol(depData_depmap), "\n")

## 检查 3 个示例细胞系的相关性（快速验证数据对齐）
cat("\n=== 示例细胞系 Pearson 相关 ===\n")
for (id in c("ACH-000001", "ACH-000007", "ACH-000227")) {
  if (id %in% common_cell) {
    ct <- cor.test(depData_depmap[id, ], depData_sanger[id, ])
    cat(id, ":", round(ct$estimate, 4), "(p =", format(ct$p.value, digits = 3), ")\n")
  }
}

#############################################################
## 第三部分：每个细胞系的 Pearson + Spearman 相关性
#############################################################

cat("\n正在计算每个细胞系的相关性...\n")

## 转置矩阵：行 = 基因，列 = 细胞系（便于按列取基因效应向量）
depData_depmap_t <- t(depData_depmap)
depData_sanger_t <- t(depData_sanger)

corData <- data.frame(ModelID = common_cell, stringsAsFactors = F)
for (i in seq_along(common_cell)) {
  if (i %% 50 == 0) cat("  进度:", i, "/", length(common_cell), "\n")

  dd_pearson  <- cor.test(depData_depmap_t[, i], depData_sanger_t[, i], method = "pearson")
  dd_spearman <- cor.test(depData_depmap_t[, i], depData_sanger_t[, i], method = "spearman")

  corData$Cor_pearson[i]  <- dd_pearson$estimate
  corData$pvalue_pearson[i]  <- dd_pearson$p.value
  corData$Cor_spearman[i] <- dd_spearman$estimate
  corData$pvalue_spearman[i] <- dd_spearman$p.value
}
cat("  完成！\n")

## 整理为长格式（Pearson vs Spearman），便于并排画图
corData1 <- corData[, c("ModelID", "Cor_pearson", "pvalue_pearson")]
colnames(corData1) <- c("ModelID", "Cor", "pvalue")
corData1$group <- "Pearson"
corData1$method_mean <- mean(corData1$Cor, na.rm = T)

corData2 <- corData[, c("ModelID", "Cor_spearman", "pvalue_spearman")]
colnames(corData2) <- c("ModelID", "Cor", "pvalue")
corData2$group <- "Spearman"
corData2$method_mean <- mean(corData2$Cor, na.rm = T)

corDataPlot <- rbind(corData1, corData2)
cat("  Pearson 均值:", round(mean(corData1$Cor, na.rm = T), 4), "\n")
cat("  Spearman 均值:", round(mean(corData2$Cor, na.rm = T), 4), "\n")

## 细胞系水平相关性密度图（Pearson + Spearman 并排）
p1 <- ggplot(corDataPlot, aes(x = Cor, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  labs(title = "Per-cell-line: DepMap vs Sanger correlation",
       x = "Correlation coefficient", fill = "Method") +
  theme(legend.position = "top")
print(p1)
ggsave(file.path(output_dir, "15_cell_correlation_pearson_spearman.png"),
       p1, width = 7, height = 5, dpi = 150)
cat("  已保存: 15_cell_correlation_pearson_spearman.png\n")

saveRDS(corData, file.path(output_dir, "15_cell_corData.rds"))

#############################################################
## 第四部分：每个基因的相关性（全部 + 高 SD 子集）
#############################################################

cat("\n正在计算每个基因的相关性...\n")

corData_gene <- data.frame(Gene = common_gene, Cor = NA_real_, pvalue = NA_real_,
                           stringsAsFactors = F)
for (i in seq_along(common_gene)) {
  if (i %% 2000 == 0) cat("  进度:", i, "/", length(common_gene), "\n")

  ## 某些基因在某一组中 NA 太多，cor.test 会报错，用 tryCatch 保护
  dd <- tryCatch(
    cor.test(depData_depmap[, i], depData_sanger[, i]),
    error = function(e) NULL
  )
  if (!is.null(dd)) {
    corData_gene$Cor[i]    <- dd$estimate
    corData_gene$pvalue[i] <- dd$p.value
  }
}
cat("  完成！有效基因:", sum(!is.na(corData_gene$Cor)), "/", nrow(corData_gene), "\n")

## 计算每个基因在 DepMap 数据中的 SD（衡量变异程度）
sdData <- apply(depData_depmap, 2, sd, na.rm = T)
corData_gene$SD <- sdData

## 高 SD 基因子集（SD > 0.3，说明效应值变异大，更值得比较）
corData_gene_sd <- corData_gene[corData_gene$SD > 0.3, ]
cat("  全部基因:", nrow(corData_gene), "  高 SD 基因:", nrow(corData_gene_sd), "\n")

## 基因水平相关性合并图（全部 vs 高 SD）
corData_gene$group <- "ALL genes"
corData_gene_sd$group <- "High SD genes"
gene_plot_data <- rbind(corData_gene, corData_gene_sd)

p2 <- ggplot(gene_plot_data, aes(x = Cor, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.3, bins = 50) +
  geom_density(alpha = 0.2) +
  theme_bw() +
  labs(title = "Per-gene: DepMap vs Sanger correlation",
       x = "Correlation coefficient", fill = "") +
  theme(legend.position = "top")
print(p2)
ggsave(file.path(output_dir, "15_gene_correlation_all_vs_highSD.png"),
       p2, width = 7, height = 5, dpi = 150)
cat("  已保存: 15_gene_correlation_all_vs_highSD.png\n")

saveRDS(corData_gene, file.path(output_dir, "15_gene_corData.rds"))

#############################################################
## 第五部分：单基因散点图（示例：RPS29）
#############################################################

## 选择一个在两套数据中都有数据的基因
gene <- "RPS29"
if (gene %in% common_gene) {
  corplotdata <- data.frame(
    depmapData = depData_depmap[, gene],
    sangerData = depData_sanger[, gene]
  )
  ct_gene <- cor.test(corplotdata$depmapData, corplotdata$sangerData)

  p3 <- ggplot(corplotdata, aes(x = depmapData, y = sangerData)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", color = "red", se = T) +
    theme_bw() +
    labs(title = paste0(gene, ": DepMap vs Sanger (r = ",
                        round(ct_gene$estimate, 3), ")"),
         x = "DepMap Broad Gene Effect", y = "Sanger KY Gene Effect")
  print(p3)
  ggsave(file.path(output_dir, paste0("15_", gene, "_scatter.png")),
         p3, width = 5, height = 5, dpi = 150)
  cat("  已保存: 15_", gene, "_scatter.png\n")
} else {
  cat("  基因", gene, "不在共有基因列表中，跳过散点图\n")
}

#############################################################
## 第六部分：Top-100 必需基因重合度
#############################################################

cat("\n正在计算 Top-100 必需基因重合度...\n")

## 对每个细胞系，分别在两套数据中找到效应值最低的 100 个基因（最必需）
genes <- common_gene

data1 <- apply(depData_depmap, 1, function(row) {
  genes[order(row)[1:100]]
})
data1 <- as.data.frame(data1)

data2 <- apply(depData_sanger, 1, function(row) {
  genes[order(row)[1:100]]
})
data2 <- as.data.frame(data2)

## 逐个细胞系计算两套数据 Top-100 的交集大小
cell_lines <- colnames(data1)
intersection_sizes <- sapply(cell_lines, function(col) {
  length(intersect(data1[, col], data2[, col]))
})

## 画密度图
overlap_df <- data.frame(overlap = intersection_sizes)
p4 <- ggplot(overlap_df, aes(x = overlap)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey70", color = "white") +
  geom_density(color = "steelblue", linewidth = 0.8) +
  theme_bw() +
  labs(title = "Top-100 essential gene overlap per cell line",
       x = "Number of overlapping genes (out of 100)", y = "Density")
print(p4)
ggsave(file.path(output_dir, "15_top100_overlap_density.png"),
       p4, width = 6, height = 5, dpi = 150)

cat("  中位数 overlap:", median(intersection_sizes), "\n")
cat("  均值 overlap:", round(mean(intersection_sizes), 1), "\n")
cat("  已保存: 15_top100_overlap_density.png\n")

cat("\n========== 15 脚本全部完成 ==========\n")
