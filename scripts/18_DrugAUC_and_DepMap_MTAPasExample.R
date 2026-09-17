################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（23Q2 → 26Q1）

### ================================================================
### 主题：药物敏感性（AMG193 AUC）与 DepMap 多组学数据的联合分析
###
### 背景：
###   MTAP、CDKN2A/B 位于染色体 9p21.3 区域
###   MTAP 是代谢底物甲硫腺苷（MTA）的降解酶，缺失后 MTA 积累
###   MTA 是 PRMT5 的天然抑制剂 → MTAP 缺失 → MTA 积累 → PRMT5 抑制
###   AMG 193 是安进公司的 PRMT5 抑制剂
###
###   通过将 AMG193 的 AUC（药物敏感性指标）与 DepMap 的多种组学数据
###   做批量 Pearson 相关性，找到预测药物敏感性的生物标志物。
###
### 四个分析模块（同一个 AUC 数据 × 不同组学数据源）：
###   1. CRISPR Gene Effect（基因依赖性）
###   2. RNAi DEMETER2（RNAi 基因依赖性）
###   3. 拷贝数变异（CNV）
###   4. 基因表达量
###
### 参考阅读:
###   https://www.nature.com/articles/s41588-025-02336-6
###   https://mp.weixin.qq.com/s/chi8OmbcYuQPko1JwGc5DQ
### ================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggrepel)

output_dir <- "TM00/DepMap_TM00/output"
auc_file <- "TM00/DepMap_TM00/resource/AMG193AUC.csv"

#############################################################
## 通用函数：加载 AUC 数据
#############################################################

load_AUC <- function(auc_path) {
  aucData <- data.table::fread(auc_path, data.table = F) %>%
    filter(!is.na(AMG193_AUC)) %>%
    distinct(model_id, .keep_all = T)
  ## 用 setNames 显式保留 names 属性（避免 $ 提取时丢失）
  return(setNames(aucData$AMG193_AUC, aucData$model_id))
}

#############################################################
## 通用函数：批量相关性 + 火山图
#############################################################

## 批量计算基因特征 vs AUC 的 Pearson 相关性
batch_cor_auc <- function(featureMat, aucVec, progress_name = "") {
  cat("  正在批量计算相关性:", progress_name, "...\n")

  corResult <- apply(featureMat, 2, function(col) {
    dd <- tryCatch(
      cor.test(col, aucVec, method = "pearson"),
      error = function(e) NULL
    )
    if (is.null(dd)) {
      c(NA_real_, NA_real_)
    } else {
      c(dd$estimate, dd$p.value)
    }
  })
  corResult <- as.data.frame(t(corResult))
  colnames(corResult) <- c("correlation", "pvalue")
  corResult$p.adjust <- p.adjust(corResult$pvalue, method = "BH")
  corResult$gene <- rownames(corResult)
  return(corResult)
}

## 火山图
plot_volcano <- function(corData, genesToshow, title) {
  df <- corData
  df$Significance <- "NS"
  df$Significance[df$p.adjust < 0.05 & df$correlation > 0] <- "Up"
  df$Significance[df$p.adjust < 0.05 & df$correlation < 0] <- "Down"

  p <- ggplot(df, aes(x = correlation, y = -log10(pvalue), color = Significance)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", panel.grid.minor = element_blank()) +
    labs(x = "Pearson correlation", y = expression(-log[10](pvalue)),
         color = "", title = title) +
    geom_text_repel(
      data = subset(df, gene %in% genesToshow),
      aes(label = gene),
      size = 4, box.padding = 0.4, max.overlaps = Inf, color = "black"
    )
  return(p)
}

#############################################################
## 模块 1：CRISPR Gene Effect × AUC
#############################################################

cat("\n========== 模块 1: CRISPR Gene Effect ==========\n")

## 读取 Gene Effect
cat("正在读取 CRISPR Gene Effect...\n")
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
rownames(geneEffect) <- geneEffect[, 1]
geneEffect <- geneEffect[, -1]
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))

## 加载 AUC
aucVec <- load_AUC(auc_file)

## 取交集
coID <- intersect(rownames(geneEffect), names(aucVec))
cat("  共有细胞系:", length(coID), "\n")
geneEffect_sub <- geneEffect[coID, ]
auc_sub <- aucVec[coID]

## 批量相关性
corData1 <- batch_cor_auc(geneEffect_sub, auc_sub, "CRISPR GeneEffect")
saveRDS(corData1, file.path(output_dir, "18_corData_CRISPR_GeneEffect.rds"))

## 火山图（PRMT5 和 WDR77 是 PRMT5 复合体的成员）
p1 <- plot_volcano(corData1, c("PRMT5", "WDR77"),
                   "CRISPR Gene Effect vs AMG193 AUC")
print(p1)
ggsave(file.path(output_dir, "18_volcano_CRISPR_GeneEffect.png"),
       p1, width = 7, height = 5, dpi = 150)
cat("  已保存: 18_volcano_CRISPR_GeneEffect.png\n")

rm(geneEffect, geneEffect_sub)

#############################################################
## 模块 2：RNAi DEMETER2 × AUC
#############################################################

cat("\n========== 模块 2: RNAi DEMETER2 ==========\n")

## 读取细胞系信息（用于 CCLEName → ModelID 映射）
cellInfor <- readRDS(file.path(output_dir, "cellinfor.rds"))
cellMap <- cellInfor[, c("ModelID", "CCLEName")]

## 读取 RNAi DEMETER2 数据
## 原始格式：行 = 基因，列 = CCLE_Name（如 "143B_BONE"）
cat("正在读取 RNAi DEMETER2...\n")
geneEffect_rnai <- data.table::fread("data/DEMETER2_Data_v6/D2_Achilles_gene_dep_scores.csv",
                                      data.table = F)
geneEffect_rnai <- geneEffect_rnai %>%
  column_to_rownames("V1") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("CCLEName") %>%
  inner_join(cellMap, ., by = "CCLEName") %>%
  select(-CCLEName) %>%
  column_to_rownames("ModelID")

## 清理列名：去掉含 "&" 的列（旧格式残留），去掉 (EntrezID) 后缀
geneEffect_rnai <- geneEffect_rnai[, !grepl("&", colnames(geneEffect_rnai))]
colnames(geneEffect_rnai) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect_rnai))
geneEffect_rnai[is.na(geneEffect_rnai)] <- 0
cat("  RNAi 矩阵:", nrow(geneEffect_rnai), "模型 ×", ncol(geneEffect_rnai), "基因\n")

## 取交集
coID <- intersect(rownames(geneEffect_rnai), names(aucVec))
cat("  共有细胞系:", length(coID), "\n")
geneEffect_rnai_sub <- geneEffect_rnai[coID, ]
auc_sub <- aucVec[coID]

## 批量相关性
corData2 <- batch_cor_auc(geneEffect_rnai_sub, auc_sub, "RNAi DEMETER2")
saveRDS(corData2, file.path(output_dir, "18_corData_RNAi_DEMETER2.rds"))

## 火山图
p2 <- plot_volcano(corData2, c("PRMT5", "WDR77"),
                   "RNAi DEMETER2 vs AMG193 AUC")
print(p2)
ggsave(file.path(output_dir, "18_volcano_RNAi_DEMETER2.png"),
       p2, width = 7, height = 5, dpi = 150)
cat("  已保存: 18_volcano_RNAi_DEMETER2.png\n")

rm(geneEffect_rnai, geneEffect_rnai_sub)

#############################################################
## 模块 3：拷贝数变异 (CNV) × AUC
#############################################################

cat("\n========== 模块 3: 拷贝数变异 (CNV) ==========\n")

## 新版 CNV 用 ModelConditionID（MC-编号），需映射为 ACH-编号
cat("正在读取 CNV 数据（MC→ACH 映射）...\n")
cnData_raw <- data.table::fread("data/OmicsCNGeneMC_WES.csv", data.table = F)
cnData_raw <- cnData_raw[cnData_raw$IsDefaultEntryForMC == "Yes", ]

mcMap <- data.table::fread("data/ModelCondition.csv", data.table = F,
                           select = c("ModelConditionID", "ModelID"))
cnData_raw$ModelID <- mcMap$ModelID[match(cnData_raw$ModelConditionID, mcMap$ModelConditionID)]
cnData_raw <- cnData_raw[!is.na(cnData_raw$ModelID), ]
cnData_raw <- cnData_raw[!duplicated(cnData_raw$ModelID), ]
rownames(cnData_raw) <- cnData_raw$ModelID

## 去掉元数据列，只保留基因列
meta_cols <- c("ModelConditionID", "IsDefaultEntryForMC", "ModelID")
cnData <- cnData_raw[, !(colnames(cnData_raw) %in% meta_cols)]
colnames(cnData) <- gsub("\\s+\\(\\d+\\)", "", colnames(cnData))
cat("  CNV 矩阵:", nrow(cnData), "模型 ×", ncol(cnData), "基因\n")
rm(cnData_raw, mcMap)

## 取交集
coID <- intersect(rownames(cnData), names(aucVec))
cat("  共有细胞系:", length(coID), "\n")
cnData_sub <- cnData[coID, ]
auc_sub <- aucVec[coID]

## 批量相关性
corData3 <- batch_cor_auc(cnData_sub, auc_sub, "CNV")
saveRDS(corData3, file.path(output_dir, "18_corData_CNV.rds"))

## 火山图（MTAP/CDKN2A/CDKN2B 位于 9p21.3，是 AMG193 敏感性的关键标志物）
p3 <- plot_volcano(corData3, c("MTAP", "CDKN2A", "CDKN2B"),
                   "Copy Number vs AMG193 AUC")
print(p3)
ggsave(file.path(output_dir, "18_volcano_CNV.png"),
       p3, width = 7, height = 5, dpi = 150)
cat("  已保存: 18_volcano_CNV.png\n")

rm(cnData, cnData_sub)

#############################################################
## 模块 4：基因表达量 × AUC
#############################################################

cat("\n========== 模块 4: 基因表达量 ==========\n")

## 新版表达量含元数据列，需过滤
cat("正在读取表达量数据...\n")
expr_raw <- data.table::fread("data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                               data.table = F)
expr_raw <- expr_raw[expr_raw$IsDefaultEntryForMC == "Yes", ]
expr_raw <- expr_raw[!duplicated(expr_raw$ModelID), ]
rownames(expr_raw) <- expr_raw$ModelID

## 去掉元数据列
meta_cols <- c("V1", "SequencingID", "ModelConditionID", "ModelID",
                "IsDefaultEntryForMC", "IsDefaultEntryForModel")
exprData <- expr_raw[, !(colnames(expr_raw) %in% meta_cols)]
rm(expr_raw)
colnames(exprData) <- gsub("\\s+\\(\\d+\\)", "", colnames(exprData))
cat("  表达量矩阵:", nrow(exprData), "模型 ×", ncol(exprData), "基因\n")

## 取交集
coID <- intersect(rownames(exprData), names(aucVec))
cat("  共有细胞系:", length(coID), "\n")
exprData_sub <- exprData[coID, ]
auc_sub <- aucVec[coID]

## 批量相关性
corData4 <- batch_cor_auc(exprData_sub, auc_sub, "Expression")
saveRDS(corData4, file.path(output_dir, "18_corData_Expression.rds"))

## 火山图
p4 <- plot_volcano(corData4, c("MTAP", "CDKN2A", "CDKN2B"),
                   "Expression vs AMG193 AUC")
print(p4)
ggsave(file.path(output_dir, "18_volcano_Expression.png"),
       p4, width = 7, height = 5, dpi = 150)
cat("  已保存: 18_volcano_Expression.png\n")

cat("\n========== 18 脚本全部完成 ==========\n")
