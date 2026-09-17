################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 05.3 通路活性 vs PRISM 药物敏感性
### ================================================================
### 思路：用 PROGENy 11 条通路活性分，预测细胞系对 PRISM 药物库的敏感性
###       AUC 越低 = 越敏感（药物效果越好）
### 生物学意义：哪些通路活性能预测药物响应？
### ================================================================

rm(list = ls())

## =================== 读取数据 ===================

## 读取 PROGENy 通路活性分
progeny_scores <- readRDS(file = "TM00/DepMap_TM00/output/progeny_pathway_scores.rds")
progeny_df <- as.data.frame(progeny_scores)

## 读取 PRISM AUC 药物敏感性矩阵（行：细胞系 ACH-，列：药物 DPC-）
cat("读取 PRISM 药物敏感性数据...\n")
prism_auc <- data.table::fread("data/PRISM_Repurposing_AUC_Matrix.csv", data.table = F)
rownames(prism_auc) <- prism_auc[, 1]
prism_auc <- prism_auc[, -1]

## 读取药物注释信息
prism_compounds <- data.table::fread("data/PRISM_Repurposing_Compound_Conditions.csv", data.table = F)

cat("PRISM 药物数：", ncol(prism_auc), "\n")
cat("PROGENy 通路数：", ncol(progeny_df), "\n")

## =================== 对齐细胞系 ===================

common_cells <- intersect(rownames(progeny_df), rownames(prism_auc))
progeny_aligned <- progeny_df[common_cells, ]
prism_aligned <- prism_auc[common_cells, ]

cat("共有的细胞系数：", length(common_cells), "\n")

## =================== 通路活性 vs 药物敏感性批量相关性 ===================

## 对每条通路，计算它与所有药物的 AUC 相关性
pathways <- colnames(progeny_aligned)
ndrugs <- ncol(prism_aligned)

## 存储所有结果
all_results <- data.frame()

cat("\n批量计算通路活性 vs 药物 AUC 相关性...\n")
pb <- txtProgressBar(min = 0, max = length(pathways), style = 3)

for (p in seq_along(pathways)) {

  pathway_name <- pathways[p]
  pathway_scores <- as.numeric(progeny_aligned[, pathway_name])

  ## 逐药物计算相关性
  for (d in 1:ndrugs) {

    drug_auc <- as.numeric(prism_aligned[, d])
    drug_name <- colnames(prism_aligned)[d]

    ## 去掉 NA
    valid <- !is.na(pathway_scores) & !is.na(drug_auc)
    if (sum(valid) < 50) next  ## 至少需要 50 个有效样本

    dd <- cor.test(pathway_scores[valid], drug_auc[valid])

    all_results <- rbind(all_results, data.frame(
      Pathway = pathway_name,
      Drug    = drug_name,
      cor     = dd$estimate,
      pvalue  = dd$p.value,
      n       = sum(valid)
    ))
  }

  setTxtProgressBar(pb, p)
}
close(pb)

cat("\n总结果行数：", nrow(all_results), "\n")

## =================== FDR 校正 ===================

all_results$FDR <- p.adjust(all_results$pvalue, method = "BH")

## 筛选显著的通路-药物关联
sig_results <- all_results[all_results$FDR < 0.05 &
                             abs(all_results$cor) > 0.3, ]

cat("显著关联（FDR < 0.05, |cor| > 0.3）：", nrow(sig_results), " 对\n")

## =================== 药物名称映射 ===================

## 将 DPC- 编号映射为药物名称
dpc_map <- setNames(prism_compounds$CompoundName, prism_compounds$CompoundID)
sig_results$DrugName <- dpc_map[sig_results$Drug]

## 按 |cor| 排序
sig_results <- sig_results[order(-abs(sig_results$cor)), ]

cat("\nTop 20 通路-药物关联：\n")
print(head(sig_results[, c("Pathway", "DrugName", "cor", "FDR")], 20))

## =================== 可视化：热图 ===================

## 构建通路 × 药物 的相关系数矩阵（只取 Top 50 最显著药物）
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

top_drugs <- unique(head(sig_results$Drug, 50))
top_results <- sig_results[sig_results$Drug %in% top_drugs, ]

## 构建矩阵
heat_matrix <- reshape2::dcast(top_results, Pathway ~ Drug, value.var = "cor")
rownames(heat_matrix) <- heat_matrix[, 1]
heat_matrix <- heat_matrix[, -1]
heat_matrix <- as.matrix(heat_matrix)

## 药物名称映射：找不到名称的保留原始 DPC 编号
drug_colnames <- colnames(heat_matrix)
mapped_names <- dpc_map[drug_colnames]
mapped_names[is.na(mapped_names)] <- drug_colnames[is.na(mapped_names)]
colnames(heat_matrix) <- mapped_names

## 清理 NA / NaN / Inf，全部替换为 0，避免聚类报错
heat_matrix[!is.finite(heat_matrix)] <- 0

## 去掉全为 0 的空行/空列
heat_matrix <- heat_matrix[rowSums(abs(heat_matrix)) > 0, , drop = FALSE]
heat_matrix <- heat_matrix[, colSums(abs(heat_matrix)) > 0, drop = FALSE]

## 只有行数 >= 2 且列数 >= 2 时才能聚类
can_cluster <- nrow(heat_matrix) >= 2 && ncol(heat_matrix) >= 2

## 绘制热图
pheatmap(heat_matrix,
         cluster_rows = can_cluster,
         cluster_cols = can_cluster,
         main = "通路活性 vs 药物敏感性（Pearson cor）",
         fontsize_col = 6,
         angle_col = 45)

## =================== 保存结果 ===================

saveRDS(all_results,
        file = "TM00/DepMap_TM00/output/pathway_drug_sensitivity.rds")

write.csv(sig_results,
          file = "TM00/DepMap_TM00/output/pathway_drug_significant.csv",
          row.names = FALSE)

cat("\n结果已保存到 output/pathway_drug_significant.csv\n")
