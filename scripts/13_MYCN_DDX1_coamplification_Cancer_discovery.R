################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（21Q1 → 26Q1）

### ================================================================
### 主题：双基因共扩增下的易感基因筛选（In silico CRISPR Screen）
### 灵感来源：https://aacrjournals.org/cancerdiscovery/article/14/3/492/734910/
###
### 分析思路：
###   MYCN 扩增的神经母细胞瘤中，DDX1 常常共扩增（passenger co-amplification）
###   问题：MYCN+DDX1 共扩增 vs MYCN 单扩增，对哪些基因依赖性不同？
###   经典发现：共扩增组更依赖 RPTOR/MTOR（mTOR 通路）
###
### 与 12 的区别（单基因 → 双基因）：
###   12：CCNE1 扩增 vs 野生型
###   13：MYCN+DDX1 共扩增 vs MYCN 单扩增（对照组也要求 MYCN 扩增）
###        → 排除 MYCN 本身的影响，只看 DDX1 共扩增的"附带"效应
###
### 新版数据格式变化（同 12 脚本）：
###   Achilles_gene_effect.csv → CRISPRGeneEffect.csv（ACH-编号）
###   CCLE_gene_cn.csv → OmicsCNGeneMC_WES.csv（MC-编号 → 映射为 ACH-编号）
### ================================================================

############################################################
## =================== 第一部分：数据加载 + 批量计算 ===================
############################################################
rm(list = ls())

## ---- 1. 读取基因效应数据（GeneEffect，ACH-编号）----
## CRISPRGeneEffect.csv：行=细胞系(ACH-编号)，列=基因
## 值 = 基因效应得分（越负 = 越依赖），Less is More
depData <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)

## 列名去掉 (EntrezID) 后缀
colnames(depData) <- gsub("\\s+\\(\\d+\\)", "", colnames(depData))

## 第一列设为行名
rownames(depData) <- depData[, 1]
depData <- depData[, -1]

## ---- 2. 读取拷贝数数据（CNV，MC-编号 → 映射为 ACH-编号）----
## 处理方式同 12 脚本
cnData_raw <- data.table::fread("data/OmicsCNGeneMC_WES.csv", data.table = F)

## 过滤默认条目
cnData_raw <- cnData_raw[cnData_raw$IsDefaultEntryForMC == "Yes", ]

## 通过 ModelCondition.csv 将 MC-编号映射为 ACH-编号
mcMap <- data.table::fread("data/ModelCondition.csv", data.table = F,
                           select = c("ModelConditionID", "ModelID"))
cnData_raw$ModelID <- mcMap$ModelID[match(cnData_raw$ModelConditionID,
                                           mcMap$ModelConditionID)]

## 去除映射失败和重复
cnData_raw <- cnData_raw[!is.na(cnData_raw$ModelID), ]
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
index <- intersect(rownames(depData), rownames(cnData))
depData <- depData[index, ]
cnData <- cnData[index, ]

cat("交集细胞系数:", length(index), "\n\n")

## ---- 4. 双基因共扩增分组 ----
## 提取 MYCN 和 DDX1 的拷贝数
cnGeneData <- data.frame(DepMap_ID = rownames(cnData),
                         MYCN = cnData[, "MYCN"],
                         DDX1 = cnData[, "DDX1"])

## 阈值 cutoff = 2（log2 ratio）
## 比单基因分析的 1.58 更严格，因为要同时满足两个基因
cutoff <- 2

## 实验组：MYCN 扩增 + DDX1 共扩增
ampCells <- cnGeneData$DepMap_ID[cnGeneData$MYCN >= cutoff & cnGeneData$DDX1 >= cutoff]

## 对照组：MYCN 扩增 + DDX1 不共扩增
## 注意：对照组也要求 MYCN 扩增，排除 MYCN 本身的影响
wtCells <- cnGeneData$DepMap_ID[cnGeneData$MYCN >= cutoff & cnGeneData$DDX1 < cutoff]

cat("MYCN+DDX1 共扩增组:", length(ampCells), "\n")
cat("MYCN 单扩增组（对照）:", length(wtCells), "\n\n")

## 提取两组的基因效应数据
ampCellsDepData <- depData[ampCells, ]
wtCellsDepData <- depData[wtCells, ]

## ---- 5. 单基因验证（RPTOR）----
## 阳性对照：RPTOR 是已知的共扩增敏感基因
cat("=== 阳性对照验证：RPTOR ===\n")
cat("共扩增组 median:", median(ampCellsDepData[, "RPTOR"], na.rm = TRUE), "\n")
cat("单扩增组 median:", median(wtCellsDepData[, "RPTOR"], na.rm = TRUE), "\n")
cat("dep.FC:", median(ampCellsDepData[, "RPTOR"], na.rm = TRUE) -
                   median(wtCellsDepData[, "RPTOR"], na.rm = TRUE), "\n")
cat("Wilcoxon p:", wilcox.test(ampCellsDepData[, "RPTOR"],
                                wtCellsDepData[, "RPTOR"],
                                alternative = "less")$p.value, "\n\n")

## ---- 6. 批量遍历所有基因 ----
## 对每个基因计算：共扩增组 vs 单扩增组的差异 + 统计检验
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
  results[i, 2] <- median(ampD, na.rm = TRUE)
  results[i, 3] <- median(wtD, na.rm = TRUE)
  results[i, 4] <- median(ampD, na.rm = TRUE) - median(wtD, na.rm = TRUE)

  ## 5-7. mean 差异
  results[i, 5] <- mean(ampD, na.rm = TRUE)
  results[i, 6] <- mean(wtD, na.rm = TRUE)
  results[i, 7] <- mean(ampD, na.rm = TRUE) - mean(wtD, na.rm = TRUE)

  ## 8. Wilcoxon 秩和检验（单侧：共扩增组更依赖）
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

## 增加 FDR 校正
results$Wilcox_p_value_fdr <- p.adjust(results$Wilcox_p_value, method = "fdr")
results$Ttest_p_value_fdr  <- p.adjust(results$Ttest_p_value, method = "fdr")

## 保存结果
saveRDS(results, file = "TM00/DepMap_TM00/output/MYCN_DDX1_amp_Cancer_discovery_fig2_batchData.rds")
cat("\n批量计算完成，结果已保存\n")


############################################################
## =================== 第二部分：火山图复现 ===================
############################################################
rm(list = ls())

## 读取批量结果
data <- readRDS("TM00/DepMap_TM00/output/MYCN_DDX1_amp_Cancer_discovery_fig2_batchData.rds")

library(ggplot2)
library(ggrepel)

## 要标注的关键基因（来自 Cancer Discovery 论文的发现）
genes_to_show <- c("RPTOR", "MTOR")

## 火山图：dep.FC vs -log10(Wilcoxon p-value)
## dep.FC < 0 = 共扩增组更依赖该基因
p1 <- ggplot(data, aes(dep.FC, -log10(Wilcox_p_value))) +
  geom_point(shape = 21, colour = "black", fill = "lightyellow") +
  geom_point(data = subset(data, GeneSymbol %in% genes_to_show),
             shape = 21, colour = "black", fill = "skyblue3", size = 3) +
  geom_text_repel(data = subset(data, GeneSymbol %in% genes_to_show), aes(label = GeneSymbol)) +
  theme_bw() +
  xlab(expression(paste(Delta, "Median (", italic("DDX1-MYCN"),
                        " coamplification vs no coamplification)"))) +
  ylab(expression(paste("-log"[10], "(Wilcoxon ", italic("P"), " value)"))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

print(p1)
ggsave("TM00/DepMap_TM00/output/13_MYCN_DDX1_volcano.png", p1, width = 8, height = 6, dpi = 300)
cat("火山图已保存\n")


############################################################
## =================== 第三部分：GO 富集分析 ===================
## 对显著依赖的基因做 GO 功能富集
############################################################
rm(list = ls())

library(dplyr)
library(clusterProfiler)   ## 富集分析工具

## 读取批量结果
data <- readRDS("TM00/DepMap_TM00/output/MYCN_DDX1_amp_Cancer_discovery_fig2_batchData.rds")

## 筛选显著依赖的基因：
##   1. Wilcoxon p < 0.05（统计显著）
##   2. dep.FC < -0.1（效应量足够大：共扩增组依赖性更强）
##   3. 取 top 300（按 dep.FC 排序）
gene <- data %>%
  filter(Wilcox_p_value < 0.05) %>%
  filter(dep.FC < -0.1) %>%
  arrange(dep.FC) %>%
  slice(1:300)

cat("筛选到", nrow(gene), "个显著依赖基因\n")

## 基因 Symbol → Entrez ID（GO 分析需要）
gene_id <- bitr(gene$GeneSymbol,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

## ---- GO CC（细胞组分）富集 ----
ego_CC <- enrichGO(gene          = gene_id$ENTREZID,
                   OrgDb         = "org.Hs.eg.db",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

cat("\nGO CC 富集到", nrow(as.data.frame(ego_CC)), "个条目\n")

p_cc <- barplot(ego_CC, title = "GO Cellular Component")
print(p_cc)
ggsave("TM00/DepMap_TM00/output/13_MYCN_DDX1_GO_CC.png", p_cc, width = 8, height = 6, dpi = 300)

## ---- GO BP（生物学过程）富集 ----
ego_BP <- enrichGO(gene          = gene_id$ENTREZID,
                   OrgDb         = "org.Hs.eg.db",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

cat("GO BP 富集到", nrow(as.data.frame(ego_BP)), "个条目\n")

p_bp <- barplot(ego_BP, title = "GO Biological Process")
print(p_bp)
ggsave("TM00/DepMap_TM00/output/13_MYCN_DDX1_GO_BP.png", p_bp, width = 8, height = 6, dpi = 300)

cat("\n=== 脚本全部完成 ===\n")
