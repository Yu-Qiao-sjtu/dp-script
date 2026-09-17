################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（22Q4/22Q2 → 26Q1）

### ================================================================
### 主题：不同肿瘤类型中特定基因（CDK4/CDK6）的依赖性分布
###
### 分析思路：
###   利用 DepMap CRISPR Gene Effect 数据，按照肿瘤亚型分组，
###   展示 CDK4 和 CDK6 在各癌种中的依赖程度（CERES score）。
###   乳腺癌额外区分 ER 阳性/阴性；前列腺癌区分原发/转移。
###
### 数据来源（26Q1）：
###   - 基因效应：CRISPRGeneEffect.csv（ACH-编号）
###   - 表达量：OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
###     （新版含元数据列，需过滤 IsDefaultEntryForMC）
###   - 细胞系信息：cellinfor.rds（07 脚本生成，基于 Model.csv）
### ================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

output_dir <- "TM00/DepMap_TM00/output"

#############################################################
## 第一部分：数据加载
#############################################################

## ---- 1. 读取 DepMap Gene Effect ----
cat("正在读取 Gene Effect...\n")
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
rownames(geneEffect) <- geneEffect[, 1]     ## 第一列 V1 = ACH-编号
geneEffect <- geneEffect[, -1]
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))
cat("  Gene Effect:", nrow(geneEffect), "模型 ×", ncol(geneEffect), "基因\n")

## ---- 2. 读取 CCLE 表达量（新版含元数据列，需过滤） ----
cat("正在读取表达量数据...\n")
exprSet_raw <- data.table::fread("data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", data.table = F)
## 新版前 6 列是元数据：V1, SequencingID, ModelConditionID, ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel
## 只保留默认条件 (IsDefaultEntryForMC == "Yes")
exprSet_raw <- exprSet_raw[exprSet_raw$IsDefaultEntryForMC == "Yes", ]
## 去掉可能的重复 ModelID（个别模型有多条默认记录）
exprSet_raw <- exprSet_raw[!duplicated(exprSet_raw$ModelID), ]
rownames(exprSet_raw) <- exprSet_raw$ModelID
## 去掉元数据列，只保留基因列
meta_cols <- c("V1", "SequencingID", "ModelConditionID", "ModelID",
                "IsDefaultEntryForMC", "IsDefaultEntryForModel")
exprSet <- exprSet_raw[, !(colnames(exprSet_raw) %in% meta_cols)]
rm(exprSet_raw)
colnames(exprSet) <- gsub("\\s+\\(\\d+\\)", "", colnames(exprSet))
cat("  表达量:", nrow(exprSet), "模型 ×", ncol(exprSet), "基因\n")

## ---- 3. 读取细胞系信息 ----
cellinfor <- readRDS(file = "TM00/DepMap_TM00/output/cellinfor.rds")
rownames(cellinfor) <- cellinfor$ModelID

#############################################################
## 第二部分：数据修剪 —— 共有细胞系和基因
#############################################################

## 三套数据共有的细胞系
commonindex <- intersect(rownames(geneEffect), rownames(exprSet))
commonindex <- intersect(commonindex, rownames(cellinfor))
cat("三套数据共有细胞系:", length(commonindex), "\n")

geneEffect <- geneEffect[commonindex, ]
exprSet    <- exprSet[commonindex, ]
cellinfor  <- cellinfor[commonindex, ]

## 共有基因
commonGenes <- intersect(colnames(geneEffect), colnames(exprSet))
cat("共有基因:", length(commonGenes), "\n")
geneEffect <- geneEffect[, commonGenes]
exprSet    <- exprSet[, commonGenes]

#############################################################
## 第三部分：肿瘤亚型分类
#############################################################

## ---- 乳腺癌：ER+ / ER- ----
## ModelSubtypeFeatures 包含 ER+/ER-/HER2+/TNBC/luminal 等信息
breast_idx <- cellinfor$OncotreeLineage == "Breast"
msf <- as.character(cellinfor$ModelSubtypeFeatures)
## ER+ 包含 "ER+" 或 "ER," (如 "luminal ER, PR+")
is_er_positive <- grepl("ER\\+", msf) | grepl("ER,", msf)
cellinfor$OncotreeSubtype[breast_idx] <-
  ifelse(is_er_positive[breast_idx],
         "ER_Positive_Breast",
         "ER_Negative_Breast")

## ---- 前列腺癌：原发 / 转移 ----
prostate_idx <- cellinfor$OncotreeLineage == "Prostate"
pom <- as.character(cellinfor$PrimaryOrMetastasis)
cellinfor$OncotreeSubtype[prostate_idx] <-
  ifelse(grepl("Metastatic", pom[prostate_idx]),
         "PRAD_Metastatic",
         "PRAD_Primary")

## ---- 筛选频率 >= 3 的亚型 ----
cellTypeFreq <- as.data.frame(table(cellinfor$OncotreeSubtype)) %>%
  filter(Freq >= 3)

## 保留的亚型顺序（按频率从高到低排列，画图时更美观）
cellType <- as.character(cellTypeFreq$Var1)
cat("频率 >= 3 的亚型数:", length(cellType), "\n")

## 过滤细胞系：只保留属于高频亚型的
cellinfor_new <- cellinfor[cellinfor$OncotreeSubtype %in% cellType, ]
cat("过滤后细胞系数:", nrow(cellinfor_new), "\n")

#############################################################
## 第四部分：提取目标基因的依赖性数据
#############################################################

## 构造函数：提取某基因在各肿瘤亚型中的 CERES 评分
getDepData <- function(gene) {
  myData <- data.frame(
    ModelID  = rownames(cellinfor_new),
    CERES    = geneEffect[rownames(cellinfor_new), gene],
    TumorType = cellinfor_new$OncotreeSubtype,
    stringsAsFactors = F
  )
  ## 按频率排序亚型 levels
  myData$TumorType <- factor(myData$TumorType, levels = cellType)
  myData$group <- gene
  return(myData)
}

## 提取 CDK4 和 CDK6
data1 <- getDepData("CDK4")
data2 <- getDepData("CDK6")
myData <- rbind(data1, data2)

#############################################################
## 第五部分：可视化 —— 按肿瘤亚型的依赖性箱线图
#############################################################

## 画图函数：CERES 评分 × 肿瘤类型
plot_gene <- function(data, gene, show_y = TRUE) {
  g <- ggplot(data[data$group == gene, ],
              aes(x = CERES, y = TumorType)) +
    geom_vline(xintercept = 0,  colour = "grey60",   linewidth = 1.3) +
    geom_vline(xintercept = -1, colour = "seagreen4", linewidth = 1.3) +
    geom_boxplot(coef = 1.5, outlier.shape = NA,
                 width = .7, fill = "orange", alpha = .10) +
    geom_jitter(width = 0, height = 0.25, size = .8, alpha = .12) +
    stat_summary(fun = median, geom = "point",
                 shape = 95, size = 4, colour = "red") +
    scale_x_reverse(breaks = c(0, -1, -2),
                    limits = c(0.25, min(-2, min(data$CERES, na.rm = TRUE))),
                    expand = expansion(add = c(.25, 0))) +
    scale_y_discrete(limits = rev(levels(data$TumorType))) +
    labs(x = paste0(gene, " CERES"), y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  ## 不需要 y 轴信息时隐藏
  if (!show_y) {
    g <- g + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  g
}

## 生成左右两图：CDK4（有 y 轴标签） | CDK6（无 y 轴标签）
p_cdk4 <- plot_gene(myData, "CDK4", show_y = TRUE)
p_cdk6 <- plot_gene(myData, "CDK6", show_y = FALSE)

## 用 patchwork 并排
library(patchwork)
p_combined <- (p_cdk4 | p_cdk6) + plot_layout(widths = c(1, 1))
print(p_combined)
ggsave(file.path(output_dir, "16_CDK4_CDK6_dependency_by_tumor_type.png"),
       p_combined, width = 12, height = 8, dpi = 150)
cat("已保存: 16_CDK4_CDK6_dependency_by_tumor_type.png\n")

## 保存数据供后续分析使用
saveRDS(myData, file.path(output_dir, "16_CDK4_CDK6_depData.rds"))

cat("\n========== 16 脚本完成 ==========\n")

### 拓展：
### 1. 一个是 gene effect，一个是 expression → 可以做两者的相关性
### 2. 两个都是 gene effect → 可以做基因间的共依赖分析
### 3. 究竟多少类，为什么要分 ER+，Prostate Metastasis → 关注亚型特异性依赖
