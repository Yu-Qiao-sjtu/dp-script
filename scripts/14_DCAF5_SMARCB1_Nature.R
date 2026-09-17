################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（21Q1 → 26Q1）

### ================================================================
### 主题：特定癌症（Rhabdoid/ATRT）特定基因（SMARCB1）突变下的易感基因筛选
### 灵感来源 Nature 文章:
###   Targeting DCAF5 suppresses SMARCB1-mutant cancer by stabilizing SWI/SNF
### 支线任务：学习如何把别人的代码用起来
###
### 分析思路：
###   SMARCB1 突变（缺失）的 Rhabdoid 肿瘤 → 对哪些基因更加依赖？
###   核心发现：DCAF5 是 SMARCB1 突变肿瘤的特异性依赖基因
###             → 抑制 DCAF5 可以稳定残余的 SWI/SNF 复合体
###
### 两个部分：
###   第一部分：论文原始代码（source 论文的 run_lm_stats_limma 函数）
###   第二部分：自己用 limma 手写差异分析（更直观易懂）
###
### 新版数据中识别 Rhabdoid 细胞系：
###   旧版需要 200 行重分类代码，现在直接用 OncotreeSubtype 匹配
###   ATRT（CNS/Brain，6株）+ MRT（Kidney，9株）+ Extrarenal Rhabdoid（7株）= 22 株
### ================================================================

############################################################
## =================== 第一部分：论文原始代码流程 ===================
## 使用论文 GitHub 仓库的 run_lm_stats_limma 函数做差异分析
############################################################
rm(list = ls())

## ---- 1. 加载 R 包和论文的辅助函数 ----
library(dplyr)
library(stringr)

## 从论文 GitHub 仓库读取 run_lm_stats_limma 函数
## 该函数封装了 limma 的线性模型 + 贝叶斯检验 + 单侧 p 值转换
source("TM00/DepMap_TM00/resource/SWISNF.DCAF5.Dependency-main/r_scripts/depmap/functions.R")

## ---- 2. 读取细胞系信息 ----
## 新版 Model.csv 字段名与旧版 sample_info.csv 不同：
##   DepMap_ID → ModelID
##   CCLE_Name → CCLEName
##   primary_disease → OncotreePrimaryDisease
##   lineage → OncotreeLineage
##   lineage_subtype → OncotreeSubtype
##   age → Age
sample_info <- data.table::fread(file = "data/Model.csv", data.table = F)

## 构建细胞系信息表（新版字段映射到旧版变量名，保持与论文代码兼容）
mf <- sample_info %>%
  filter(!is.na(ModelID) & ModelID != '') %>%
  mutate(EngineeredModel = as.character(EngineeredModel)) %>%
  ## 过滤工程化细胞系和非肿瘤细胞系
  filter(EngineeredModel != "True") %>%
  filter(!grepl('MATCHED_NORMAL', CCLEName)) %>%
  dplyr::select(DepMap_ID = ModelID,
                CCLE_name = CCLEName,
                Type = OncotreePrimaryDisease,
                Age) %>%
  distinct()

## ---- 3. 直接识别 Rhabdoid 细胞系（替代旧版 200 行重分类代码）----
## 新版 OncotreeSubtype / OncotreeCode 已经精确标注了 Rhabdoid 肿瘤：
##   OncotreeCode == "ATRT" → Atypical Teratoid/Rhabdoid Tumor（CNS/Brain，6 株）
##   OncotreeCode == "MRT"  → Malignant Rhabdoid Tumor（Kidney，9 株）
##   OncotreeSubtype 含 "Rhabdoid" 且非 "Rhabdomyo" → Extrarenal Rhabdoid（7 株）
rhabdoid_info <- sample_info %>%
  filter(grepl("Rhabdoid", as.character(OncotreeSubtype))) %>%
  pull(ModelID)

mf <- mf %>%
  mutate(Type = ifelse(DepMap_ID %in% rhabdoid_info, "Rhabdoid", Type))

cat("Rhabdoid 细胞系数:", sum(mf$Type == "Rhabdoid"), "\n")
cat("总细胞系数:", nrow(mf), "\n\n")

## 保存处理后的细胞系信息
saveRDS(mf, file = "TM00/DepMap_TM00/output/cell_line_infor_trimmed.rds")

## ---- 4. 读取基因效应数据 ----
## 旧版 Achilles_gene_effect.csv → 新版 CRISPRGeneEffect.csv
gene_effect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
colnames(gene_effect) <- gsub("\\s+\\(\\d+\\)", "", colnames(gene_effect))

## 过滤到只含 mf 中的细胞系
gene_effect <- dplyr::filter(gene_effect, V1 %in% mf$DepMap_ID)
rownames_temp <- gene_effect$V1
gene_effect <- as.matrix(gene_effect[, -1])

## 将行名从 DepMap_ID 映射为 CCLE_name（论文代码使用 CCLE_name 做行名）
DepMap_to_CCLE <- mf$CCLE_name
names(DepMap_to_CCLE) <- mf$DepMap_ID
rownames(gene_effect) <- DepMap_to_CCLE[rownames_temp]

## ---- 5. 构建 Rhabdoid 分组向量 ----
## context_matrix: 1 = Rhabdoid, 0 = 其他
Rhabdoid_lines <- mf %>% dplyr::filter(Type == "Rhabdoid") %>% pull(CCLE_name)
cat("Rhabdoid CCLE names:", length(Rhabdoid_lines), "\n")

context_vector <- rep(0, nrow(gene_effect))
names(context_vector) <- rownames(gene_effect)
context_vector[rownames(gene_effect) %in% Rhabdoid_lines] <- 1
cat("Rhabdoid=1:", sum(context_vector == 1), " Rhabdoid=0:", sum(context_vector == 0), "\n\n")

## ---- 6. 差异分析（使用论文的 run_lm_stats_limma 函数）----
## 对每个基因做线性模型：gene_effect ~ Rhabdoid vs all others
ttest_Rhabdoid_v_all <- run_lm_stats_limma(gene_effect, context_vector)

## 基因名去掉后缀（Entrez ID）
ttest_Rhabdoid_v_all$Gene <- gsub("\\ .*", "", ttest_Rhabdoid_v_all$Gene)

cat("差异分析完成，基因数:", nrow(ttest_Rhabdoid_v_all), "\n")
cat("DCAF5 q.value:", ttest_Rhabdoid_v_all$q.value[ttest_Rhabdoid_v_all$Gene == "DCAF5"], "\n\n")

## 保存结果
saveRDS(ttest_Rhabdoid_v_all,
        file = "TM00/DepMap_TM00/output/DCAF5_Rhabdoid_vs_all_limma.rds")

## ---- 7. 火山图（论文 Figure 风格）----
library(ggplot2)
library(ggrepel)

## 标注的关键基因（论文发现的 DCAF5 + 对照基因）
highlight_genes <- c("DCAF5", "TP53", "NABP2", "PUM3",
                     "CUL4A", "DCAF15", "DCAF13")
plot_data <- ttest_Rhabdoid_v_all
plot_data$highlight <- plot_data$Gene %in% highlight_genes

p1 <- ggplot(plot_data, aes(x = EffectSize, y = -log10(q.value))) +
  geom_point(size = 2, pch = 21, stroke = 0.25, fill = "black") +
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        text = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 10, color = "black")) +
  xlab("Effect Size") + ylab("-log10(q-value)") +
  ## DCAF5 红色标注
  geom_point(data = subset(plot_data, Gene == "DCAF5"),
             aes(x = EffectSize, y = -log10(q.value)),
             fill = "red", size = 3, pch = 21, stroke = 0.25) +
  geom_text_repel(data = subset(plot_data, Gene == "DCAF5"),
                  aes(label = Gene), size = 3, nudge_x = 0.1, nudge_y = 0.1) +
  ## TP53 蓝色标注
  geom_point(data = subset(plot_data, Gene == "TP53"),
             aes(x = EffectSize, y = -log10(q.value)),
             fill = "#00aeef", size = 2, pch = 21, stroke = 0.25) +
  geom_text_repel(data = subset(plot_data, Gene == "TP53"),
                  aes(label = Gene), size = 2, nudge_x = -0.1, nudge_y = 0.1) +
  ## 其他 DCAF 家族成员
  geom_point(data = subset(plot_data, Gene %in% c("DCAF15", "DCAF13", "CUL4A")),
             aes(x = EffectSize, y = -log10(q.value)),
             fill = "magenta", size = 2, pch = 21, stroke = 0.25) +
  geom_text_repel(data = subset(plot_data, Gene %in% c("DCAF15", "DCAF13", "CUL4A")),
                  aes(label = Gene), size = 2) +
  ggtitle(paste0("26Q1: Rhabdoid (n=", sum(context_vector == 1),
                 ") v all others (n=", sum(context_vector == 0), ")"))

print(p1)
ggsave("TM00/DepMap_TM00/output/14_Rhabdoid_vs_all_volcano.png", p1, width = 8, height = 6, dpi = 300)
cat("火山图已保存\n\n")


############################################################
## =================== 第二部分：自己用 limma 手写 ===================
## 不依赖论文的 run_lm_stats_limma，直接用 limma 的四步法
## lmFit → eBayes → topTable（果子老师"四步法"）
############################################################
rm(list = ls())

library(limma)

## ---- 1. 读取数据 ----
## 细胞系信息（第一部分已保存）
cell_line <- readRDS(file = "TM00/DepMap_TM00/output/cell_line_infor_trimmed.rds")

## 基因效应数据
gene_effect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
rownames(gene_effect) <- gene_effect[, 1]
gene_effect <- gene_effect[, -1]
colnames(gene_effect) <- gsub("\\s+\\(\\d+\\)", "", colnames(gene_effect))

## 取交集细胞系
coID <- intersect(cell_line$DepMap_ID, rownames(gene_effect))

## 确认 Rhabdoid 细胞系
## 注意：论文直接提供了 14 个突变细胞系名单
## 实际分析时应该参照前面的教程，从突变矩阵确认 SMARCB1 突变状态
cell_line_common <- cell_line[cell_line$DepMap_ID %in% coID, ]
cellRT <- cell_line_common[cell_line_common$Type == "Rhabdoid", ]
cat("Rhabdoid 细胞系数:", nrow(cellRT), "\n")

cellRTname <- cellRT$DepMap_ID
cellNotRTname <- setdiff(cell_line_common$DepMap_ID, cellRTname)

## ---- 2. 调整矩阵：行=基因，列=细胞系 ----
gene_effect_t <- t(gene_effect)
## 列顺序：先 Rhabdoid，后其他（方便分组）
gene_effect_t <- gene_effect_t[, c(cellRTname, cellNotRTname)]

exprSet <- gene_effect_t

## ---- 3. limma 四步法 ----

## 步骤 1：创建分组向量
nRT <- length(cellRTname)
nOther <- length(cellNotRTname)
group <- c(rep("Mut", nRT), rep("Wt", nOther))
## levels 中对照组放前面（limma 惯例）
group <- factor(group, levels = c("Wt", "Mut"))

## 步骤 2：构建设计矩阵
design <- model.matrix(~group)
colnames(design) <- levels(group)

## 步骤 3：线性模型拟合 + 贝叶斯检验
fit <- lmFit(exprSet, design)
fit2 <- eBayes(fit)

## 步骤 4：提取差异分析结果
## coef=2 表示第二列（Mut）与第一列（Wt）的比较
allDiff <- topTable(fit2, adjust = 'BH', coef = 2, number = Inf)

## 重命名列（与论文格式一致）
colnames(allDiff)[c(1, 5)] <- c("EffectSize", "q.value")
allDiff$Gene <- rownames(allDiff)

cat("\nTop 10 差异基因:\n")
print(head(allDiff[, c("Gene", "EffectSize", "q.value")], 10))
cat("\nDCAF5 排名:\n")
dcaf5_row <- which(allDiff$Gene == "DCAF5")
cat("DCAF5 EffectSize:", allDiff$EffectSize[dcaf5_row], "\n")
cat("DCAF5 q.value:", allDiff$q.value[dcaf5_row], "\n")

## 保存结果
ttest_Rhabdoid_v_all <- allDiff
saveRDS(ttest_Rhabdoid_v_all,
        file = "TM00/DepMap_TM00/output/DCAF5_Rhabdoid_vs_all_limma_manual.rds")

cat("\n=== 脚本全部完成 ===\n")

### 重要提醒:
### 1. 不是每个突变都能进行这样的分析，突变多就分癌种，突变少就整体
### 2. 论文直接提供了突变细胞系名单，自己做时要参照 07-11 脚本从突变矩阵确认
### 3. 新版 OncotreeSubtype 已精确标注 Rhabdoid，无需旧版 200 行重分类代码
