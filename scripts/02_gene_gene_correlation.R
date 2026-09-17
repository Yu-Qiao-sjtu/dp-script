################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（22Q2 → 26Q1）

### 复习批量操作
### 基因和基因的相关性

rm(list = ls())
library(ggplot2)

## 读取 01 脚本生成的表达量 RDS
## 数据结构：行 = 细胞系（ModelID，三套数据交集后的公共细胞系，约 1000 个），
##           列 = 基因（纯基因名，与 geneEffect 取过交集），
##           值 = log2(TPM + 1) 表达量（0 = 不表达，越大表达越高）
## readRDS 是 saveRDS 的逆操作，秒级恢复 01 脚本清洗好的矩阵，免重读 35GB 原始 CSV
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")

## ================================================================
## 本脚本核心概念：基因和基因的相关性（共表达分析）
##
## 相关性发生在【基因对】之间，细胞系是【样本/观测维度】：
##   基因 A 的表达向量 = A 列在 ~1000 个细胞系上的取值
##   基因 B 的表达向量 = B 列在 ~1000 个细胞系上的取值
##   cor(A 列, B 列) → A、B 的共表达相关系数
##
## 为什么需要这么多细胞系：只有在表达量有变化的群体里才能看出
## “A 高 B 也高”的共变模式；单一/少数样本无从谈相关。
##
## 生物学含义：高共表达通常暗示同通路/同蛋白复合物成员、
## 共同上游转录因子调控；注意相关 ≠ 直接相互作用。
## ================================================================

## 单个基因对的相关性验证：FOXA1 vs ESR1
## cor.test：Pearson 相关检验（t 检验框架），一次输出两个回答不同问题的统计量：
##   estimate（r）→ 共表达的【强度和方向】：-1 ~ 1，越接近 1 共表达越强
##   p.value      → 这个 r 的【可靠性】：排除“纯随机波动也能冒出这个 r”的可能。
##                  H0：真实相关系数 ρ=0（不共表达）；p 越小越难用随机解释
## p 值公式：t = r*sqrt(n-2)/sqrt(1-r^2)，服从 df=n-2 的 t 分布。
## 注意：n≈1000 的大样本下 r>0.062 就 p<0.05（万物皆显著），
## 实际筛选应主要看 |r| 大小（如 >0.3），p 值仅用于排除完全无信号者
cor.test(exprSet[, "FOXA1"], exprSet[, "ESR1"])

## 把检验结果存入 dd，分别提取两个关键统计量：
##   dd$estimate → 相关系数 r（-1 ~ 1，越接近 1 共表达越强）
##   dd$p.value  → 显著性（基于 ~1000 个细胞系的配对观测）
dd <- cor.test(exprSet[, "FOXA1"], exprSet[, "ESR1"])
dd$p.value
dd$estimate

## 批量计算 ESR1 与所有基因的相关性
## 思路：固定基因1 = ESR1，循环遍历全部 ~18000 个基因做 gene2，
## 每次都是“ESR1 列 vs 当前基因列”的 cor.test
gene1 <- "ESR1"
genedata <- exprSet[, gene1]     ## 提取 ESR1 的表达向量（长度 = 细胞系数）
corData <- data.frame()          ## 空表逐步填入结果

for (i in seq_len(ncol(exprSet))) {
  if (i %% 2000 == 0) cat("进度:", i, "/", ncol(exprSet), "\n")   ## 每 2000 个基因打印一次进度

  gene2 <- colnames(exprSet)[i]           ## 当前循环到的基因名
  dd <- cor.test(genedata, exprSet[, gene2])  ## ESR1 vs 当前基因的共表达检验

  corData[i, 1] <- gene1          ## 固定为 ESR1
  corData[i, 2] <- gene2          ## 当前基因
  corData[i, 3] <- dd$estimate    ## 相关系数 r
  corData[i, 4] <- dd$p.value     ## p 值
}

colnames(corData) <- c("Gene1", "Gene2", "cor", "pvalue")
## 下游建议筛选：|cor| > 0.3 且 p < 0.05（大样本下 p 普遍显著，
## 主要看 r 的绝对值；ESR1 的强共表达伙伴通常含 FOXA1/GATA3/XBP1 等雌激素受体通路成员）

## 保存结果
saveRDS(corData, file = "TM00/DepMap_TM00/output/02_ESR1_gene_correlation.rds")

### 推荐阅读: GZ07,批量技能
