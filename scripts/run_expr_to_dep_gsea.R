################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 本脚本由 05_from_gene_to_dependency.R 泛化而来：
###   任意目标基因的表达量 vs 全基因依赖性相关性 + Hallmark GSEA
###   用法：Rscript run_expr_to_dep_gsea.R PDE4B
###
### ================================================================
### GSEA 原理与方向解读（为什么取负号、怎么读上下富集）
###
### 1. 富集对象不是 Top 基因，而是完整排序列表：
###    输入是全部通过方差过滤的基因（每个基因携带一个得分
###    = 它与目标基因表达量的相关系数取负），排成一条长队列。
###    GSEA 把每条 Hallmark 通路的成员基因"撒"回队列里，
###    看它们扎堆在哪一段——评估的是通路【整体】在目标基因
###    表达梯度上的偏移，而非个别基因的高低。
###    （与 ORA 卡阈值取清单的做法相比，GSEA 无阈值、不丢弃
###    未过线基因的信息；弱但一致的整通路偏移也能检出。）
###
### 2. 取负号 = 翻转数轴，把"依赖更强"翻译成 GSEA 的"上调方向"：
###    Effect 生物学方向是"负 = 死 = 重要"，相关系数继承了这套方向。
###    GSEA 惯例是"排在队列前面 = 活性升高"。
###    取负后：cor 最负（依赖更强/敏化）的基因被推到队列【顶部】。
###
### 3. 方向速查表（勿搞反，换脚本/换基因时先确认有无取负）：
###    队列顶部  ← 原始 cor 为【负】→ 目标基因高表达时该基因
###               依赖更强（敏化方向，如 UQCRH vs PDE4B r=-0.91）
###    队列底部  ← 原始 cor 为【正】→ 目标基因高表达时该基因
###               依赖减弱（被保护方向，如 SAA1 vs PDE4B r=+0.69）
###    NES > 0（上富集）→ 通路基因扎堆顶部 → 目标基因高表达时
###               该通路【依赖更强】（敏化）
###    NES < 0（下富集）→ 通路基因扎堆底部 → 目标基因高表达时
###               该通路【依赖被削弱】（被保护，同 IL-6 保护
###               chromatin 基因的场景）
###    气泡图（split=".sign"）左侧面板 = NES>0 = 敏化；
###               右侧面板 = NES<0 = 保护。
###
### 4. 若不取负号直接排序，结果完全对称（NES 全部反号、
###    上下/左右面板互换），取负纯粹是为了让"左侧 = 敏化"
###    与直觉（左/上 = 激活）对齐。
### ================================================================

args <- commandArgs(trailingOnly = TRUE)
target_gene <- if (length(args) >= 1) args[1] else "PDE4B"
out_dir <- "TM00/DepMap_TM00/output"
gmt_file <- "TM00/DepMap_TM00/resource/geneSets/h.all.v2023.2.Hs.symbols.gmt"

cat("目标基因:", target_gene, "\n")

## ---- 数据加载 ----
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

if (!target_gene %in% colnames(exprSet)) stop("表达矩阵中找不到 ", target_gene)

## 行对齐防御：cor(x, y) 按位置配对，行序错乱不会报错只会全错
## （当前 01 脚本的 RDS 已验证对齐，此为换数据源时的保险）
if (!identical(rownames(exprSet), rownames(geneEffect))) {
  common_idx <- intersect(rownames(exprSet), rownames(geneEffect))
  exprSet <- exprSet[common_idx, , drop = FALSE]
  geneEffect <- geneEffect[common_idx, , drop = FALSE]
  cat("行未对齐，已按交集重排：", length(common_idx), "个细胞系\n")
}

## 表达方差检查（低方差 → 相关结果不可信）
expr_var <- var(exprSet[, target_gene], na.rm = TRUE)
cat(target_gene, "表达方差:", round(expr_var, 4), "\n")

## ---- 批量相关：PDE4B 表达 vs 全基因 Effect ----
## 用向量化 cor() 代替逐列 cor.test（快几个量级）
genes <- colnames(geneEffect)
x <- exprSet[, target_gene]
mat <- as.matrix(geneEffect)

## 逐列计算（含 NA 处理）
cor_vec <- apply(mat, 2, function(y) cor(x, y, use = "pairwise.complete.obs"))
corData <- data.frame(
  exp = target_gene,
  Dependency = genes,
  cor = cor_vec,
  stringsAsFactors = FALSE
)
## p 值（n = 配对完整样本数）
n_ok <- sum(!is.na(x))
corData$pvalue <- 2 * pt(-abs(corData$cor) * sqrt(n_ok - 2) / sqrt(1 - corData$cor^2), df = n_ok - 2)

cat("\n== Top 15 负相关（", target_gene, "高表达时依赖更强，敏化方向）==\n")
print(head(corData[order(corData$cor), ], 15))
cat("\n== Top 15 正相关（", target_gene, "高表达时依赖减弱，保护方向）==\n")
print(head(corData[order(-corData$cor), ], 15))

## 方差质检：Top 基因的 Effect 方差（低方差相关是假象）
top_bind <- rbind(head(corData[order(corData$cor), ], 15),
                  head(corData[order(-corData$cor), ], 15))
eff_var <- apply(mat[, top_bind$Dependency, drop = FALSE], 2, var, na.rm = TRUE)
cat("\n== Top 基因的 Effect 方差（接近 0 的为低方差假象）==\n")
print(round(sort(eff_var), 4))

## ---- 低方差基因过滤（进入 GSEA 前必做）----
## Effect 方差过小的基因（平线+少数离群点）与任何向量都能冒出 |r|>0.9 的假相关，
## 若不过滤会占据排序列表两端、扭曲 GSEA 的 NES 和 leading edge。
## 阈值 0.05：实测 OR5K2=0.023、NOMO3=0.003 被剔；真信号如 UQCRH=0.89 保留
min_eff_var <- 0.05
all_eff_var <- apply(mat, 2, var, na.rm = TRUE)
keep_genes <- names(all_eff_var)[all_eff_var >= min_eff_var]
cat("\nEffect 方差过滤：", ncol(geneEffect), "→", length(keep_genes),
    "个基因保留（方差阈值", min_eff_var, "）\n")
corData_gsea <- corData[corData$Dependency %in% keep_genes, ]
corData_gsea <- corData_gsea[!is.na(corData_gsea$cor), ]  ## 剔除 cor 为 NA 的基因

## ---- GSEA（Hallmark）----
library(clusterProfiler)
geneSet <- read.gmt(gmt_file)

mygeneList <- -corData_gsea$cor   ## 取负：cor 最负（敏化）→ 排序列表顶部（见头部注释"方向速查表"）
names(mygeneList) <- corData_gsea$Dependency
mygeneList <- sort(mygeneList, decreasing = TRUE)
mygeneList <- mygeneList[!is.na(mygeneList)]  ## 双保险：剔除残余 NA

mygsea <- GSEA(geneList = mygeneList, TERM2GENE = geneSet, verbose = FALSE)
res <- as.data.frame(mygsea)

cat("\n== GSEA 结果（按 NES 排序）==\n")
## 提醒：NES>0 = 敏化方向通路；NES<0 = 保护方向通路（详见头部注释）
print(res[, c("ID", "NES", "pvalue", "p.adjust", "setSize")], row.names = FALSE)

## ---- 落盘 ----
saveRDS(corData, file.path(out_dir, paste0("exprToDep_corData_", target_gene, ".rds")))
saveRDS(res, file.path(out_dir, paste0("exprToDep_gsea_", target_gene, ".rds")))
write.csv(res[, c("ID", "NES", "pvalue", "p.adjust", "setSize")],
          file.path(out_dir, paste0("exprToDep_gsea_", target_gene, ".csv")),
          row.names = FALSE)
cat("\n结果已保存:", file.path(out_dir, paste0("exprToDep_gsea_", target_gene, ".csv")), "\n")
