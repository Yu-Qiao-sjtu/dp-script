################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 05.4 WGCNA 共表达网络
### ================================================================
### 思路：用 WGCNA 将 ~1.7 万个基因聚类为共表达模块，
###       然后计算各模块与目标基因（ESR1）依赖性的关联
### 生物学意义：发现协同调控的基因功能模块
### ================================================================

rm(list = ls())

## =================== 自动安装缺失依赖 ===================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("WGCNA", "flashClust")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "WGCNA") {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## =================== 读取数据 ===================

## 读取表达矩阵（行：细胞系，列：基因）
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")

## 读取基因效应矩阵（用于后续模块-性状关联）
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## WGCNA 要求行=样本（细胞系），列=基因
## exprSet 本身就是这个格式，无需转置
datExpr <- as.matrix(exprSet)
cat("原始维度：", nrow(datExpr), "细胞系 ×", ncol(datExpr), "基因\n")

## 去掉低表达基因：只保留方差最高的 5000 个基因
## WGCNA 官方推荐 5000-12000 个基因，TOM 计算是 O(n²)，基因数减半速度快 4 倍
gene_var <- apply(datExpr, 2, var)

## 去掉方差为 NA 或 0 的基因
gene_var <- gene_var[!is.na(gene_var) & gene_var > 0]

## 取方差 top 5000
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, length(gene_var))]
datExpr <- datExpr[, top_genes]
cat("方差 Top5000 筛选后：", nrow(datExpr), "细胞系 ×", ncol(datExpr), "基因\n")

## =================== 选择软阈值 ===================

cat("\n计算软阈值...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

## 可视化软阈值选择（保存为 PNG）
## 先清理所有残留的图形设备，避免冲突
while (!is.null(dev.list())) dev.off()
png("TM00/DepMap_TM00/output/wgcna_soft_threshold.png", width = 2400, height = 1200, res = 300, type = "cairo")
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.9, col = "red")
dev.off()
cat("软阈值图已保存：wgcna_soft_threshold.png\n")

## 自动选择软阈值（无标度拟合 R^2 > 0.9 的最小 power）
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 12  ## 备用默认值
cat("选择的软阈值 power =", softPower, "\n")

## =================== 一步法构建共表达网络 ===================

cat("\n构建共表达网络（blockwiseModules）...\n")
net <- blockwiseModules(
  datExpr,
  power            = softPower,
  TOMType          = "unsigned",
  maxBlockSize     = 20000,  ## 确保一次性处理所有基因，不分 block
  minModuleSize    = 30,
  mergeCutHeight   = 0.25,
  numericLabels    = TRUE,
  saveTOMs         = TRUE,
  saveTOMFileBase  = "TM00/DepMap_TM00/output/WGCNA_TOM",
  verbose          = 3
)

## 模块颜色分配
moduleColors <- labels2colors(net$colors)
table(moduleColors)

## 模块数量
nModules <- length(unique(net$colors)) - 1  ## 减去灰色（未分配）模块
cat("检测到的模块数（不含灰色）：", nModules, "\n")

## 可视化：聚类树 + 模块颜色（保存为 PNG）
while (!is.null(dev.list())) dev.off()
png("TM00/DepMap_TM00/output/wgcna_dendrogram.png", width = 3600, height = 1800, res = 300, type = "cairo")
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
cat("聚类树图已保存：wgcna_dendrogram.png\n")

## =================== 计算模块特征值（ME） ===================

## 自动检测 datExpr 与 net 是否匹配（交互逐段运行时 net 可能是旧变量）
## datExpr 行=样本（细胞系），列=基因；net$colors 长度应等于基因数（列数）
if (ncol(datExpr) != length(net$colors)) {
  cat("\n[自动修复] datExpr 基因数 (", ncol(datExpr), ") 与 net$colors (", length(net$colors),
      ") 维度不一致，自动重新构建网络...\n", sep = "")
  net <- blockwiseModules(
    datExpr,
    power            = softPower,
    TOMType          = "unsigned",
    maxBlockSize     = 20000,
    minModuleSize    = 30,
    mergeCutHeight   = 0.25,
    numericLabels    = TRUE,
    saveTOMs         = TRUE,
    saveTOMFileBase  = "TM00/DepMap_TM00/output/WGCNA_TOM",
    verbose          = 3
  )
  moduleColors <- labels2colors(net$colors)
  cat("[自动修复] 完成！net$colors 长度 =", length(net$colors), "\n")
} else {
  moduleColors <- labels2colors(net$colors)
}
cat("维度检查通过：datExpr 基因数 =", ncol(datExpr),
    ", net$colors 长度 =", length(net$colors), "\n")

## moduleEigengenes 要求行=样本（细胞系），列=基因，colors 长度=基因数
## datExpr 已经是这个格式，无需转置
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
## 行=细胞系，列=各模块的特征值

## =================== 模块与 ESR1 依赖性关联 ===================

## 取共有细胞系
common_cells <- intersect(rownames(MEs), rownames(geneEffect))
MEs_aligned <- MEs[common_cells, ]
esr1_dep <- geneEffect[common_cells, "ESR1"]

## 计算每个模块特征值与 ESR1 依赖性的相关性
moduleTraitCor <- cor(MEs_aligned, esr1_dep, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs_aligned))

## 整理结果
moduleTrait <- data.frame(
  Module    = colnames(MEs_aligned),
  cor       = moduleTraitCor[, 1],
  pvalue    = moduleTraitPvalue[, 1],
  FDR       = p.adjust(moduleTraitPvalue[, 1], method = "BH")
)
moduleTrait <- moduleTrait[order(-abs(moduleTrait$cor)), ]

cat("\n模块与 ESR1 依赖性关联：\n")
print(moduleTrait)

## =================== 可视化（保存为 PNG）===================

## 模块-性状热图
while (!is.null(dev.list())) dev.off()
png("TM00/DepMap_TM00/output/wgcna_module_trait_heatmap.png", width = 1500, height = 2400, res = 300, type = "cairo")
par(mar = c(6, 8, 3, 3))
labeledHeatmap(
  Matrix         = moduleTraitCor,
  xLabels        = "ESR1 Dependency",
  yLabels        = colnames(MEs_aligned),
  ySymbols       = colnames(MEs_aligned),
  colorLabels    = FALSE,
  colors         = blueWhiteRed(50),
  textMatrix     = paste(signif(moduleTraitCor, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = ""),
  setStdMargins  = FALSE,
  cex.text       = 0.5,
  zlim           = c(-1, 1),
  main           = "Module-Trait Relationships (ESR1)"
)
dev.off()
cat("模块-性状热图已保存：wgcna_module_trait_heatmap.png\n")

## =================== 提取关键模块的基因 ===================

## 找与 ESR1 依赖性关联最强的模块
top_module <- gsub("ME", "", moduleTrait$Module[1])
top_genes <- names(net$colors)[net$colors == as.numeric(top_module)]
cat("\n与 ESR1 依赖性关联最强的模块：", moduleTrait$Module[1],
    "（", top_module, "），基因数：", length(top_genes), "\n")

## =================== 保存结果 ===================

saveRDS(list(
  net           = net,
  moduleColors  = moduleColors,
  MEs           = MEs,
  moduleTrait   = moduleTrait,
  softPower     = softPower,
  top_module    = moduleTrait$Module[1],
  top_genes     = top_genes
), file = "TM00/DepMap_TM00/output/wgcna_result.rds")

cat("\n结果已保存到 output/wgcna_result.rds\n")
