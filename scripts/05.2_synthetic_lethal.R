################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 05.2 协同致死基因对筛选
### ================================================================
### 思路：在 CRISPR 基因效应相关矩阵中，找强负相关的基因对
###       即：敲掉 A 能活、敲掉 B 能活，但同时失去 A+B 就致死的对
### 生物学意义：药物联用靶点发现（如 PARP + BRCA1）
### ================================================================

rm(list = ls())

## =================== 读取共依赖性矩阵 ===================

## 读取 03 脚本计算的全基因两两相关矩阵
data <- readRDS(file = "TM00/DepMap_TM00/output/co_dependency_matrix.rds")
cat("相关矩阵维度：", nrow(data), "x", ncol(data), "\n")

## =================== 筛选协同致死候选对 ===================

## 定义筛选阈值
cor_threshold <- -0.3   ## 相关系数阈值（负相关）
pval_threshold <- 0.01   ## p 值阈值（需要额外计算）

## 提取下三角（避免重复：A-B 和 B-A 只保留一个）
gene_names <- rownames(data)
n <- nrow(data)

## 逐行扫描下三角，收集负相关基因对
synthetic_pairs <- data.frame()

cat("扫描相关矩阵下三角...\n")
pb <- txtProgressBar(min = 0, max = n - 1, style = 3)

for (i in 1:(n - 1)) {

  ## 取第 i 行中第 1~(i-1) 列的值（下三角）
  cor_vals <- data[i, 1:(i - 1)]

  ## 筛选低于阈值的（强负相关）
  hits <- which(cor_vals < cor_threshold)

  if (length(hits) > 0) {
    for (j in hits) {
      synthetic_pairs <- rbind(synthetic_pairs, data.frame(
        GeneA = gene_names[i],
        GeneB = gene_names[j],
        cor   = cor_vals[j]
      ))
    }
  }

  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n筛出协同致死候选对：", nrow(synthetic_pairs), " 对\n")

## 按 |cor| 排序（相关性越负越优先）
synthetic_pairs <- synthetic_pairs[order(synthetic_pairs$cor), ]

## 查看 Top20
cat("\nTop 20 协同致死候选对：\n")
print(head(synthetic_pairs, 20))

## =================== 生物学验证：已知合成致死对 ===================

## 检查经典的合成致死对是否被筛出
known_pairs <- list(
  c("PARP1", "BRCA1"),
  c("PARP1", "BRCA2"),
  c("ATM", "PARP1"),
  c("ARID1A", "ARID1B"),
  c("STAG1", "STAG2")
)

cat("\n已知合成致死对的验证：\n")
for (pair in known_pairs) {
  if (pair[1] %in% gene_names && pair[2] %in% gene_names) {
    cor_val <- data[pair[1], pair[2]]
    cat(sprintf("  %s - %s : cor = %.4f\n", pair[1], pair[2], cor_val))
  }
}

## =================== 指定基因的合成致死伙伴 ===================

## 输入你感兴趣的基因，找它的负相关伙伴
target_gene <- "ESR1"
if (target_gene %in% gene_names) {
  target_cors <- data[target_gene, ]
  target_neg <- sort(target_cors[target_cors < cor_threshold])
  cat(sprintf("\n%s 的协同致死候选伙伴（cor < %.1f）：\n", target_gene, cor_threshold))
  print(head(target_neg, 20))
}

## =================== 保存结果 ===================

saveRDS(synthetic_pairs,
        file = "TM00/DepMap_TM00/output/synthetic_lethal_pairs.rds")

write.csv(synthetic_pairs,
          file = "TM00/DepMap_TM00/output/synthetic_lethal_pairs.csv",
          row.names = FALSE)

cat("\n结果已保存到 output/synthetic_lethal_pairs.csv\n")
