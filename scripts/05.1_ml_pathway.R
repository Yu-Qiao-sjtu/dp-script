################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 05.1 基于通路活性的机器学习预测模型
### ================================================================
### 思路：用 PROGENy 11 条通路活性分（而非单基因表达）作为特征，
###       预测细胞系对 ESR1（或任意目标基因）的依赖性
### 优势：通路活性比单基因更稳健，泛化能力更强
### ================================================================

rm(list = ls())

## =================== 读取数据 ===================

## 读取 PROGENy 通路活性分（05 脚本第三部分的输出）
progeny_scores <- readRDS(file = "TM00/DepMap_TM00/output/progeny_pathway_scores.rds")
progeny_df <- as.data.frame(progeny_scores)

## 读取基因效应矩阵
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## 设置目标基因
target_gene <- "ESR1"

## 对齐细胞系（取 PROGENy 和 geneEffect 的交集）
common_cells <- intersect(rownames(progeny_df), rownames(geneEffect))
X <- progeny_df[common_cells, ]        ## 特征：11 条通路活性分
y <- geneEffect[common_cells, target_gene]  ## 响应：目标基因依赖性

## 去掉 NA
valid <- complete.cases(X) & !is.na(y)
X <- X[valid, ]
y <- y[valid]

cat("样本数：", nrow(X), "，特征数：", ncol(X), "\n")

## =================== LASSO 回归（特征选择） ===================

## ---- 自动安装缺失依赖 ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in c("glmnet", "randomForest", "pROC", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(glmnet)

## 转为矩阵
X_mat <- as.matrix(X)

## LASSO 10 折交叉验证
set.seed(42)
cv_fit <- cv.glmnet(X_mat, y, alpha = 1, nfolds = 10)

## 提取非零系数的通路（被 LASSO 保留的预测特征）
lasso_coefs <- coef(cv_fit, s = "lambda.min")
selected_pathways <- rownames(lasso_coefs)[which(lasso_coefs[, 1] != 0)]
selected_pathways <- setdiff(selected_pathways, "(Intercept)")
cat("\nLASSO 选出的通路（", length(selected_pathways), " 条）：\n")
print(selected_pathways)

## =================== 随机森林（预测评估） ===================

library(randomForest)

## 用全部 11 条通路训练随机森林
set.seed(42)
rf_model <- randomForest(X_mat, y, ntree = 500, importance = TRUE)

## 查看模型性能（% Var explained 越高越好）
print(rf_model)

## 通路重要性排序
rf_imp <- importance(rf_model)
rf_imp <- rf_imp[order(-rf_imp[, "%IncMSE"]), ]
cat("\n通路重要性排序（%IncMSE）：\n")
print(round(rf_imp, 4))

## =================== 可视化 ===================

library(ggplot2)

## LASSO 系数路径图
plot(cv_fit)

## 随机森林重要性条形图
imp_df <- data.frame(
  Pathway = rownames(rf_imp),
  IncMSE = rf_imp[, "%IncMSE"]
)
p <- ggplot(imp_df, aes(x = reorder(Pathway, IncMSE), y = IncMSE)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = paste("PROGENy 通路对", target_gene, "依赖性的预测重要性"),
       x = "通路", y = "%IncMSE") +
  theme_bw()
print(p)

## =================== 保存结果 ===================
saveRDS(list(
  cv_fit       = cv_fit,
  selected     = selected_pathways,
  rf_model     = rf_model,
  rf_importance = rf_imp
), file = "TM00/DepMap_TM00/output/ml_pathway_model.rds")

