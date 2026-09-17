

################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 基因A的Dependency受到哪些基因调控
### 或者哪些基因的表达能够预测抑制剂疗效
### ================================================================
### 分析思路（从单变量筛选到多变量预测模型）：
###
### 第一步（本脚本）：单变量筛选
###   用 cor.test 逐个计算 ESR1 依赖性 vs 每个基因表达的相关性
###   目的：从 ~1.7 万个基因中初筛候选 biomarker
###   输出：corData 表（Gene1, Gene2, cor, pvalue）
###
### 第二步（进阶）：LASSO 特征选择 —— 从候选基因中选出最精简的预测组合
###   library(glmnet)
###   top_genes <- head(corData[order(-abs(corData$cor)), "Gene2"], 50)  # 取 Top50 候选
###   x <- exprSet[, top_genes]          # 候选基因表达矩阵
###   y <- geneEffect[, "ESR1"]          # ESR1 依赖性
###   model <- cv.glmnet(x, y, alpha = 1)  # LASSO 回归
###   目的：可能发现只需要 3~5 个基因就能很好预测 ESR1 依赖性
###
### 第三步（进阶）：随机森林 —— 评估预测能力 + 基因重要性排序
###   library(randomForest)
###   rf <- randomForest(x, y)
###   importance(rf)  # 各基因对预测的贡献排名
###   目的：验证预测准确性，捕捉基因间交互效应
###
### 生物学意义：
###   哪些基因的表达特征可以预测一个细胞系对 ESR1 抑制剂
###   （如氟维司群）的敏感性？
###   单变量 cor → 初筛候选
###   LASSO → 精简到可用的 biomarker panel（如 5 基因 qPCR panel）
###   随机森林 → 验证预测准确性
###
### ================================================================
### 结果解读：相关方向怎么翻译成生物学结论（速查）
###
### 本脚本固定的是 ESR1 的 Effect、循环全基因表达量：
###   每个数据点 = "这个细胞系里基因 X 表达多高" 配上
###               "这个细胞系里敲掉 ESR1 后死得多惨"
###
### Effect 坐标方向：越负 = 死得越多（强依赖）；0 附近 = 敲不敲无所谓；
###   死亡编码在"负"方向，活着/无影响编码在"零"方向。
###
### 负相关（r < 0）= 敏化因子 / 敏感人群 biomarker：
###   X 表达高 ↔ ESR1 Effect 更负 → X 高表达的患者用 ESR1 抑制剂疗效好
###   例：ESR1 自身表达 vs 自身 Effect 为负相关 ——
###   "免疫组化 ER 三加号的病人比一加号的疗效好"（临床已证实）；
###   PGR 同理 → ER+/PR+ 双指标把病人分成敏感人群
###
### 正相关（r > 0）= 抵抗因子 / 联合用药靶点：
###   X 表达高 ↔ ESR1 Effect 趋于 0 → X 的存在让 ESR1 敲除失效
###   （保护伞）→ 抑制 X 可恢复 ESR1 抑制剂的致死效果
###
### 【正相关的下游"富集 → 通路 → 协同致死"完整逻辑链见 05 脚本】
###   05 固定表达量、循环 Effect 并做 GSEA，正是 IL-6 那篇 Nature
###   （Hong2022, cGAS-STING/IL6/CIN）的做法；两种问法互为镜像：
###   04 问"谁能预测靶基因抑制剂的疗效"，05 问"基因 A 的表达
###   影响了哪些基因的依赖、保护了哪条通路"
###
### 警告：相关 ≠ 因果（可能只是癌种混杂，需控制 lineage / 实验闭环）；
###   n≈1000 大样本下 r>0.062 即 p<0.05（万物皆显著），
###   筛选主要看 |r| 大小（如 >0.3），p 值仅排除无信号者
### ================================================================

rm(list = ls())

## =================== 第一步：单变量筛选 ===================
## 读取 CRISPR 基因效应矩阵（行：细胞系，列：基因，值：Chronos 效应分）
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## 读取基因表达矩阵（行：细胞系，列：基因，值：log2(TPM+1)）
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")

## 快速演示：ESR1 依赖性 vs ESR1 表达的相关性
## 预期为负相关：ESR1 表达量高的细胞系，敲减 ESR1 的杀伤越强
## （临床对应：ER 免疫组化强阳性的患者，抗 ER 治疗疗效好）
genedata <- geneEffect[,"ESR1"]
cor.test(genedata, exprSet[,"ESR1"], method = "pearson")

## 批量计算：ESR1 依赖性 vs 全部基因表达的相关性
corData <- data.frame()
gene1 = "ESR1"
genedata <- geneEffect[,"ESR1"]

for (i in 1:ncol(exprSet)) {
  
  ## 1. 打印进度
  print(i)
  
  ## 2. 计算当前基因表达与 ESR1 依赖性的 Pearson 相关
  gene2 = colnames(exprSet)[i]
  dd = cor.test(genedata, exprSet[,gene2])
  
  ## 3. 存储结果
  corData[i,1] = gene1      ## Gene1 = ESR1（依赖性来源）
  corData[i,2] = gene2      ## Gene2 = 当前基因（表达来源）
  corData[i,3] = dd$estimate  ## 相关系数
  corData[i,4] = dd$p.value   ## p 值
}

colnames(corData) <- c("Gene1","Gene2","cor","pvalue")

## =================== 第二步：LASSO 特征选择 ===================
## 从单变量筛选的 Top 候选基因中，用 LASSO 选出最精简的预测组合
## 注意：机器学习（LASSO/随机森林）与本脚本的单变量相关分析本质相通
## —— 都是"表达矩阵当 X、ESR1 Effect 当 y"找预测特征，
## 只是机器学习能同时考虑多个基因的组合效应（交互项）

## ---- 自动安装缺失依赖 ----
for (pkg in c("glmnet", "randomForest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(glmnet)

## 取 |cor| 排名前 50 的候选基因作为 LASSO 输入
corData <- corData[order(-abs(corData$cor)), ]
top_genes <- head(corData$Gene2, 50)

## 构建特征矩阵 x（表达）和响应变量 y（ESR1 依赖性）
## 去掉含 NA 的细胞系
valid <- complete.cases(exprSet[, top_genes]) & complete.cases(geneEffect[, "ESR1"])
x <- as.matrix(exprSet[valid, top_genes])
y <- geneEffect[valid, "ESR1"]

## LASSO 交叉验证，选择最优 lambda
cv_fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10)

## 提取被 LASSO 保留的非零系数基因（即精简后的 biomarker panel）
lasso_coefs <- coef(cv_fit, s = "lambda.min")
selected_genes <- rownames(lasso_coefs)[which(lasso_coefs[, 1] != 0)]
selected_genes <- setdiff(selected_genes, "(Intercept)")  # 去掉截距
cat("\nLASSO 选出的 biomarker panel（", length(selected_genes), " 个基因）:\n")
print(selected_genes)

## =================== 第三步：随机森林评估 ===================
## 用选出的基因构建随机森林，评估预测能力 + 基因重要性排序

library(randomForest)

## 用 LASSO 选出的基因构建特征矩阵
rf_x <- as.matrix(exprSet[valid, selected_genes])
rf_y <- geneEffect[valid, "ESR1"]

## 训练随机森林回归模型
set.seed(42)
rf_model <- randomForest(rf_x, rf_y, ntree = 500, importance = TRUE)

## 查看模型预测性能（% Var explained 越高越好）
print(rf_model)

## 基因重要性排序（IncNodePurity 越大 = 对预测贡献越大）
rf_importance <- importance(rf_model)
rf_importance <- rf_importance[order(-rf_importance[, "%IncMSE"]), ]
cat("\n基因重要性排序（%IncMSE）:\n")
print(round(rf_importance, 4))

## =================== 保存结果 ===================
## 保存单变量筛选结果
saveRDS(corData, file = "TM00/DepMap_TM00/output/esr1_predictive_biomarkers_corData.rds")

## 保存 LASSO + 随机森林结果
saveRDS(list(
  lasso_model    = cv_fit,
  selected_genes = selected_genes,
  rf_model       = rf_model,
  rf_importance  = rf_importance
), file = "TM00/DepMap_TM00/output/esr1_predictive_biomarkers_model.rds")

## =================== 推荐阅读 ===================
### 预测模型构建,lasso,随机森林等机器学习算法
### https://bookdown.org/gongchangzhaojie/TranslationalBioinformaticsWithR/predictive-biomarkers.html