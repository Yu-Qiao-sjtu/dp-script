################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08

### 基因A的表达量调控哪些基因的Dependency
### ================================================================
### 分析思路（表达量 → 依赖性 → 通路富集）：
###
### 核心问题：基因A（如 ESR1）的表达量高低，会影响哪些基因的依赖性？
###          即：ESR1 高表达的细胞系，同时依赖哪些基因？
###
### 第一步：计算 ESR1 表达量 vs 全基因依赖性的相关性
### 第二步：将相关系数排序后做 GSEA 富集分析，看富集到哪些通路
###
### ================================================================
### 结果解读：相关方向怎么翻译成生物学结论（完整逻辑链）
###
### 本脚本固定的是【基因A的表达量】、循环【全基因的 Effect】：
###   每个数据点 = "这个细胞系里 A 表达多高" 配上
###               "这个细胞系里敲掉基因 Y 后死得多惨"
###
### Effect 坐标方向：越负 = 死得越多（强依赖）；0 附近 = 敲不敲无所谓；
###   死亡编码在"负"方向，活着/无影响编码在"零"方向。
###
### 负相关（r < 0）= A 高表达时 Y 更重要（敏化）：
###   A 表达越高 ↔ Y 的 Effect 越负 → A 高表达的肿瘤依赖 Y
###
### 正相关（r > 0）= A 高表达时 Y 被"保护"（抵抗）：
###   A 表达越高 ↔ Y 的 Effect 趋于 0 → 敲掉 Y 本来该致死，
###   A 养着它不死（保护伞作用）
###   → 对正相关基因列表做 GSEA/富集，反推 A 保护的是哪条通路
###     （Guilt by Association）
###
### 临床闭环（Hong et al. 2022 Nature, cGAS-STING/IL-6/CIN）：
###   CIN 肿瘤 → cGAS-STING 激活 → IL-6 自分泌
###   → IL-6 高表达使染色质相关基因 Effect 趋零（即与 IL6R 表达
###     正相关的基因富集到 chromatin 通路）
###   → 这群患者肿瘤不死、预后差
###   → 用 IL-6 通路抑制剂拆掉"保护伞"，恢复这群基因的致死效果
###   → 只给 IL-6 高/通路特征阳性的患者 → 协同致死（synthetic
###     lethality），把不可治疗变成可治疗；
###     范式同 ARID1A 失活肿瘤用 PARP 抑制剂：限定人群+靶向
###
### 与 04 脚本的区别（方向相反的两种问法，互为镜像）：
###   04 固定 ESR1 的 Effect、循环表达量 → 问"谁能预测 ESR1 抑制剂
###     疗效"（单变量 cor + LASSO + 随机森林，找 biomarker panel）
###   05 固定 ESR1/IL6R 的表达量、循环 Effect → 问"A 的表达影响了
###     哪些依赖、保护了哪条通路"（相关性排序 + GSEA → 联合用药靶点）
###   注：单个基因对的 cor(a,b)=cor(b,a) 是对称的，
###   区别在于固定谁、给谁做富集、回答什么临床问题
###
### 警告：相关 ≠ 因果（可能只是癌种混杂，需控制 lineage / 实验闭环）；
###   n≈1000 大样本下 r>0.062 即 p<0.05（万物皆显著），
###   筛选主要看 |r| 大小（如 >0.3）
### ================================================================

## =================== 第一部分：ESR1 表达量 → 全基因依赖性 + GSEA ===================

## 清空环境
rm(list = ls())

## 读取基因表达矩阵（行：细胞系，列：基因，值：log2(TPM+1)）
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")

## 读取 CRISPR 基因效应矩阵（行：细胞系，列：基因，值：Chronos 效应分）
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## 创建空 data.frame，用于存储结果
corData <- data.frame()

## 固定目标基因为 ESR1
gene1 = "ESR1"

## 提取 ESR1 在所有细胞系中的表达向量
genedata <- exprSet[,"ESR1"]

## 遍历 geneEffect 的每一列（每一个基因），计算其依赖性与 ESR1 表达的相关性
for (i in 1:ncol(geneEffect)) {
  
  ## 1. 打印进度
  print(i)
  
  ## 2. 获取当前基因名，计算其依赖性与 ESR1 表达的 Pearson 相关
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata, geneEffect[,gene2])
  
  ## 3. 存储结果：Gene1=ESR1(表达来源), Gene2=当前基因(依赖性来源)
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate  ## 相关系数
  corData[i,4] = dd$p.value   ## p 值
}

## 给结果表添加列名
colnames(corData) <- c("Gene1","Gene2","cor","pvalue")

## =================== GSEA 富集分析 ===================

library(clusterProfiler)

## 构建 GSEA 所需的排序基因列表
## 取负号（注意方向，勿搞反）：
##   cor 为负 = ESR1 高表达时该基因 Effect 更负（依赖强、被敏化），
##             取负后值最大，排在排序列表【顶部】
##   cor 为正 = ESR1 高表达时该基因 Effect 趋于 0（依赖弱、被保护），
##             取负后值最小，排在排序列表【底部】
## 这样 GSEA 上富集 = ESR1 高表达时【依赖更强】的通路（敏化方向）
##        下富集 = ESR1 高表达时【依赖减弱】的通路（保护/拮抗方向，
##        对应 IL-6 高表达保护 chromatin 相关基因的场景）
mygeneList <- -corData$cor
names(mygeneList) <- corData$Gene2
mygeneList <- sort(mygeneList, decreasing = T)
head(mygeneList)

## 读取 Hallmark 基因集（50 个经典通路）
geneSet <- read.gmt("TM00/DepMap_TM00/resource/geneSets/h.all.v2023.2.Hs.symbols.gmt")

## 运行 GSEA
mygsea <- GSEA(geneList = mygeneList, TERM2GENE = geneSet)

## 转为 data.frame 查看完整结果表
data <- as.data.frame(mygsea)

## 气泡图：按正/负富集分面展示 Top30 通路
library(ggplot2)
dotplot(mygsea, showCategory = 30,
        split = ".sign",
        font.size = 8,
        label_format = 60) + facet_grid(~.sign)

## 单个通路的 GSEA 富集曲线图（这里看 MYC 靶基因通路）
library(enrichplot)
gseaplot2(mygsea, "HALLMARK_MYC_TARGETS_V2", color = "red", pvalue_table = T)

## =================== 第二部分：IL6R 表达量 → 全基因依赖性 + GSEA ===================
## 换一个目标基因 IL6R，用 Reactome 通路库做富集
## IL6R 与染色体不稳定（CIN）相关，预期富集到 DNA 损伤修复等通路

## 清空环境
rm(list = ls())

## 读取表达矩阵
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")

## 读取基因效应矩阵
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## 创建空 data.frame
corData <- data.frame()

## 固定目标基因为 IL6R
gene1 = "IL6R"

## 提取 IL6R 表达向量
genedata <- exprSet[,"IL6R"]

## 遍历 geneEffect 每一列，计算 IL6R 表达与各基因依赖性的相关性
for (i in 1:ncol(geneEffect)) {
  
  ## 1. 打印进度
  print(i)
  
  ## 2. 计算当前基因依赖性与 IL6R 表达的 Pearson 相关
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata, geneEffect[,gene2])
  
  ## 3. 存储结果
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate
  corData[i,4] = dd$p.value
}

## 列名与第一部分不同：exp=表达来源基因, Dependency=依赖性来源基因
colnames(corData) <- c("exp","Dependency","cor","pvalue")

## =================== GSEA 富集分析 ===================

library(clusterProfiler)

## 构建排序基因列表（取负号逻辑同上）
mygeneList <- -corData$cor
names(mygeneList) <- corData$Dependency
mygeneList <- sort(mygeneList, decreasing = T)
head(mygeneList)

## 这次用 Reactome 通路库（比 Hallmark 更详细，上千条通路）
geneSet <- read.gmt("TM00/DepMap_TM00/resource/geneSets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

## 运行 GSEA
mygsea <- GSEA(geneList = mygeneList, TERM2GENE = geneSet)

## 查看完整结果
data <- as.data.frame(mygsea)

## 气泡图
library(ggplot2)
dotplot(mygsea, showCategory = 30,
        split = ".sign",
        font.size = 8,
        label_format = 60) + facet_grid(~.sign)

## 单个通路的富集曲线（染色体维护通路，与 IL6/CIN 文献呼应）
library(enrichplot)
gseaplot2(mygsea, "REACTOME_CHROMOSOME_MAINTENANCE", color = "red", pvalue_table = T)

## =================== 第三部分：PROGENy 通路活性量化 ===================
## 升级思路：从单基因表达升级为 11 条核心信号通路的活性打分
## PROGENy 用下游响应基因表达反推上游通路活性，比单基因更稳健
## 生物学依据（Hong et al. 2022 Nature）：
##   CIN（染色体不稳定）→ cGAS-STING 激活 → IL-6 自分泌 → 癌细胞依赖 IL6R 生存
##   用 PROGENy 的 IL6_STAT3 通路活性替代单个 IL6R 基因表达，更接近生物学真相

## ---- 确保 BiocManager 可用（必须先于下面的安装循环）----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## ---- 自动安装缺失依赖 ----
for (pkg in c("progeny", "decoupleR", "dorothea", "dplyr", "tibble", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("progeny", "decoupleR", "dorothea")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(progeny)
library(decoupleR)
library(dorothea)
library(dplyr)
library(tibble)
library(tidyr)

## 重新读取表达矩阵（前面 rm 清掉了）
exprSet <- readRDS(file = "TM00/DepMap_TM00/output/ccle_exprSet.rds")
geneEffect <- readRDS(file = "TM00/DepMap_TM00/output/depmap_geneEffect.rds")

## ---- PROGENy：计算 11 条信号通路活性分 ----
## PROGENy 需要 ExpressSet 或 matrix，且基因名需为 HGNC symbol
## 将表达矩阵转置（PROGENy 要求行=基因，列=样本）
exprMat <- t(as.matrix(exprSet))

## 计算 11 条通路活性分（Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT）
progeny_scores <- progeny(exprMat, scale = TRUE, organism = "Human")

## 转为 data.frame，行=细胞系，列=通路活性分
progeny_df <- as.data.frame(progeny_scores)

## 查看实际列名（确认 PROGENy 返回的通路名称）
colnames(progeny_df)
head(progeny_df)

## ---- 验证文献假说：JAK-STAT 通路活性 vs 全基因依赖性 ----
## 用 JAK-STAT 通路活性替代单个 IL6R 基因表达
stat3_activity <- progeny_df[, "JAK-STAT"]

## 批量计算 JAK-STAT 通路活性 vs 全基因依赖性的相关性
corData_progeny <- data.frame()
gene1 = "JAK-STAT_activity"
genedata <- stat3_activity

for (i in 1:ncol(geneEffect)) {
  print(i)
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata, geneEffect[,gene2])
  corData_progeny[i,1] = gene1
  corData_progeny[i,2] = gene2
  corData_progeny[i,3] = dd$estimate
  corData_progeny[i,4] = dd$p.value
}
colnames(corData_progeny) <- c("Pathway","Gene","cor","pvalue")

## 构建排序基因列表，做 GSEA
library(clusterProfiler)
mygeneList_progeny <- -corData_progeny$cor
names(mygeneList_progeny) <- corData_progeny$Gene
mygeneList_progeny <- sort(mygeneList_progeny, decreasing = T)

## GSEA（用 Reactome 通路库）
geneSet_reactome <- read.gmt("TM00/DepMap_TM00/resource/geneSets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
mygsea_progeny <- GSEA(geneList = mygeneList_progeny, TERM2GENE = geneSet_reactome)

## 气泡图
library(ggplot2)
dotplot(mygsea_progeny, showCategory = 30,
        split = ".sign",
        font.size = 8,
        label_format = 60) + facet_grid(~.sign)

## =================== 第四部分：decoupleR 通路活性推断 ===================
## 升级思路：用 decoupleR 推断每个细胞系的下游通路活性
## 与 PROGENy 互补：decoupleR 支持多种基因集库，覆盖面更广

## 用 decoupleR 的 enrich基因集库（基于 Footprint）
## 这里使用 decoupleR 内置的 Dorothea 转录因子活性推断作为示例

## 获取 Dorothea 转录因子调控网络（置信度 A-C，高质量 TF-Target 关系）
## 直接用 dorothea 包内置数据，无需 OmnipathR
library(dorothea)
data("dorothea_hs", package = "dorothea")
net <- dorothea_hs
## 仅保留置信度 A、B、C 的调控关系
net <- net[net$confidence %in% c("A", "B", "C"), ]
cat("转录因子调控网络：", nrow(net), "条 TF-Target 关系，",
    length(unique(net$tf)), "个转录因子\n")

## 计算每个细胞系的转录因子活性分（用 ulm 方法：无权重最小二乘）
## decoupleR 2.x API：参数 net→network，.mor 通过 args 传递
net <- net[, c("tf", "target", "mor")]  ## 确保只保留需要的列
tf_scores <- decouple(mat = exprMat,
                      network = net,
                      .source = "tf",
                      .target = "target",
                      statistics = "ulm",
                      args = list(.mor = "mor"),
                      minsize = 5)

## 提取 ulm 分数，转为矩阵形式（行=细胞系，列=转录因子）
tf_mat <- tf_scores |>
  filter(statistic == "ulm") |>
  pivot_wider(id_cols = condition, names_from = source, values_from = score) |>
  as.data.frame()
rownames(tf_mat) <- tf_mat[,1]
tf_mat <- tf_mat[,-1]
tf_mat <- as.matrix(tf_mat)

## 查看结果：每个细胞系在各转录因子上的活性分
dim(tf_mat)
head(tf_mat[, 1:6])

## ---- 验证假说：STAT3 转录因子活性 vs 全基因依赖性 ----
## 用 STAT3 转录因子活性分替代单个 IL6R 表达
tf_target <- "STAT3"
if (tf_target %in% colnames(tf_mat)) {
  tf_target_data <- tf_mat[, tf_target]
  corData_tf <- data.frame()
  genedata <- tf_target_data
  gene1 <- tf_target

  for (i in 1:ncol(geneEffect)) {
    print(i)
    gene2 = colnames(geneEffect)[i]
    dd = cor.test(genedata, geneEffect[,gene2])
    corData_tf[i,1] = gene1
    corData_tf[i,2] = gene2
    corData_tf[i,3] = dd$estimate
    corData_tf[i,4] = dd$p.value
  }
  colnames(corData_tf) <- c("TF","Gene","cor","pvalue")

  ## GSEA
  mygeneList_tf <- -corData_tf$cor
  names(mygeneList_tf) <- corData_tf$Gene
  mygeneList_tf <- sort(mygeneList_tf, decreasing = T)

  mygsea_tf <- GSEA(geneList = mygeneList_tf, TERM2GENE = geneSet_reactome)

  dotplot(mygsea_tf, showCategory = 30,
          split = ".sign",
          font.size = 8,
          label_format = 60) + facet_grid(~.sign)
}

## =================== 保存结果 ===================
saveRDS(progeny_scores, file = "TM00/DepMap_TM00/output/progeny_pathway_scores.rds")
saveRDS(tf_scores, file = "TM00/DepMap_TM00/output/decoupleR_tf_scores.rds")

## =================== 推荐阅读 ===================
### cGAS–STING drives the IL-6-dependent survival of chromosomally instable cancers
### Hong et al. 2022 Nature
### https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/pw_bk.html
### https://saezlab.github.io/progeny/articles/progeny.html
###
### 升级逻辑：
### 单基因表达（IL6R）→ PROGENy 通路活性（JAK-STAT3）→ decoupleR 转录因子活性（STAT3）
### 越往上越接近生物学机制，越不受单个基因波动影响
###
### 拓展方向
### 机器学习（基于通路活性的预测模型）
### 协同致死（通路活性相关的依赖性基因对）
### 药物敏感性（通路活性 vs PRISM/GDSC 药物筛选）
### 基因功能模块（WGCNA 共表达网络）