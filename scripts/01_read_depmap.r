################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（22Q2 → 26Q1）

### ================================================================
### DepMap 数据读取与预处理（入门版）
###
### 本脚本演示如何加载 DepMap 26Q1 的 5 类核心数据，
### 对齐细胞系和基因，保存为 RDS 供后续脚本使用。
###
### 新版数据的主要变化（22Q2 → 26Q1）：
###   - CRISPR_gene_effect.csv → CRISPRGeneEffect.csv（第一列 V1 = ACH-编号）
###   - CRISPR_gene_dependency.csv → CRISPRGeneDependency.csv
###   - CCLE_expression.csv → OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
###     （新增 6 列元数据，需过滤 IsDefaultEntryForMC）
###   - sample_info.csv → Model.csv（字段名：DepMap_ID → ModelID）
###   - CCLE_mutations.csv → OmicsSomaticMutationsMatrixDamaging.csv
###     （格式从长表变为宽矩阵；损伤性判定字段同步更替：
###       旧长表 isDeleterious → 中间版 CCLEDeleterious →
###       现行宽矩阵由 LikelyLoF 聚合而来，详见脚本 07/08 头部注释）
###
### 注意：完整的突变数据处理见脚本 07/08
### ================================================================

rm(list = ls())
library(data.table)
output_dir <- "TM00/DepMap_TM00/output"

## ---- 1. Gene Effect（CRISPR 基因效应）----
## 数据结构：行 = 细胞系（约 1000+ 个，以 ACH 编号/ModelID 标识，如 ACH-000001）
##           列 = 基因（约 18000+ 个，列名形如 "TP53 (7157)"，即 基因名 (Entrez ID)）
## 矩阵值 = Gene Effect 分数：基于 Chronos 算法校正拷贝数偏倚后的 CRISPR 敲除效应值
##   - 越负（如 -1.5）→ 该细胞系对该基因依赖性越强（敲除后细胞存活受损严重）
##   - 接近 0        → 敲除该基因对该细胞系无显著影响
##   - 通常 < -0.5 左右视为"依赖"（对应第 2 步 Gene Dependency 的概率化阈值）
## fread 多线程读取大矩阵（远快于 read.csv）；data.table = F 使返回值为普通 data.frame，
## 便于后续直接用 data.frame 语法操作
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
## 第一列说明：原文件第一列是细胞系编号（DepMap Model ID，如 ACH-000001），
## 但其表头为空，fread 读入后自动命名为 V1，值形如：
##   V1,          A1BG (1),  A1CF (29974), ...
##   ACH-000001,  -0.05,     0.12,        ...
##   ACH-000002,  -0.89,     0.03,        ...
## 若不处理，V1 会以字符型数据列混在数值矩阵里，干扰后续 cor()/scale() 等运算
rownames(geneEffect) <- geneEffect[, 1]     ## 把 V1 列（ACH 编号）设为行名，作为细胞系索引
geneEffect <- geneEffect[, -1]              ## 删掉 V1 列（编号已存进行名，避免字符列混入数值矩阵）
## 去掉列名中的 Entrez ID 后缀："TP53 (7157)" → "TP53"
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))
cat("Gene Effect:", nrow(geneEffect), "模型 ×", ncol(geneEffect), "基因\n")

## ---- 2. Gene Dependency（基因依赖性概率）----
## 与第 1 步 Gene Effect 同结构，但矩阵值是 0~1 的依赖概率（基于贝叶斯推断），
## 通常 > 0.5 视为该细胞系依赖该基因
geneDependency <- data.table::fread("data/CRISPRGeneDependency.csv", data.table = F)
## 第一列说明：与第 1 步相同，原文件第一列是细胞系编号（ACH-000001 等），
## 表头为空，fread 读入后自动命名为 V1
rownames(geneDependency) <- geneDependency[, 1]   ## 把 V1 列（ACH 编号）设为行名，作为细胞系索引
geneDependency <- geneDependency[, -1]            ## 删掉 V1 列（编号已存进行名，避免字符列混入数值矩阵）
## 去掉列名中的 Entrez ID 后缀："TP53 (7157)" → "TP53"，便于按基因名操作和跨数据集对齐
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneDependency))
cat("Gene Dependency:", nrow(geneDependency), "模型 ×", ncol(geneDependency), "基因\n")

## ---- 3. CCLE 细胞系表达量（新版含元数据列，需过滤）----
exprRaw <- data.table::fread("data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                              data.table = F)
## 新版前 6 列是元数据（与旧版 CCLE_expression.csv 不同，旧版首列直接是 ACH 编号）：
##   V1                  无表头索引列（fread 自动命名）
##   SequencingID        测序唯一标识（同一细胞系可能有多次测序）
##   ModelConditionID    模型培养条件 ID
##   ModelID             细胞系编号（ACH-xxx），与 Gene Effect 行名同源，是对齐钥匙
##   IsDefaultEntryForMC           该测序是否为该模型条件的默认代表（Yes/No）
##   IsDefaultEntryForModel        该测序是否为该模型的默认代表（Yes/No）
## ★ 最重要的是前两者：
##   ModelID —— 整个 DepMap 体系的主键，表达量/依赖性/突变/细胞系信息全靠它对齐
##   IsDefaultEntryForMC —— 决定每个细胞系用哪次测序的过滤开关，不过滤则同一细胞系多行重复
##   （其余 4 列为辅助信息：V1 是冗余行号、IsDefaultEntryForModel 突变数据用、
##    SequencingID/ModelConditionID 仅追溯原始测序/培养条件时才用）
## 第 7 列起为基因列（约 19000+ 个），列名形如 "TP53 (7157)"，值 = log2(TPM + 1)：
##   0 表示不表达，数值越大表达越高；log2+1 变换使分布近似正态，便于相关/差异分析
## 同一细胞系可能有多次测序/多个培养条件（行数 > 细胞系数），
## 故需：① 只保留每个模型条件的默认测序；② 去重（过滤后仍可能重复，
## 不去重直接设行名会报 duplicate 'row.names' 错误）
## 两行均为“行筛选”（自我覆盖：用筛选后的子集覆盖原变量，不产生新名字）：
##   df[条件, ]      → 按行保留满足条件的行；
##   df$列 == "Yes"   → 逻辑向量，该列值为 Yes 的行为 TRUE；
##   duplicated(向量) → 向量中重复出现的后续行为 TRUE，首次出现为 FALSE，
##                       取反 ！后即“每个值只留首次出现的那行”
exprRaw <- exprRaw[exprRaw$IsDefaultEntryForMC == "Yes", ]   ## ① 只留 IsDefaultEntryForMC == "Yes" 的行
exprRaw <- exprRaw[!duplicated(exprRaw$ModelID), ]          ## ② 按 ModelID 去重，保证行名唯一
## 为什么②仍需要：①只保证“每个模型条件一个默认测序”，但同一 ModelID 可以有多个模型条件
## （如 2D 培养和器官筹培养、不同代次/培养基），每个条件各有一条默认测序，
## 过滤后同一 ModelID 仍可能出现多行；若直接设行名会报错：
##   Error in `.rowNamesDF<-`: duplicate 'row.names' are not allowed
## 故必须先去重。duplicated() 返回重复出现的后续行，取反即每个 ModelID 只留首次出现的那行
rownames(exprRaw) <- exprRaw$ModelID
## 去掉元数据列，只保留基因列（反向筛选：保留列名不在 meta_cols 里的列，
## 按列名匹配而非位置删除，对版本列变动更稳健）
meta_cols <- c("V1", "SequencingID", "ModelConditionID", "ModelID",
                "IsDefaultEntryForMC", "IsDefaultEntryForModel")
exprSet <- exprRaw[, !(colnames(exprRaw) %in% meta_cols)]   
## colnames %in% meta_cols 逐列判断是否在待删清单，取反后保留基因列
rm(exprRaw)                                                  
## 删除中间变量，释放内存（此时信息已全部转移到 exprSet）
## 去掉列名中的 Entrez ID 后缀："TP53 (7157)" → "TP53"，与第 1/2 步保持一致
## 便于后续与 Gene Effect 按基因名取交集对齐
colnames(exprSet) <- gsub("\\s+\\(\\d+\\)", "", colnames(exprSet))
cat("Expression:", nrow(exprSet), "模型 ×", ncol(exprSet), "基因\n")

## ---- 4. 细胞系信息（Model.csv）----
## 数据结构：一行 = 一个细胞系模型（约 2000+ 个，含少量未做 CRISPR/测序的模型），
## 一列 = 一项注释字段（共 47 列），是全部 DepMap 分析的“注释字典”
## 最常用的字段（后续脚本取子集/分组全靠它们）：
##   ModelID               细胞系编号（ACH-xxx），主键，与前三份数据的行名同源
##   CellLineName          常用细胞系名（如 MCF-7、NIH:OVCAR-3）
##   OncotreeLineage       组织谱系（如 Breast、Myeloid）——按大谱系分组用
##   OncotreePrimaryDisease 原发疾病（如 Breast Cancer、Lung Adenocarcinoma）——按癌种分组用
##   OncotreeSubtype       分型（如 Luminal、High-Grade Serous Ovarian Cancer）
##   Sex / Age / SampleCollectionSite  患者性别/年龄/取样部位
##   GrowthPattern         生长方式（Adherent 贴壁 / Suspension 悬浮）
## 其余字段为溯源/培养/治疗信息（RRID、SourceType、CatalogNumber、培养基、分期分级等）
## 注意：ModelID 列不删除（不像前三份数据）——后续需用 cellinfor$xxx 取注释列，保留完整表更好用
## 全部 47 列的逐列说明（含示例与使用频率）见：docs/Model_csv_列名说明文档.md
##
## ⚠ 模型质控提醒：模型库中混有四类可能污染分析的“问题模型”，
## 全部可通过本表字段识别（识别方法与已验证清单见：docs/Model_质控说明文档.md）：
##   ① 体外诱导耐药系（8 株，如 A-375_DAB_R）—— EngineeredModelDetails="drug adapted" + CulturedResistanceDrug
##   ② 工程改造克隆（如 HCC-827-GR5，肺癌，MET 扩增耐药）—— EngineeredModelDetails 非空
##   ③ 非癌细胞系（几十株，如 MCF 10A、HA1E、成纤维细胞）—— OncotreePrimaryDisease="Non-Cancerous"
##   ④ 病毒转化/身份存疑系（EBV/HPV 转化、衍生物、疑似错鉴）—— ModelTreatment / PublicComments
## 其中 7 株耐药系已验证在本项目 CRISPR 矩阵内（经 CRISPRScreenMap.csv 核实）。
## 做“癌种 vs 其他”对比时背景组尤其要清洗（!＝ 会把非癌系全部收进背景组）；
## 做耐药机制研究时耐药系应保留并与亲本配对。质控代码见文档“综合质控代码”一节。
cellinfor <- data.table::fread("data/Model.csv", data.table = F)
rownames(cellinfor) <- cellinfor$ModelID      ## 设 ModelID 为行名；旧版字段 DepMap_ID → ModelID
cat("Model info:", nrow(cellinfor), "细胞系\n")

## ---- 5. 突变数据（宽矩阵格式，旧版为长表 CCLE_mutations.csv）----
## 损伤性判定字段演进：isDeleterious（旧，已废弃）→ CCLEDeleterious
##   → LikelyLoF（现行：OncoKB LoF 注释 OR VEP impact HIGH）。
## 本宽矩阵即对 MAF 长表 LikelyLoF==True 聚合的结果，值 0/1/2。
## 完整说明见 08 脚本头部与 docs/26Q1_Mutation_Pipeline_Documentation.pdf
mutData <- data.table::fread("data/OmicsSomaticMutationsMatrixDamaging.csv",
                              data.table = F)
## 新版前 6 列是元数据，需过滤 IsDefaultEntryForModel == "Yes"
mutData <- mutData[mutData$IsDefaultEntryForModel == "Yes", ]
rownames(mutData) <- mutData$ModelID
mutData <- mutData[, !(colnames(mutData) %in%
                         c("V1", "ModelID", "ModelConditionID",
                           "IsDefaultEntryForMC", "IsDefaultEntryForModel",
                           "GenomeVersion"))]
colnames(mutData) <- gsub("\\s+\\(\\d+\\)", "", colnames(mutData))
cat("Mutation matrix:", nrow(mutData), "模型 ×", ncol(mutData), "基因\n")

################################
## 修剪：三套数据取交集
## 目的：五份数据的细胞系/基因覆盖范围不同（如 Model.csv 有 2155 株，
## CRISPR 篛库只有 ~1000+ 株），不取交集则交叉分析时会出现查不到/NA。
## 交集后三套数据不仅行数相同，行顺序也完全一致（第 i 行 = 同一株细胞），
## 这是后续按位置对应操作的前提。注意 mutData 不参与本轮交集（见脚本 07/08）
################################

## 第一次交集：CRISPR（geneEffect） ∩ 表达量（exprSet）的细胞系
## intersect() 返回两向量共同元素；rownames = ModelID（ACH-xxx）
commonindex <- intersect(rownames(geneEffect), rownames(exprSet))
## 第二次交集：再 ∩ 细胞系注释表（cellinfor）
commonindex <- intersect(commonindex, rownames(cellinfor))
cat("三套数据共有细胞系:", length(commonindex), "\n")

## 三套数据按同一份 commonindex 同步修剪（行对齐）
geneEffect <- geneEffect[commonindex, ]
exprSet    <- exprSet[commonindex, ]
cellinfor  <- cellinfor[commonindex, ]

## 共有基因：geneEffect（~18000+ 基因） ∩ exprSet（仅蛋白编码 ~19000）
## 两矩阵基因范围不同，取交集后列也完全对齐
commonGenes <- intersect(colnames(geneEffect), colnames(exprSet))
cat("共有基因:", length(commonGenes), "\n")
geneEffect <- geneEffect[, commonGenes]
exprSet    <- exprSet[, commonGenes]

################################
## 保存 RDS（供后续脚本使用）
## RDS 是 R 专属二进制序列化格式：读取快（行名列名直接恢复）、体积小；
## 目的：一次清洗处处复用——后续脚本只需 readRDS() 秒级加载，
## 不必重跑本脚本（读 35GB 原始 CSV + 清洗需几分钟）
################################

saveRDS(geneEffect, file = file.path(output_dir, "depmap_geneEffect.rds"))   ## CRISPR 基因效应矩阵
saveRDS(exprSet, file = file.path(output_dir, "ccle_exprSet.rds"))          ## 表达量矩阵
## ⚠ cellinfor 存别名 cellinfor_exprset.rds，不再覆盖共享的 cellinfor.rds：
##   此处 cellinfor 已被裁剪到 commonindex（geneEffect∩exprSet 共有，~1140 株）；
##   而共享的 cellinfor.rds 应保持全量版（Model.csv 全部 ~2154 株，由 07/08 脚本
##   写入），供 09-18 突变×依赖类脚本做三重交集（否则 coID 会被压到 1140，
##   与 dep∩mut 双交集 1208 口径不一致——2026-09-09 已修复过一次）
saveRDS(cellinfor, file = file.path(output_dir, "cellinfor_exprset.rds"))    ## 细胞系注释表（裁剪版，与 exprSet 同行同序）
cat("已保存 3 个 RDS 到", output_dir, "\n")
