# tm00-script — DepMap 主力 R 脚本集

本仓库是 DepMap 数据工作台的分析源头，下游所有应用（tauri_app 的 Parquet、depmap_agent_v2 的计算模块）均派生自这里。
数据本体（DepMap 26Q1 CSV、rds 中间产物）不入库，只跟踪脚本。

## 主线一：基因 → 依赖性

| 脚本 | 内容 |
|---|---|
| `01_read_depmap.r` | 数据入口，读入 26Q1 各层矩阵，产出 `ccle_exprSet.rds`（约 1000 细胞系 × 18000 基因 log2(TPM+1)） |
| `02_gene_gene_correlation.R` | 基因共表达相关性：单基因对 cor.test 示例 + 固定基因 vs 全基因批量（ESR1 示例） |
| `03_co_dependency.R` | 共依赖分析：基因效应（CERES）层面基因间相关性 |
| `04_predivtive_biomarkers.R` | 预测性生物标志物 |
| `05_from_gene_to_dependency.R` | 从单基因出发的依赖性全景 |
| `05.1_ml_pathway.R` | ML 视角通路分析 |
| `05.2_synthetic_lethal.R` | 合成致死筛选 |
| `05.3_drug_sensitivity.R` | 药物敏感性关联 |
| `05.4_wgcna.R` | WGCNA 加权共表达网络 |
| `12_CCNE1_AMP_PKMYT1.R` | 案例复现：CCNE1 扩增 → PKMYT1 |
| `13_MYCN_DDX1_coamplification_Cancer_discovery.R` | 案例复现：MYCN-DDX1 共扩增（Cancer Discovery） |
| `14_DCAF5_SMARCB1_Nature.R` | 案例复现：DCAF5-SMARCB1（Nature） |
| `15_Sanger_CRISPR.R` | Sanger CRISPR 筛选数据对照 |
| `16_Dependency_nagative_correlation.R` | 不同肿瘤类型中基因依赖性分布（CDK4/CDK6，负相关视角） |
| `17_bipolar_dependency_ASB7_as_example.R` | 双向依赖基因对（ASB7 示例：H3K9me3 负调控） |
| `18_DrugAUC_and_DepMap_MTAPasExample.R` | 药敏 AUC 与多层组学（CNV/GeneEffect/表达/RNAi）相关性，MTAP 示例 |

## 主线二：突变 → 靶点（26Q1）

| 脚本 | 内容 |
|---|---|
| `06_mut_anchor_gene_selection_26Q1.R` | 突变锚定基因选择（四维框架：统计功效/功能性质/临床价值/可验证性） |
| `07_mutant_dependency_26Q1.R` | 突变组 vs 野生型依赖性差异（26Q1 口径） |
| `08_mutData_updata_26Q1.R` | 突变矩阵更新（Damaging/LoF、Hotspot 分套） |
| `09_batch_from_mut_to_target_26Q1.R` | 批量"突变后最依赖什么"（Score = Sensitivity × Difference，ARID1A 示例） |
| `10_batch_from_mut_to_target_26Q1_add_celltype.R` | 同上，按细胞类型分层 |
| `11_batch_from_gene_to_mut_26Q1_add_celltype.R` | 反向批量：靶点 → 突变背景（含细胞类型） |

## 主线三：互为 top-k 共依赖网络 → 模块发现（26Q1，新）

| 脚本 | 内容 |
|---|---|
| `20_mutual_topk_modules_26Q1.R` | TM00-13 课堂方法（CoffeeTime）的完整实现：全基因共依赖 Pearson 矩阵（一次算好存 RDS）→ mutual rank（S = rank+rankᵀ）长表 → k 网格（1~100）× 正/负相关建图 → Louvain 社区 → 跨尺度 Jaccard 稳定性 → 核心模块基因表（接 09~11 做实验候选）。阳性对照：已知复合体/通路应被捞回 |

## 工具

- `19_find_TLG_all(1).R`：全库 TLG（true love gene）查找
- `run_expr_to_dep_gsea.R`：表达→依赖性 GSEA

## 数据版本

- 当前统一口径 **26Q1**（2026-08-12 迁移完成）；02~05、12~18 均已升级。
- `20` 号脚本兼容 25Q3（课堂原型）与 26Q1，仅数据路径不同。

## 约定

- 每次更新脚本后在本仓库提交，commit message 写明改了哪个编号的脚本和数据版本（如 `06-08: 升级到 26Q1`）。
- 原型脚本验证后，按 `build_depmap_*` 参数化规范移植到生产仓库 Yu-Qiao-sjtu/wisp-science-depmap。
- 大矩阵类中间产物（相关矩阵 RDS 等）放 `TM00/DepMap_TM00/output/`，只算一次、下游复用。
