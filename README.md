# tm00-script — DepMap 主力 R 脚本集

本仓库是 DepMap 数据工作台的分析源头，下游所有应用（tauri_app 的 Parquet、depmap_agent_v2 的计算模块）均派生自这里。
数据本体（DepMap 26Q1 CSV、rds 中间产物）不入库，只跟踪脚本。

## 主线一：基因 → 依赖性

| 脚本 | 内容 |
|---|---|
| `01_read_depmap.r` | 数据入口，读入 26Q1 各层矩阵 |
| `02~04` | 基因共相关、共依赖、预测生物标志物 |
| `05*` | 基因出发的多视角分析：ML 通路 / 合成致死 / 药敏 / WGCNA |
| `12~18` | 已发表案例复现（CCNE1-PKMYT1、MYCN-DDX1、DCAF5-SMARCB1、Sanger、MTAP 等） |

## 主线二：突变 → 靶点（26Q1）

| 脚本 | 内容 |
|---|---|
| `06` | 突变锚定基因选择 |
| `07` | 突变依赖性分析 |
| `08` | 突变数据更新（26Q1） |
| `09~11` | 突变↔靶点双向批量版（含细胞类型） |

## 工具

- `19_find_TLG_all(1).R`：全库 TLG（true love gene）查找
- `run_expr_to_dep_gsea.R`：表达→依赖性 GSEA

## 约定

- 每次更新脚本后在本仓库提交，commit message 写明改了哪个编号的脚本和数据版本（如 `06-08: 升级到 26Q1`）。
- 原型脚本验证后，按 `build_depmap_*` 参数化规范移植到生产仓库 Yu-Qiao-sjtu/wisp-science-depmap。
