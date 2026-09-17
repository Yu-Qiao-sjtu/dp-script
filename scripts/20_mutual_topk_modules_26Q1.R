################################################
################################################
### 作者：果子（课堂原型 12_CO_IC/CK_modifier_rank.R、05_KD_modules_graph.R）
### 改编与扩展：乔宇（上海交通大学，免疫药理PhD）
### 创建时间：2026-09-17
### 数据口径：DepMap 26Q1（课堂原型为 25Q3）

### ================================================================
### 主题：互为 top-k 共依赖网络（mutual top-k co-fitness）→ 跨尺度稳定模块发现
###
### 科学思想（对应 TM00-13 课程"CoffeeTime / 从小到大发文章"）：
###   1. 两个基因的"关系"有多种强度口径。k=1（互为排名第一）排他性太强，
###      只能得到孤立的基因对；把 k 放大，基因对逐渐连成社区（模块）。
###   2. mutual rank（互为排名和 S = rank_A|B + rank_B|A）是对称化的
###      top-k 口径：一个基因很强但另一个很弱的关系会被 S 惩罚掉，
###      避免超级枢纽（如必需基因）把整张图吸成一颗星。
###   3. 跨尺度稳定性 = 生物学真实性的过滤器：真实的功能单元（复合体、
###      通路、共适应回路）在 k 从小到大的过程中模块组成基本不动；
###      只在大 k 才出现、随 k 剧烈变形的，多是稀释出来的噪声。
###   4. 结果必须自带阳性对照：已知模块（血红素、蛋白酶体、MCM 复合体
###      等）能被捞回来 = 方法正确；跨尺度稳定但注释未知/弱的模块 =
###      发现空间，接 mutation→target 流水线（09-11 脚本）做实验候选。
###
### 方法本质（图论视角，先说清楚"且"为什么不能改成"或"）：
###   把每个基因 |rho| 排名前 k 的伙伴记为它的 k 近邻集合 N(A)，
###   无向边 (A,B) 存在当且仅当 B ∈ N(A) 且 A ∈ N(B)——这就是图论
###   标准的 mutual k-NN 图（互近邻图）构造。"且"是定义本身，不是
###   实现选择：
###     · "或"（union k-NN 图）允许单向边：A 单方面把 B 排进前 k 就
###       连边。k 小时无妨，k 放大后几乎每对中度相关基因都会落入
###       某一方的名单，弱关系大量涌入，整图塌成一个连通块，
###       Louvain 只能给出一个巨大的"万能模块"，什么也分不出来。
###     · "且"（mutual 图）每条边经过两端一致确认，天然无向、稀疏、
###       无 Hub 病——是社区发现（Louvain/谱聚类）建图的文献标准。
###   ⚠ 课堂口述（23:03）说"这个地方要是一个或"，但控制台代码
###   （23:54）写的是 filter(rank.x < k, rank.y < k)，逗号 = AND；
###   网站页面"伙伴排名是否彼此靠前"也是双向判定。手比嘴诚实，
###   以代码为准：本脚本用且，勿对照课程回放把这里"纠正"回或。
###
###   在 SNA 骨架上有三层方法论增量：
###   · 骨架 = SNA 社区发现：基因是节点，互为 top-k 的关系是边，
###     |rho| 是边权，cluster_louvain() 找模块度高的社区——Louvain
###     本就是 2008 年为大规模社交网络设计的，换个节点类型就能用。
###   · 增量一：建图用 mutual k-NN（源自机器学习）。直接卡 |rho|
###     阈值会得到以必需基因为中心的"名人图"——人人都连着它，
###     整图糊成一团；互近邻强制留下互相选择的关系，图才散得开。
###   · 增量二：关系带符号。正相关 = 共同必需，负相关 = 合成致死/
###     互斥，是两种生物学，分开建两张图分别聚类；负网络里的模块
###     恰是找靶点最有价值的。
###   · 增量三：跨尺度稳定性检验 = 把聚类稳健性分析用成生物学发现
###     标准。Louvain 对分辨率敏感是社区发现的著名缺陷；"k 从小到
###     大滑、找一直存在的模块"把"这是真功能单元还是算法伪影"变成
###     可计算判据（Jaccard ≥ 0.7 across k），同时充当方法校验
###     （已知复合体被捞回）和发现过滤器（稳定但未知的模块才值得
###     做实验）。
###
###   mutual rank（S = rank + rankᵀ）与 mutual 图是同一思想的软硬
###   两态：S 是软化形态（不设门槛，按两边排名之和全排序，适合输出
###   长表 browsable）；k 截断建图是硬化形态（严格互近邻）。本脚本
###   第 3 部分是软的，第 4 部分是硬的。
###
### 这个脚本的科学思想，剥到底就是一句话：
###   基因间关系的强弱不是一条阈值线，而是一个可滑动的尺度 k；
###   在每个尺度下用互近邻图问一次"谁和谁抱团"，那么在所有尺度
###   下都稳定抱团的模块，大概率是真实的生物学功能单元——已知的
###   是方法正确的阳性对照，未知的是值得做实验的发现空间。
###   其余一切（SNA、Louvain、Jaccard）都是把这句话变成可计算
###   流程的工具。
###
### 流程（对应课堂两个脚本，k 网格循环 + 稳定性为本脚本新增）：
###   第 1 部分：读 26Q1 GeneEffect，清洗列名（同 07/09 口径）
###   第 2 部分：全基因 Pearson 共依赖矩阵（课堂约 4 min，存 RDS 复用）
###   第 3 部分：mutual rank 上三角长表（课堂原样照搬，保证可对齐）
###   第 4 部分：k 网格 × 正/负相关 → igraph 图 → Louvain 社区
###   第 5 部分：跨尺度 Jaccard 稳定性 → 核心模块识别
###   第 6 部分：模块汇总表与图
###
### 运行环境：
###   建议 ≥16 GB 内存（相关矩阵 + 排名矩阵 + S 各约 2.7 GB，峰值 ~8 GB）
###   R 包：data.table、igraph（首次使用需 install.packages("igraph")）
###   工作目录：仓库根目录（与 07/09 相同，相对路径 data/... TM00/...）
### ================================================================

rm(list = ls())
library(data.table)
library(igraph)

output_dir <- "TM00/DepMap_TM00/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "20_modules"), showWarnings = FALSE, recursive = TRUE)

### 可调参数（课堂演示值：k=20、rho_use=0.2）
k_grid    <- c(1, 2, 3, 5, 10, 20, 50, 100)  # k 尺度网格：从"互为第一"到宽松
rho_min   <- 0.2                             # 边的最低 |rho| 门槛（课堂 rho_use=0.2）
ref_k     <- 10                              # 稳定性参考尺度（模块已成型又不至于过稀释）
stab_thre <- 0.7                             # Jaccard ≥ 0.7 视为"该尺度下模块保持"
q_var     <- 0                               # 低方差基因过滤分位数（0 = 不过滤，与课堂全量对齐）

#############################################################
## 第 1 部分：读取 26Q1 CRISPR Gene Effect
#############################################################
## 数据结构（同 07 脚本口径）：
##   行 = 细胞系（ACH- 编号），列 = 基因（"SYMBOL (ENTREZ)" 需清洗）
##   值 = CERES/Chronos 效应评分，越负 = 越依赖
## 这里行是样本、列是基因，cor() 直接在列（基因）之间算，
## 即"两个基因在哪些细胞系里共同必需/共同不必需"的共依赖相关。
cat("读取 26Q1 Gene Effect...\n")
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = FALSE)
rownames(geneEffect) <- geneEffect[[1]]
geneEffect[[1]] <- NULL
clean_names <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))
colnames(geneEffect) <- clean_names
cat("维度：", nrow(geneEffect), "细胞系 ×", ncol(geneEffect), "基因\n")

## 可选：去掉几乎所有细胞系里取值恒定的基因（方差极低 → 相关无意义）
## 课堂全量跑（q_var=0）；内存紧张时可设 q_var=0.1 砍掉最不活跃的 10%
if (q_var > 0) {
  v <- apply(geneEffect, 2, var, na.rm = TRUE)
  geneEffect <- geneEffect[, v > quantile(v, q_var, na.rm = TRUE)]
  cat("低方差过滤后：", ncol(geneEffect), "基因\n")
}

#############################################################
## 第 2 部分：基因×基因 共依赖 Pearson 相关矩阵
#############################################################
## 课堂原型（25Q3，1186×18435）约 4 分钟；26Q1 细胞系更多，5-10 分钟级。
## 结果存 RDS：这一步是全流程最贵的计算，后续任何 k/阈值调整都秒级复用。
cor_file <- file.path(output_dir, "20_geneEffect26Q1_cor_matrix.rds")
if (file.exists(cor_file)) {
  cat("发现已算好的相关矩阵，直接复用：", cor_file, "\n")
  M <- readRDS(cor_file)
} else {
  cat("计算全基因 Pearson 相关（一次矩阵运算，勿用循环）...\n")
  M <- cor(as.matrix(geneEffect), use = "pairwise.complete.obs")
  saveRDS(M, cor_file)
}
diag(M) <- 0   # 基因与自身相关 = 1，必须清掉，否则排名永远是自己第一

#############################################################
## 第 3 部分：mutual rank 长表（课堂 CK_modifier_rank.R 原样逻辑）
#############################################################
## Mrank[i,j] = 基因 j 在基因 i 眼中的 |rho| 排名（1 = 最相关）
## S = Mrank + t(Mrank)：S[i,j] = i 在 j 眼中的名次 + j 在 i 眼中的名次
##   S=2（1+1）即互为第一；k=1 的互为 top-1 筛选就是 S 最小的那些对
##   S 小 ≠ rho 大：一个 rho=0.9 排第 8 的关系可能输给两个 rho=0.6 互排第 2 的
cat("计算 mutual rank 并提取上三角长表...\n")
Mrank <- apply(M, 2, rank)          # 列方向排名（对每个基因给所有伙伴排序）
S <- Mrank + t(Mrank)
diag(S) <- 0
idx <- which(upper.tri(S), arr.ind = TRUE)   # 只取上三角，避免 (A,B)/(B,A) 重复

result <- data.frame(
  geneA = rownames(S)[idx[, "row"]],
  geneB = colnames(S)[idx[, "col"]],
  rankA = Mrank[cbind(idx[, "col"], idx[, "row"])],   # A 在 B 眼中的名次
  rankB = Mrank[cbind(idx[, "row"], idx[, "col"])],   # B 在 A 眼中的名次
  value = S[upper.tri(S)],                            # 互为排名和（mutual rank）
  rho   = M[cbind(idx[, "row"], idx[, "col"]])
)
rm(Mrank, S, idx); invisible(gc())

## 保存完整长表（下游 k/阈值任意切，不必重算）
saveRDS(result, file.path(output_dir, "20_mutualrank_pairs_full.rds"))

## 正/负相关分开导 top-k 榜（课堂 bipolar 口径：正负相关是两类不同生物学，
## 共同必需 → 正相关；合成致死/互斥 → 负相关，后者是找靶点的富矿）
topk <- 100
for (sgn in c("positive", "negative")) {
  sub <- if (sgn == "positive") result[result$rho >= rho_min, ] else result[result$rho <= -rho_min, ]
  sub <- sub[order(sub$value), ]
  fwrite(head(sub, topk * 5),
         file.path(output_dir, "20_modules", paste0("bipolar_", sgn, "_top", topk, ".csv")))
}

#############################################################
## 第 4 部分：k 网格 × 正/负相关 → 图 → Louvain 社区
#############################################################
## 每个尺度 k 下：
##   边 = 互为前 k（rankA<=k 且 rankB<=k）且 |rho|>=rho_min 且同号
##   注意"且"比课堂口述的"或"更严格：rankA<=k | rankB<=k 会把单向抱大腿
##   的边也放进图，k 大时整张图连成一片；mutual 的本意是"且"。
## 每个基因在每个 (方向, k) 组合下有唯一模块编号，汇总成 membership 长表。
cat("遍历 k 网格构建图并做 Louvain 社区发现...\n")
genes_all <- rownames(M)
membership_list <- list()   # [[paste(sgn, k)]] = data.frame(gene, module, k, sign)

for (sgn in c("positive", "negative")) {
  sub <- if (sgn == "positive") result[result$rho >= rho_min, ] else result[result$rho <= -rho_min, ]

  for (k in k_grid) {
    edges <- sub[sub$rankA <= k & sub$rankB <= k, c("geneA", "geneB", "rho")]
    if (nrow(edges) < 10) {
      membership_list[[paste(sgn, k)]] <-
        data.frame(gene = character(0), module = integer(0), k = k, sign = sgn)
      next
    }
    g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = genes_all))
    E(g)$weight <- abs(E(g)$rho)             # 边权 = |rho|，供 Louvain 使用
    comm <- cluster_louvain(g, weights = E(g)$weight)
    mem <- membership(comm)   # 孤立点各自成模块（负的编号），Louvain 模块为正整数
    membership_list[[paste(sgn, k)]] <-
      data.frame(gene = names(mem), module = as.integer(mem), k = k, sign = sgn)
    cat(sprintf("  %s k=%-3d 边数=%-7d 模块数=%d\n",
                sgn, k, nrow(edges), length(unique(mem[mem > 0]))))
  }
}
membership_all <- rbindlist(membership_list)
fwrite(membership_all, file.path(output_dir, "20_modules", "membership_all_k.csv"))

#############################################################
## 第 5 部分：跨尺度 Jaccard 稳定性 → 核心模块
#############################################################
## 以 ref_k（默认 10）的划分为基准：
##   对 ref_k 的每个模块 M0，在更大 k' 下找 Jaccard 最高的模块 M'，
##   stable(k') = Jaccard(M0, M')。
## 核心模块 = 在所有 k' >= ref_k 的尺度上 stable 都 >= stab_thre。
## 生物学解读：k 放大 = 允许更多外围基因加入；核心模块"吸人但不散架"，
## 说明它是一个真实的、自洽的功能单元，而不是阈值的人为产物。
cat("计算跨尺度 Jaccard 稳定性...\n")

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

ref_key <- paste("positive", ref_k)
ref_mem <- membership_list[[ref_key]]
ref_mods <- split(ref_mem$gene, ref_mem$module)
ref_mods <- ref_mods[ lengths(ref_mods) >= 5 ]   # 太小的模块无稳定性可言

bigger_keys <- paste("positive", k_grid[k_grid > ref_k])
stab_rows <- list()
for (mname in names(ref_mods)) {
  genes0 <- ref_mods[[mname]]
  row <- data.frame(module_ref = mname, size = length(genes0))
  for (key in bigger_keys) {
    mem2 <- membership_list[[key]]
    mem2 <- mem2[mem2$module > 0, ]
    mods2 <- split(mem2$gene, mem2$module)
    row[[paste0("stab_", sub("positive ", "k", key))]] <-
      max(vapply(mods2, function(x) jaccard(genes0, x), numeric(1)))
  }
  stab_rows[[mname]] <- row
}
stab <- rbindlist(stab_rows)
stab_cols <- grep("^stab_", colnames(stab), value = TRUE)
stab$core <- apply(stab[, ..stab_cols], 1, function(x) all(x >= stab_thre))
stab <- stab[order(-stab$core, -stab$size), ]
fwrite(stab, file.path(output_dir, "20_modules", "module_stability.csv"))

## 核心模块成员基因表（接 09-11 mutation→target 流水线的入口）
core_ids <- stab$core
core_genes <- do.call(rbind, lapply(names(core_ids)[core_ids], function(mname) {
  data.frame(module = mname, gene = ref_mods[[mname]])
}))
fwrite(core_genes, file.path(output_dir, "20_modules", "core_module_genes.csv"))
cat("核心（跨尺度稳定）模块数：", sum(core_ids), "/", nrow(stab), "\n")

#############################################################
## 第 6 部分：汇总图
#############################################################
## 图 1：k 放大过程中模块数量与边数的变化（"从小到大"的全景）
## 图 2：核心模块 vs 非核心模块的稳定性分布
pdf(file.path(output_dir, "20_modules", "20_summary_plots.pdf"), width = 9, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

n_mod <- sapply(k_grid, function(k) {
  m <- membership_list[[paste("positive", k)]]
  length(unique(m$module[m$module > 0]))
})
plot(k_grid, n_mod, type = "b", pch = 19, log = "x",
     xlab = "k（互为前 k）", ylab = "Louvain 模块数（正相关网络）",
     main = "k 放大：模块如何涌现")

if (sum(core_ids) >= 1) {
  boxplot(as.numeric(unlist(stab[core  == TRUE,  ..stab_cols])),
          as.numeric(unlist(stab[core  == FALSE, ..stab_cols])),
          names = c("core", "non-core"), ylim = c(0, 1),
          ylab = sprintf("Jaccard 稳定性（相对 k=%d）", ref_k),
          main = "跨尺度稳定 = 核心模块")
  abline(h = stab_thre, lty = 2)
}
dev.off()

cat("完成。输出目录：", file.path(output_dir, "20_modules"), "\n")
