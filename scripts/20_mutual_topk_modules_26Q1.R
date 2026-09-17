################################################
################################################
### 作者：果子（课堂原型 12_CO_IC/CK_modifier_rank.R、05_KD_modules_graph.R）
### 改编与扩展：乔宇（上海交通大学，免疫药理PhD）
### 创建时间：2026-09-17（v1）｜2026-09-17（v2 修订）
### 数据口径：DepMap 26Q1（课堂原型为 25Q3）

### ================================================================
### 主题：互为 top-k 共依赖网络（mutual top-k co-fitness）→ 跨尺度稳定模块发现
###
### v2 修订清单（对照"一句话科学思想"逐条兑现）：
###   [修1] mutual rank 排名方向：改为 |rho| 降序 rank(-abs(x))。
###         v1 照搬课堂 apply(M,2,rank) 是裸升序 = "最负排第一"，
###         正相关网络根本建不起来（有 7:01 截图 rho=+0.84 占 rank1-2
###         为证，他的真实口径是强度降序）。
###   [修2] 建边矩阵化：不再物化 1.7 亿行长表（20GB 级），直接
###         mutual <- (Mrank<=k) & t(Mrank<=k) 出边，k 网格才滑得动。
###   [修3] Louvain 随机性控制：每个 (符号,k) 跑 n_louvain_runs 次取
###         模块度最高的划分。不加这个，"模块不稳定"可能只是抽签
###         运气，Jaccard 判据被算法噪声污染。
###   [修4] 稳定性同时报告模块大小比 |M'|/|M0|：暴增=被稀释吞并、
###         暴减=碎裂，都算"散架"，只看 max(Jaccard) 太宽松。
###   [修5] 正、负相关网络各出一张稳定性表（bipolar 思想贯彻到模块层）。
###   [修6] 阳性对照自动化：对核心模块做 GMT 基因集超几何富集，
###         有显著注释=方法正确；无注释但跨尺度稳定=发现空间排序表
###         （本方法论的最终产出，接 09~11 做实验候选）。
###   [修7] 健壮性：重复 SYMBOL 去重、缓存文件名带基因数、NA 比例
###         检查决定 cor 的 use 策略。
###   [新]  节点中心性（SNA 指标，见下"中心性模块"说明）。

### ================================================================
### 科学思想（对应 TM00-13 课程"CoffeeTime / 从小到大发文章"）：
###   1. 两个基因的"关系"有多种强度口径。k=1（互为排名第一）排他性太强，
###      只能得到孤立的基因对；把 k 放大，基因对逐渐连成社区（模块）。
###   2. mutual rank（互为排名和 S = rank_A|B + rank_B|A）是对称化的
###      top-k 口径：一个基因很强但另一个很弱的关系会被 S 惩罚掉，
###      避免超级枢纽（如必需基因）把整张图吸成一颗星。
###   3. 跨尺度稳定性 = 生物学真实性的过滤器：真实的功能单元（复合体、
###      通路、共适应回路）在 k 从小到大的过程中模块组成基本不动；
###      只在大 k 才出现、随 k 剧烈变形的，多是稀释出来的噪声。
###   4. 结果必须自带阳性对照：已知模块能被捞回来 = 方法正确；
###      跨尺度稳定但注释未知/弱的模块 = 发现空间，接 mutation→target
###      流水线（09-11 脚本）做实验候选。

### ================================================================
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
###   两态：S 是软化形态（不设门槛，按两边排名之和全排序）；k 截断
###   建图是硬化形态（严格互近邻）。本脚本第 3 部分是软的，
###   第 4 部分是硬的。

### ================================================================
### 中心性模块（v2 新增）：把 SNA 节点指标引入基因网络，怎么选、怎么读
###   · 强度 strength（加权度）：基因"共依赖广度"的直接度量。高但
###     与 Common Essential 重叠者多为必需性假象，先排除再解读。
###   · 中介中心性 betweenness：桥接基因 = 连接两个模块的路径必经
###     点。若一头连已知模块、一头连未知核心模块，它就是"未知模块
###     如何接进已知生物学"的最佳入口，也是药理上"一点撬两面"的
###     候选靶点。⚠ 必须在分模块后的稀疏图上算；边权取 1/|rho| 作
###     距离（相似度转长度，不转会 prefer 弱关系，图论经典坑）。
###   · 特征向量中心性 eigen：近似"模块核心度"，找和重要基因扎堆
###     的基因。
###   计算口径：在 ref_k 尺度的正、负两张图上分别计算（稀疏、且
###   模块结构已成型），结果按模块归组输出。

### 这个脚本的科学思想，剥到底就是一句话：
###   基因间关系的强弱不是一条阈值线，而是一个可滑动的尺度 k；
###   在每个尺度下用互近邻图问一次"谁和谁抱团"，那么在所有尺度
###   下都稳定抱团的模块，大概率是真实的生物学功能单元——已知的
###   是方法正确的阳性对照，未知的是值得做实验的发现空间。
###   其余一切（SNA、Louvain、Jaccard、中心性）都是把这句话变成
###   可计算流程的工具。

### ================================================================
### 流程：
###   第 1 部分：读 26Q1 GeneEffect，清洗列名、去重、查 NA（同 07/09 口径）
###   第 2 部分：全基因 Pearson 共依赖矩阵（最贵一步，存 RDS 复用）
###   第 3 部分：mutual rank（|rho| 降序）+ bipolar top-k 软排序输出
###   第 4 部分：k 网格 × 正/负相关 → mutual k-NN 图 → Louvain（多次取优）
###   第 5 部分：跨尺度 Jaccard + 大小比稳定性 → 核心/碎裂/吞并分类
###   第 6 部分：GMT 基因集超几何富集 → 阳性对照 / 发现空间注释表
###   第 7 部分：节点中心性（strength / betweenness / eigen）@ ref_k
###   第 8 部分：汇总图
###
### 运行环境：
###   建议 ≥16 GB 内存（M + Mrank + S 各约 2.7 GB，峰值 ~8 GB）
###   R 包：data.table、igraph（首次使用需 install.packages("igraph")）
###   工作目录：仓库根目录（与 07/09 相同，相对路径 data/... TM00/...）
### ================================================================

rm(list = ls())
library(data.table)
library(igraph)

set.seed(20260917)   # Louvain 是随机算法，固定种子保证结果可复现

output_dir <- "TM00/DepMap_TM00/output"
mod_dir    <- file.path(output_dir, "20_modules")
dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)

### 可调参数（课堂演示值：k=20、rho_use=0.2）
k_grid          <- c(1, 2, 3, 5, 10, 20, 50, 100)  # k 尺度网格：从"互为第一"到宽松
rho_min         <- 0.2            # 边的最低 |rho| 门槛（课堂 rho_use=0.2）
ref_k           <- 10             # 稳定性/中心性参考尺度（模块已成型又不至于过稀释）
stab_thre       <- 0.7            # Jaccard ≥ 0.7 视为"该尺度下模块保持"
size_ratio_thre <- c(0.5, 2)      # 匹配模块大小比在此区间外 = 碎裂/吞并，不算稳定
min_module_size <- 5              # 参与稳定性检验的最小模块规模
n_louvain_runs  <- 3              # 每个 (符号,k) 的 Louvain 重复次数，取模块度最高者
q_var           <- 0              # 低方差基因过滤分位数（0 = 不过滤，与课堂全量对齐）
export_full_pairs <- FALSE        # 是否导出完整互排名长表（1.7 亿行，慎开）
gmt_file        <- NULL           # MSigDB/CRISPR 通路 GMT 路径；NULL = 跳过富集
                                  # 例："data/msigdb_hallmark.gmt"（gene Set2 基因\t基因...）

#############################################################
## 第 1 部分：读取 26Q1 CRISPR Gene Effect
#############################################################
## 数据结构（同 07 脚本口径）：
##   行 = 细胞系（ACH- 编号），列 = 基因（"SYMBOL (ENTREZ)" 需清洗）
##   值 = CERES/Chronos 效应评分，越负 = 越依赖
cat("读取 26Q1 Gene Effect...\n")
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = FALSE)
rownames(geneEffect) <- geneEffect[[1]]
geneEffect[[1]] <- NULL
clean_names <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))

## [修7a] 重复 SYMBOL 去重：ENTREZ 不同但 SYMBOL 相同的列会让 cor() 产出
## 重名列，下游按名字取值会静默出错。保留第一个，记录被丢弃者。
dup_flag <- duplicated(clean_names)
if (any(dup_flag)) {
  cat("发现重复基因名，丢弃", sum(dup_flag), "列：",
      paste(clean_names[dup_flag], collapse = ", "), "\n")
  geneEffect <- geneEffect[, !dup_flag, drop = FALSE]
  clean_names <- clean_names[!dup_flag]
}
colnames(geneEffect) <- clean_names
cat("维度：", nrow(geneEffect), "细胞系 ×", ncol(geneEffect), "基因\n")

## [修7c] NA 比例检查：26Q1 GeneEffect 通常干净；干净就用 everything
## （pairwise.complete.obs 在 18000² 规模上慢且可能产出病态相关）
na_frac <- mean(is.na(geneEffect))
cor_use <- if (na_frac < 1e-4) "everything" else "pairwise.complete.obs"
cat(sprintf("NA 比例：%.4f%%，cor use = %s\n", 100 * na_frac, cor_use))

## 可选：低方差基因过滤（方差极低 → 相关无意义）
if (q_var > 0) {
  v <- apply(geneEffect, 2, var, na.rm = TRUE)
  geneEffect <- geneEffect[, v > quantile(v, q_var, na.rm = TRUE)]
  cat("低方差过滤后：", ncol(geneEffect), "基因\n")
}

#############################################################
## 第 2 部分：基因×基因 共依赖 Pearson 相关矩阵
#############################################################
## 课堂原型（25Q3，1186×18435）约 4 分钟；26Q1 细胞系更多，5-10 分钟级。
## [修7b] 缓存文件名带基因数：改 q_var / 基因集后不会误读旧缓存。
n_genes <- ncol(geneEffect)
cor_file <- file.path(output_dir, sprintf("20_geneEffect26Q1_cor_matrix_g%d.rds", n_genes))
if (file.exists(cor_file)) {
  cat("发现已算好的相关矩阵，直接复用：", cor_file, "\n")
  M <- readRDS(cor_file)
} else {
  cat("计算全基因 Pearson 相关（一次矩阵运算，勿用循环）...\n")
  M <- cor(as.matrix(geneEffect), use = cor_use)
  saveRDS(M, cor_file)
  cat("相关矩阵已缓存：", cor_file, "\n")
}
rm(geneEffect); invisible(gc())
diag(M) <- 0   # 基因与自身相关 = 1，必须清掉，否则排名永远是自己第一

#############################################################
## 第 3 部分：mutual rank（|rho| 降序）+ bipolar top-k
#############################################################
## [修1] 排名方向：rank(-abs(x)) → |rho| 越大名次越靠前，正负同权。
## 课堂 apply(M,2,rank) 裸升序会把"最负"排第一，与 bipolar 正负分网的
## 设计自相矛盾（正相关网络会一张边都不剩）。
## 符号不参与排名，只留给第 4 部分分图用。
cat("计算 mutual rank（|rho| 降序）...\n")
Mrank <- apply(M, 2, function(x) rank(-abs(x)))
S <- Mrank + t(Mrank)      # S[i,j] = i 在 j 眼中名次 + j 在 i 眼中名次，越小越互为靠前
diag(S) <- 0

## 软排序输出（bipolar top-k）：正/负相关分开，按 S 升序取前列。
## 注意这里只物化 top 候选，不做全量长表（[修2]，全量 = 1.7 亿行）。
topk <- 100
for (sgn in c("positive", "negative")) {
  sign_mask <- if (sgn == "positive") M >= rho_min else M <= -rho_min
  cand <- S * sign_mask          # 不满足符号条件的 S 置 0（S 恒 >= 2，不会混入）
  cand[lower.tri(cand)] <- 0     # 只取上三角
  n_take <- min(topk * 5, sum(cand > 0))
  take_idx <- which(cand > 0, arr.ind = TRUE)[seq_len(n_take), , drop = FALSE]
  take_ord <- order(cand[take_idx])
  take_idx <- take_idx[take_ord, , drop = FALSE]
  out <- data.frame(
    geneA = rownames(M)[take_idx[, "row"]],
    geneB = colnames(M)[take_idx[, "col"]],
    rankA = Mrank[cbind(take_idx[, "row"], take_idx[, "col"])],  # A 在 B 眼中名次
    rankB = Mrank[cbind(take_idx[, "col"], take_idx[, "row"])],  # B 在 A 眼中名次
    value = cand[take_idx],
    rho   = M[take_idx]
  )
  fwrite(out, file.path(mod_dir, paste0("bipolar_", sgn, "_top", topk, ".csv")))
}
rm(S); invisible(gc())

#############################################################
## 第 4 部分：k 网格 × 正/负相关 → mutual k-NN 图 → Louvain
#############################################################
## [修2] 全程不物化长表：互为前 k 直接用对称逻辑矩阵判定，
## (Mrank<=k) & t(Mrank<=k) 就是 mutual k-NN 邻接矩阵的定义式。
## [修3] Louvain 跑 n_louvain_runs 次取模块度最高划分，把算法
## 随机性和尺度扰动解耦——第 5 部分的 Jaccard 只反映尺度效应。
cat("遍历 k 网格构建图并做 Louvain 社区发现...\n")
membership_list <- list()   # [[paste(sgn, k)]] = data.frame(gene, module)

for (sgn in c("positive", "negative")) {
  for (k in k_grid) {
    key <- paste(sgn, k)
    topk_flag <- Mrank <= k
    mutual <- topk_flag & t(topk_flag)
    diag(mutual) <- FALSE
    idx <- which(mutual & upper.tri(mutual), arr.ind = TRUE)
    if (nrow(idx) < 10) {
      membership_list[[key]] <- data.frame(gene = character(0), module = integer(0))
      cat(sprintf("  %s k=%-3d 边数不足，跳过\n", sgn, k)); next
    }
    rho_e <- M[idx]
    keep <- if (sgn == "positive") rho_e >= rho_min else rho_e <= -rho_min
    edges <- data.frame(
      geneA = rownames(M)[idx[keep, 1]],
      geneB = colnames(M)[idx[keep, 2]],
      rho   = rho_e[keep]
    )
    g <- graph_from_data_frame(edges, directed = FALSE)
    E(g)$weight <- abs(E(g)$rho)
    best <- NULL
    for (r in seq_len(n_louvain_runs)) {
      comm <- cluster_louvain(g, weights = E(g)$weight)
      if (is.null(best) || modularity(comm) > best$mod) {
        best <- list(mod = modularity(comm),
                     mem = as.integer(membership(comm)),
                     names = names(membership(comm)))
      }
    }
    membership_list[[key]] <- data.frame(gene = best$names, module = best$mem)
    cat(sprintf("  %s k=%-3d 边数=%-7d 模块数=%d（Q=%.3f，%d 次取优）\n",
                sgn, k, nrow(edges),
                length(unique(best$mem)), best$mod, n_louvain_runs))
  }
}
membership_all <- data.table::rbindlist(lapply(names(membership_list), function(key) {
  df <- membership_list[[key]]
  if (nrow(df) == 0) return(df)
  sgn_k <- strsplit(key, " ")[[1]]
  cbind(df, sign = sgn_k[1], k = as.integer(sgn_k[2]))
}))
fwrite(membership_all, file.path(mod_dir, "membership_all_k.csv"))

#############################################################
## 第 5 部分：跨尺度稳定性（Jaccard + 大小比）→ 核心模块
#############################################################
## [修4][修5] 对每个符号的网络，以 ref_k 划分为基准：
##   stable(k') = max_Jaccard(M0, M')，size_ratio = |M'|/|M0|
##   核心 = 所有 k' > ref_k 上 Jaccard ≥ thre 且大小比在区间内。
##   暴增 = 被稀释吞并；暴减 = 碎裂。两者都算"散架"。
cat("计算跨尺度 Jaccard 稳定性...\n")
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

stability_all <- list()
for (sgn in c("positive", "negative")) {
  ref_mem <- membership_list[[paste(sgn, ref_k)]]
  if (nrow(ref_mem) == 0) next
  ref_mods <- split(ref_mem$gene, ref_mem$module)
  ref_mods <- ref_mods[lengths(ref_mods) >= min_module_size]
  bigger_keys <- paste(sgn, k_grid[k_grid > ref_k])
  rows <- list()
  for (mname in names(ref_mods)) {
    genes0 <- ref_mods[[mname]]
    row <- data.frame(module_ref = mname, size = length(genes0))
    for (key in bigger_keys) {
      kk <- strsplit(key, " ")[[1]][2]
      mem2 <- membership_list[[key]]
      if (nrow(mem2) == 0) { row[[paste0("stab_k", kk)]] <- NA; next }
      mods2 <- split(mem2$gene, mem2$module)
      js <- vapply(mods2, function(x) jaccard(genes0, x), numeric(1))
      best_i <- which.max(js)
      row[[paste0("stab_k", kk)]] <- js[best_i]
      row[[paste0("ratio_k", kk)]] <- length(mods2[[best_i]]) / length(genes0)
    }
    rows[[mname]] <- row
  }
  stab <- do.call(rbind, rows)
  if (is.null(stab)) next
  stab_cols <- grep("^stab_", colnames(stab), value = TRUE)
  ratio_cols <- grep("^ratio_", colnames(stab), value = TRUE)
  stab$core <- apply(stab[, stab_cols, drop = FALSE], 1, function(x)
    all(!is.na(x) & x >= stab_thre)) &&
    apply(stab[, ratio_cols, drop = FALSE], 1, function(x)
      all(!is.na(x) & x >= size_ratio_thre[1] & x <= size_ratio_thre[2]))
  stab <- stab[order(-stab$core, -stab$size), ]
  stab$sign <- sgn
  stability_all[[sgn]] <- stab

  ## 核心模块成员表（接 09-11 mutation→target 流水线的入口）
  core_genes <- do.call(rbind, lapply(stab$module_ref[stab$core], function(mname) {
    data.frame(module = mname, gene = ref_mods[[mname]])
  }))
  if (!is.null(core_genes)) {
    fwrite(core_genes, file.path(mod_dir, paste0("core_module_genes_", sgn, ".csv")))
  }
  cat(sprintf("[%s] 参考模块 %d 个，核心（跨尺度稳定）%d 个\n",
              sgn, nrow(stab), sum(stab$core)))
}
saveRDS(stability_all, file.path(mod_dir, "stability_all.rds"))

#############################################################
## 第 6 部分：GMT 超几何富集 → 阳性对照 / 发现空间注释
#############################################################
## 思想的最后一块拼图：核心模块里"有显著注释的"= 阳性对照（方法正确），
## "无注释但跨尺度稳定的"= 发现空间排序表。gmt_file 为 NULL 时跳过。
## GMT 格式：每行 "set_name\t来源说明\t基因1\t基因2..."（MSigDB 标准）。
for (sgn in names(stability_all)) {
  if (is.null(gmt_file) || !file.exists(gmt_file)) {
    cat("[", sgn, "] 未提供 GMT，跳过富集（阳性对照未检验）\n"); next
  }
  lines <- readLines(gmt_file)
  sets <- lapply(strsplit(lines, "\t"), function(x) unique(x[-c(1, 2)]))
  names(sets) <- vapply(strsplit(lines, "\t"), function(x) x[1], character(1))

  ref_mem <- membership_list[[paste(sgn, ref_k)]]
  ref_mods <- split(ref_mem$gene, ref_mem$module)
  ref_mods <- ref_mods[lengths(ref_mods) >= min_module_size]
  stab_df <- stability_all[[sgn]]
  universe <- unique(unlist(ref_mods))
  ann_rows <- list()
  for (mname in names(ref_mods)) {
    genes_m <- ref_mods[[mname]]
    is_core <- stab_df$core[match(mname, stab_df$module_ref)]
    for (s in names(sets)) {
      hit_m <- sum(genes_m %in% sets[[s]])
      if (hit_m < 2) next
      p <- phyper(hit_m - 1, length(intersect(sets[[s]], universe)),
                  length(universe) - length(intersect(sets[[s]], universe)),
                  length(genes_m), lower.tail = FALSE)
      ann_rows[[length(ann_rows) + 1]] <- data.frame(
        module = mname, core = is_core, set = s,
        overlap = hit_m, module_size = length(genes_m),
        set_in_universe = length(intersect(sets[[s]], universe)), p = p)
    }
  }
  if (length(ann_rows)) {
    ann <- do.call(rbind, ann_rows)
    ann$p_BH <- p.adjust(ann$p, "BH")
    ann <- ann[order(ann$p), ]
    fwrite(ann, file.path(mod_dir, paste0("module_annotation_", sgn, ".csv")))
    n_core_ann <- length(unique(ann$module[ann$core & ann$p_BH < 0.05]))
    cat(sprintf("[%s] 核心模块中有显著注释（BH<0.05）的：%d 个 = 阳性对照\n",
                sgn, n_core_ann))
  }
}

#############################################################
## 第 7 部分：节点中心性（SNA 指标）@ ref_k
#############################################################
## 在正、负两张 ref_k 图上分别计算（稀疏、模块结构已成型）：
##   strength    = 加权度，共依赖广度（先与 Common Essential 对照再解读）
##   betweenness = 桥接基因，边权取 1/|rho| 作距离（相似度转长度！）
##   eigen       = 特征向量中心性 ≈ 模块核心度
## 高 betweenness 且横跨"已知注释模块 ↔ 无注释核心模块"的基因，
## 是"未知模块如何接进已知生物学"的最佳入口候选。
cat("计算节点中心性 @ ref_k =", ref_k, "...\n")
if (requireNamespace("igraph", quietly = TRUE)) {
  for (sgn in c("positive", "negative")) {
    k <- ref_k
    topk_flag <- Mrank <= k
    mutual <- topk_flag & t(topk_flag); diag(mutual) <- FALSE
    idx <- which(mutual & upper.tri(mutual), arr.ind = TRUE)
    if (nrow(idx) < 10) next
    rho_e <- M[idx]
    keep <- if (sgn == "positive") rho_e >= rho_min else rho_e <= -rho_min
    edges <- data.frame(geneA = rownames(M)[idx[keep, 1]],
                        geneB = colnames(M)[idx[keep, 2]],
                        rho   = rho_e[keep])
    g <- graph_from_data_frame(edges, directed = FALSE)
    E(g)$weight <- abs(E(g)$rho)
    mem <- membership_list[[paste(sgn, k)]]

    node_metrics <- data.frame(
      gene        = V(g)$name,
      sign        = sgn,
      module      = mem$module[match(V(g)$name, mem$gene)],
      degree      = degree(g),
      strength    = strength(g, weights = E(g)$weight),
      betweenness = betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE),
      eigen       = eigen_centrality(g, weights = E(g)$weight)$vector
    )
    node_metrics <- node_metrics[order(-node_metrics$betweenness), ]
    fwrite(node_metrics, file.path(mod_dir, paste0("node_centrality_", sgn, "_k", k, ".csv")))
    cat(sprintf("[%s] 节点 %d，betweenness top5：%s\n", sgn, nrow(node_metrics),
                paste(head(node_metrics$gene, 5), collapse = ", ")))
  }
}

#############################################################
## 第 8 部分：汇总图
#############################################################
pdf(file.path(mod_dir, "20_summary_plots.pdf"), width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

## 图 1：k 放大过程中模块数量的变化（"从小到大"的全景，正负并排）
for (sgn in c("positive", "negative")) {
  n_mod <- sapply(k_grid, function(k) {
    m <- membership_list[[paste(sgn, k)]]
    if (is.null(m) || nrow(m) == 0) 0 else length(unique(m$module))
  })
  if (sgn == "positive") {
    plot(k_grid, n_mod, type = "b", pch = 19, log = "x", col = "firebrick",
         xlab = "k（互为前 k）", ylab = "Louvain 模块数",
         main = "k 放大：模块如何涌现")
  } else {
    points(k_grid, n_mod, type = "b", pch = 17, col = "steelblue")
  }
}
legend("topright", c("正相关", "负相关"), col = c("firebrick", "steelblue"),
       pch = c(19, 17), lty = 1)

## 图 2：核心 vs 非核心稳定性（正相关网络）
if (!is.null(stability_all$positive)) {
  stab <- stability_all$positive
  stab_cols <- grep("^stab_", colnames(stab), value = TRUE)
  boxplot(as.numeric(unlist(stab[stab$core == TRUE,  stab_cols])),
          as.numeric(unlist(stab[stab$core == FALSE, stab_cols])),
          names = c("core", "non-core"), ylim = c(0, 1), col = c("gold", "grey"),
          ylab = sprintf("Jaccard 稳定性（相对 k=%d）", ref_k),
          main = "跨尺度稳定 = 核心模块")
  abline(h = stab_thre, lty = 2)
}

## 图 3：桥接基因（正相关网络 betweenness top20）
bw_file <- file.path(mod_dir, paste0("node_centrality_positive_k", ref_k, ".csv"))
if (file.exists(bw_file)) {
  nm <- fread(bw_file)
  top <- head(nm[order(-betweenness)], 20)
  barplot(rev(setNames(top$betweenness, top$gene)), horiz = TRUE, las = 1,
          col = "darkorange", xlab = "normalized betweenness",
          main = "桥接基因 top20（正相关）")
}
dev.off()

cat("完成。输出目录：", mod_dir, "\n")
