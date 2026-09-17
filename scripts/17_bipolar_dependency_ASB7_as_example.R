################################################
################################################
### 作者：果子
### 维护者：乔宇（上海交通大学，免疫药理PhD）
### 更新时间：2026-09-08
### 数据迁移：2026-08-12（23Q2 → 26Q1）

### ================================================================
### 主题：双极依赖（Bipolar Dependency）分析
### 复现文章: Science: ASB7 is a negative regulator of H3K9me3 homeostasis
###
### 核心概念：
###   "双极依赖" = 两个基因的依赖性呈反向相关
###   例如 ASB7 和 SUV39H1：
###     - ASB7 高依赖的细胞系，SUV39H1 反而不依赖（负相关）
###     - 这种反向关系提示两者处于同一调控通路的两端
###
### 分析流程：
###   1. 加载 Gene Effect 数据
###   2. ASB7 vs SUV39H1 散点图验证（ggscatterstats）
###   3. 批量计算每个基因与 ASB7/SUV39H1 的相关性
###   4. 复现论文"条带图"（strip plot）—— Top100 相关基因的分布
###   5. 扩展分析：ASB7 负相关基因的递归搜索
###   6. TP53/MDM2 基因对验证
### ================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

output_dir <- "TM00/DepMap_TM00/output"

#############################################################
## 第一部分：数据加载
#############################################################

## 读取 Gene Effect（用全量 1208 细胞系，不做表达量交集过滤）
cat("正在读取 Gene Effect...\n")
geneEffect <- data.table::fread("data/CRISPRGeneEffect.csv", data.table = F)
rownames(geneEffect) <- geneEffect[, 1]     ## 第一列 V1 = ACH-编号
geneEffect <- geneEffect[, -1]
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)", "", colnames(geneEffect))
geneEffect <- as.matrix(geneEffect)
cat("  Gene Effect:", nrow(geneEffect), "模型 ×", ncol(geneEffect), "基因\n")

#############################################################
## 第二部分：ASB7 vs SUV39H1 散点图验证
#############################################################

cat("\n=== ASB7 vs SUV39H1 相关性 ===\n")
ct <- cor.test(geneEffect[, "ASB7"], geneEffect[, "SUV39H1"])
cat("  r =", round(ct$estimate, 4), " p =", format(ct$p.value, digits = 3), "\n")

## 用 ggstatsplot 画带统计信息的散点图
library(ggstatsplot)
p_scatter <- ggscatterstats(
  data  = as.data.frame(geneEffect),
  x     = SUV39H1,
  y     = ASB7,
  xlab  = "SUV39H1 Gene Effect",
  ylab  = "ASB7 Gene Effect",
  digits = 4L
)
print(p_scatter)
ggsave(file.path(output_dir, "17_ASB7_SUV39H1_scatterstats.png"),
       p_scatter, width = 8, height = 6, dpi = 150)
cat("  已保存: 17_ASB7_SUV39H1_scatterstats.png\n")

#############################################################
## 第三部分：批量相关性计算函数
#############################################################

## 封装函数：计算某基因与所有其他基因的相关性
## 返回 data.frame: geneA, geneB, cor, pvalue
batch_cor <- function(geneA, geneEffect_mat, show_progress = TRUE) {
  geneAdata <- as.numeric(geneEffect_mat[, geneA])
  geneNames <- colnames(geneEffect_mat)

  if (show_progress) {
    library(pbapply)
    pbo <- pboptions(type = "timer", char = "=")
  }

  result <- do.call(rbind, lapply(geneNames, function(i) {
    dd <- tryCatch(
      cor.test(geneAdata, geneEffect_mat[, i]),
      error = function(e) NULL
    )
    if (is.null(dd)) {
      data.frame(geneA = geneA, geneB = i, cor = NA_real_, pvalue = NA_real_)
    } else {
      data.frame(geneA = geneA, geneB = i, cor = dd$estimate, pvalue = dd$p.value)
    }
  }))

  if (show_progress) pboptions(pbo)
  return(result)
}

#############################################################
## 第四部分：ASB7 和 SUV39H1 的批量相关性
#############################################################

cat("\n=== 批量计算 SUV39H1 相关性 ===\n")
corData1 <- batch_cor("SUV39H1", geneEffect)
saveRDS(corData1, file.path(output_dir, "17_corData_SUV39H1.rds"))

cat("\n=== 批量计算 ASB7 相关性 ===\n")
corData2 <- batch_cor("ASB7", geneEffect)
saveRDS(corData2, file.path(output_dir, "17_corData_ASB7.rds"))

#############################################################
## 第五部分：复现论文"条带图"（strip plot）
#############################################################

library(cowplot)
library(grid)

## 条带图绘图函数（支持自适应坐标 & 右侧标签）
strip_plot <- function(df, target_gene,
                       bg_fill, pt_fill, right_label,
                       xlim_vec, brks, labs, r = .12) {

  targ_x <- df$cor[df$geneB == target_gene]
  ## 如果 target_gene 不在 Top100 中，取全量相关性中该基因的值
  if (length(targ_x) == 0) {
    cat("  注意:", target_gene, "不在 Top100 中，跳过箭头标注\n")
    targ_x <- NA
  }

  p <- ggplot(df, aes(cor, 0)) +
    ## 背景圆角矩形
    annotation_custom(
      roundrectGrob(x = .5, y = .5, width = 1, height = 1,
                    r = unit(r, "snpc"),
                    gp = gpar(fill = bg_fill, col = NA)),
      xmin = xlim_vec[1] - .02, xmax = xlim_vec[2] + .02,
      ymin = -.5, ymax =  .5) +
    ## 主轴（用 annotate 避免函数作用域问题）
    annotate("segment", x = xlim_vec[1], xend = xlim_vec[2], y = 0, yend = 0,
             linewidth = .9) +
    ## 刻度线
    annotate("segment", x = brks, xend = brks, y = 0, yend = -.08,
             linewidth = .8) +
    annotate("text", x = brks, y = -.18, label = labs, size = 4) +
    ## 数据点
    geom_point(shape = 21, size = 5, stroke = .35,
               fill = pt_fill, colour = "black")

  ## 仅当 targ_x 有效时添加箭头和标签
  if (!is.na(targ_x[1])) {
    p <- p +
      annotate("segment", x = targ_x, xend = targ_x, y = .28, yend = .06,
               arrow = arrow(length = unit(.24, "cm")), linewidth = 1) +
      annotate("text", x = targ_x, y = .36, label = target_gene,
               fontface = "italic")
  }

  p <- p +
    ## 右侧纵向标签
    annotate("text", x = Inf, y = 0, label = right_label,
             angle = -90, hjust = 0.5, vjust = 0.2,
             fontface = "italic") +
    scale_x_continuous(limits = xlim_vec,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(-.5, .5), clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      plot.title  = element_text(face = "italic", size = 16,
                                 hjust = 0, margin = margin(b = 2)),
      plot.margin = margin(2, 20, 2, 24)
    )
  return(p)
}

## ---- ASB7 / SUV39H1 条带图 ----

## 取 |cor| Top 100（去掉自身 = 第 1 行）
geneA_df <- arrange(corData1, desc(abs(cor)))[2:101, ]
geneB_df <- arrange(corData2, desc(abs(cor)))[2:101, ]

## 自动计算 x 轴范围
xl <- floor(min(c(geneA_df$cor, geneB_df$cor)) * 10) / 10
xu <- ceiling(max(c(geneA_df$cor, geneB_df$cor)) * 10) / 10
axis_breaks <- seq(xl, xu, by = 0.1)
axis_labels <- sprintf("%.1f", axis_breaks)

## 色彩参数
geneA_col <- "#E76F51"; geneB_col <- "#457B9D"
bg_red    <- "#F9D9D3"; bg_blue   <- "#D6EAF8"

## 绘制上下两张条带
p1 <- strip_plot(geneA_df, target_gene = "ASB7",
                 bg_red, geneA_col, right_label = "SUV39H1",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)

p2 <- strip_plot(geneB_df, target_gene = "SUV39H1",
                 bg_blue, geneB_col, right_label = "ASB7",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)

## 主标题栏
title_bar <- ggdraw() +
  draw_label("Correlation", fontface="bold", size=19, x=.5, y=.8, hjust=.5) +
  draw_line(x=c(.43,.17), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Inverse", x=.13, y=.55, size=11) +
  draw_line(x=c(.57,.83), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Direct",  x=.87, y=.55, size=11)

## 拼图输出
panel <- plot_grid(p1, p2, ncol=1, align="v")
p_asb7 <- plot_grid(title_bar, panel, ncol=1, rel_heights=c(.16,1))
print(p_asb7)
ggsave(file.path(output_dir, "17_ASB7_SUV39H1_strip_plot.png"),
       p_asb7, width = 8, height = 5, dpi = 150)
cat("已保存: 17_ASB7_SUV39H1_strip_plot.png\n")

#############################################################
## 第六部分：扩展 —— ASB7 负相关基因的递归搜索
#############################################################

## 在 ASB7 最负相关的 Top10 基因中，寻找是否还有类似 ASB7 的"双极依赖"模式
cat("\n=== ASB7 负相关基因递归搜索 ===\n")
corData2_sorted <- dplyr::arrange(corData2, cor)
genesTocheck <- corData2_sorted$geneB[1:10]
cat("  待检查基因:", paste(genesTocheck, collapse = ", "), "\n")

corData3 <- do.call(rbind, lapply(genesTocheck, function(myGene) {
  cat("  正在计算:", myGene, "...\n")
  geneAdata <- as.numeric(geneEffect[, myGene])

  mycorData <- do.call(rbind, lapply(colnames(geneEffect), function(i) {
    dd <- tryCatch(
      cor.test(geneAdata, geneEffect[, i]),
      error = function(e) NULL
    )
    if (is.null(dd)) {
      data.frame(geneA = myGene, geneB = i, cor = NA_real_, pvalue = NA_real_)
    } else {
      data.frame(geneA = myGene, geneB = i, cor = dd$estimate, pvalue = dd$p.value)
    }
  }))

  mycorData <- dplyr::arrange(mycorData, cor)
  mycorData <- mycorData[1:10, ]                  ## Top10 最负相关
  mycorData$isGeneAin <- mycorData$geneB %in% "ASB7"  ## ASB7 是否在其中
  return(mycorData)
}))
saveRDS(corData3, file.path(output_dir, "17_corData3_ASB7_recursive.rds"))
cat("  已保存: 17_corData3_ASB7_recursive.rds\n")

## SUV39H1 负相关基因递归搜索
cat("\n=== SUV39H1 负相关基因递归搜索 ===\n")
corData1_sorted <- dplyr::arrange(corData1, cor)
genesTocheck <- corData1_sorted$geneB[1:10]
cat("  待检查基因:", paste(genesTocheck, collapse = ", "), "\n")

corData4 <- do.call(rbind, lapply(genesTocheck, function(myGene) {
  cat("  正在计算:", myGene, "...\n")
  geneAdata <- as.numeric(geneEffect[, myGene])

  mycorData <- do.call(rbind, lapply(colnames(geneEffect), function(i) {
    dd <- tryCatch(
      cor.test(geneAdata, geneEffect[, i]),
      error = function(e) NULL
    )
    if (is.null(dd)) {
      data.frame(geneA = myGene, geneB = i, cor = NA_real_, pvalue = NA_real_)
    } else {
      data.frame(geneA = myGene, geneB = i, cor = dd$estimate, pvalue = dd$p.value)
    }
  }))

  mycorData <- dplyr::arrange(mycorData, cor)
  mycorData <- mycorData[1:10, ]
  mycorData$isGeneAin <- mycorData$geneB %in% "SUV39H1"
  return(mycorData)
}))
saveRDS(corData4, file.path(output_dir, "17_corData4_SUV39H1_recursive.rds"))
cat("  已保存: 17_corData4_SUV39H1_recursive.rds\n")

#############################################################
## 第七部分：TP53 / MDM2 基因对验证
#############################################################

cat("\n=== TP53 批量相关性 ===\n")
corData5 <- batch_cor("TP53", geneEffect)
saveRDS(corData5, file.path(output_dir, "17_corData_TP53.rds"))

cat("\n=== MDM2 批量相关性 ===\n")
corData6 <- batch_cor("MDM2", geneEffect)
saveRDS(corData6, file.path(output_dir, "17_corData_MDM2.rds"))

## ---- TP53 / MDM2 条带图 ----

geneA_df <- arrange(corData5, desc(abs(cor)))[2:101, ]
geneB_df <- arrange(corData6, desc(abs(cor)))[2:101, ]

## 自动计算 x 轴范围
xl <- floor(min(c(geneA_df$cor, geneB_df$cor)) * 10) / 10
xu <- ceiling(max(c(geneA_df$cor, geneB_df$cor)) * 10) / 10
axis_breaks <- seq(xl, xu, by = 0.1)
axis_labels <- sprintf("%.1f", axis_breaks)

p1 <- strip_plot(geneA_df, target_gene = "MDM2",
                 bg_red, geneA_col, right_label = "TP53",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)

p2 <- strip_plot(geneB_df, target_gene = "TP53",
                 bg_blue, geneB_col, right_label = "MDM2",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)

## 主标题栏
title_bar <- ggdraw() +
  draw_label("Correlation", fontface="bold", size=19, x=.5, y=.8, hjust=.5) +
  draw_line(x=c(.43,.17), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Inverse", x=.13, y=.55, size=11) +
  draw_line(x=c(.57,.83), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Direct",  x=.87, y=.55, size=11)

panel <- plot_grid(p1, p2, ncol=1, align="v")
p_tp53 <- plot_grid(title_bar, panel, ncol=1, rel_heights=c(.16,1))
print(p_tp53)
ggsave(file.path(output_dir, "17_TP53_MDM2_strip_plot.png"),
       p_tp53, width = 8, height = 5, dpi = 150)
cat("已保存: 17_TP53_MDM2_strip_plot.png\n")

cat("\n========== 17 脚本全部完成 ==========\n")

### 拓展：
### Depmap 功能相关性 → 共依赖分析
### Depmap In silico CRISPR/screen → 合成致死筛选
