################################################
################################################
### 作者：果子
### 更新时间：2026-09-09
### 微信公众号:果子学生信
### 私人微信：guotosky
### 复现文章: Scince: ASB7 is a negative regulator of H3K9me3 homeostasis

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
### 把Depmap数据读取进来
## Gene effect
geneEffect <- data.table::fread("data/DepMap_Public_25Q3/CRISPRGeneEffect.csv",data.table = F)
test <- geneEffect[1:10,1:10]
rownames(geneEffect) <- geneEffect[,1]
geneEffect <- geneEffect[,-1]
test <- geneEffect[1:10,1:10]
colnames(geneEffect) <- gsub("\\s+\\(\\d+\\)","",colnames(geneEffect))
### 1186 18435
saveRDS(geneEffect,file = "output/depmap_geneEffect25Q3.rds")

### 4min
M = cor(geneEffect)
saveRDS(M,file = "output/depmap_geneEffect25Q3_cor_matrix_20251005.rds")

M <- readRDS(file = "output/depmap_geneEffect25Q3_cor_matrix_20251005.rds")

### 提取负相关,最小的rank会相加,总和越小越好 2min
Mrank = apply(M,2,rank)
S <- Mrank + t(Mrank)          # S[i,j] = M[i,j] + M[j,i]
diag(S) <- 0
# 上三角索引 (不含对角线) 的行列坐标
idx <- which(upper.tri(S), arr.ind = TRUE)   # matrix of (row, col)
# 构建结果 data.frame：geneA, geneB, value
result <- data.frame(
  geneA = rownames(S)[idx[,"row"]],
  geneB = colnames(S)[idx[,"col"]],
  rankA = Mrank[cbind(idx[, "col"], idx[, "row"])],
  rankB = Mrank[cbind(idx[, "row"], idx[, "col"])],
  value = S[upper.tri(S)],
  rho   = M[cbind(idx[, "row"], idx[, "col"])],
  row.names = NULL,
  stringsAsFactors = FALSE
)

result2 = result[result$rankA<=20&result$rankB <=20,]
saveRDS(result2,file = "output/depmap_geneEffect25Q3_biopolar_result_negative_top20.rds")

######### 提取正相关
rm(list = ls())
M <- readRDS(file = "output/depmap_geneEffect25Q3_cor_matrix_20251005.rds")
### 排除自相关
diag(M) <- NA
### 提取正相关,矩阵乘以-1,把相关性系数颠倒,最小的rank会相加,总和越小越好
Mrank = apply(-M,2,rank)
S <- Mrank + t(Mrank)          # S[i,j] = M[i,j] + M[j,i]
diag(S) <- 0
# 上三角索引 (不含对角线) 的行列坐标
idx <- which(upper.tri(S), arr.ind = TRUE)   # matrix of (row, col)
# 构建结果 data.frame：geneA, geneB, value
result <- data.frame(
  geneA = rownames(S)[idx[,"row"]],
  geneB = colnames(S)[idx[,"col"]],
  rankA = Mrank[cbind(idx[, "col"], idx[, "row"])],
  rankB = Mrank[cbind(idx[, "row"], idx[, "col"])],
  value = S[upper.tri(S)],
  rho   = M[cbind(idx[, "row"], idx[, "col"])],
  row.names = NULL,
  stringsAsFactors = FALSE
)

result2 = result[result$rankA<=20&result$rankB <=20,]
saveRDS(result2,file = "output/depmap_geneEffect25Q3_biopolar_result_positive_top20.rds")

result3 = result[abs(result$rankA)>=0.1,]
################################################################
rm(list = ls())
geneEffect = readRDS(file = "output/depmap_geneEffect25Q3.rds")
biopolar_neg_top20 = readRDS(file = "output/depmap_geneEffect25Q3_biopolar_result_negative_top20.rds")
biopolar_pos_top20 = readRDS(file = "output/depmap_geneEffect25Q3_biopolar_result_positive_top20.rds")

data = dplyr::arrange(biopolar_neg_top20,value,rho)
data = data[data$value==2& data$rho <= -0.3,]

dput(paste(data$geneA,data$geneB,sep=":"))


library(pbapply)

### 画图
geneA = "TP53"
geneAdata = as.numeric(geneEffect[,geneA])
corData5 <- do.call(rbind, pblapply(colnames(geneEffect),function(i){
  dd = cor.test(geneAdata ,geneEffect[,i])
  data.frame(geneA = geneA,geneB =i,cor = dd$estimate,pvalue = dd$p.value)
  
}))

#########
geneA = "MDM2"
geneAdata = as.numeric(geneEffect[,geneA])
corData6 <- do.call(rbind, pblapply(colnames(geneEffect),function(i){
  dd = cor.test(geneAdata ,geneEffect[,i])
  data.frame(geneA = geneA,geneB =i,cor = dd$estimate,pvalue = dd$p.value)
  
}))


### 画图
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)

## ── 1. 数据准备 ────────────────────────────────────────────
geneA_df  <- arrange(corData5, desc(abs(cor)))[2:101, ]
geneB_df  <- arrange(corData6, desc(abs(cor)))[2:101, ]


## ── 2. 根据两组数据范围自动计算 x 轴 ──────────────────────
xl <- floor(min(c(geneA_df$cor, geneB_df$cor)) * 10) / 10  # 向下 0.1
xu <- ceiling(max(c(geneA_df$cor, geneB_df$cor)) * 10) / 10 # 向上 0.1
axis_breaks <- seq(xl, xu, by = 0.1)
axis_labels <- sprintf("%.1f", axis_breaks)

## ── 3. 色彩参数 ────────────────────────────────────────────
geneA_col  <- "#E76F51" ; geneB_col  <- "#457B9D"
bg_red   <- "#F9D9D3" ; bg_blue  <- "#D6EAF8"

## ── 4. 单条带绘图函数（支持自适应坐标 & 右侧标签）────────
strip_plot <- function(df, target_gene,
                       bg_fill, pt_fill, right_label,
                       xlim_vec, brks, labs, r = .12) {
  
  targ_x <- df$cor[df$geneB == target_gene]
  
  ggplot(df, aes(cor, 0)) +
    # 背景圆角矩形
    annotation_custom(
      roundrectGrob(x = .5, y = .5, width = 1, height = 1,
                    r = unit(r, "snpc"),
                    gp = gpar(fill = bg_fill, col = NA)),
      xmin = xlim_vec[1] - .02, xmax = xlim_vec[2] + .02,
      ymin = -.5, ymax =  .5) +
    
    # 主轴与刻度
    geom_segment(aes(x = xlim_vec[1], xend = xlim_vec[2], y = 0, yend = 0),
                 linewidth = .9) +
    geom_segment(data = data.frame(x = brks),
                 aes(x = x, xend = x, y = 0, yend = -.08),
                 linewidth = .8) +
    geom_text(data = data.frame(x = brks, y = -.18, lab = labs),
              aes(x, y, label = lab), size = 4) +
    
    # 数据点 & 箭头
    geom_point(shape = 21, size = 5, stroke = .35,
               fill = pt_fill, colour = "black") +
    geom_segment(x = targ_x, xend = targ_x, y = .28, yend = .06,
                 arrow = arrow(length = unit(.24, "cm")), linewidth = 1) +
    annotate("text", x = targ_x, y = .36, label = target_gene,
             fontface = "italic") +
    
    # ── ① 右侧纵向标签放在 x = Inf
    geom_text(x = Inf, y = 0, label = right_label,
              angle = -90, hjust = 0.5, vjust = 0.2,
              fontface = "italic") +
    
    # ── ② x 轴右侧扩展 
    scale_x_continuous(limits = xlim_vec,
                       expand = expansion(mult = c(0, 0.05))) +
    
    # ── ③ clip = "off" 以免被裁剪
    coord_cartesian(ylim = c(-.5, .5), clip = "off") +
    
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      plot.title  = element_text(face = "italic", size = 16,
                                 hjust = 0, margin = margin(b = 2)),
      plot.margin = margin(2, 20, 2, 24)   # 上 右 下 左
    ) 
}


## ── 5. 绘制上下两张图 ──────────────────────────────────────
p1 <- strip_plot(geneA_df,"MDM2",
                 bg_red, geneA_col,  right_label = "TP53",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)

p2 <- strip_plot(geneB_df,"TP53",
                 bg_blue, geneB_col,  right_label = "MDM2",
                 xlim_vec = c(xl, xu), brks = axis_breaks, labs = axis_labels)


## ── 6. 主标题 (保持不变) ───────────────────────────────────
title_bar <- ggdraw() +
  draw_label("Correlation", fontface="bold", size=19, x=.5, y=.8, hjust=.5) +
  draw_line(x=c(.43,.17), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Inverse", x=.13, y=.55, size=11) +
  draw_line(x=c(.57,.83), y=c(.55,.55),
            arrow=arrow(type="closed", length=unit(.14,"cm"))) +
  draw_label("Direct",  x=.87, y=.55, size=11)

## ── 7. 拼图输出 ───────────────────────────────────────────
panel <- plot_grid(p1, p2, ncol=1, align="v")
plot_grid(title_bar, panel, ncol=1, rel_heights=c(.16,1))

### 推荐两个我们开发的AI工具
### wispterm：https://xuzhougeng.github.io/wispterm/
### wisp science：https://xuzhougeng.github.io/wisp-science/ 
### wisp figure: 
### wisp depmap: 


