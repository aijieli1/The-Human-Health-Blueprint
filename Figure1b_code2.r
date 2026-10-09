library(dplyr)
library(ggsci)
library(scales)
library(pheatmap)

b      <- read.table("rg_square.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
group  <- read.table("group.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
p_FDR  <- read.table("rg_square_p_FDR.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

b[is.na(b)] <- 0
p_FDR[is.na(p_FDR)] <- 1

fix_names <- function(x) gsub("_adjbmi", "*", x)
rownames(b)     <- fix_names(rownames(b));     colnames(b)     <- fix_names(colnames(b))
rownames(p_FDR) <- fix_names(rownames(p_FDR)); colnames(p_FDR) <- fix_names(colnames(p_FDR))
rownames(group) <- fix_names(rownames(group))

row_group <- group %>% filter(!Group %in% c("one","two","three","four"))
col_group <- group %>% filter( Group %in% c("one","two","three","four"))

choose_row <- intersect(rownames(row_group), rownames(b))
choose_col <- intersect(rownames(col_group), colnames(b))

b_choose     <- b[choose_row, choose_col, drop = FALSE]
p_FDR_choose <- p_FDR[choose_row, choose_col, drop = FALSE]

ordered_col <- intersect(rownames(col_group), colnames(b_choose))
b_choose     <- b_choose[, ordered_col, drop = FALSE]
p_FDR_choose <- p_FDR_choose[, ordered_col, drop = FALSE]

n_group <- length(unique(group$Group))
base_pal   <- pal_npg("nrc")(10)
npg_colors <- colorRampPalette(base_pal)(n_group)
ann_colors <- list(Group = setNames(alpha(npg_colors, 0.6), sort(unique(group$Group))))

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
p_FDR_choose[is.na(p_FDR_choose)] <- 1

display_numbers <- matrix(sprintf("%.2f", as.numeric(b_choose)),
                          nrow = nrow(b_choose), ncol = ncol(b_choose),
                          dimnames = dimnames(b_choose))
display_numbers[p_FDR_choose > 0.05] <- " "

# 绘图
pheatmap(
  b_choose,
  color = c(
    colorRampPalette(c("#032f6c","white"))(length(bk)/2),
    colorRampPalette(c("white","red4"))(length(bk)/2)
  ),
  breaks = bk,
  legend_breaks = seq(-1, 1, 0.5),
  annotation_col = col_group[ordered_col, , drop = FALSE],
  annotation_row = row_group[choose_row, , drop = FALSE],
  annotation_colors = ann_colors,
  border_color = "grey80",
  number_color = "black",
  cellwidth = 16, cellheight = 15,
  angle_col = 45, fontsize_number = 6/1.2,
  show_colnames = TRUE,
  display_numbers = display_numbers,
  fontsize = 9,
  cluster_rows = FALSE, cluster_cols = FALSE,
  number_format = "%.2f",
  filename = "output.pdf",
  main = "Genetic correlation"
)
