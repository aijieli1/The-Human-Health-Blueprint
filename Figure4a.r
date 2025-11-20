# ------------------------------------------------------------
# Circular heatmap of SNP counts per chromosome and subgroup
# Sectors: CHR 1..22
# Rings:   prefix (7 subgroups)
# Color:   SNP_CHR_COUNT
# Label:   Gene in each ring/sector cell
# Changes:
#   - Add black border to cells where Gene is duplicated across the grid.
#   - Change the in-cell second line to (count/SNP_CHR_COUNT); if `count` column
#     is not present in df, it falls back to SNP_CHR_COUNT.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(viridisLite)
  library(ComplexHeatmap)
  library(grid)
})

plot_snp_circular_heatmap <- function(df,
                                      prefix_levels = c("A","B","N","C","D","E","F"),
                                      chr_levels = as.character(1:22),
                                      palette = viridisLite::viridis(200, direction = -1),
                                      na_color = "#EEEEEE",
                                      track_height = 0.075,
                                      gap_degree = 2,
                                      start_degree = 90,
                                      label_cex = 0.35,
                                      label_color = "black",
                                      truncate_labels = TRUE,
                                      max_label_char = 12,
                                      show_chr_labels = TRUE,
                                      chr_label_cex = 0.6,
                                      legend_title = "SNP_CHR_COUNT",
                                      legend_side = c("right","left","bottom","top")[1]) {
  
  # Basic checks
  required_cols <- c("prefix","Gene","CHR","SNP_CHR_COUNT")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Handle optional 'count' column (for label count/SNP_CHR_COUNT).
  has_count_col <- "count" %in% names(df)
  
  # Normalize columns
  df <- df %>%
    mutate(
      prefix = as.character(prefix),
      Gene   = as.character(Gene),
      CHR    = as.character(CHR),
      SNP_CHR_COUNT = as.numeric(SNP_CHR_COUNT),
      count = if (has_count_col) as.numeric(.data$count) else as.numeric(SNP_CHR_COUNT)
    ) %>%
    filter(CHR %in% chr_levels)
  
  # Keep prefix order that actually exists in data
  prefix_levels <- prefix_levels[prefix_levels %in% unique(df$prefix)]
  if (length(prefix_levels) == 0) stop("No valid prefixes found in data.")
  
  # If duplicate (prefix, CHR), keep the one with the largest SNP_CHR_COUNT
  df_dedup <- df %>%
    group_by(prefix, CHR) %>%
    arrange(desc(SNP_CHR_COUNT), .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup()
  
  # Complete 7 x 22 grid (allow NA)
  grid_df <- tidyr::expand_grid(prefix = prefix_levels, CHR = chr_levels) %>%
    left_join(df_dedup, by = c("prefix","CHR")) %>%
    mutate(Gene = ifelse(is.na(Gene), "", Gene))
  
  # Identify duplicated genes across the entire grid (ignore empty "")
  dup_genes <- grid_df %>%
    filter(Gene != "") %>%
    count(Gene, name = "n") %>%
    filter(n > 1) %>%
    pull(Gene)
  
  # Color mapping
  vmin <- min(grid_df$SNP_CHR_COUNT, na.rm = TRUE)
  vmax <- max(grid_df$SNP_CHR_COUNT, na.rm = TRUE)
  if (!is.finite(vmin) || !is.finite(vmax)) {
    stop("SNP_CHR_COUNT appears to be all NA or non-finite.")
  }
  col_fun <- circlize::colorRamp2(c(vmin, vmax), c(palette[1], palette[length(palette)]))
  
  # Circos setup
  circos.clear()
  circos.par(
    start.degree = start_degree,
    gap.degree   = gap_degree,
    track.margin = c(0.002, 0.002),
    cell.padding = c(0,0,0,0),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.30, 1.45),
    canvas.ylim = c(-1.20, 1.20)
  )
  # Initialize sectors (equal width for 22 chromosomes)
  chr_fct_levels <- factor(chr_levels, levels = chr_levels)
  xlim_mat <- cbind(rep(0, length(chr_levels)), rep(1, length(chr_levels)))
  circos.initialize(factors = chr_fct_levels, xlim = xlim_mat)
  
  # Optional outer track for chromosome labels
  if (show_chr_labels) {
    circos.trackPlotRegion(
      ylim = c(0,1),
      track.height = 0.10,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector = get.cell.meta.data("sector.index")
        circos.text(
          x = 0.5, y = 0.5,
          labels = paste0("Chr", sector),
          cex = chr_label_cex, 
          font = 2,            
          facing = "inside",
          niceFacing = FALSE
        )
      }
    )
  }
  
  # Small integer formatter
  fmt_int <- function(x) formatC(x, digits = 0, format = "f", big.mark = ",")
  
  # Draw one ring per prefix
  for (pref in prefix_levels) {
    d <- grid_df %>% filter(prefix == pref)
    value_by_chr <- setNames(d$SNP_CHR_COUNT, d$CHR)
    gene_by_chr  <- setNames(d$Gene,          d$CHR)
    count_by_chr <- setNames(d$count,          d$CHR)
    
    circos.trackPlotRegion(
      ylim = c(0,1),
      track.height = track_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector = get.cell.meta.data("sector.index")
        val  <- value_by_chr[[sector]]      # SNP_CHR_COUNT
        gene <- gene_by_chr[[sector]]
        cnt  <- count_by_chr[[sector]]      # count (or fallback)
        
        fill_col <- if (is.na(val)) na_color else col_fun(val)
        is_dup_gene <- (!is.na(gene) && nzchar(gene) && gene %in% dup_genes)
        
        # Fill cell with optional black border if Gene duplicated
        circos.rect(
          xleft = 0, ybottom = 0, xright = 1, ytop = 1,
          col = fill_col,
          border = if (is_dup_gene) "black" else NA,
          lwd = if (is_dup_gene) 1.2 else 0.5
        )
        
        # Labels: line 1 Gene, line 2 (count/SNP_CHR_COUNT)
        if (!is.na(gene) && nzchar(gene)) {
          lab <- gene
          if (truncate_labels && nchar(lab) > max_label_char) {
            lab <- paste0(substr(lab, 1, max_label_char - 1), "…")
          }
          # Line 1: Gene (upper)
          circos.text(
            x = 0.5, y = 0.72, labels = lab,
            cex = label_cex,
            col = label_color,
            facing = "inside",
            niceFacing = TRUE
          )
          # Line 2: (count/SNP_CHR_COUNT), only if val is not NA
          if (!is.na(val)) {
            # if cnt is NA (shouldn't when val exists), fallback to val
            disp_cnt <- if (!is.na(cnt)) cnt else val
            cnt_lab <- paste0("(", fmt_int(disp_cnt), "/", fmt_int(val), ")")
            circos.text(
              x = 0.5, y = 0.28, labels = cnt_lab,
              cex = label_cex * 0.9,
              col = label_color,
              facing = "inside",
              niceFacing = TRUE
            )
          }
        }
      }
    )
  }
  
  # Legend for color scale
  lg <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    title = legend_title,
    at = pretty(c(vmin, vmax), n = 5),
    labels = pretty(c(vmin, vmax), n = 5)
  )
  
  # Place legend
  if (legend_side %in% c("right","left","top","bottom")) {
    draw(lg, just = c(switch(legend_side,
                             right  = "right",
                             left   = "left",
                             top    = "center",
                             bottom = "center"),
                      switch(legend_side,
                             right  = "center",
                             left   = "center",
                             top    = "top",
                             bottom = "bottom")),
         x = switch(legend_side,
                    right  = unit(1, "npc") - unit(5, "mm"),
                    left   = unit(5, "mm"),
                    top    = unit(0.5, "npc"),
                    bottom = unit(0.5, "npc")),
         y = switch(legend_side,
                    right  = unit(0.5, "npc"),
                    left   = unit(0.5, "npc"),
                    top    = unit(1, "npc") - unit(5, "mm"),
                    bottom = unit(5, "mm"))
    )
  } else {
    draw(lg, just = c("right","center"), x = unit(1, "npc") - unit(5, "mm"), y = unit(0.5, "npc"))
  }
  
  invisible(list(col_fun = col_fun, vmin = vmin, vmax = vmax))
  
  # Prefix order legend (left side)
  prefix_lg <- ComplexHeatmap::Legend(
    title = "Subgroups (outer -> inner)",
    labels = paste0(seq_along(prefix_levels), ". ", prefix_levels),
    type = "labels",
    labels_gp = grid::gpar(fontsize = 10)
  )
  draw(prefix_lg, x = unit(6, "mm"), y = unit(0.5, "npc"), just = c("left", "center"))
}

# ----------------------------

# ----------------------------
# Example usage:
# ----------------------------
# df <- readr::read_tsv("D://snp_count_gene.tsv", show_col_types = FALSE)
# pdf("circos_snp_heatmap.pdf", width = 8, height = 8)
# plot_snp_circular_heatmap(df)
# dev.off()
# ----------------------------
# Usage examples
# ----------------------------

