library(ggplot2)
library(scatterpie)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggnewscale)
library(scales)
library(patchwork)

rm(list=ls())

df <- read.table("D://result.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

df <- df %>% filter(verify != "no")

gene_levels <- unique(df$casual_SNP_Gene)
files_levels <- c("A","B","C","D","E","F","G")

scale_by_count <- TRUE  
r_min <- 0.30          
r_max <- 0.60          
r_fixed <- 0.40        
use_fixed_ratio <- FALSE 
ratio_y_per_x <- 0.7   
show_y_group_strip <- TRUE
group_by <- "show_atc3_name" 
strip_width <- 0.8      
strip_x_center <- 0.4   

pal_files <- scales::alpha(c(
  "A" = "#E64B35FF",
  "B"  = "#4DBBD5FF",
  "C"  = "#00A087FF",
  "D"= "#3C5488FF",
  "E"  = "#F39B7FFF",
  "F"  = "#8491B4FF",
  "G" = "#E6B800"
), 0.8)

for (nm in intersect(c("J","Q","K","A","2"), names(df))) {
  df[[nm]] <- ifelse(is.na(df[[nm]]), NA, trimws(df[[nm]]))
}

required_cols <- c("J","Q","K")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf(" %s", paste(missing_cols, collapse = ", ")))
}

dat <- df %>%
  transmute(
    gene = as.character(.data$casual_SNP_Gene),
    drug = as.character(.data$generic_name),
    file = as.character(.data$file)
  ) %>%
  filter(!is.na(gene), !is.na(drug), !is.na(file),
         file %in% files_levels) %>% distinct() %>%
  mutate(file = factor(file, levels = files_levels))

agg <- dat %>%
  group_by(gene, drug, file) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(gene, drug, file = files_levels, fill = list(n = 0))

wide <- agg %>%
  tidyr::pivot_wider(names_from = file, values_from = n, values_fill = 0)

drug_levels <- unique(df$generic_name[df$generic_name %in% unique(wide$drug)])

plot_df <- wide %>%
  mutate(
    x = as.numeric(factor(gene, levels = gene_levels)),
    y = as.numeric(factor(drug, levels = drug_levels)),
    total_n = rowSums(dplyr::across(all_of(files_levels)))
  )

if (scale_by_count) {
  rng <- range(plot_df$total_n, na.rm = TRUE)
  plot_df <- plot_df %>%
    mutate(r = ifelse(total_n > 0, rescale(total_n, to = c(r_min+0.1, r_max+0)), 0))
} else {
  plot_df <- plot_df %>% mutate(r = ifelse(total_n > 0, r_fixed, 0))
}

p <- ggplot()

if (show_y_group_strip) {
  if (!group_by %in% c("show_atc3_name","show_atc3_name")) {
    stop("group_by  'show_atc3_name'  'show_atc3_name'")
  }
  if (!group_by %in% names(df)) {
    stop(sprintf("%s", group_by))
  }
  
  atc_col <- df[[group_by]]
  atc_levels <- unique(atc_col[!is.na(atc_col)])
  
  drug_meta <- df %>%
    select(generic_name, show_atc1_name, show_atc3_name) %>%
    distinct(generic_name, .keep_all = TRUE)
  
  group_df <- tibble(drug = drug_levels) %>%
    left_join(drug_meta, by = c("drug" = "generic_name")) %>%
    mutate(
      group = factor(.data[[group_by]], levels = atc_levels),
      y = as.numeric(factor(drug, levels = drug_levels))
    )
  
  pal_groups <- setNames(
    c("#1F77B480", "#FF7F0E80"),
    atc_levels
  )
  
  strip_df <- group_df %>%
    transmute(
      x = strip_x_center,
      y = y,
      width = strip_width,
      height = 0.9,
      group = group
    )
  
  p <- p +
    geom_tile(
      data = strip_df,
      aes(x = x, y = y, width = width, height = height, fill = group),
      color = NA
    ) +
    scale_fill_manual(
      values = pal_groups,
      limits = atc_levels,
      drop = FALSE,
      name = paste0(group_by, " group")
    ) +
    ggnewscale::new_scale_fill()
}

p <- p +
  scatterpie::geom_scatterpie(
    data = plot_df,
    aes(x = x, y = y, r = r),
    cols = files_levels,
    color = "white",
    size = 0
  ) +
  scale_fill_manual(values = pal_files, name = "File")

x_min <- if (show_y_group_strip) min(0, strip_x_center - strip_width/2) else 0.5
x_max <- length(gene_levels) + 0.6

if (use_fixed_ratio) {
  p <- p + coord_fixed(ratio = ratio_y_per_x,
                       xlim = c(x_min, x_max),
                       ylim = c(0.5, length(drug_levels) + 0.5),
                       expand = FALSE)
} else {
  p <- p + coord_cartesian(xlim = c(x_min, x_max),
                           ylim = c(0.5, length(drug_levels) + 0.5),
                           expand = FALSE)
}

p <- p +
  scale_x_continuous(
    breaks = seq_along(gene_levels),
    labels = gene_levels
  ) +
  scale_y_continuous(
    breaks = seq_along(drug_levels),
    labels = drug_levels,
    trans = "reverse"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 6.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  labs(
    x = "x",
    y = "y",
    title = "title"
  )

p <- p + coord_fixed(ratio = 1,
                     xlim = c(x_min, x_max),
                     expand = FALSE)

df2 <- df %>%
  select(casual_SNP_Gene, GETx.eQTL.MR.Direction, file) %>% unique()

heatmap_data <- df2 %>%
  select(casual_SNP_Gene, GETx.eQTL.MR.Direction, file) %>%
  distinct() %>%
  group_by(casual_SNP_Gene, file, GETx.eQTL.MR.Direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(casual_SNP_Gene, file, GETx.eQTL.MR.Direction, fill = list(count = 0)) %>%
  filter(count > 0)

gene_order <- gene_levels
file_order <- rev(c("A","B","C","D","E","F","G"))

red <- "red"
blue <- "blue"
shape_mapping <- c("up" = 24, "down" = 25)

p2 <- ggplot(heatmap_data, 
             aes(x = factor(casual_SNP_Gene, levels = gene_order), 
                 y = factor(file, levels = file_order))) +
  scale_x_discrete(
    breaks = gene_levels,
    labels = gene_levels
  ) +
  geom_vline(xintercept = seq(1, length(gene_order), by = 1), color = "gray90", linewidth = 0.2) +
  geom_hline(yintercept = seq(1, length(file_order), by = 1), color = "gray90", linewidth = 0.2) +
  geom_point(aes(shape = GETx.eQTL.MR.Direction, fill = GETx.eQTL.MR.Direction),
             size = 3, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c("up" = red, "down" = blue), name = "Direction") +
  scale_shape_manual(values = shape_mapping, name = "Direction") +
  labs(title = "title", x = "x", y = "y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "right",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
  ) +
  coord_fixed(ratio = 0.5) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE)

p2 <- p2 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    plot.title = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  coord_fixed(ratio = 1)

final_plot <- p / p2 + plot_layout(heights = c(2.5, 1), guides = "collect")

print(final_plot)

ggsave("D://result.pdf", final_plot, width = 12.5/1.5, height = 30/1.5, dpi = 300)