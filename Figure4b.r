library(tidyverse)
library(patchwork)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

npg_colors <- pal_npg("nrc")(7)
npg_colors <- alpha(npg_colors, 0.6)

setwd("D://")

data <- read.table("data.tsv", sep = "\t", fill = TRUE, header = TRUE) %>%
  mutate(pvalue = as.numeric(pvalue),
         log_pvalue = -log10(pvalue))

data$Description <- factor(data$Description, levels = rev(unique(data$Description)))
data$ONTOLOGY   <- factor(data$ONTOLOGY,   levels = unique(data$ONTOLOGY))

data <- data %>%
  mutate(gene_count = str_count(geneID, "/") + 1) %>%
  select(Description, pvalue, ONTOLOGY, Ontology, category, gene_count, geneID) %>%
  mutate(log_pvalue = -log10(as.numeric(pvalue)))

data <- within(data, geneID <- sapply(strsplit(geneID, "/"),
                                      function(x) paste(head(trimws(x), 2), collapse = "/")))

onto_levels <- levels(factor(data$Ontology, levels = unique(data$Ontology)))
onto_colors <- setNames(
  ggsci::pal_nejm()(max(3, length(onto_levels)))[seq_along(onto_levels)],
  onto_levels
)

n_groups <- length(unique(data$ONTOLOGY))
group_colors <- npg_colors
names(group_colors) <- unique(data$ONTOLOGY)

desc_strip <- data %>%
  distinct(Description, Ontology) %>%
  mutate(
    Description = factor(Description, levels = levels(data$Description)),
    Ontology    = factor(Ontology,    levels = onto_levels)
  )

p0 <- ggplot(desc_strip, aes(x = Description, y = 1, fill = Ontology)) +
  geom_tile(height = 1) +
  scale_x_discrete(position = "top", limits = rev(levels(data$Description))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = scales::alpha(onto_colors, 0.6), name = "name") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.margin = margin(b = 0)
  )

p1 <- ggplot(data, aes(x = Description, y = ONTOLOGY)) +
  geom_vline(xintercept = seq(0.5, length(unique(data$Description)) + 0.5, by = 1),
             color = "gray88", linewidth = 0.2) +
  geom_hline(yintercept = seq(0.5, length(unique(data$ONTOLOGY)) + 0.5, by = 1),
             color = "gray88", linewidth = 0.2) +
  geom_point(aes(size = gene_count,
                 fill = ONTOLOGY,
                 color = log_pvalue,
                 shape = category),
             stroke = 0.8) +
  scale_fill_manual(values = scales::alpha(group_colors, 0.7), guide = "none") +
  scale_color_gradientn(colors = c("white", "#1E4565"), name = "-log10(pvalue)",
                        guide = guide_colorbar(reverse = FALSE)) +
  scale_shape_manual(values = c(GO = 24, KEGG = 21), name = "name") +
  scale_size_continuous(range = c(0.5, 3), name = "name") +
  labs(title = "", x = "x", y = "y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 0, b = 0)
  ) +
  scale_x_discrete(position = "top", limits = rev(levels(data$Description))) +
  scale_y_discrete(limits = rev(levels(data$ONTOLOGY)))

gene_data <- data %>%
  separate_rows(geneID, sep = "/") %>%
  select(geneID, ONTOLOGY, Description, log_pvalue) %>%
  mutate(
    geneID = factor(geneID, levels = unique(geneID)),
    Description = factor(Description, levels = levels(data$Description))
  )

plot_data <- gene_data %>%
  group_by(geneID, Description) %>%
  mutate(
    position = seq_len(n()),
    total = n(),
    ymin = as.numeric(Description) - 0.4 + (position - 1) * (0.8/total),
    ymax = as.numeric(Description) - 0.4 + position * (0.8/total)
  ) %>%
  ungroup()

plot_data$Description <- factor(plot_data$Description, levels = rev(levels(plot_data$Description)))
plot_data$geneID      <- factor(plot_data$geneID,      levels = rev(unique(plot_data$geneID)))

p2 <- ggplot(plot_data, aes(x = Description, y = geneID)) +
  geom_tile(width = 0.9, height = 0.9, fill = "grey95", color = NA, show.legend = FALSE) +
  geom_vline(xintercept = seq(0.5, length(unique(plot_data$Description)) + 0.5, by = 1),
             color = "gray88", linewidth = 0.1) +
  geom_hline(yintercept = seq(0.5, length(unique(plot_data$geneID)) + 0.5, by = 1),
             color = "gray88", linewidth = 0.2) +
  geom_rect(
    aes(ymin = as.numeric(geneID) - 0.4,
        ymax = as.numeric(geneID) + 0.4,
        xmin = as.numeric(Description) - 0.4 + (position - 1) * (0.8/total),
        xmax = as.numeric(Description) - 0.4 + position * (0.8/total),
        fill = ONTOLOGY),
    color = NA,
    alpha = 0.6
  ) +
  scale_y_discrete(limits = levels(plot_data$geneID), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_fill_manual(values = group_colors, name = "Subtype",
                    guide = guide_legend(override.aes = list(alpha = 0.7))) +
  labs(title = NULL, y = "Gene", x = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

p2 <- p2 +
  guides(
    fill = guide_legend(
      keywidth  = grid::unit(5, "mm"),
      keyheight = grid::unit(5, "mm")
    )
  )

p0 <- p0 +
  guides(
    fill = guide_legend(
      keywidth  = grid::unit(5, "mm"),
      keyheight = grid::unit(5, "mm")
    )
  )

p1 <- p1 +
  guides(
    shape = guide_legend(
      override.aes = list(
        size  = c(2.2, 3.6),
        stroke = 0.8
      ),
      keywidth  = grid::unit(5, "mm"),
      keyheight = grid::unit(5, "mm")
    )
  )

combined_plot <- p1 / p0 / p2 +
  plot_layout(heights = c(1.3, 0.15, 7), guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

ggsave(
  combined_plot,
  filename = "heatmap.png",
  width = 8, height = 10, dpi = 300
)

ggsave(
  combined_plot,
  filename = "heatmap.pdf",
  width = 8 * 1.2, height = 10 * 1.1
)