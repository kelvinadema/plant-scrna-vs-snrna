.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(scattermore)
library(ggrepel)

out_dir <- file.path(fig_main_dir, "F5")
species_levels <- c("medicago", "arabidopsis")
species_strip  <- labeller(species = as_labeller(species_expr, label_parsed))

comp <- read_tsv_q(section_table("05_composition", "fig_composition.tsv")) |>
  filter(primary_comparison) |>
  mutate(species   = factor(species, levels = species_levels),
         modality  = factor(modality, levels = c("nuclei", "protoplast")),
         cell_type = factor(cell_type_labels[cell_type], levels = cell_type_labels),
         treatment_shape = case_when(species == "arabidopsis" ~ "Untreated (A. thaliana)",
                                     treatment == "mock" ~ "Mock",
                                     treatment %in% c("NF", "nod_factor_3h") ~ "Nod factor",
                                     TRUE ~ "BAP"))

comp_mean <- comp |>
  group_by(species, modality, cell_type) |>
  summarise(pct = 100 * mean(proportion), .groups = "drop")

dodge <- position_dodge(width = 0.8)

pA <- ggplot(comp, aes(cell_type, 100 * proportion, fill = modality)) +
  geom_col(data = comp_mean, aes(y = pct), position = dodge, width = 0.75, alpha = 0.35) +
  geom_point(aes(colour = modality, shape = treatment_shape, group = modality),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1),
             size = 2.0, stroke = 0.4, fill = "white") +
  facet_grid(~species, scales = "free_x", space = "free_x", labeller = species_strip) +
  scale_fill_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_shape_manual(values = c("Mock" = 21, "Nod factor" = 24, "BAP" = 22,
                                "Untreated (A. thaliana)" = 23), name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Retained cells (%)") +
  guides(fill = guide_legend(order = 1), colour = "none",
         shape = guide_legend(order = 2, override.aes = list(colour = ink_secondary))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, colour = ink_primary),
        axis.ticks.x = element_blank(),
        legend.position = "right", legend.key.size = unit(2.5, "mm"),
        legend.spacing.y = unit(1, "mm"))

umap <- read_tsv_q(section_table("05_composition", "fig_umap_joint.tsv")) |>
  filter(primary_comparison) |>
  arrange(plot_order) |>
  mutate(cell_type = factor(cell_type_labels[cell_type], levels = names(col_cell_type)),
         modality  = factor(modality, levels = c("protoplast", "nuclei")))

theme_umap <- theme_void(base_size = base_size) +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = base_size, hjust = 0, margin = margin(b = 2)),
        plot.tag = element_text(size = base_size + 3, face = "bold"),
        axis.title.x = element_text(size = base_size - 1, colour = ink_secondary, hjust = 0),
        axis.title.y = element_text(size = base_size - 1, colour = ink_secondary, hjust = 0,
                                    angle = 90),
        legend.position = "bottom",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = base_size - 1),
        plot.margin = margin(3, 3, 3, 3, "mm"))

umap_title <- function(sp, what, d) {
  bquote(italic(.(lab_species[[sp]])) ~ "\u00b7" ~ .(what) ~
           "(" * .(scales::comma(nrow(d))) ~ "cells)")
}

umap_celltype <- function(sp) {
  d <- umap |> filter(species == sp)
  lab <- d |>
    filter(cell_type != "Unassigned") |>
    group_by(cell_type) |>
    filter(n() >= 30) |>
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
  ggplot(d, aes(umap_1, umap_2, colour = cell_type)) +
    geom_scattermore(pointsize = 2.2, pixels = c(1400, 1400)) +
    geom_label_repel(data = lab, aes(label = cell_type), colour = ink_primary,
                     size = (base_size - 1.5) / .pt, label.size = 0,
                     label.padding = unit(0.06, "lines"), fill = alpha("white", 0.75),
                     min.segment.length = 0.2, segment.size = 0.2, seed = 1,
                     max.overlaps = Inf) +
    scale_colour_manual(values = col_cell_type, limits = names(col_cell_type), drop = FALSE,
                        name = NULL) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2))) +
    labs(title = umap_title(sp, "cell type", d), x = "UMAP 1", y = "UMAP 2") +
    theme_umap
}

umap_modality <- function(sp) {
  d <- umap |> filter(species == sp)
  ggplot(d, aes(umap_1, umap_2, colour = modality)) +
    geom_scattermore(pointsize = 2.2, pixels = c(1400, 1400), alpha = 0.7) +
    scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL,
                        breaks = c("nuclei", "protoplast")) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = umap_title(sp, "modality", d), x = "UMAP 1", y = "UMAP 2") +
    theme_umap
}

legend_right <- theme(legend.position = "right")
pB <- umap_modality("medicago")
pC <- umap_modality("arabidopsis") + legend_right
pD <- umap_celltype("medicago")
pE <- umap_celltype("arabidopsis") + legend_right +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2)))
no_leg <- theme(legend.position = "none")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "p5A.png"), pA, width = 180, height = 65, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "p5B.png"), pB + no_leg, width = 80, height = 80, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "p5C.png"), pC, width = 100, height = 80, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "p5D.png"), pD + no_leg, width = 80, height = 80, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "p5E.png"), pE, width = 110, height = 80, units = "mm", dpi = 600, bg = "white")

mod_legend <- wrap_elements(full = cowplot::get_plot_component(pC, "guide-box-right"))
ct_legend  <- wrap_elements(full = cowplot::get_plot_component(pE, "guide-box-right"))

layout <- "
AAA
BCF
DEG
"
fig5 <- pA + (pB + no_leg) + (pC + no_leg) + (pD + no_leg) + (pE + no_leg) + mod_legend + ct_legend +
  plot_layout(design = layout, widths = c(1, 1, 0.38), heights = c(0.8, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "", "")))

save_fig(fig5, "Fig5_composition_umap", width_mm = 180, height_mm = 215, dir = out_dir)
