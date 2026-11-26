.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
out_dir <- file.path(fig_main_dir, "F4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s4 <- "04_read_yield"
min_cells_curve <- 50

species_levels <- c("medicago", "arabidopsis")
species_strip  <- labeller(species = as_labeller(species_expr, label_parsed))

lib_order <- function(df) {
  df |>
    mutate(species = factor(species, levels = species_levels),
           modality = factor(modality, levels = c("nuclei", "protoplast"))) |>
    arrange(species, modality, sample_label) |>
    mutate(sample_label = factor(sample_label, levels = rev(unique(sample_label))))
}

yield <- read_tsv_q(section_table(s4, "fig_yield_perlib.tsv")) |>
  filter(primary_comparison) |>
  lib_order()

fate_labels <- c(unmapped                        = "Unmapped",
                 mapped_outside_permit_barcodes  = "Invalid barcode",
                 umi_duplicates                  = "Duplicate",
                 removed_by_cellbender           = "Ambient",
                 molecules_not_in_retained_cells = "Not in retained cells",
                 retained                        = "Informative")
loss_cols <- setNames(c("#e6e5e1", "#c4c3bd", "#9d9c96", "#6f6e69", "#3f3e3a"), names(fate_labels)[1:5])

readfate <- read_tsv_q(section_table(s4, "fig_readfate_stack.tsv")) |>
  filter(primary_comparison) |>
  lib_order() |>
  mutate(fill = if_else(fate == "retained", as.character(modality), fate),
         fill = factor(fill, levels = c(names(loss_cols), "nuclei", "protoplast")))

pA <- ggplot(readfate, aes(100 * fraction_of_input_reads, sample_label, fill = fill,
                            group = factor(fate, levels = names(fate_labels)))) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2, position = position_stack(reverse = TRUE)) +
  geom_text(data = yield, aes(x = 101.5, y = sample_label, label = sprintf("%.1f%%", 100 * frac_retained),
                              colour = modality),
            inherit.aes = FALSE, hjust = 0, size = (base_size - 1) / .pt, fontface = "bold") +
  facet_grid(species ~ ., scales = "free_y", space = "free_y", labeller = species_strip) +
  scale_fill_manual(values = c(loss_cols, col_modality),
                    labels = c(fate_labels[1:5], nuclei = "Informative (nuclei)",
                               protoplast = "Informative (protoplast)"),
                    name = NULL) +
  scale_colour_manual(values = col_modality, guide = "none") +
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(0, 122), expand = expansion(0)) +
  labs(x = "Input read pairs (%)", y = NULL) +
  guides(fill = guide_legend(nrow = 1)) + theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary),
        strip.text.y = element_text(angle = -90, hjust = 0.5),
        legend.position = "top", legend.justification = "center",
        legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = base_size - 1.5))

curve <- read_tsv_q(section_table(s4, "fig_yield_curve.tsv")) |>
  filter(primary_comparison, valid) |>
  mutate(species = factor(species, levels = species_levels))

pB <- ggplot(curve, aes(reads_per_fixed_cell, 100 * yield, colour = modality, group = sample)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.5, shape = 21, fill = "white", stroke = 0.4) +
  facet_wrap(~species, nrow = 1, labeller = species_strip) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5), labels = c("1k", "10k", "100k")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Reads per cell", y = "Informative reads (%)") +
  theme(legend.position = "top", legend.justification = "center",
        legend.margin = margin(0, 0, -3, 0))

matched <- read_tsv_q(section_table(s4, "fig_detection_matched.tsv")) |>
  lib_order()

rare <- read_tsv_q(section_table(s4, "fig_rarefaction_curves.tsv")) |>
  filter(definition == "corrected", primary_comparison, n_cells_eligible >= min_cells_curve) |>
  mutate(species = factor(species, levels = species_levels))

matched_depth <- matched |> distinct(species, matched_depth)

pC <- ggplot(rare, aes(depth, genes_detected, colour = modality, group = sample)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.5, shape = 21, fill = "white", stroke = 0.4) +
  facet_wrap(~species, nrow = 1, scales = "free_x", labeller = species_strip) +
  theme_bw() +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "UMIs per cell", y = "Genes per cell") +
  theme(legend.position = "none", panel.spacing.x = unit(5, "mm"))

ggsave(filename = 'p4A.png', pA, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p4B.png', pB, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p4C.png', pC, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)

pC <- pC +
  facet_wrap(~species, nrow = 1, scales = "free_x", labeller = species_strip) +
  scale_x_continuous(breaks = c(1000, 5000, 10000), labels = c("1k", "5k", "10k"),
                     expand = expansion(mult = c(0, 0.03)))

modality_legend <- wrap_elements(full = cowplot::get_plot_component(
  pB + theme(legend.position = "bottom", legend.direction = "horizontal"),
  "guide-box-bottom"))
pB <- pB + theme(legend.position = "none")
pC <- pC + theme(legend.position = "none")

layout <- "
AB
AD
AC
"

fig4 <- pA + pB + pC + modality_legend +
  plot_layout(design = layout, widths = c(1.5, 1), heights = c(1, 0.08, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "")))

save_fig(fig4, "Fig4_read_yield", width_mm = 180, height_mm = 150, dir = out_dir)
