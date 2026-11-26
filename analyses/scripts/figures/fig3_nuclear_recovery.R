.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
out_dir <- file.path(fig_main_dir, "F3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
library(ggbeeswarm)

s3 <- "03_recovery_attrition"
floor_cells <- 500

libs <- read_tsv_q(section_table(s3, "fig_mapping_readfate_perlib.tsv")) |>
  filter(is_used_library(sample, species, modality)) |>
  mutate(group   = factor(prep_groups[arm], levels = prep_groups),
         species = factor(species, levels = c("medicago", "arabidopsis")),
         passed  = cells_retained >= floor_cells)

group_counts <- libs |>
  group_by(species, group, modality) |>
  summarise(label = sprintf("%d/%d", sum(passed), n()), .groups = "drop")

pA <- ggplot(libs, aes(group, cells_retained)) +
  geom_hline(yintercept = floor_cells, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_quasirandom(aes(colour = modality, fill = if_else(passed, modality, "below")),
                   shape = 21, size = 1.8, stroke = 0.5, width = 0.22) +
  geom_text(data = group_counts, aes(group, 3e4, label = label, colour = modality),
            size = base_size / .pt, fontface = "bold") +
  geom_text(data = tibble(species = factor("medicago", levels = levels(libs$species)),
                          group = factor(prep_groups[[3]], levels = prep_groups)),
            aes(group, floor_cells), label = "500-cell floor", hjust = 0.5, vjust = -0.6,
            size = (base_size - 1.5) / .pt, colour = ink_muted) +
  facet_grid(~species, scales = "free_x", space = "free_x",
             labeller = labeller(species = as_labeller(species_expr, label_parsed))) +
  scale_colour_manual(values = col_modality, guide = "none") +
  scale_fill_manual(values = c(col_modality, below = "white"), guide = "none") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 10, 100, 1000, 10000), labels = scales::label_comma(),
                     limits = c(0, 4e4), expand = expansion(mult = c(0.03, 0))) +
  labs(x = NULL, y = "Retained cells / nuclei") + theme_bw(base_size = base_size) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x  = element_text(colour = ink_primary, lineheight = 0.95))

fail <- read_failure_modes() |>
  select(sample, failure_mode)

nuc <- libs |>
  filter(modality == "nuclei") |>
  inner_join(fail, by = "sample") |>
  mutate(failure_mode = factor(failure_mode, levels = names(failure_labels)),
         row_label = if_else(species == "arabidopsis", sample_label, sample)) |>
  arrange(failure_mode, species, plot_order) |>
  mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

qc_medians <- read_tsv_q(section_table(s3, "fig_qc_percell_summary.tsv")) |>
  filter(stage == "count_gate_candidate", metric %in% c("mito_fraction", "rrna_fraction")) |>
  select(sample, metric, median) |>
  pivot_wider(names_from = metric, values_from = median, names_prefix = "candidate_")

libs <- libs |> left_join(qc_medians, by = "sample")
nuc  <- nuc  |> left_join(qc_medians, by = "sample")

prot_range <- function(col) range(libs[[col]][libs$modality == "protoplast"])

row_layout <- list(
  facet_grid(failure_mode ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(failure_mode = failure_labels)),
  labs(y = NULL),
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
        panel.spacing.y    = unit(1.2, "mm"),
        strip.text.y       = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.y        = element_blank(),
        axis.title.x       = element_text(size = base_size - 0.5, lineheight = 0.95),
        plot.margin        = margin(3, 1.5, 3, 1.5, "mm"))
)

band <- function(col, scale = 1) {
  r <- prot_range(col) * scale
  annotate("rect", xmin = r[1], xmax = r[2], ymin = -Inf, ymax = Inf,
           fill = col_modality[["protoplast"]], alpha = 0.14)
}

dot <- function(x) {
  geom_point(aes(x = {{ x }}), shape = 21, size = 1.7, stroke = 0.4,
             colour = "white", fill = col_modality[["nuclei"]])
}

pB <- ggplot(nuc, aes(y = row_label)) +
  band("frac_mapping", 100) +
  dot(100 * frac_mapping) +
  scale_x_continuous(limits = c(70, 100), breaks = c(70, 80, 90, 100)) +
  labs(x = "Reads mapped\n(%)") +
  row_layout + theme_bw(base_size = base_size) +
  theme(axis.text.y = element_text(colour = ink_primary, size = base_size - 1), strip.text = element_blank())

fate_labels <- c(outside_cb_model           = "Empty droplets",
                 removed_within_cb_model    = "Removed as ambient",
                 corrected_below_count_gate = "Below cell-calling cutoff",
                 corrected_later_qc_removed = "Removed by QC",
                 retained                   = "Retained")
fate_cols <- c(outside_cb_model = "#dcdbd6", removed_within_cb_model = "#b3b2ab",
               corrected_below_count_gate = "#86857f", corrected_later_qc_removed = "#3f3e3a",
               retained = col_modality[["nuclei"]])

fates <- read_tsv_q(section_table(s3, "fig_readfate_stack.tsv")) |>
  filter(partition == "raw_umis") |>
  inner_join(nuc |> select(sample, row_label, failure_mode), by = "sample") |>
  mutate(fate = factor(fate, levels = names(fate_labels)))

pC <- ggplot(fates, aes(100 * fraction, row_label, fill = fate)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = fate_cols, labels = fate_labels, name = NULL) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Deduplicated molecules\n(%)") +
  row_layout +
  guides(fill = guide_legend(nrow = 1, title = "Deduplicated molecules:")) +
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = base_size - 1, colour = ink_secondary),
        legend.margin = margin(0, 0, 0, 0), strip.text = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

pD <- ggplot(nuc, aes(y = row_label)) +
  band("raw_rank_1000_to_10000_ratio") +
  dot(raw_rank_1000_to_10000_ratio) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50), limits = c(1.5, 90)) +
  labs(x = "Barcode-rank knee\n(UMIs rank 1k / 10k)") + theme_bw(base_size = base_size) +
  row_layout

pE <- ggplot(nuc, aes(y = row_label)) +
  band("candidate_mito_fraction", 100) +
  geom_vline(xintercept = 2, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  dot(100 * candidate_mito_fraction) +
  scale_x_continuous(limits = c(0, 16), breaks = c(0, 5, 10, 15)) +
  labs(x = "Median mitochondrial\nfraction (%)") + theme_bw(base_size = base_size) +
  row_layout

pF <- ggplot(nuc, aes(y = row_label)) +
  band("candidate_rrna_fraction", 100) +
  geom_vline(xintercept = 2, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  dot(100 * candidate_rrna_fraction) +
  scale_x_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "Median rRNA\nfraction (%)") +
  row_layout + theme_bw(base_size = base_size) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = base_size - 1,
                                    colour = ink_primary, lineheight = 0.9),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(filename = 'p3A.png', pA, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p3B.png', pB, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p3C.png', pC, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p3D.png', pD, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p3E.png', pE, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p3F.png', pF, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)

fate_legend <- cowplot::get_plot_component(pC, "guide-box-bottom")
pC <- pC + theme(legend.position = "none")
key_text <- grid::textGrob(
  "Blue band: protoplast range     Dashed line: 2% cap  nuclei",
  gp = grid::gpar(fontsize = base_size - 1, col = ink_secondary))
key_row <- wrap_elements(full = cowplot::plot_grid(fate_legend, key_text, ncol = 1,
                                                   rel_heights = c(1, 0.8)))

rows <- (pB | pC | pD | pE | pF) +
  plot_layout(widths = c(0.9, 1.5, 1, 1, 1))

fig3 <- (pA / rows / key_row) +
  plot_layout(heights = c(0.62, 1, 0.08)) +
  plot_annotation(tag_levels = list(c(LETTERS[1:6], "")))

save_fig(fig3, "Fig3_nuclear_recovery", width_mm = 180, height_mm = 200, dir = out_dir)
