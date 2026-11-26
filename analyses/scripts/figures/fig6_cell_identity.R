.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))

out_dir <- file.path(fig_main_dir, "F6")
species_levels <- c("medicago", "arabidopsis")
species_strip  <- labeller(species = as_labeller(species_expr, label_parsed))

test_labels <- c(own = "Same library", within_modality = "Same modality",
                 cross_modality = "Other modality")
set_labels  <- c(single_gene = "Single gene", programme = "Programme")

markers <- read_tsv_q(section_table("06_cell_identity", "fig_markers_auroc_summary.tsv")) |>
  filter(primary, test_type %in% names(test_labels)) |>
  mutate(species    = factor(species, levels = species_levels),
         test_type  = factor(test_type, levels = names(test_labels), labels = test_labels),
         marker_set = factor(marker_set, levels = names(set_labels), labels = set_labels))

pC <- ggplot(markers, aes(test_type, auroc_mean, colour = source_modality)) +
  geom_hline(yintercept = 0.5, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6, seed = 1),
             size = 2, alpha = 0.85) +
  facet_grid(species ~ marker_set, labeller = labeller(species = as_labeller(species_expr, label_parsed))) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = "Markers from") +
  guides(colour = guide_legend(title.position = "top")) +
  scale_y_continuous(limits = c(0.45, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  labs(x = "Markers tested on", y = "Held-out marker AUROC") +
  theme_bw(base_size = base_size) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = ink_primary, angle = 40, hjust = 1,
                                   size = base_size - 1.5),
        strip.text.x = element_text(hjust = 0.5),
        panel.spacing = unit(2, "mm"),
        legend.position = "top", legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        legend.margin = margin(0, 0, -3, 0))

cen_mat <- read_tsv_q(section_table("06_cell_identity", "fig_centroid_matrix.tsv"), na = "") |>
  filter(primary_arm, kind == "modality")
cen_type <- read_tsv_q(section_table("06_cell_identity", "fig_centroid_per_type.tsv"), na = "") |>
  filter(primary_arm, kind == "modality") |>
  mutate(matched = centroid_matched(top1_reciprocal, at_null))

cen_avg <- function(sp, min_pairs = 3) {
  m <- cen_mat |> filter(species == sp)
  keep <- m |> distinct(pair_id, type_row) |> count(type_row) |>
    filter(n >= min(min_pairs, n_distinct(m$pair_id))) |> pull(type_row)
  matched <- cen_type |> filter(species == sp, cell_type %in% keep) |>
    group_by(cell_type) |> summarise(label = sprintf("%d/%d", sum(as.logical(matched)), n()), .groups = "drop")
  m |>
    filter(type_row %in% keep, type_col %in% keep) |>
    group_by(type_row, type_col) |>
    summarise(rho = mean(rho), n_pairs = n(), .groups = "drop") |>
    left_join(matched, by = c("type_row" = "cell_type")) |>
    mutate(label = if_else(type_row == type_col, label, NA_character_),
           type_row = factor(cell_type_labels[type_row], levels = rev(cell_type_labels)),
           type_col = factor(cell_type_labels[type_col], levels = cell_type_labels))
}
cen_mt <- cen_avg("medicago")
cen_at <- cen_avg("arabidopsis")
rho_lim <- max(abs(c(cen_mt$rho, cen_at$rho)))

centroid_panel <- function(d, sp) {
  d <- d |> mutate(species_lab = species_expr(sp))
  ggplot(d, aes(type_col, type_row, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = d |> filter(!is.na(label)), fill = NA, colour = ink_primary, linewidth = 0.35,
              width = 0.94, height = 0.94) +
    geom_text(aes(label = label), size = (base_size - 2) / .pt, colour = ink_primary, na.rm = TRUE) +
    scale_fill_gradient2(low = col_diverging[["low"]], mid = col_diverging[["mid"]],
                         high = col_diverging[["high"]], midpoint = 0,
                         limits = c(-rho_lim, rho_lim), breaks = c(-0.4, 0, 0.4),
                         name = "Spearman ρ") +
    facet_wrap(~species_lab, labeller = label_parsed) +
    labs(x = "Protoplast", y = "Nuclei") +
    theme(aspect.ratio = 1,
          strip.text = element_text(size = base_size, face = "plain", hjust = 0.5,
                                    margin = margin(b = 2)),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = ink_primary, size = base_size - 1.5),
          axis.text.y = element_text(colour = ink_primary, size = base_size - 1.5),
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.key.width = unit(2.5, "mm"), legend.key.height = unit(4, "mm"))
}

p_cen_mt <- centroid_panel(cen_mt, "medicago") +
  theme(legend.position = "none")
p_cen_at <- centroid_panel(cen_at, "arabidopsis") +
  theme(legend.position = "bottom", legend.key.width = unit(7, "mm"), legend.justification = "center",
        legend.key.height = unit(2.5, "mm"), legend.title = element_text(vjust = 0.8))

dir_labels <- c(nuclei_to_protoplast = "Nuclei \u2192 protoplast",
                protoplast_to_nuclei = "Protoplast \u2192 nuclei")

knn <- read_tsv_q(section_table("06_cell_identity", "fig_knn_purity_by_type.tsv")) |>
  filter(pair %in% knn_modality_pairs, condition == "base", shared_both_sides) |>
  mutate(direction = paste0(source_modality, "_to_", target_modality)) |>
  group_by(species, cell_type, direction, source_modality) |>
  summarise(purity_mean = mean(purity), purity_min = min(purity), purity_max = max(purity),
            null_mean = mean(null_purity), n_pairs = n(), .groups = "drop") |>
  mutate(species   = factor(species, levels = species_levels),
         cell_type = factor(cell_type_labels[cell_type], levels = rev(cell_type_labels)),
         direction = factor(direction, levels = names(dir_labels)))
stopifnot(!anyNA(knn$cell_type))

knn_dodge <- position_dodge(width = 0.6)

pB <- ggplot(knn, aes(y = cell_type, colour = source_modality, group = direction)) +
  geom_linerange(aes(xmin = purity_min, xmax = purity_max), position = knn_dodge,
                 linewidth = 0.5, alpha = 0.6) +
  geom_point(aes(x = null_mean), position = knn_dodge, shape = 21, fill = "white",
             size = 1.2, stroke = 0.4) +
  geom_point(aes(x = purity_mean), position = knn_dodge, size = 1.4) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y", labeller = species_strip) +
  scale_colour_manual(values = col_modality,
                      labels = c(nuclei = dir_labels[["nuclei_to_protoplast"]],
                                 protoplast = dir_labels[["protoplast_to_nuclei"]]),
                      name = "Query → neighbours") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(mult = 0.02)) +
  labs(x = "Same-label neighbours, mean (range) over pairs\n(open: shuffled null)", y = NULL) +
  guides(colour = guide_legend(ncol = 1, title.position = "top")) +
  theme_bw(base_size = base_size) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary),
        strip.text.y = element_text(angle = 0, hjust = 0),
        legend.position = "top", legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        legend.margin = margin(0, 0, -3, 0),
        legend.key.spacing.y = unit(0, "mm"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pA <- (p_cen_mt / p_cen_at) + plot_layout(heights = c(2.7, 8.5))
panel_files <- list(p6A = list(pA, 75, 150), p6B = list(pB, 100, 90), p6C = list(pC, 100, 80))
for (nm in names(panel_files)) {
  x <- panel_files[[nm]]
  ggsave(file.path(out_dir, paste0(nm, ".png")), x[[1]], width = x[[2]], height = x[[3]],
         units = "mm", dpi = 600, bg = "white")
}

left  <- (p_cen_mt / p_cen_at) + plot_layout(heights = c(2.7, 8.5))
right <- (pB / pC) + plot_layout(heights = c(1.15, 1))

fig6 <- (left | right) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "C")))

save_fig(fig6, "Fig6_cell_identity", width_mm = 180, height_mm = 175, dir = out_dir)
