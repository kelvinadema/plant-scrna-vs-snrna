.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

pair_labels <- c(
  ls__pMock__nMOCK = "Mt \u00b7 mock",
  ls__pNOD__nNOD   = "Mt \u00b7 Nod factor",
  ls__pBAP__nBAP   = "Mt \u00b7 BAP (nBAP flagged)",
  `ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep2` = "At \u00b7 Ryu rep1 \u00d7 Farmer rep2",
  `ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep4` = "At \u00b7 Ryu rep1 \u00d7 Farmer rep4",
  `ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep5` = "At \u00b7 Ryu rep1 \u00d7 Farmer rep5",
  `ls__GSE123013_WT_rep2__GSE155304_sNucRNA-seq_rep2` = "At \u00b7 Ryu rep2 \u00d7 Farmer rep2\n(not quoted)"
)
dir_labels <- c(nuclei = "Nuclei \u2192 protoplast", protoplast = "Protoplast \u2192 nuclei")

knn <- read_tsv_q(section_table("06_cell_identity", "fig_knn_purity_by_type.tsv")) |>
  filter(condition == "base", kind == "modality", shared_both_sides) |>
  mutate(pair      = factor(pair, levels = names(pair_labels), labels = pair_labels),
         cell_type = factor(cell_type_labels[cell_type], levels = rev(cell_type_labels)),
         above     = delta_purity > 0 & z_purity >= 2)

stopifnot(!anyNA(knn$pair), !anyNA(knn$cell_type), n_distinct(knn$pair) == 7)

p <- ggplot(knn, aes(delta_purity, cell_type, colour = source_modality)) +
  geom_vline(xintercept = 0, colour = ink_muted, linewidth = 0.3) +
  geom_point(aes(fill = if_else(above, source_modality, "white")), shape = 21, size = 1.8,
             stroke = 0.6, position = position_dodge(width = 0.55)) +
  facet_wrap(~pair, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = col_modality, labels = dir_labels, name = "Query \u2192 neighbours") +
  scale_fill_manual(values = c(col_modality, white = "white"), guide = "none") +
  guides(colour = guide_legend(override.aes = list(fill = unname(col_modality[c("nuclei", "protoplast")])))) +
  scale_x_continuous(limits = c(-0.15, 0.85), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = "Same-label neighbours above shuffled null (open: z < 2)", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary),
        strip.text = element_text(face = "plain", size = base_size - 1),
        legend.position = "top", legend.justification = "center")

pooled <- read_tsv_q(section_table("06_cell_identity", "fig_knn_pooled_by_direction.tsv")) |>
  filter(condition == "base") |>
  rename(pair_kind = kind) |>
  group_by(species, pair, pair_kind, source_modality) |>
  summarise(delta = mean(delta_purity), .groups = "drop") |>
  group_by(species, pair, pair_kind) |>
  summarise(delta = mean(delta),
            kind = case_when(first(pair_kind) == "modality" ~ "Across modalities",
                             all(source_modality == "protoplast") ~ "Within protoplast",
                             TRUE ~ "Within nuclei"),
            .groups = "drop") |>
  mutate(species  = factor(species, levels = c("medicago", "arabidopsis")),
         kind     = factor(kind, levels = c("Across modalities", "Within protoplast", "Within nuclei")),
         kind_col = c("Across modalities" = "across", "Within protoplast" = "protoplast",
                      "Within nuclei" = "nuclei")[as.character(kind)],
         quoted   = pair %in% knn_modality_pairs |
                    (pair_kind == "yardstick" & str_detect(pair, "GSE123013_WT_rep1__GSE123013_WT_rep2|pMock__pNOD")))
stopifnot(nrow(pooled) == 12, sum(pooled$quoted) == 8)

p_b <- ggplot(pooled, aes(kind, delta, colour = kind_col)) +
  geom_point(aes(fill = if_else(quoted, kind_col, "open")), shape = 21, size = 2, stroke = 0.6,
             position = position_jitter(width = 0.12, seed = 1)) +
  facet_wrap(~species, nrow = 1, labeller = labeller(species = as_labeller(species_expr, label_parsed))) +
  scale_colour_manual(values = c(across = ink_primary, col_modality), guide = "none") +
  scale_fill_manual(values = c(across = ink_primary, col_modality, open = "white"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Agreement above null\n(pooled over cell types)") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

figS6 <- (p / p_b) + plot_layout(heights = c(1.7, 1)) + plot_annotation(tag_levels = "A")

save_fig(figS6, "FigS6_knn_all_pairs", width_mm = 180, height_mm = 170, dir = fig_supp_dir)
