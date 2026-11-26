.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s6 <- "06_cell_identity"
species_levels <- c("medicago", "arabidopsis")
pair_title <- function(x) str_replace(x, " x ", "\n\u00d7 ")

yardstick_roles <- c("yardstick_replicate_protoplast", "yardstick_replicate_nuclei", "yardstick_treatment")
keep_pairs <- function(df) df |> filter(
  primary_arm,
  kind == "modality" | (kind == "yardstick" & pair_role %in% yardstick_roles))
mt_across_levels <- c("nMOCK x pMock", "nNOD x pNOD", "nBAP x pBAP")

comparison_of <- function(kind, pair_id) case_when(
  kind == "modality" ~ "across",
  str_detect(pair_id, "Ryu|pMock") ~ "protoplast",
  TRUE ~ "nuclei")
cmp_labels <- c(across = "Nuclei vs protoplast", protoplast = "Protoplast vs protoplast",
                nuclei = "Nuclei vs nuclei")

per_type <- read_tsv_q(section_table(s6, "fig_centroid_per_type.tsv"), na = "") |>
  keep_pairs() |>
  mutate(species    = factor(species, levels = species_levels),
         matched    = centroid_matched(top1_reciprocal, at_null),
         comparison = factor(comparison_of(kind, pair_id), levels = names(cmp_labels)),
         cell_type  = factor(cell_type_labels[cell_type], levels = cell_type_labels))

mat <- read_tsv_q(section_table(s6, "fig_centroid_matrix.tsv"), na = "") |>
  keep_pairs() |>
  mutate(comparison = comparison_of(kind, pair_id),
         diagonal = type_row == type_col,   # V2 matrix has no `diagonal` column
         type_row = factor(cell_type_labels[type_row], levels = rev(cell_type_labels)),
         type_col = factor(cell_type_labels[type_col], levels = cell_type_labels))

stopifnot(n_distinct(mat$pair_id[mat$kind == "modality" & mat$species == "arabidopsis"]) == 9,
          n_distinct(mat$pair_id[mat$kind == "yardstick"]) == 7,
          !anyNA(mat$type_row), !anyNA(per_type$cell_type))

rho_lim <- max(abs(mat$rho))

centroid_heatmap <- function(m, d, ncol, xlab, ylab) {
  ggplot(m, aes(type_col, type_row, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_tile(data = d, aes(linetype = matched), fill = NA, colour = ink_primary,
              linewidth = 0.4, width = 0.9, height = 0.9) +
    facet_wrap(~panel, ncol = ncol, scales = "free") +
    scale_fill_gradient2(low = col_diverging[["low"]], mid = col_diverging[["mid"]],
                         high = col_diverging[["high"]], midpoint = 0,
                         limits = c(-rho_lim, rho_lim), name = "Spearman \u03c1") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                          labels = c(`TRUE` = "matched", `FALSE` = "not matched"),
                          name = "Same label", drop = FALSE) +
    labs(x = xlab, y = ylab) +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 2),
          axis.text.y = element_text(size = base_size - 2),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          strip.text = element_text(hjust = 0.5, size = base_size - 1.5, lineheight = 0.9))
}

diag_of <- function(m) m |>
  filter(diagonal) |>
  left_join(per_type |> transmute(pair_id, type_col = cell_type, matched), by = c("pair_id", "type_col"))

stopifnot(n_distinct(mat$pair_id[mat$kind == "modality" & mat$species == "medicago"]) == 3)

across <- mat |> filter(comparison == "across") |>
  mutate(panel = factor(pair_title(pair_id),
                        levels = pair_title(c(mt_across_levels, sort(unique(pair_id[species == "arabidopsis"]))))))
is_mt <- across$species == "medicago"
d_across <- diag_of(across)
pA_mt <- centroid_heatmap(droplevels(across[is_mt, ]), d_across[d_across$species == "medicago", ], 1,
                          "Protoplast cell type", "Nuclear cell type") +
  theme(legend.position = "none")
pA_at <- centroid_heatmap(droplevels(across[!is_mt, ]), d_across[d_across$species == "arabidopsis", ], 3,
                          "Protoplast cell type", NULL) +
  theme(legend.position = "right")

ryu <- mat |> filter(str_detect(pair_id, "^Ryu"))
ryu_mean <- ryu |>
  group_by(species, type_row, type_col, diagonal) |>
  summarise(rho = mean(rho), .groups = "drop") |>
  mutate(pair_id = "Ryu mean", panel = "Ryu rep1\u20133\n(mean of 3 pairs)")
ryu_diag <- ryu_mean |> filter(diagonal) |>
  left_join(per_type |> filter(str_detect(pair_id, "^Ryu")) |>
              group_by(type_col = cell_type) |> summarise(matched = all(matched), .groups = "drop"),
            by = "type_col")

within_levels <- c("pMock x pNOD", "nMOCK x nNOD", "Ryu rep1\u20133\n(mean of 3 pairs)",
                   "Farmer rep2 x Farmer rep5", "Farmer rep4 x Farmer rep5")
within <- bind_rows(
  mat |> filter(comparison != "across", !str_detect(pair_id, "^Ryu")) |> mutate(panel = pair_id),
  ryu_mean) |>
  mutate(panel = factor(pair_title(panel), levels = pair_title(within_levels)))
d_within <- bind_rows(diag_of(mat |> filter(comparison != "across", !str_detect(pair_id, "^Ryu")) |>
                                mutate(panel = pair_id)),
                      ryu_diag) |>
  mutate(panel = factor(pair_title(panel), levels = pair_title(within_levels)))
stopifnot(!anyNA(within$panel))

pB <- centroid_heatmap(within, d_within, 5, "Cell type, second library", "Cell type, first library") +
  theme(legend.position = "none")

pC <- ggplot(per_type, aes(cell_type, margin, colour = comparison)) +
  geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.3) +
  geom_point(aes(fill = if_else(matched, as.character(comparison), "open")),
             shape = 21, size = 1.6, stroke = 0.5,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75, seed = 1)) +
  facet_grid(~species, scales = "free_x", space = "free_x",
             labeller = labeller(species = as_labeller(species_expr, label_parsed))) +
  scale_colour_manual(values = c(across = ink_primary, col_modality), labels = cmp_labels, name = NULL) +
  scale_fill_manual(values = c(across = ink_primary, col_modality, open = "white"), guide = "none") +
  guides(colour = guide_legend(override.aes = list(fill = unname(c(ink_primary, col_modality[["protoplast"]],
                                                                   col_modality[["nuclei"]]))))) +
  labs(x = NULL, y = "Margin (same-label \u03c1 \u2212\nbest other-label \u03c1; open: not matched)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "top", legend.justification = "center")

layout <- "
AB
CC
DD
"
figS5 <- pA_mt + pA_at + pB + pC +
  plot_layout(design = layout, widths = c(0.3, 1), heights = c(2.6, 0.9, 1.1)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "C")))

save_fig(figS5, "FigS5_centroids", width_mm = 180, height_mm = 250, dir = fig_supp_dir)
