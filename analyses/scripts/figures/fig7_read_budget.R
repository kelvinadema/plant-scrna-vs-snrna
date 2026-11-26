.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))

out_dir <- file.path(fig_main_dir, "F7")
species_levels <- c("medicago", "arabidopsis")
species_strip  <- labeller(species = as_labeller(species_expr, label_parsed))
target <- 0.70

sample_label_of <- function(s) {
  s |>
    str_replace("^GSE123013_WT_", "Ryu ") |>
    str_replace("^GSE155304_sNucRNA-seq_", "Farmer ")
}

curve <- read_tsv_q(section_table("07_read_budget", "fig_auroc_vs_reads.tsv")) |>
  filter(layer == "corrected", !technical_draw, !is.na(auroc)) |>
  mutate(species = factor(species, levels = species_levels),
         sample_label = sample_label_of(sample))

crossing <- read_tsv_q(section_table("07_read_budget", "fig_target_crossings.tsv")) |>
  filter(layer == "corrected", target == "absolute_0.70") |>
  distinct(species, sample, modality, reads_per_cell, reached, status) |>
  mutate(species = factor(species, levels = species_levels))

stopifnot(n_distinct(curve$sample) == 11,
          setequal(setdiff(curve$sample, crossing$sample), "pNOD"),
          all(!is.na(crossing$reads_per_cell[crossing$reached])))

ends <- curve |>
  filter(used_for_crossing) |>
  group_by(species, sample, sample_label, modality) |>
  slice_max(reads_per_cell, n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(crossing |> select(sample, reached), by = "sample") |>
  mutate(label = case_when(is.na(reached) ~ paste0(sample_label, " (no pair)"),
                           reached        ~ sample_label,
                           TRUE           ~ paste0(sample_label, " (n.r.)")))

pA <- ggplot(curve, aes(reads_per_cell, auroc, colour = modality)) +
  geom_hline(yintercept = target, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_line(data = curve |> filter(used_for_crossing), aes(group = sample_label),
            linewidth = 0.45, alpha = 0.9) +
  geom_linerange(aes(ymin = auroc - auroc_split_sd, ymax = auroc + auroc_split_sd),
                 linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(fill = if_else(used_for_crossing, modality, "excluded")),
             shape = 21, size = 1.2, stroke = 0.4) +
  geom_point(data = crossing |> filter(reached), aes(x = reads_per_cell, y = target),
             shape = 124, size = 2.5) +
  ggrepel::geom_text_repel(data = ends, aes(label = label), size = (base_size - 2) / .pt,
                           direction = "y", hjust = 0, nudge_x = 0.08, segment.size = 0.2,
                           segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.1,
                           show.legend = FALSE, seed = 1) +
  facet_wrap(~species, nrow = 1, labeller = species_strip) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL,
                      breaks = c("nuclei", "protoplast")) +
  scale_fill_manual(values = c(col_modality, excluded = "white"), guide = "none") +
  scale_x_log10(breaks = c(1e3, 3e3, 1e4, 3e4, 1e5, 3e5),
                labels = c("1k", "3k", "10k", "30k", "100k", "300k"),
                expand = expansion(mult = c(0.03, 0.22))) +
  scale_y_continuous(limits = c(0.55, 0.87), breaks = seq(0.55, 0.85, 0.05)) +
  labs(x = "Reads per cell (subsampled)", y = "Held-out marker AUROC") +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.6))) +
  theme_bw() +
  theme(panel.spacing = unit(5, "mm"),
        strip.text = element_text(face = "plain", hjust = 0),
        legend.position = "top", legend.justification = "center",
        legend.margin = margin(0, 0, -3, 0))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "p7A.png"), pA, width = 180, height = 80, units = "mm", dpi = 600, bg = "white")
save_fig(pA, "Fig7_read_budget", width_mm = 180, height_mm = 80, dir = out_dir)
