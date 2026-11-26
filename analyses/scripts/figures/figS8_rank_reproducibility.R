.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

rule_labels <- c(top100 = "Top 100", top200 = "Top 200", top500 = "Top 500",
                 `2fold_reproducible` = "All 2-fold reproducible")
cmp_labels <- c(rank_both_rho_within_nuclei_quarter     = "Within nuclei",
                rank_both_rho_within_protoplast_quarter = "Within protoplast",
                rank_both_rho_cross_quarter             = "Nuclei vs protoplast")

rho <- read_tsv_q(section_table("08_nf_response", "fig_refree_overlap_perdraw.tsv")) |>
  filter(scope == "without_root_hair", exclusion == "either_up", direction == "induced") |>
  pivot_longer(all_of(names(cmp_labels)), names_to = "comparison", values_to = "rho") |>
  group_by(rule, contrast, comparison) |>
  summarise(mean = mean(rho, na.rm = TRUE), sd = sd(rho, na.rm = TRUE), n = sum(!is.na(rho)),
            .groups = "drop") |>
  mutate(rule       = factor(rule, levels = names(rule_labels), labels = rule_labels),
         comparison = factor(comparison, levels = names(cmp_labels), labels = cmp_labels),
         col        = c("Within nuclei" = "nuclei", "Within protoplast" = "protoplast",
                        "Nuclei vs protoplast" = "across")[as.character(comparison)],
         contrast   = factor(contrast, levels = c("NOD", "BAP"), labels = c("Nod factor", "BAP")))
stopifnot(nrow(rho) == 24)
rho <- rho |> filter(n >= 5)

dodge <- position_dodge(width = 0.5)
p <- ggplot(rho, aes(comparison, mean, colour = col, shape = contrast)) +
  geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.3) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), position = dodge, linewidth = 0.4) +
  geom_point(position = dodge, size = 2, fill = "white", stroke = 0.6) +
  facet_wrap(~rule, nrow = 1) +
  scale_colour_manual(values = c(col_modality, across = ink_primary), guide = "none") +
  scale_shape_manual(values = c("Nod factor" = 16, "BAP" = 21), name = NULL) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Rank correlation of induced genes\n(Spearman \u03c1)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "top", legend.justification = "center")

save_fig(p, "FigS8_rank_reproducibility", width_mm = 180, height_mm = 80, dir = fig_supp_dir)
