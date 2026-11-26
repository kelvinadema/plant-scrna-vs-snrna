.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))

s8 <- "08_nf_response"
out_dir <- file.path(fig_main_dir, "F8")
mod_levels <- c("nuclei", "protoplast")
dodge <- position_dodge(width = 0.5)

panel_theme <- theme(legend.position = "none",
                     axis.ticks.x = element_blank(),
                     strip.text = element_text(face = "plain", hjust = 0.5))

recall <- read_tsv_q(section_table(s8, "fig_recall_auroc_direction_summary.tsv")) |>
  filter(arm == "matched_cells_UMIs", exclusion == "either_up", universe == "per_modality")

metric_labels <- c(auroc_up = "AUROC, induced genes", dir_agree = "Direction agreement")
recall_long <- bind_rows(
  recall |> transmute(scope, modality, metric = "auroc_up",  value = auroc_up,  sd = auroc_up_draw_sd),
  recall |> transmute(scope, modality, metric = "dir_agree", value = dir_agree, sd = dir_agree_draw_sd)) |>
  mutate(metric   = factor(metric, levels = names(metric_labels), labels = metric_labels),
         scope    = factor(scope, levels = c("without_root_hair", "all_cells"),
                           labels = c("Without\nroot hairs", "All cells")),
         modality = factor(modality, levels = mod_levels))

pA <- ggplot(recall_long, aes(scope, value, colour = modality)) +
  geom_hline(yintercept = 0.5, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_linerange(aes(ymin = value - sd, ymax = value + sd), position = dodge, linewidth = 0.4) +
  geom_point(position = dodge, size = 2) +
  facet_wrap(~metric) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_y_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = "Recovery of bulk NF response") +
  theme_bw() +
  theme(legend.position = "top", legend.justification = "center",
        legend.margin = margin(0, 0, -3, 0),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "plain", hjust = 0.5))

rule_labels <- c(top100 = "Top\n100", top200 = "Top\n200", top500 = "Top\n500",
                 `2fold_reproducible` = "All 2-fold\nreproducible")
overlap <- read_tsv_q(section_table(s8, "fig_refree_overlap_summary.tsv")) |>
  filter(scope == "without_root_hair", exclusion == "either_up", direction == "induced") |>
  mutate(rule = factor(rule, levels = names(rule_labels), labels = rule_labels),
         contrast = factor(contrast, levels = c("NOD", "BAP"), labels = c("Nod factor", "BAP")))

pB <- ggplot(overlap, aes(rule, enrichment, shape = contrast, group = contrast)) +
  geom_hline(yintercept = 1, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_line(colour = ink_secondary, linewidth = 0.4) +
  geom_point(size = 2, colour = ink_primary, fill = "white", stroke = 0.6) +
  scale_shape_manual(values = c("Nod factor" = 16, "BAP" = 21), name = NULL) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20), limits = c(0.9, 30)) +
  labs(x = "Induced genes compared", y = "Shared induced genes\n(fold over chance)") +
  theme_bw() +
  theme(legend.position = "top", legend.justification = "center",
        legend.margin = margin(0, 0, -3, 0),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = ink_primary, lineheight = 0.9))

genes <- read_tsv_q(section_table(s8, "fig_gene_lfc_single_cell_vs_bulk.tsv")) |>
  filter(a17_set == "up", scope == "without_root_hair", frac_draws_in_universe >= 0.5) |>
  mutate(modality = factor(modality, levels = mod_levels))

gene_n <- genes |> count(modality)
lim <- range(c(genes$bulk_log2fc, genes$sc_log2fc_mean))

show_cor <- "spearman"
rho_tab <- read_tsv_q(section_table(s8, "fig_within_set_rho_summary.tsv")) |>
  filter(arm == "matched_cells_UMIs", universe == "per_modality", scope == "without_root_hair",
         gene_set == "induced_311_signature") |>
  select(modality, rho)

cor_lab <- genes |>
  group_by(modality) |>
  summarise(r = cor(bulk_log2fc, sc_log2fc_mean, method = "pearson"), .groups = "drop") |>
  left_join(rho_tab, by = "modality") |>
  mutate(label = paste(c(if ("spearman" %in% show_cor) sprintf("Mean Spearman ρ = %.2f", rho),
                         if ("pearson" %in% show_cor) sprintf("Pearson r = %.2f", r)),
                       collapse = "\n"),
         .by = modality)

pC <- ggplot(genes, aes(bulk_log2fc, sc_log2fc_mean, colour = modality)) +
  geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_text(data = cor_lab, aes(x = lim[1] + 0.2, y = lim[2] - 0.2, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, lineheight = 0.95,
            size = (base_size - 1) / .pt, colour = ink_primary) +
  facet_wrap(~modality, labeller = as_labeller(function(x)
    sprintf("%s (%d plotted)", lab_modality[x], gene_n$n[match(x, gene_n$modality)]))) +
  scale_colour_manual(values = col_modality, guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("Bulk NF response (log"[2] * "FC)"),
       y = expression("Single-cell NF response (log"[2] * "FC)")) +
  theme_bw() +
  theme(strip.text = element_text(face = "plain", hjust = 0.5),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
panel_files <- list(p8A = list(pA, 90, 75), p8B = list(pB, 70, 75), p8C = list(pC, 110, 65))
for (nm in names(panel_files)) {
  x <- panel_files[[nm]]
  ggsave(file.path(out_dir, paste0(nm, ".png")), x[[1]], width = x[[2]], height = x[[3]],
         units = "mm", dpi = 600, bg = "white")
}

row2 <- (pB | pC) + plot_layout(widths = c(1, 1.6))

fig8 <- (pA / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

save_fig(fig8, "Fig8_nf_response", width_mm = 180, height_mm = 150, dir = out_dir)
