.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s3  <- "03_recovery_attrition"

group_labels <- c(protoplast = "Protoplasts", failure_labels)

nuclear <- read_failure_modes() |>
  select(sample, sample_label, species, plot_order, group = failure_mode)
protoplast <- read_tsv_q(section_table(s3, "fig_mapping_readfate_perlib.tsv")) |>
  filter(modality == "protoplast", is_used_library(sample, species, modality)) |>
  transmute(sample, sample_label, species, plot_order, group = "protoplast")

libs <- bind_rows(protoplast, nuclear) |>
  mutate(group     = factor(group, levels = names(group_labels)),
         species   = factor(species, levels = c("medicago", "arabidopsis")),
         row_label = if_else(species == "arabidopsis", sample_label, sample)) |>
  arrange(group, species, plot_order) |>
  mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

stopifnot(nrow(libs) == 29, !anyNA(libs$group))

metric_labels <- c(mito_fraction = "Mitochondrial", chloroplast_fraction = "Chloroplast",
                   rrna_fraction = "rRNA")

cells <- read_tsv_q(section_table(s3, "fig_qc_percell.tsv.gz"),
                    col_select = c(sample, stage, mito_fraction, chloroplast_fraction,
                                   rrna_fraction, mito_cap)) |>
  filter(stage == "count_gate_candidate") |>
  inner_join(libs, by = "sample") |>
  pivot_longer(names(metric_labels), names_to = "metric", values_to = "fraction") |>
  filter(!is.na(fraction)) |>
  mutate(metric   = factor(metric, levels = names(metric_labels)),
         modality = if_else(group == "protoplast", "protoplast", "nuclei"))

stopifnot(all(cells$mito_cap[cells$modality == "nuclei"] == 0.02),
          all(cells$mito_cap[cells$modality == "protoplast"] == 0.05))
caps <- cells |> distinct(group, cap = mito_cap)

above_cap <- cells |>
  group_by(row_label, group, metric) |>
  summarise(n = n(), frac_above = mean(fraction > mito_cap), .groups = "drop") |>
  mutate(label = sprintf("%.0f%%", 100 * frac_above))

n_cands <- cells |> filter(metric == "mito_fraction") |> count(row_label)

p <- ggplot(cells, aes(100 * fraction, row_label)) +
  geom_vline(data = caps, aes(xintercept = 100 * cap), linetype = "22", colour = ink_muted,
             linewidth = 0.3) +
  geom_boxplot(aes(fill = modality), outlier.shape = NA, width = 0.62, linewidth = 0.3, coef = 1.5,
               colour = ink_secondary, alpha = 0.55) +
  geom_text(data = above_cap, aes(x = Inf, y = row_label, label = label),
            hjust = 1.1, size = (base_size - 1.5) / .pt, colour = ink_secondary) +
  facet_grid(group ~ metric, scales = "free_y", space = "free_y",
             labeller = labeller(group = group_labels, metric = metric_labels)) +
  scale_fill_manual(values = col_modality, guide = "none") +
  scale_x_sqrt(breaks = c(0, 2, 5, 10, 25, 50), expand = expansion(mult = c(0.01, 0.12))) +
  coord_cartesian(xlim = c(0, 60)) +   # zoom only; box statistics use all cells
  scale_y_discrete(labels = function(x) {
    n <- n_cands$n[match(x, n_cands$row_label)]
    sprintf("%s (%s)", x, scales::comma(n))
  }) +
  labs(x = "Fraction of corrected counts per candidate barcode (%)",
       y = "Library (number of candidate barcodes)") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
        panel.spacing.y = unit(1.2, "mm"),
        panel.spacing.x = unit(3, "mm"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary),
        strip.text.x = element_text(hjust = 0.5),
        strip.text.y = element_text(angle = 0, hjust = 0, size = base_size - 1, lineheight = 0.9))

save_fig(p, "FigS2_qc_fractions", width_mm = 180, height_mm = 170, dir = fig_supp_dir)
