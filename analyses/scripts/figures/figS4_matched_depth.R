.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s4 <- "04_read_yield"
species_levels <- c("medicago", "arabidopsis")
species_strip  <- labeller(species = as_labeller(species_expr, label_parsed))
short_lib <- function(x) x |>
  str_replace("GSE155304_sNucRNA-seq_rep", "Farmer rep") |>
  str_replace("GSE123013_WT_rep", "Ryu rep")

pairs <- read_tsv_q(section_table(s4, "fig_matched_yield_pairs.tsv")) |>
  mutate(species = factor(species, levels = species_levels),
         pair = paste(short_lib(nuclei_sample), "\u00d7", short_lib(protoplast_sample)),
         pair = factor(pair, levels = rev(unique(pair))))
stopifnot(nrow(pairs) == 12)

pairs_long <- pairs |>
  select(species, pair, `Matched reads per cell` = yield_ratio_matched, `Full depth` = yield_ratio_full) |>
  pivot_longer(-c(species, pair), names_to = "depth", values_to = "ratio") |>
  mutate(depth = factor(depth, levels = c("Matched reads per cell", "Full depth")))

pA <- ggplot(pairs, aes(y = pair)) +
  geom_vline(xintercept = 1, colour = ink_muted, linewidth = 0.3) +
  geom_segment(aes(x = yield_ratio_matched, xend = yield_ratio_full, yend = pair),
               colour = ink_muted, linewidth = 0.4) +
  geom_point(data = pairs_long, aes(x = ratio, fill = depth), shape = 21, size = 2,
             stroke = 0.5, colour = ink_primary) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y", labeller = species_strip) +
  scale_fill_manual(values = c("Matched reads per cell" = ink_primary, "Full depth" = "white"),
                    name = NULL) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  labs(x = "Protoplast / nuclei informative yield", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = -90),
        legend.position = "top", legend.justification = "center")

matched <- read_tsv_q(section_table(s4, "fig_detection_matched.tsv")) |>
  mutate(species  = factor(species, levels = species_levels),
         modality = factor(modality, levels = c("nuclei", "protoplast"))) |>
  arrange(species, modality, sample_label) |>
  mutate(sample_label = factor(sample_label, levels = rev(unique(sample_label))))
stopifnot(nrow(matched) == 12)

depth_lab <- matched |> distinct(species, matched_depth)
pB <- ggplot(matched, aes(genes_mean, sample_label, colour = modality)) +
  geom_linerange(aes(xmin = genes_q25, xmax = genes_q75), linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 2) +
  geom_text(aes(x = Inf, label = sprintf("%.0f%%", 100 * frac_eligible)), hjust = 1.1,
            size = (base_size - 1.5) / .pt, colour = ink_secondary) +
  facet_wrap(~species, ncol = 1, scales = "free",
             labeller = labeller(species = as_labeller(function(x) {
               d <- depth_lab$matched_depth[match(x, depth_lab$species)]
               sprintf("%s ~ '(%s UMIs per cell)'", species_expr(x), scales::comma(d))
             }, label_parsed))) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL,
                      breaks = c("nuclei", "protoplast")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  labs(x = "Genes detected per cell at matched UMI depth", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top", legend.justification = "center")

figS4 <- (pA | pB) + plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = "A")

save_fig(figS4, "FigS4_matched_depth", width_mm = 180, height_mm = 110, dir = fig_supp_dir)
