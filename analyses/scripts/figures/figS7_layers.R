.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s8 <- "08_nf_response"
mod_levels <- c("nuclei", "protoplast")

layers <- read_tsv_q(section_table(s8, "fig_layers_auroc_summary.tsv")) |>
  filter(arm == "B_matched_resolved", scope == "without_root_hair", contrast == "NOD",
         component %in% c("res_s", "res_u")) |>
  mutate(component = factor(component, levels = c("res_s", "res_u"), labels = c("Spliced", "Unspliced")),
         modality  = factor(modality, levels = mod_levels))
stopifnot(nrow(layers) == 4)

p <- ggplot(layers, aes(component, auroc_up, colour = modality, group = modality)) +
  geom_hline(yintercept = 0.5, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  geom_linerange(aes(ymin = auroc_up - auroc_up_draw_sd, ymax = auroc_up + auroc_up_draw_sd),
                 linewidth = 0.4) +
  geom_point(size = 2) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_y_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = "AUROC, bulk NF-induced genes") +
  theme(legend.position = "top", legend.justification = "center",
        axis.ticks.x = element_blank())

save_fig(p, "FigS7_layers", width_mm = 70, height_mm = 75, dir = fig_supp_dir)
