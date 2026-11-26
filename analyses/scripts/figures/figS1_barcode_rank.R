.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s3 <- "03_recovery_attrition"

prep_short <- c("unsorted", "tip", "FANS", "SZ", "tip", "Farmer", "Ryu")

libs <- read_tsv_q(section_table(s3, "fig_mapping_readfate_perlib.tsv")) |>
  filter(is_used_library(sample, species, modality)) |>
  left_join(read_failure_modes() |> select(sample, failure_mode),
            by = "sample") |>
  mutate(group = factor(prep_groups[arm], levels = prep_groups),
         panel = sprintf("%s %s  · %s",
                         if_else(species == "medicago", "Mt", "At"),
                         if_else(species == "arabidopsis", sample_label, sample),
                         prep_short[as.integer(group)])) |>
  arrange(group, plot_order) |>
  mutate(panel = factor(panel, levels = unique(panel)))

stopifnot(nrow(libs) == 29)

curves <- read_tsv_q(section_table(s3, "fig_barcode_rank_curves.tsv")) |>
  filter(in_raw_downsample, raw_umis > 0) |>
  inner_join(libs |> select(sample, panel), by = "sample")

annot <- libs |>
  mutate(counts = paste(scales::comma(barcodes_count_gate), "\u2192",
                        scales::comma(cells_retained), "cells"),
         mode   = if_else(modality == "nuclei",
                          str_replace_all(failure_labels[failure_mode], "\n", " "), NA_character_))

outcome_levels <- c(protoplast = "Protoplasts",
                    retained   = "Nuclei, retained",
                    no_sep     = "Nuclei, no droplet separation",
                    qc_loss    = "Nuclei, lost at QC\norganellar RNA, rRNA, low yield")
outcomes <- libs |>
  transmute(sample, species,
            outcome = case_when(modality == "protoplast" ~ "protoplast",
                                failure_mode == "retained" ~ "retained",
                                failure_mode == "no_droplet_separation" ~ "no_sep",
                                TRUE ~ "qc_loss"),
            outcome = factor(outcome_levels[outcome], levels = outcome_levels))
stopifnot(!anyNA(outcomes$outcome))
outcome_n <- outcomes |> count(outcome)

pA <- curves |>
  inner_join(outcomes |> select(sample, outcome), by = "sample") |>
  mutate(species = factor(lab_species[species], levels = lab_species)) |>
  ggplot(aes(rank, raw_umis, group = sample)) +
  geom_hline(yintercept = 500, linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_line(aes(colour = modality, linetype = species), linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~outcome, nrow = 1, labeller = as_labeller(function(x)
    str_replace(x, "^[^\n]*", function(first) sprintf("%s (n = %d)", first,
                                                     outcome_n$n[match(x, outcome_n$outcome)])))) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL, guide = "none") +
  scale_linetype_manual(values = c("M. truncatula" = "solid", "A. thaliana" = "42"), name = NULL) +
  scale_x_log10(breaks = c(1, 1e2, 1e4), labels = scales::label_log(), limits = c(1, 5e5)) +
  scale_y_log10(breaks = c(1e1, 1e3, 1e5), labels = scales::label_log(), limits = c(1, 5e5)) +
  labs(x = "Barcode rank", y = "Raw UMIs per barcode") +
  guides(linetype = guide_legend(override.aes = list(colour = ink_secondary, linewidth = 0.5))) +
  theme(strip.text = element_text(size = base_size - 1, hjust = 0.5, lineheight = 0.9),
        legend.position = "top", legend.justification = "center",
        legend.text = element_text(face = "italic"),
        panel.spacing = unit(2, "mm"))

p <- ggplot(curves, aes(rank, raw_umis)) +
  geom_hline(data = libs, aes(yintercept = count_gate_threshold),
             linetype = "22", colour = ink_muted, linewidth = 0.3) +
  geom_line(aes(colour = modality), linewidth = 0.5) +
  geom_text(data = annot, aes(x = 1.3, y = 1.3, label = counts), hjust = 0, vjust = 0,
            size = (base_size - 1.5) / .pt, colour = ink_secondary) +
  geom_text(data = annot |> filter(!is.na(mode)), aes(x = 1.3, y = 6, label = mode),
            hjust = 0, vjust = 0, size = (base_size - 1.5) / .pt, colour = ink_primary,
            fontface = "bold") +
  facet_wrap(~panel, ncol = 5) +
  scale_colour_manual(values = col_modality, labels = lab_modality, name = NULL) +
  scale_x_log10(breaks = c(1, 1e2, 1e4), labels = scales::label_log(), limits = c(1, 5e5)) +
  scale_y_log10(breaks = c(1e1, 1e3, 1e5), labels = scales::label_log(), limits = c(1, 5e5)) +
  labs(x = "Barcode rank", y = "Raw UMIs per barcode") +
  theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25),
        strip.text = element_text(size = base_size - 1.5, face = "bold", hjust = 0),
        legend.position = "top", legend.justification = "left",
        panel.spacing = unit(2, "mm"))

figS1 <- (pA / p) + plot_layout(heights = c(0.22, 1)) + plot_annotation(tag_levels = "A")

save_fig(figS1, "FigS1_barcode_rank", width_mm = 180, height_mm = 250, dir = fig_supp_dir)
