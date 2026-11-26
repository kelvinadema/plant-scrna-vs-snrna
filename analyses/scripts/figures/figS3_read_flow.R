.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
theme_set(theme_ch3_bw())

s3 <- "03_recovery_attrition"
prep_short <- c("Unsorted\nSZ nuclei", "Tip nuclei", "FANS\nSZ nuclei", "SZ\nprotoplasts",
                "Tip\nprotoplasts", "Farmer\nnuclei", "Ryu\nprotoplasts")

greys <- c("#e6e5e1", "#c4c3bd", "#9d9c96", "#6f6e69", "#3f3e3a")

tidy_libs <- function(df) {
  df |>
    filter(is_used_library(sample, species, modality)) |>
    mutate(group = factor(prep_groups[arm], levels = prep_groups),
           row_label = paste(if_else(species == "medicago", "Mt", "At"),
                             if_else(species == "arabidopsis", sample_label, sample))) |>
    arrange(group, plot_order) |>
    mutate(row_label = factor(row_label, levels = rev(unique(row_label))))
}

readfate <- read_tsv_q(section_table(s3, "fig_readfate_stack.tsv")) |> tidy_libs()
stages   <- read_tsv_q(section_table(s3, "fig_qc_stage_losses.tsv")) |> tidy_libs()
stopifnot(n_distinct(readfate$sample) == 29, n_distinct(stages$sample) == 29)

flow_panel <- function(df, value, levels, labels, kept, x_title) {
  loss_levels <- setdiff(levels, kept)
  fill_vals <- c(setNames(greys[seq_along(loss_levels)], loss_levels), col_modality)
  fill_labs <- c(labels[loss_levels], protoplast = paste(labels[[kept]], "(protoplast)"),
                 nuclei = paste(labels[[kept]], "(nuclei)"))
  df |>
    mutate(stage_f = factor(.data[[value]], levels = levels),
           fill = if_else(.data[[value]] == kept, modality, .data[[value]]),
           fill = factor(fill, levels = c(loss_levels, "nuclei", "protoplast"))) |>
    ggplot(aes(100 * share, row_label, fill = fill, group = stage_f)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2,
             position = position_stack(reverse = TRUE)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(group = setNames(prep_short, prep_groups))) +
    scale_fill_manual(values = fill_vals, labels = fill_labs, name = NULL, drop = TRUE) +
    scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0.04))) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = x_title, y = NULL) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
          panel.spacing.y = unit(1.2, "mm"),
          axis.ticks.y = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "bottom", legend.justification = "left",
          legend.key.size = unit(2.5, "mm"),
          plot.margin = margin(3, 2, 3, 2, "mm"))
}

p_a <- readfate |>
  filter(partition == "input_reads") |>
  rename(share = fraction) |>
  flow_panel("fate",
             levels = c("unmapped", "barcode_loss", "duplicate_reads", "deduplicated_molecules"),
             labels = c(unmapped = "Unmapped", barcode_loss = "Barcode not recovered",
                        duplicate_reads = "Duplicate reads",
                        deduplicated_molecules = "Deduplicated molecules"),
             kept = "deduplicated_molecules", x_title = "Input read pairs (%)") +
  theme(axis.text.y = element_text(colour = ink_primary))

p_b <- readfate |>
  filter(partition == "raw_umis") |>
  rename(share = fraction) |>
  flow_panel("fate",
             levels = c("outside_cb_model", "removed_within_cb_model", "corrected_below_count_gate",
                        "corrected_later_qc_removed", "retained"),
             labels = c(outside_cb_model = "Empty droplets (outside CellBender window)",
                        removed_within_cb_model = "Removed within CellBender window",
                        corrected_below_count_gate = "Below cell-calling cutoff",
                        corrected_later_qc_removed = "In barcodes removed by QC",
                        retained = "Retained"),
             kept = "retained", x_title = "Molecules in all barcodes (%)") +
  theme(axis.text.y = element_blank())

p_c <- stages |>
  rename(share = fraction_of_candidates) |>
  flow_panel("stage",
             levels = c("MAD_or_miQC", "absolute_caps", "regression", "other", "manual", "retained"),
             labels = c(MAD_or_miQC = "Adaptive (MAD / miQC)", absolute_caps = "Absolute caps",
                        regression = "Regression", other = "Other filters",
                        manual = "Manual / doublets", retained = "Retained"),
             kept = "retained", x_title = "Candidate barcodes (%)") +
  theme(axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = base_size - 1, lineheight = 0.9))

fig <- (p_a | p_b | p_c) + plot_annotation(tag_levels = "A")

save_fig(fig, "FigS3_read_flow", width_mm = 180, height_mm = 185, dir = fig_supp_dir)
