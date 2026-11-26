.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

s6 <- "06_cell_identity"
role_label <- c(primary = "quoted", robustness = "quoted",
                robustness_flagged = "quoted (flagged library)",
                yardstick_replicate_protoplast = "yardstick: protoplast replicates",
                yardstick_replicate_nuclei = "yardstick: nuclear replicates",
                yardstick_treatment = "yardstick: same modality, other treatment",
                yardstick_treatment_flagged = "yardstick: same modality, other treatment (flagged)")

per_type <- read_tsv_q(section_table(s6, "fig_centroid_per_type.tsv"), na = "") |>
  filter(primary_arm, kind == "modality") |>
  mutate(matched = centroid_matched(top1_reciprocal, at_null)) |>
  arrange(desc(species), pair_role != "primary", pair_id,
          factor(cell_type, levels = names(cell_type_labels))) |>
  transmute(Species = if_else(species == "medicago", "M. truncatula", "A. thaliana"),
            Pair = str_replace(pair_id, " x ", " \u00d7 "), Role = role_label[pair_role],
            `Cell type` = cell_type_labels[cell_type],
            `Same-label rho` = rho_match, Margin = margin, `Margin draw SD` = margin_draw_sd,
            `Reciprocal top-1 (fraction of draws)` = top1_reciprocal,
            `Null margin mean` = null_mean, `Null margin SD` = null_sd, `z vs null` = z_null,
            `At null (|z| < 2)` = if_else(as.logical(at_null), "yes", "no"),
            `Yardstick margin` = margin_yardstick, `Margin / yardstick` = ratio_to_yardstick,
            Matched = if_else(as.logical(matched), "yes", "no"),
            `Cells (min of both sides)` = n_cells, `Cell types in pair (K)` = K)

pair_level <- read_tsv_q(section_table(s6, "fig_centroid_pair_level.tsv"), na = "") |>
  filter(primary_arm) |>
  arrange(desc(species), kind != "modality", pair_role != "primary", pair_id) |>
  transmute(Species = if_else(species == "medicago", "M. truncatula", "A. thaliana"),
            Pair = str_replace(pair_id, " x ", " \u00d7 "),
            Kind = if_else(kind == "modality", "nuclei \u00d7 protoplast", "within modality"),
            Role = role_label[pair_role], `Cell types (K)` = K,
            `Median margin` = median_margin, `Fraction of types with margin > 0` = frac_margin_pos,
            `Reciprocal top-1 accuracy` = reciprocal_top1_acc,
            `Null median margin` = null_median_margin, `z vs null` = z_null,
            `Yardstick median margin` = yardstick_median_margin, `Margin / yardstick` = ratio_to_yardstick)

stopifnot(!anyNA(per_type$`Cell type`), !anyNA(per_type$Role), !anyNA(pair_level$Role))

readme <- tibble(
  Field = c("Title", "Method", "Margin", "Matched", "Yardstick", "Quoted pairs", "Caveat"),
  Description = c(
    "Supplementary Table S5. Cross-modal correspondence of cell-type centroids.",
    "Per library, pseudobulk of each cell type from corrected counts, double-centred (genes and cell types); Spearman correlation between every nuclear and every protoplast cell-type centroid; medians over 10 draws.",
    "Same-label rho minus the highest rho of that cell type with any other label (row and column).",
    "Reciprocal top-1 match in >= 50% of draws and margin not at the label-shuffled null (|z| < 2).",
    "The same statistic between replicate libraries of one modality (Arabidopsis: Ryu protoplast replicates); Margin / yardstick expresses the cross-modal margin relative to it.",
    "All nuclei x protoplast pairs: M. truncatula nMOCK x pMock, nNOD x pNOD, nBAP x pBAP (nBAP flagged); A. thaliana Farmer rep2/4/5 x Ryu rep1-3 (rep4 flagged). The text quotes mean (range) over them.",
    "V2 labels: M. truncatula pairs share 5 cell types (double centring still constrains small panels); A. thaliana pairs 7-9."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(16, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
for (nm in c("Per cell type", "Per pair")) {
  df <- if (nm == "Per cell type") per_type else pair_level
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = hdr)
  freezePane(wb, nm, firstRow = TRUE)
  num <- which(vapply(df, is.double, logical(1)))
  addStyle(wb, nm, createStyle(numFmt = "0.000"), rows = 2:(nrow(df) + 1), cols = num, gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_along(df), widths = 14)
  setColWidths(wb, nm, cols = which(names(df) %in% c("Pair", "Role")), widths = 30)
  setRowHeights(wb, nm, rows = 1, heights = 42)
}

out_file <- file.path(table_dir, "TableS5_centroids.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
