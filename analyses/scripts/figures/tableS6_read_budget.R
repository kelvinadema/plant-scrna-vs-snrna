.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

s7 <- "07_read_budget"
sp_lab    <- c(medicago = "M. truncatula", arabidopsis = "A. thaliana")
layer_lab <- c(corrected = "Corrected (CellBender)", raw = "Raw")
lib_label <- function(s) s |>
  str_replace("^GSE123013_WT_", "Ryu ") |>
  str_replace("^GSE155304_sNucRNA-seq_", "Farmer ")
target_lab <- function(t) if_else(t == "absolute_0.70", "AUROC 0.70", "Protoplast full-depth AUROC - 0.01")

saturation <- read_tsv_q(section_table("04_read_yield", "fig_yield_perlib.tsv")) |>
  filter(primary_comparison) |>
  transmute(Species = sp_lab[species], Library = sample_label, Modality = lab_modality[modality],
            `Input read pairs` = total_reads, `Retained cells` = n_cells,
            `Reads per retained cell` = reads_per_retained_cell,
            `Reads per retained UMI` = reads_per_retained_umi,
            `Sequencing saturation` = sequencing_saturation)

# V2 stores crossings per nuclei x protoplast pair (the second target depends on the pair).
crossing <- read_tsv_q(section_table(s7, "supplementary_target_crossings.tsv")) |>
  arrange(desc(species), layer, target, pair, modality) |>
  transmute(Species = sp_lab[species], Pair = str_replace_all(pair, "_", " "),
            Library = lib_label(sample), Modality = lab_modality[modality],
            Counts = layer_lab[layer], Target = target_lab(target), `Target AUROC` = target_auc,
            `Reads per cell at target` = reads_per_cell, Reached = if_else(reached, "yes", "no"),
            Status = status, `Full-depth reads per cell` = reads_full_depth,
            `Full-depth AUROC` = auc_full_depth)

ratios <- read_tsv_q(section_table(s7, "supplementary_exchange_ratio.tsv")) |>
  mutate(quoted = layer == "corrected" & target == "absolute_0.70") |>
  arrange(desc(quoted), desc(species), layer, target, nuclei_library, protoplast_library) |>
  transmute(Species = sp_lab[species], Nuclei = lib_label(nuclei_library),
            Protoplast = lib_label(protoplast_library), Counts = layer_lab[layer],
            Target = target_lab(target), `Target AUROC` = target_auc,
            `Nuclei reads per cell` = reads_n, `Protoplast reads per cell` = reads_p,
            `Nuclei status` = status_n,
            `Ratio (nuclei / protoplast)` = ratio_measured, `Ratio kind` = ratio_kind,
            `Ratio, second seed` = ratio_measured_seed2,
            `Nuclei full-depth AUROC` = auc_full_depth_n, `Protoplast full-depth AUROC` = auc_full_depth_p,
            `Nuclei saturation` = saturation_n, `Protoplast saturation` = saturation_p,
            `Quoted in text` = if_else(quoted, "yes", "no"))

stopifnot(nrow(saturation) == 12, nrow(crossing) == 88, nrow(ratios) == 44,
          !anyNA(crossing$Counts), !anyNA(ratios$Counts))

readme <- tibble(
  Field = c("Title", "Design", "AUROC", "Crossing point", "Ratio kind", "Exclusions", "Confound"),
  Description = c(
    "Supplementary Table S6. Sequencing read budget for cell-identity discrimination.",
    "Each library's reads were subsampled to 5, 10, 25, 50 and 75% (plus full depth), remapped and put through the same QC; the held-out marker test (V2: 25 training + 25 test cells per class) was re-run on the same cells at every fraction. nNOD cannot supply 25 + 25 cells per class, so the Nod factor pair is not assessable.",
    "Held-out single-gene marker AUROC: median over 50 markers, averaged over cell types and splits.",
    "Reads per (fixed) cell at which the AUROC curve first reaches the target, interpolated linearly in log10(reads); no extrapolation beyond the tested range.",
    "interpolated = both libraries cross within the tested range; not_reached = the nuclear library never reaches the target at any tested depth (no ratio; do not extrapolate). Farmer rep4 only just reaches AUROC 0.70, so its ratio is unstable.",
    "Corrected counts: pMock 5% and 10% (CellBender analysed < 60% of the fixed cells) and nBAP 75% (CellBender anomaly) are excluded from interpolation. Raw counts use every subsample.",
    "A. thaliana nuclei (Farmer et al., 10x v3/v2) and protoplasts (Ryu et al., 10x v2) come from different studies and chemistries; ratios are library comparisons, not a pure modality constant."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(16, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
for (nm in c("Saturation", "Crossing points", "Exchange ratios")) {
  df <- switch(nm, "Saturation" = saturation, "Crossing points" = crossing, "Exchange ratios" = ratios)
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = hdr)
  freezePane(wb, nm, firstRow = TRUE)
  num <- which(vapply(df, is.double, logical(1)))
  addStyle(wb, nm, createStyle(numFmt = "0.000"), rows = 2:(nrow(df) + 1), cols = num, gridExpand = TRUE)
  big <- which(names(df) %in% c("Input read pairs", "Reads per retained cell", "Reads per cell at target",
                                "Full-depth reads per cell", "Nuclei reads per cell", "Protoplast reads per cell"))
  addStyle(wb, nm, createStyle(numFmt = "#,##0"), rows = 2:(nrow(df) + 1), cols = big, gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_along(df), widths = 15)
  setRowHeights(wb, nm, rows = 1, heights = 42)
}

out_file <- file.path(table_dir, "TableS6_read_budget.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
