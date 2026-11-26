.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

comp <- read_tsv_q(section_table("05_composition", "fig_composition.tsv")) |>
  mutate(cell_type_f = factor(cell_type_labels[cell_type], levels = cell_type_labels))
stopifnot(!anyNA(comp$cell_type_f))

library_label <- function(sample) {
  sample |>
    str_replace("GSE155304_sNucRNA-seq_rep", "Farmer rep") |>
    str_replace("GSE123013_WT_rep", "Ryu rep")
}

long <- comp |>
  arrange(desc(species), modality, sample, cell_type_f) |>
  transmute(Species = if_else(species == "medicago", "M. truncatula", "A. thaliana"),
            Library = library_label(sample), Modality = lab_modality[modality],
            Treatment = treatment, `Used in cross-modality comparisons` = if_else(primary_comparison, "yes", "no"),
            `Cell type` = as.character(cell_type_f), Cells = n_cells,
            `Retained cells in library` = n_total, `Share (%)` = 100 * proportion)

wide <- comp |>
  arrange(desc(species), modality, sample) |>
  mutate(Library = library_label(sample)) |>
  select(species, cell_type_f, Library, proportion) |>
  mutate(proportion = 100 * proportion) |>
  pivot_wider(names_from = Library, values_from = proportion) |>
  arrange(desc(species), cell_type_f) |>
  mutate(Species = if_else(species == "medicago", "M. truncatula", "A. thaliana"),
         `Cell type` = as.character(cell_type_f), .before = 1) |>
  select(-species, -cell_type_f)

readme <- tibble(
  Field = c("Title", "Libraries", "Denominator", "Labels", "Sheets"),
  Description = c(
    "Supplementary Table S4. Cell-type composition of retained cells per library.",
    "M. truncatula susceptible-zone nuclei (nMOCK, nNOD, nBAP; FANS) and protoplasts (pMock, pNOD, pBAP); A. thaliana root-tip nuclei (Farmer et al. rep1-5) and protoplasts (Ryu et al. rep1-3). Farmer rep1 and rep3 fall below the 500-cell floor and are not used in cross-modality comparisons.",
    "All retained cells of the library, including cells left Unassigned.",
    "Shared per-species label space (integration within species across modality, annotation once); 'pericycle cell' merged into 'pericycle'. Cell types with zero cells are listed explicitly.",
    "Long: one row per library x cell type with counts. Wide: share (%) per cell type (rows) and library (columns)."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(16, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
for (nm in c("Long", "Wide")) {
  df <- if (nm == "Long") long else wide
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = hdr)
  freezePane(wb, nm, firstRow = TRUE)
  pct <- if (nm == "Long") which(names(df) == "Share (%)") else 3:ncol(df)
  addStyle(wb, nm, createStyle(numFmt = "0.00"), rows = 2:(nrow(df) + 1), cols = pct, gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_along(df), widths = "auto")
}

out_file <- file.path(table_dir, "TableS4_composition.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
