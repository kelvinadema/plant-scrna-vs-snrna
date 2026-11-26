.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

s2 <- "02_nuclear_recovery"
s3 <- "03_recovery_attrition"

recovery <- read_tsv_q(section_table(s2, "fig_library_recovery.tsv")) |>
  select(sample, study, technology, prep, tissue, cells_barcode_rank, cells_automatic,
         retained_cells, main_analysis_selected, median_umis_retained, reads_per_retained_cell)

mapping <- read_tsv_q(section_table(s3, "fig_mapping_readfate_perlib.tsv")) |>
  filter(is_used_library(sample, species, modality))

fail <- read_failure_modes() |> select(sample, failure_mode)

qc_medians <- read_tsv_q(section_table(s3, "fig_qc_percell_summary.tsv")) |>
  filter(stage == "count_gate_candidate",
         metric %in% c("mito_fraction", "chloroplast_fraction", "rrna_fraction")) |>
  select(sample, metric, median) |>
  pivot_wider(names_from = metric, values_from = median)

study_label <- c(in_house = "This study", GSE155304 = "Farmer et al. (GSE155304)",
                 GSE123013 = "Ryu et al. (GSE123013)")
prep_label  <- c(facs_sorted = "FANS-sorted", unsorted = "Unsorted",
                 tip_unfixed_or_dtt = "Unfixed or DTT-fixed",
                 standard = "As published (not recorded)")
treatment_label <- c(mock = "Mock", NF = "Nod factor", BAP = "BAP (cytokinin)",
                     S_meliloti = "S. meliloti", none = "Untreated")

inv <- mapping |>
  left_join(recovery |> select(sample, !any_of(names(mapping))), by = "sample") |>
  left_join(fail, by = "sample") |>
  left_join(qc_medians, by = "sample") |>
  mutate(group = factor(prep_groups[arm], levels = prep_groups)) |>
  arrange(group, plot_order) |>
  transmute(
    Library               = if_else(species == "arabidopsis", sample_label, sample),
    `Sample ID`           = sample,
    Species               = if_else(species == "medicago", "M. truncatula", "A. thaliana"),
    Modality              = lab_modality[modality],
    Preparation           = str_replace_all(prep_groups[arm], "\n", " "),
    `Nuclear prep`        = if_else(modality == "nuclei", prep_label[prep], NA_character_),
    Study                 = study_label[study],
    Chemistry             = technology,
    Tissue                = if_else(species == "arabidopsis", "root tip", tissue),
    Treatment             = coalesce(treatment_label[treatment], "Untreated"),
    `Input read pairs`    = total_reads,
    `Reads mapped (%)`    = 100 * frac_mapping,
    `Sequencing saturation` = sequencing_saturation,
    `Barcode-rank knee (UMIs rank 1k / 10k)` = raw_rank_1000_to_10000_ratio,
    `Cell-calling cutoff (corrected UMIs)` = count_gate_threshold,
    `Candidate barcodes`  = barcodes_count_gate,
    `Cells after automatic QC` = cells_automatic,
    `Retained cells`      = cells_retained,
    `Passes 500-cell floor` = if_else(cells_retained >= 500, "yes", "no"),
    `Used in cross-modality comparisons` = if_else(main_analysis_selected, "yes", "no"),
    `Median UMIs per retained cell` = median_umis_retained,
    `Input reads per retained cell` = reads_per_retained_cell,
    `Reads to retained UMIs (%)` = 100 * frac_read_to_retained_umi,
    `Molecules in empty droplets, outside CellBender window (%)` = 100 * raw_umi_fraction_outside_cb_model,
    `Molecules removed within CellBender window (%)` = 100 * raw_umi_fraction_removed_within_cb_model,
    `Molecules below cell-calling cutoff (%)` = 100 * raw_umi_fraction_corrected_below_count_gate,
    `Molecules in QC-removed barcodes (%)` = 100 * raw_umi_fraction_corrected_later_qc_removed,
    `Molecules retained (%)` = 100 * raw_umi_fraction_retained,
    `Median mitochondrial fraction, candidates (%)` = 100 * mito_fraction,
    `Median chloroplast fraction, candidates (%)`   = 100 * chloroplast_fraction,
    `Median rRNA fraction, candidates (%)`          = 100 * rrna_fraction,
    `Organellar / rRNA cap (%)` = if_else(modality == "nuclei", 2, 5),
    `Failure mode` = if_else(modality == "nuclei",
                             str_replace_all(failure_labels[failure_mode], "\n", " "), NA_character_)
  )

mt_nuc <- inv |> filter(Species == "M. truncatula", Modality == "Nuclei")
stopifnot(nrow(inv) == 29, nrow(mt_nuc) == 14, sum(mt_nuc$`Retained cells` >= 500) == 3,
          all(inv$`Retained cells`[inv$Modality == "Protoplast"] >= 500))

readme <- tibble(
  Field = c("Title", "Libraries", "Candidate barcodes", "Cells after automatic QC", "Retained cells",
            "Passes 500-cell floor", "Barcode-rank knee", "Molecules (%)", "Median fractions",
            "Failure mode", "Excluded libraries"),
  Description = c(
    "Supplementary Table S3. Library inventory, cell recovery and quality control of nuclear and protoplast libraries.",
    "All M. truncatula nuclear libraries (14), the M. truncatula protoplast libraries used in this chapter (7), and the A. thaliana root-tip nuclear (Farmer et al., 5) and protoplast (Ryu et al., 3) libraries, all 10x Genomics.",
    "Barcodes above the cell-calling cutoff (knee of the barcode-rank curve on CellBender-corrected UMIs).",
    "Candidates remaining after the automatic QC cascade (adaptive MAD/miQC filters, absolute organellar/rRNA caps, regression filter).",
    "Cells or nuclei retained after the final (manual) QC step.",
    "Libraries with fewer than 500 retained cells were not used for downstream analyses.",
    "Raw UMIs of the barcode at rank 1,000 divided by rank 10,000; low values indicate a linearly degrading curve without a knee.",
    "Fate of all deduplicated molecules (raw UMIs) in all permit-list barcodes, including empty droplets; the five columns sum to 100%. Barcodes outside the CellBender analysis window keep no corrected counts (empty droplets); the first two columns together are all molecules CellBender assigned to background.",
    "Medians over candidate barcodes, before adaptive QC, as fractions of CellBender-corrected counts.",
    "Rule-based, nuclear libraries only: retained (>= 500 cells); no droplet separation (< 500 candidates); organellar RNA (candidate median mitochondrial fraction above cap); rRNA (candidate median rRNA above cap); otherwise adaptive QC attrition.",
    "BD Rhapsody protoplast libraries, M. truncatula 2 mm root-tip protoplasts (p1-p3) and the S. meliloti protoplast timecourse (B24H-B96H) are not used in this chapter."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(26, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)

addWorksheet(wb, "Libraries")
writeData(wb, "Libraries", inv, headerStyle = hdr)
freezePane(wb, "Libraries", firstRow = TRUE, firstCol = TRUE)
pct_cols <- which(str_detect(names(inv), "\\(%\\)|knee"))
int_cols <- which(names(inv) %in% c("Input read pairs", "Candidate barcodes", "Cells after automatic QC",
                                    "Retained cells", "Median UMIs per retained cell",
                                    "Input reads per retained cell", "Cell-calling cutoff (corrected UMIs)"))
addStyle(wb, "Libraries", createStyle(numFmt = "0.00"), rows = 2:(nrow(inv) + 1),
         cols = pct_cols, gridExpand = TRUE)
addStyle(wb, "Libraries", createStyle(numFmt = "#,##0"), rows = 2:(nrow(inv) + 1),
         cols = int_cols, gridExpand = TRUE)
addStyle(wb, "Libraries", createStyle(numFmt = "0.000"), rows = 2:(nrow(inv) + 1),
         cols = which(names(inv) == "Sequencing saturation"), gridExpand = TRUE)
setColWidths(wb, "Libraries", cols = seq_along(inv), widths = 14)
setColWidths(wb, "Libraries", cols = which(names(inv) %in% c("Sample ID", "Preparation", "Study")),
             widths = 26)
setRowHeights(wb, "Libraries", rows = 1, heights = 45)

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(table_dir, "TableS3_library_inventory.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
