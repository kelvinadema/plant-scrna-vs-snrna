.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

s1 <- "01_stress_response"

per_gene <- read_tsv_q(section_table(s1, "fig_lfc_per_gene.tsv")) |>
  filter(is_deg)

mt_overlap <- read_tsv_q(section_table(s1, "fig_lfc_mt_protoplast_vs_nuclei.tsv")) |>
  transmute(gene,
            overlap = recode(category,
              shared_up   = "shared (induced in both)",
              shared_down = "shared (repressed in both)",
              protoplast_only_up = "protoplast only", protoplast_only_down = "protoplast only",
              nuclei_only_up     = "nuclei only",     nuclei_only_down     = "nuclei only",
              opposite_protoplast_up_nuclei_down = "opposite direction",
              opposite_protoplast_down_nuclei_up = "opposite direction"))

mt_annot <- read_tsv_q(mt_annot_file, quote = "", col_types = cols(.default = col_character())) |>
  transmute(gene = gene_id, gene_name = na_if(gene_name, ""), product = na_if(product, ""))

deg_sheet <- function(sp, mod) {
  out <- per_gene |>
    filter(species == sp, modality == mod) |>
    arrange(desc(direction == "up"), padj_lfc_threshold) |>
    transmute(gene_id = gene, direction = if_else(direction == "up", "induced", "repressed"),
              log2FC_shrunk = log2fc, log2FC = log2fc_unshrunk, lfcSE = lfc_se,
              baseMean = base_mean, pvalue, padj, padj_lfcThreshold1 = padj_lfc_threshold)
  if (sp == "medicago") {
    out <- out |>
      left_join(mt_annot, by = c("gene_id" = "gene")) |>
      left_join(mt_overlap, by = c("gene_id" = "gene")) |>
      relocate(gene_name, product, .after = gene_id) |>
      relocate(overlap, .after = direction)
  }
  out
}

sheets <- list(
  Mt_protoplast = deg_sheet("medicago", "protoplast"),
  Mt_nuclei     = deg_sheet("medicago", "nuclei"),
  At_protoplast = deg_sheet("arabidopsis", "protoplast")
)

stopifnot(nrow(sheets$Mt_protoplast) == 3562, nrow(sheets$Mt_nuclei) == 1184,
          nrow(sheets$At_protoplast) == 2453)

readme <- tibble(
  Field = c("Title", "Contrasts", "Gene-calling gate", "Direction", "log2FC_shrunk", "log2FC",
            "padj", "padj_lfcThreshold1", "overlap (Mt sheets)", "gene_name / product",
            "Medicago data", "Arabidopsis data", "Software"),
  Description = c(
    "Supplementary Table S1. Genes differentially expressed after protoplasting or nuclei extraction, relative to intact root tissue.",
    "Mt_protoplast: M. truncatula protoplasts (P1-P3) vs tissue (T1-T3). Mt_nuclei: M. truncatula isolated nuclei (N6-N8) vs tissue. At_protoplast: A. thaliana protoplasts (RP1-RP2) vs tissue (RU1-RU2).",
    "DESeq2 Wald test against |log2FC| > 1 (lfcThreshold = 1), adjusted P < 0.01. Nuclear-chromosome genes only (organellar contigs removed).",
    "Sign of the unshrunk log2 fold change (preparation / tissue).",
    "apeglm-shrunken log2 fold change.",
    "Unshrunk (MLE) log2 fold change.",
    "Benjamini-Hochberg adjusted P of the standard test against log2FC = 0 (for reference).",
    "Benjamini-Hochberg adjusted P of the lfcThreshold = 1 test; the gate used for this table.",
    "Membership in Figure 2B: shared = called in both Medicago preparations in the same direction; opposite direction = called in both with opposite sign.",
    "Annotations retained from the supplied Table S1; empty where unannotated.",
    "Bulk RNA-seq, kallisto (index Mt5_r1.9).",
    "Denyer et al. 2019, reprocessed with kallisto (TAIR10 / Araport11).",
    "R 4.3.2, DESeq2 1.42.0, apeglm."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(24, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)

num_cols <- c("log2FC_shrunk", "log2FC", "lfcSE", "baseMean")
p_cols   <- c("pvalue", "padj", "padj_lfcThreshold1")
for (nm in names(sheets)) {
  df <- sheets[[nm]]
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = hdr)
  freezePane(wb, nm, firstRow = TRUE)
  addStyle(wb, nm, createStyle(numFmt = "0.00"), rows = 2:(nrow(df) + 1),
           cols = which(names(df) %in% num_cols), gridExpand = TRUE)
  addStyle(wb, nm, createStyle(numFmt = "0.00E+00"), rows = 2:(nrow(df) + 1),
           cols = which(names(df) %in% p_cols), gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_along(df), widths = "auto")
  if ("product" %in% names(df)) setColWidths(wb, nm, cols = which(names(df) == "product"), widths = 50)
}

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(table_dir, "TableS1_DEG_lists.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
