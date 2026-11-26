.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
library(openxlsx)

go_file <- function(...) file.path(bulk_dir, ...)
ont_lab <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")

tidy_go <- function(df) {
  df |>
    arrange(list_name, padj) |>
    transmute(`GO ID` = GO, Ontology = ont_lab[ontology], Term = term,
              `Genes in list` = n_in_list, `Genes in universe` = n_in_universe,
              Expected = expected, `Fold enrichment` = fold_enrich,
              `P value` = pvalue, `Adjusted P` = padj,
              `Significant (adj. P < 0.05)` = if_else(padj < 0.05, "yes", "no"))
}

list_lab <- c(shared_up = "Shared, induced", shared_down = "Shared, repressed",
              protoplast_up_unique = "Protoplast only, induced", protoplast_down_unique = "Protoplast only, repressed",
              nuclei_up_unique = "Nuclei only, induced", nuclei_down_unique = "Nuclei only, repressed")
shared <- read_tsv_q(go_file("isolation-effects-go-strict/go_enrichment_shared_unique.tsv")) |>
  mutate(list_name = factor(list_name, levels = names(list_lab))) |>
  group_by(list_name) |>
  group_modify(~ tidy_go(.x |> mutate(list_name = .y$list_name))) |>
  ungroup() |>
  transmute(`Gene list` = list_lab[as.character(list_name)], across(-list_name))

per_prep_lab <- c(protoplast_vs_tissue_up = "induced", protoplast_vs_tissue_down = "repressed",
                  nuclei_vs_tissue_up = "induced", nuclei_vs_tissue_down = "repressed")
read_prep <- function(path, species) {
  read_tsv_q(path) |>
    mutate(prep = if_else(str_starts(list_name, "protoplast"), "Protoplast vs tissue", "Nuclei vs tissue"),
           direction = per_prep_lab[list_name]) |>
    group_by(prep, direction) |>
    group_modify(~ tidy_go(.x |> mutate(list_name = "x"))) |>
    ungroup() |>
    transmute(Species = species, Contrast = prep, Direction = direction, across(-c(prep, direction)))
}
per_prep <- bind_rows(
  read_prep(go_file("isolation-effects-go-strict/go_enrichment_all_lists.tsv"), "M. truncatula"),
  read_prep(go_file("arabidopsis-isolation-effects-go-strict/go_enrichment_all_lists.tsv"), "A. thaliana"))

agree_lab <- c(both = "Both species", Medicago_only = "M. truncatula only",
               Arabidopsis_only = "A. thaliana only", neither = "Neither")
cross <- read_tsv_q(go_file("cross-species-protoplast-strict/go_concordance_terms.tsv")) |>
  arrange(desc(direction), factor(agreement, levels = names(agree_lab)), padj_Mt) |>
  transmute(Direction = if_else(direction == "up", "induced", "repressed"),
            `GO ID` = GO, Ontology = ont_lab[ontology], Term = term,
            `Mt genes in list` = n_in_list_Mt, `Mt genes in universe` = n_in_universe_Mt,
            `Mt fold enrichment` = fold_enrich_Mt, `Mt adjusted P` = padj_Mt,
            `At genes in list` = n_in_list_At, `At genes in universe` = n_in_universe_At,
            `At fold enrichment` = fold_enrich_At, `At adjusted P` = padj_At,
            `Enriched in (adj. P < 0.05)` = agree_lab[agreement],
            `In Figure 2C` = if_else(direction == "up" & ontology == "BP" & agreement == "both", "yes", "no"))

concord <- read_tsv_q(go_file("cross-species-protoplast-strict/go_concordance_summary.tsv")) |>
  transmute(Direction = if_else(direction == "up", "induced", "repressed"),
            Ontology = if_else(ontology == "all", "All", ont_lab[ontology]),
            `Terms tested in both` = co_tested, `Enriched in both` = enriched_both,
            `M. truncatula only` = medicago_only, `A. thaliana only` = arabidopsis_only,
            Neither = neither, `Odds ratio` = odds_ratio, `Fisher P` = fisher_p,
            `Spearman rho (fold enrichment)` = spearman_rho_fold)

discord <- read_tsv_q(go_file("cross-species-protoplast-strict/go_direction_discordant.tsv")) |>
  transmute(`GO ID` = GO, Ontology = ont_lab[ontology], Term = term,
            `Mt induced` = sig_Mt_up, `At induced` = sig_At_up,
            `Mt repressed` = sig_Mt_down, `At repressed` = sig_At_down)

stopifnot(sum(cross$`In Figure 2C` == "yes") == 27,
          concord$`Enriched in both`[concord$Direction == "induced" & concord$Ontology == "Biological process"] == 27,
          !anyNA(shared$`Gene list`), !anyNA(per_prep$Direction))

readme <- tibble(
  Field = c("Title", "Gene lists", "Test", "Shared and unique", "Per preparation", "Cross-species",
            "Concordance", "Direction-discordant"),
  Description = c(
    "Supplementary Table S2. GO enrichment of genes responsive to protoplasting and nuclei extraction.",
    "Genes differentially expressed vs intact root tissue (DESeq2 test against |log2FC| > 1, adjusted P < 0.01; Table S1), split into induced and repressed.",
    "One-sided hypergeometric test per GO term against all tested genes (universe), terms with 10-2,000 annotated genes; Benjamini-Hochberg adjusted P. M. truncatula annotation MtrunA17r5.0; A. thaliana direct gene-GO assignments.",
    "M. truncatula: shared = changed by both preparations in the same direction; only = changed by one preparation. All tested terms.",
    "M. truncatula protoplast and nuclei, and A. thaliana protoplast, each against tissue. All tested terms.",
    "GO terms tested in both species for the protoplast lists; enrichment call per species at adjusted P < 0.05. 'In Figure 2C' marks the 27 biological-process terms induced in both.",
    "Overlap of enriched terms between species (Fisher's exact test on the 2 x 2 counts) and Spearman correlation of fold enrichment over co-tested terms, per direction and ontology.",
    "Terms enriched among induced genes in one species and among repressed genes in the other."
  ))

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:2, widths = c(20, 120))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
sheets <- list(`Shared and unique` = shared, `Per preparation` = per_prep, `Cross-species` = cross,
               Concordance = concord, `Direction-discordant` = discord)
for (nm in names(sheets)) {
  df <- sheets[[nm]]
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = hdr)
  freezePane(wb, nm, firstRow = TRUE)
  pc <- which(str_detect(names(df), "P$|P value|Adjusted P|Fisher P"))
  fc <- which(str_detect(names(df), "fold enrichment|Fold enrichment|Expected|Odds ratio|rho"))
  addStyle(wb, nm, createStyle(numFmt = "0.00E+00"), rows = 2:(nrow(df) + 1), cols = pc, gridExpand = TRUE)
  addStyle(wb, nm, createStyle(numFmt = "0.00"), rows = 2:(nrow(df) + 1), cols = fc, gridExpand = TRUE)
  setColWidths(wb, nm, cols = seq_along(df), widths = 14)
  setColWidths(wb, nm, cols = which(names(df) == "Term"), widths = 50)
  setRowHeights(wb, nm, rows = 1, heights = 30)
}

out_file <- file.path(table_dir, "TableS2_GO_enrichment.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Wrote ", out_file)
