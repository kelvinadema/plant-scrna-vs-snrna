.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
if (!file.exists(file.path(chapter_dir, "input/bulk/medicago_gene_annotation.tsv"))) {
  stop("Supply the original Medicago gAnnotation.tsv as input/bulk/medicago_gene_annotation.tsv before rebuilding GO enrichment.")
}
scripts <- list.files(dirname(normalizePath(.script_file)), pattern = "^[0-9].*\\.R$", full.names = TRUE)
for (script in scripts) {
  message(basename(script))
  source(script, local = new.env(parent = globalenv()))
}
