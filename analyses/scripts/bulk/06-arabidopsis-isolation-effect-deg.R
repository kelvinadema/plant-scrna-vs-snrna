suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

DROP_ORGANELLAR <- TRUE

.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
counts_file <- file.path(chapter_dir,
  "input/bulk/arabidopsis/gene_counts.tsv")
out_dir     <- file.path(chapter_dir,
                         if (DROP_ORGANELLAR) "output/bulk/arabidopsis-isolation-effects-deg-nuclear-only"
                         else                 "output/bulk/arabidopsis-isolation-effects-deg")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ALPHA      <- 0.01   # adjusted p-value cutoff
LFC_CUTOFF <- 1      # |log2FC| cutoff (2-fold)

sample_table <- data.frame(
  sample = c("RU1", "RU2", "RP1", "RP2"),
  method = factor(
    c("tissue", "tissue", "protoplast", "protoplast"),
    levels = c("tissue", "protoplast")   # first level = reference
  ),
  stringsAsFactors = FALSE
)
rownames(sample_table) <- sample_table$sample

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
stopifnot(all(sample_table$sample %in% colnames(counts)))
counts <- round(as.matrix(counts[, sample_table$sample, drop = FALSE]))
mode(counts) <- "integer"

message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))

markers <- c(LOX3 = "AT1G17420", BHLH92 = "AT5G43650",
             HSP101 = "AT1G74310", WRKY33 = "AT2G38470")
cpm <- t(t(counts) / colSums(counts)) * 1e6
present <- markers[markers %in% rownames(counts)]
stopifnot(length(present) >= 2)
marker_tab <- data.frame(
  marker         = names(present),
  gene           = unname(present),
  cpm_tissue     = round(rowMeans(cpm[present, c("RU1", "RU2"), drop = FALSE]), 1),
  cpm_protoplast = round(rowMeans(cpm[present, c("RP1", "RP2"), drop = FALSE]), 1),
  row.names      = NULL
)
marker_tab$log2FC <- round(log2((marker_tab$cpm_protoplast + 1) /
                                (marker_tab$cpm_tissue + 1)), 2)
cat("\n=== orientation check: protoplasting-induced markers (Denyer 2019) ===\n")
print(marker_tab, row.names = FALSE)
if (median(marker_tab$log2FC) <= 0) {
  stop("Marker genes are NOT up in the arm labelled 'protoplast'. The ",
       "sample-to-arm assignment has changed -- see caveat 2 in this ",
       "script's header and the note in design.tsv before going further.")
}
cat("Markers are up in the protoplast arm; labelling matches the paper.\n")
write.table(marker_tab, file.path(out_dir, "orientation_check_markers.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (DROP_ORGANELLAR) {
  organellar <- grepl("^AT[MC]G", rownames(counts))
  message(sprintf("Dropping %d organellar genes (%d MT, %d CP) before normalisation",
                  sum(organellar),
                  sum(grepl("^ATMG", rownames(counts))),
                  sum(grepl("^ATCG", rownames(counts)))))
  counts <- counts[!organellar, , drop = FALSE]
}

dds <- DESeqDataSetFromMatrix(counts, colData = sample_table, design = ~ method)

keep <- rowSums(counts(dds) >= 10) >= 2
message(sprintf("Kept %d / %d genes after pre-filtering", sum(keep), nrow(dds)))
dds <- dds[keep, ]

dds <- DESeq(dds)

sf <- data.frame(sample = colnames(dds), method = dds$method,
                 sizeFactor = round(sizeFactors(dds), 4))
write.table(sf, file.path(out_dir, "size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = "method", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = method, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = sprintf("PC1: %d%% variance", pct[1]),
       y = sprintf("PC2: %d%% variance", pct[2]),
       title = "Arabidopsis root: protoplasting separates on PC1 (VST, blind)") +
  coord_fixed() +
  theme_bw()
ggsave(file.path(out_dir, "pca_method.png"), p, width = 6, height = 5, dpi = 300)

sample_dist <- as.matrix(dist(t(assay(vsd))))
write.table(sample_dist, file.path(out_dir, "sample_distances.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

contrasts <- list(
  protoplast_vs_tissue = c("method", "protoplast", "tissue")
)

summarise_contrast <- function(name, contrast) {
  res <- results(dds, contrast = contrast, alpha = ALPHA)

  res_thr <- results(dds, contrast = contrast, alpha = ALPHA,
                     lfcThreshold = LFC_CUTOFF, altHypothesis = "greaterAbs")

  coef_name <- paste0("method_", contrast[2], "_vs_", contrast[3])
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)

  out <- as.data.frame(res_shrunk)
  out$gene <- rownames(out)
  out$log2FoldChange_unshrunk <- res$log2FoldChange[match(out$gene, rownames(res))]
  out$is_DEG <- !is.na(out$padj) & out$padj < ALPHA &
    abs(out$log2FoldChange_unshrunk) > LFC_CUTOFF

  out$padj_lfcThreshold <- res_thr$padj[match(out$gene, rownames(res_thr))]
  out$is_DEG_lfcThreshold <- !is.na(out$padj_lfcThreshold) &
    out$padj_lfcThreshold < ALPHA

  out <- out[order(out$padj, -abs(out$log2FoldChange)), ]
  out <- out[, c("gene", "baseMean", "log2FoldChange", "log2FoldChange_unshrunk",
                 "lfcSE", "pvalue", "padj", "is_DEG",
                 "padj_lfcThreshold", "is_DEG_lfcThreshold")]
  write.table(out, file.path(out_dir, sprintf("deseq2_%s.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  png(file.path(out_dir, sprintf("ma_%s.png", name)), width = 1800, height = 1500, res = 300)
  plotMA(res_shrunk, ylim = c(-8, 8), main = paste("Arabidopsis:", gsub("_", " ", name)))
  dev.off()

  padj_ok <- !is.na(res$padj) & res$padj < ALPHA
  deg     <- padj_ok & abs(res$log2FoldChange) > LFC_CUTOFF
  deg_thr <- !is.na(res_thr$padj) & res_thr$padj < ALPHA

  data.frame(
    contrast         = name,
    genes_tested     = sum(!is.na(res$padj)),
    DEGs             = sum(deg),
    DEGs_up          = sum(deg & res$log2FoldChange > 0),
    DEGs_down        = sum(deg & res$log2FoldChange < 0),
    pct_DEG          = round(100 * sum(deg) / sum(!is.na(res$padj)), 2),
    median_abs_lfc_of_DEGs = round(median(abs(res$log2FoldChange[deg])), 3),
    padj_only        = sum(padj_ok),
    DEGs_lfcThreshold = sum(deg_thr)
  )
}

summary_tab <- do.call(rbind, Map(summarise_contrast, names(contrasts), contrasts))
rownames(summary_tab) <- NULL

write.table(summary_tab, file.path(out_dir, "deg_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n==========================================================\n")
cat(sprintf("Arabidopsis root: DEGs vs intact tissue  (padj < %g AND |log2FC| > %g)\n",
            ALPHA, LFC_CUTOFF))
cat(sprintf("Gene set: %s\n",
            if (DROP_ORGANELLAR) "NUCLEAR ONLY (MT + CP genes dropped pre-normalisation)"
            else                 "all genes (nuclear + organellar)"))
cat("==========================================================\n")
print(summary_tab, row.names = FALSE)
cat("\nn = 2 per group: treat the DEG count as a floor, and compare gene\n")
cat("classes and direction with Medicago rather than raw counts.\n\n")
cat("Outputs written to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
