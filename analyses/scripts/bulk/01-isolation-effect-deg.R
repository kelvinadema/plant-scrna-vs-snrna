suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

DROP_ORGANELLAR <- TRUE

.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
counts_file <- file.path(chapter_dir, "input/bulk/medicago/gene_counts.tsv")
out_dir     <- file.path(chapter_dir,
                         if (DROP_ORGANELLAR) "output/bulk/isolation-effects-deg-nuclear-only"
                         else                 "output/bulk/isolation-effects-deg")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ALPHA      <- 0.01   # adjusted p-value cutoff
LFC_CUTOFF <- 1      # |log2FC| cutoff (2-fold)

sample_table <- data.frame(
  sample = c("T1", "T2", "T3", "P1", "P2", "P3", "N6", "N7", "N8"),
  method = factor(
    c("tissue", "tissue", "tissue",
      "protoplast", "protoplast", "protoplast",
      "nuclei", "nuclei", "nuclei"),
    levels = c("tissue", "protoplast", "nuclei")   # first level = reference
  ),
  row.names = NULL,
  stringsAsFactors = FALSE
)
rownames(sample_table) <- sample_table$sample

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
stopifnot(all(sample_table$sample %in% colnames(counts)))
counts <- round(as.matrix(counts[, sample_table$sample, drop = FALSE]))
mode(counts) <- "integer"

message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))

if (DROP_ORGANELLAR) {
  organellar <- grepl("^MtrunA17_(MTg|CPg)", rownames(counts))
  message(sprintf("Dropping %d organellar genes (%d MT, %d CP) before normalisation",
                  sum(organellar),
                  sum(grepl("^MtrunA17_MTg", rownames(counts))),
                  sum(grepl("^MtrunA17_CPg", rownames(counts)))))
  counts <- counts[!organellar, , drop = FALSE]
}

dds <- DESeqDataSetFromMatrix(counts, colData = sample_table, design = ~ method)

keep <- rowSums(counts(dds) >= 10) >= 3
message(sprintf("Kept %d / %d genes after pre-filtering", sum(keep), nrow(dds)))
dds <- dds[keep, ]

dds <- DESeq(dds)

sf <- data.frame(sample = colnames(dds), method = dds$method,
                 sizeFactor = round(sizeFactors(dds), 4))
write.table(sf, file.path(out_dir, "size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = "method", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = method, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = sprintf("PC1: %d%% variance", pct[1]),
       y = sprintf("PC2: %d%% variance", pct[2]),
       title = "Isolation method separates on PC1/PC2 (VST, blind)") +
  coord_fixed() +
  theme_bw()
ggsave(file.path(out_dir, "pca_method.png"), p, width = 6, height = 5, dpi = 300)

sample_dist <- as.matrix(dist(t(assay(vsd))))
write.table(sample_dist, file.path(out_dir, "sample_distances.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

contrasts <- list(
  protoplast_vs_tissue = c("method", "protoplast", "tissue"),
  nuclei_vs_tissue     = c("method", "nuclei",     "tissue")
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
  plotMA(res_shrunk, ylim = c(-8, 8), main = gsub("_", " ", name))
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
cat(sprintf("DEGs vs intact tissue  (padj < %g AND |log2FC| > %g)\n",
            ALPHA, LFC_CUTOFF))
cat(sprintf("Gene set: %s\n",
            if (DROP_ORGANELLAR) "NUCLEAR ONLY (MT + CP genes dropped pre-normalisation)"
            else                 "all genes (nuclear + organellar)"))
cat("==========================================================\n")
print(summary_tab, row.names = FALSE)
cat("\nLarger DEG count = larger effect of the isolation method on the transcriptome.\n")
cat("Reminder: the nuclei arm is confounded with sequencing batch/flowcell,\n")
cat("so its count is an upper bound. See the header of this script.\n\n")
cat("Outputs written to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
