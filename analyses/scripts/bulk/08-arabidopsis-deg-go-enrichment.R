suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(org.At.tair.db)
})

DEG_FLAG <- "is_DEG_lfcThreshold"
go_suffix <- if (DEG_FLAG == "is_DEG_lfcThreshold") "-strict" else ""

.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
deg_dir     <- file.path(chapter_dir, "output/bulk/arabidopsis-isolation-effects-deg-nuclear-only")
out_dir     <- file.path(chapter_dir,
                         paste0("output/bulk/arabidopsis-isolation-effects-go", go_suffix))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ALPHA   <- 0.01   # must match script 06
MIN_SET <- 10     # smallest GO term to test (genes in universe)
MAX_SET <- 2000   # largest  GO term to test
TOP_N   <- 20     # terms printed per list

g2go <- as.list(org.At.tairGO)
g2go <- g2go[!is.na(g2go)]
g2go <- lapply(g2go, names)   # element names are the GO ids
message(sprintf("Annotation: %d genes with >= 1 GO term", length(g2go)))

contrasts <- c("protoplast_vs_tissue")

res_list <- lapply(contrasts, function(cn)
  read.delim(file.path(deg_dir, sprintf("deseq2_%s.tsv", cn))))
names(res_list) <- contrasts

universe <- res_list[[1]]$gene[!is.na(res_list[[1]]$padj)]
universe <- intersect(universe, names(g2go))
message(sprintf("Universe: %d tested + annotated genes", length(universe)))

g2go <- g2go[universe]

anc <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
expand <- function(terms) {
  terms <- terms[terms %in% names(anc)]
  if (!length(terms)) return(character(0))
  out <- unique(c(terms, unlist(anc[terms], use.names = FALSE)))
  out[!is.na(out) & out != "all"]
}
g2go_full <- lapply(g2go, expand)

pairs <- data.frame(
  gene = rep(names(g2go_full), lengths(g2go_full)),
  go   = unlist(g2go_full, use.names = FALSE),
  stringsAsFactors = FALSE
)
term2gene <- split(pairs$gene, pairs$go)
term2gene <- term2gene[lengths(term2gene) >= MIN_SET &
                       lengths(term2gene) <= MAX_SET]
message(sprintf("Testing %d GO terms (size %d-%d)", length(term2gene), MIN_SET, MAX_SET))

go_ont  <- AnnotationDbi::Ontology(names(term2gene))
go_term <- AnnotationDbi::Term(names(term2gene))

enrich <- function(sig_genes, label) {
  sig_genes <- intersect(sig_genes, universe)
  N <- length(universe); K <- length(sig_genes)

  hits <- vapply(term2gene, function(g) sum(g %in% sig_genes), integer(1))
  size <- lengths(term2gene)

  p <- phyper(hits - 1, size, N - size, K, lower.tail = FALSE)

  tab <- data.frame(
    list_name     = label,
    GO            = names(term2gene),
    ontology      = unname(go_ont[names(term2gene)]),
    term          = unname(go_term[names(term2gene)]),
    n_in_list     = hits,
    n_in_universe = size,
    expected      = round(size * K / N, 1),
    fold_enrich   = round((hits / K) / (size / N), 2),
    pvalue        = p,
    stringsAsFactors = FALSE
  )
  tab <- tab[!is.na(tab$ontology), ]
  tab$padj <- ave(tab$pvalue, tab$ontology, FUN = function(x) p.adjust(x, "BH"))
  tab <- tab[order(tab$padj, -tab$fold_enrich), ]
  tab[tab$n_in_list >= 5, ]
}

all_tabs <- list()
for (cn in contrasts) {
  r <- res_list[[cn]]
  r <- r[r[[DEG_FLAG]] %in% TRUE, ]
  for (dir in c("up", "down")) {
    genes <- if (dir == "up") r$gene[r$log2FoldChange_unshrunk > 0]
             else             r$gene[r$log2FoldChange_unshrunk < 0]
    label <- sprintf("%s_%s", cn, dir)
    tab <- enrich(genes, label)
    all_tabs[[label]] <- tab

    cat(sprintf("\n=== %s  (%d DEGs, %d annotated) ===\n",
                label, length(genes), length(intersect(genes, universe))))
    top <- head(tab[tab$padj < 0.05, ], TOP_N)
    if (!nrow(top)) {
      cat("  no GO term enriched at BH < 0.05\n")
    } else {
      print(data.frame(
        ont  = top$ontology,
        term = substr(top$term, 1, 52),
        n    = sprintf("%d/%d", top$n_in_list, top$n_in_universe),
        fold = top$fold_enrich,
        padj = signif(top$padj, 2)
      ), row.names = FALSE)
    }
  }
}

out <- do.call(rbind, all_tabs)
write.table(out, file.path(out_dir, "go_enrichment_all_lists.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nFull tables written to:", out_dir, "\n")
