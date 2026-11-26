suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
})

DEG_FLAG <- "is_DEG_lfcThreshold"
.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
annot_file  <- file.path(chapter_dir, "input/bulk/medicago_gene_annotation.tsv")
deg_dir     <- file.path(chapter_dir, "output/bulk/isolation-effects-deg-nuclear-only")
go_suffix   <- if (DEG_FLAG == "is_DEG_lfcThreshold") "-strict" else ""
out_dir     <- file.path(chapter_dir,
                         paste0("output/bulk/isolation-effects-go", go_suffix))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MIN_SET <- 10
MAX_SET <- 2000
TOP_N   <- 20

if (!file.exists(annot_file)) stop("Missing original Medicago GO annotation: input/bulk/medicago_gene_annotation.tsv (LOCUS_TAG, GO, PRODUCT).")
annot <- read.delim(annot_file, quote = "", comment.char = "",
                    colClasses = "character")
annot <- annot[annot$GO != "", c("LOCUS_TAG", "GO", "PRODUCT")]

g2go <- strsplit(annot$GO, ";", fixed = TRUE)
names(g2go) <- annot$LOCUS_TAG

contrasts <- c("protoplast_vs_tissue", "nuclei_vs_tissue")
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

split_dir <- function(df) {
  sig <- df[df[[DEG_FLAG]] %in% TRUE, ]
  list(up   = sig$gene[sig$log2FoldChange_unshrunk > 0],
       down = sig$gene[sig$log2FoldChange_unshrunk < 0])
}
P <- split_dir(res_list[["protoplast_vs_tissue"]])
N <- split_dir(res_list[["nuclei_vs_tissue"]])

lists <- list()
for (dir in c("up", "down")) {
  lists[[sprintf("shared_%s", dir)]]            <- intersect(N[[dir]], P[[dir]])
  lists[[sprintf("nuclei_%s_unique", dir)]]     <- setdiff(N[[dir]], P[[dir]])
  lists[[sprintf("protoplast_%s_unique", dir)]] <- setdiff(P[[dir]], N[[dir]])
}

file_for <- c(shared_up              = "nuclei_up_shared.txt",
              nuclei_up_unique       = "nuclei_up_unique.txt",
              protoplast_up_unique   = "protoplast_up_unique.txt",
              shared_down            = "nuclei_down_shared.txt",
              nuclei_down_unique     = "nuclei_down_unique.txt",
              protoplast_down_unique = "protoplast_down_unique.txt")
for (nm in names(lists)) {
  f <- file.path(deg_dir, file_for[[nm]])
  if (file.exists(f)) {
    old <- readLines(f)
    if (!setequal(old, lists[[nm]]))
      warning(sprintf("%s on disk (%d genes) differs from the recomputed %s (%d) -- NOT overwritten",
                      basename(f), length(old), nm, length(lists[[nm]])))
  } else {
    writeLines(sort(lists[[nm]]), f)
    message(sprintf("wrote %s (%d genes)", basename(f), length(lists[[nm]])))
  }
}

enrich <- function(sig_genes, label) {
  sig_genes <- intersect(sig_genes, universe)
  N_u <- length(universe); K <- length(sig_genes)

  hits <- vapply(term2gene, function(g) sum(g %in% sig_genes), integer(1))
  size <- lengths(term2gene)

  p <- phyper(hits - 1, size, N_u - size, K, lower.tail = FALSE)

  tab <- data.frame(
    list_name     = label,
    GO            = names(term2gene),
    ontology      = unname(go_ont[names(term2gene)]),
    term          = unname(go_term[names(term2gene)]),
    n_in_list     = hits,
    n_in_universe = size,
    expected      = round(size * K / N_u, 1),
    fold_enrich   = round((hits / K) / (size / N_u), 2),
    pvalue        = p,
    stringsAsFactors = FALSE
  )
  tab <- tab[!is.na(tab$ontology), ]
  tab$padj <- ave(tab$pvalue, tab$ontology, FUN = function(x) p.adjust(x, "BH"))
  tab <- tab[order(tab$padj, -tab$fold_enrich), ]
  tab[tab$n_in_list >= 5, ]
}

all_tabs <- list()
for (nm in names(lists)) {
  genes <- lists[[nm]]
  tab <- enrich(genes, nm)
  all_tabs[[nm]] <- tab

  cat(sprintf("\n=== %s  (%d genes, %d annotated) ===\n",
              nm, length(genes), length(intersect(genes, universe))))
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

out <- do.call(rbind, all_tabs)
write.table(out, file.path(out_dir, "go_enrichment_shared_unique.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nFull tables written to:",
    file.path(out_dir, "go_enrichment_shared_unique.tsv"), "\n")
