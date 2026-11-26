DEG_FLAG  <- "is_DEG_lfcThreshold"
go_suffix <- if (DEG_FLAG == "is_DEG_lfcThreshold") "-strict" else ""

.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_paths.R"))
mt_file <- file.path(chapter_dir,
  paste0("output/bulk/isolation-effects-go", go_suffix, "/go_enrichment_all_lists.tsv"))
at_file <- file.path(chapter_dir,
  paste0("output/bulk/arabidopsis-isolation-effects-go", go_suffix,
         "/go_enrichment_all_lists.tsv"))
out_dir <- file.path(chapter_dir,
  paste0("output/bulk/cross-species-protoplast", go_suffix))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SIG    <- 0.05   # BH cutoff for calling a GO term enriched
TOP_N  <- 25

mt <- read.delim(mt_file, quote = "", comment.char = "")
at <- read.delim(at_file, quote = "", comment.char = "")

mt_deg <- read.delim(file.path(chapter_dir,
  "output/bulk/isolation-effects-deg-nuclear-only/deg_summary.tsv"))
at_deg <- read.delim(file.path(chapter_dir,
  "output/bulk/arabidopsis-isolation-effects-deg-nuclear-only/deg_summary.tsv"))

headline <- rbind(
  data.frame(species = "Medicago",    n_per_group = 3,
             mt_deg[mt_deg$contrast == "protoplast_vs_tissue", ]),
  data.frame(species = "Arabidopsis", n_per_group = 2,
             at_deg[at_deg$contrast == "protoplast_vs_tissue", ])
)
cat("\n=== protoplast vs tissue, nuclear-only, both species ===\n")
print(headline, row.names = FALSE)
cat("\nn differs (3 vs 2 per group), so pct_DEG is the comparable column,\n")
cat("not DEGs. Arabidopsis is the lower-powered of the two.\n")
write.table(headline, file.path(out_dir, "deg_summary_both_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

compare_direction <- function(dir, ont = NULL) {
  lab <- sprintf("protoplast_vs_tissue_%s", dir)
  m <- mt[mt$list_name == lab, ]
  a <- at[at$list_name == lab, ]
  if (!is.null(ont)) m <- m[m$ontology %in% ont, ]

  j <- merge(m[, c("GO", "ontology", "term", "n_in_list", "n_in_universe",
                   "fold_enrich", "padj")],
             a[, c("GO", "n_in_list", "n_in_universe", "fold_enrich", "padj")],
             by = "GO", suffixes = c("_Mt", "_At"))

  j$sig_Mt <- j$padj_Mt < SIG
  j$sig_At <- j$padj_At < SIG
  j$agreement <- ifelse(j$sig_Mt & j$sig_At, "both",
                 ifelse(j$sig_Mt, "Medicago_only",
                 ifelse(j$sig_At, "Arabidopsis_only", "neither")))
  j$direction <- dir

  tab <- table(factor(j$agreement,
                      levels = c("both", "Medicago_only", "Arabidopsis_only", "neither")))

  ft <- fisher.test(table(j$sig_Mt, j$sig_At))

  rho <- suppressWarnings(
    cor(j$fold_enrich_Mt, j$fold_enrich_At, method = "spearman",
        use = "complete.obs"))

  cat(sprintf("\n=== protoplast_vs_tissue_%s (%s) : GO terms co-tested in both species ===\n",
              dir, if (is.null(ont)) "all ontologies" else paste(ont, collapse = "+")))
  cat(sprintf("co-tested terms      : %d\n", nrow(j)))
  cat(sprintf("enriched in both     : %d\n", tab[["both"]]))
  cat(sprintf("Medicago only        : %d\n", tab[["Medicago_only"]]))
  cat(sprintf("Arabidopsis only     : %d\n", tab[["Arabidopsis_only"]]))
  cat(sprintf("neither              : %d\n", tab[["neither"]]))
  cat(sprintf("co-enrichment odds ratio : %.1f (Fisher p = %s)\n",
              unname(ft$estimate), format.pval(ft$p.value, digits = 3)))
  cat(sprintf("Spearman rho of fold-enrichment, all co-tested terms : %.2f\n", rho))

  both <- j[j$agreement == "both", ]
  both <- both[order(-(both$fold_enrich_Mt + both$fold_enrich_At)), ]
  cat(sprintf("\n  top %d terms enriched in BOTH species (by summed fold-enrichment):\n", TOP_N))
  if (!nrow(both)) {
    cat("  none\n")
  } else {
    print(data.frame(
      ont     = both$ontology,
      term    = substr(both$term, 1, 46),
      fold_Mt = both$fold_enrich_Mt,
      fold_At = both$fold_enrich_At,
      padj_Mt = signif(both$padj_Mt, 2),
      padj_At = signif(both$padj_At, 2)
    )[seq_len(min(TOP_N, nrow(both))), ], row.names = FALSE)
  }

  list(join = j,
       stats = data.frame(
         direction        = dir,
         ontology         = if (is.null(ont)) "all" else paste(ont, collapse = "+"),
         co_tested        = nrow(j),
         enriched_both    = as.integer(tab[["both"]]),
         medicago_only    = as.integer(tab[["Medicago_only"]]),
         arabidopsis_only = as.integer(tab[["Arabidopsis_only"]]),
         neither          = as.integer(tab[["neither"]]),
         odds_ratio       = round(unname(ft$estimate), 2),
         fisher_p         = signif(ft$p.value, 3),
         spearman_rho_fold = round(rho, 3)
       ))
}

res <- lapply(c("up", "down"), compare_direction)
names(res) <- c("up", "down")

res_bp <- lapply(c("up", "down"), compare_direction, ont = "BP")
names(res_bp) <- c("up", "down")

stats <- do.call(rbind, c(lapply(res, `[[`, "stats"),
                          lapply(res_bp, `[[`, "stats")))
rownames(stats) <- NULL
write.table(stats, file.path(out_dir, "go_concordance_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

joined <- do.call(rbind, lapply(res, function(x) {
  d <- x$join
  d[, c("direction", "GO", "ontology", "term",
        "n_in_list_Mt", "n_in_universe_Mt", "fold_enrich_Mt", "padj_Mt",
        "n_in_list_At", "n_in_universe_At", "fold_enrich_At", "padj_At",
        "agreement")]
}))
joined <- joined[order(joined$direction, joined$agreement, joined$padj_Mt), ]
write.table(joined, file.path(out_dir, "go_concordance_terms.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

up   <- res$up$join
down <- res$down$join
all_go <- union(up$GO, down$GO)
flag <- function(df, col) {
  v <- df[[col]][match(all_go, df$GO)]
  !is.na(v) & v
}
cross <- data.frame(
  GO         = all_go,
  ontology   = up$ontology[match(all_go, up$GO)],
  term       = up$term[match(all_go, up$GO)],
  sig_Mt_up  = flag(up, "sig_Mt"),   sig_At_up   = flag(up, "sig_At"),
  sig_Mt_down = flag(down, "sig_Mt"), sig_At_down = flag(down, "sig_At"),
  stringsAsFactors = FALSE
)
missing_term <- is.na(cross$term)
cross$ontology[missing_term] <- down$ontology[match(cross$GO[missing_term], down$GO)]
cross$term[missing_term]     <- down$term[match(cross$GO[missing_term], down$GO)]

cross$discordant <- (cross$sig_Mt_up & cross$sig_At_down) |
                    (cross$sig_Mt_down & cross$sig_At_up)
cross$concordant <- (cross$sig_Mt_up & cross$sig_At_up) |
                    (cross$sig_Mt_down & cross$sig_At_down)
n_disc <- sum(cross$discordant)
n_conc <- sum(cross$concordant)
cat("\n=== direction check ===\n")
cat("Counted over terms called enriched in BOTH species, in either direction.\n")
cat("A broad term can be both concordant and discordant -- enriched up AND\n")
cat("down within one species, on different gene subsets -- so the two counts\n")
cat("are not a partition and the discordant tail is not pure disagreement.\n")
cat(sprintf("same direction in both species     : %d\n", n_conc))
cat(sprintf("opposite direction between species : %d\n", n_disc))
cat(sprintf("share of cross-species calls that are opposite : %.0f%%\n",
            100 * n_disc / (n_conc + n_disc)))
cat("A majority-concordant result means the two protoplasting responses sit\n")
cat("on the same axis; the discordant tail is listed in the output file.\n")
write.table(cross[cross$discordant,
                  c("GO", "ontology", "term",
                    "sig_Mt_up", "sig_At_up", "sig_Mt_down", "sig_At_down")],
            file.path(out_dir, "go_direction_discordant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nOutputs written to:", out_dir, "\n")
