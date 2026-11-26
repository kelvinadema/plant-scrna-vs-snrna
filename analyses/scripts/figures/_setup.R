# Shared paths, labels and plotting style.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
})

script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
chapter_dir <- normalizePath(file.path(dirname(script_file), "../.."))
analyses_dir <- file.path(chapter_dir, "output/tables")
bulk_dir <- file.path(chapter_dir, "output/bulk")
fig_main_dir <- file.path(chapter_dir, "output/figures/main")
fig_supp_dir <- file.path(chapter_dir, "output/figures/supplementary")
table_dir <- file.path(chapter_dir, "output/tables/supplementary")
mt_annot_file <- file.path(chapter_dir, "input/bulk/medicago_deg_annotations.tsv")
section_table <- function(section, file) file.path(analyses_dir, section, file)

read_tsv_q    <- function(path, ...) read_tsv(path, show_col_types = FALSE, progress = FALSE, ...)

col_modality <- c(protoplast = "#2a78d6", nuclei = "#eb6834")
col_third    <- "#1baf7a"   # reserved for a third category (e.g. tissue / bulk)

seq_blue <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

col_diverging <- c(low = "#762a83", mid = "#f0efec", high = "#1b7837")

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#8a8984"
grid_col      <- "#e6e5e1"

lab_modality <- c(protoplast = "Protoplast", nuclei = "Nuclei")
lab_species  <- c(medicago = "M. truncatula", arabidopsis = "A. thaliana")

mt_protoplast_used <- c("rep1", "rep2", "rep3", "rep4", "pMock", "pNOD", "pBAP")

prep_groups <- c(
  "Medicago SZ unsorted nuclei"         = "Unsorted\nSZ nuclei",
  "Medicago tip nuclei, unfixed or DTT" = "Tip nuclei\n(unfixed / DTT)",
  "Medicago SZ FACS nuclei"             = "FANS\nSZ nuclei",
  "Medicago SZ protoplasts"             = "SZ\nprotoplasts",
  "Medicago 10 mm tip protoplasts"      = "Tip\nprotoplasts",
  "Arabidopsis nuclei"                  = "Nuclei\n(Farmer et al.)",
  "Arabidopsis protoplast"              = "Protoplasts\n(Ryu et al.)"
)

is_used_library <- function(sample, species, modality) {
  species == "arabidopsis" | modality == "nuclei" | sample %in% mt_protoplast_used
}

knn_modality_pairs <- c("ls__pMock__nMOCK", "ls__pNOD__nNOD", "ls__pBAP__nBAP",
                        "ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep2",
                        "ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep4",
                        "ls__GSE123013_WT_rep1__GSE155304_sNucRNA-seq_rep5")

centroid_matched <- function(top1_reciprocal, at_null) {
  top1_reciprocal >= 0.5 & !as.logical(at_null)
}

failure_labels <- c(
  retained              = "Retained",
  organellar            = "Organellar\nRNA",
  rrna                  = "rRNA\n(DTT)",
  no_droplet_separation = "No droplet\nseparation",
  insufficient_complexity    = "Insufficient\ncomplexity"
)

read_failure_modes <- function() {
  read_tsv_q(file.path(analyses_dir, "03_recovery_attrition", "fig_failure_modes.tsv")) |>
    mutate(failure_mode = recode(failure_mode, adaptive_attrition = "insufficient_complexity")) |>
    (\(d) { stopifnot(all(d$failure_mode %in% names(failure_labels))); d })()
}

cell_type_labels <- c(
  "root hair cell"               = "Root hair",
  "non-hair root epidermal cell" = "Non-hair epidermis",
  "lateral root cap"             = "Lateral root cap",
  "cortex"                       = "Cortex",
  "endodermis"                   = "Endodermis",
  "pericycle"                    = "Pericycle",
  "procambium"                   = "Procambium",
  "xylem"                        = "Xylem",
  "sieve element"                = "Sieve element",
  "companion cell"               = "Companion cell",
  "phloem"                       = "Phloem",
  "root meristem"                = "Root meristem",
  "Unassigned"                   = "Unassigned"
)

col_cell_type <- c(
  "Root hair"          = "#c51b7d",
  "Non-hair epidermis" = "#e992c1",
  "Lateral root cap"   = "#a6761d",
  "Cortex"             = "#4d9221",
  "Endodermis"         = "#a6d96a",
  "Pericycle"          = "#762a83",
  "Procambium"         = "#b491c8",
  "Xylem"              = "#5b3a1a",
  "Sieve element"      = "#01665e",
  "Companion cell"     = "#5ab4ac",
  "Phloem"             = "#35978f",
  "Root meristem"      = "#e6ab02",
  "Unassigned"         = "#cfcfcc"
)

base_size <- 7

theme_ch3 <- function(base = base_size) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      text              = element_text(colour = ink_primary),
      axis.text         = element_text(colour = ink_secondary, size = base - 1),
      axis.title        = element_text(colour = ink_primary, size = base),
      axis.line         = element_line(colour = ink_muted, linewidth = 0.3),
      axis.ticks        = element_line(colour = ink_muted, linewidth = 0.3),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base, hjust = 0),
      legend.title      = element_text(size = base - 0.5),
      legend.text       = element_text(size = base - 1),
      legend.key.size   = unit(3, "mm"),
      plot.title        = element_text(size = base + 0.5, face = "bold"),
      plot.tag          = element_text(size = base + 3, face = "bold"),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}
theme_set(theme_ch3())

theme_ch3_bw <- function(base = base_size) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(
      text               = element_text(colour = ink_primary),
      axis.text          = element_text(colour = ink_primary, size = base - 1),
      axis.title         = element_text(colour = ink_primary, size = base),
      panel.border       = element_rect(colour = ink_muted, linewidth = 0.4),
      panel.grid.major   = element_line(colour = grid_col, linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      strip.text         = element_text(face = "bold", size = base, colour = ink_primary),
      legend.title       = element_text(size = base - 0.5),
      legend.text        = element_text(size = base - 1),
      legend.key.size    = unit(3, "mm"),
      plot.tag           = element_text(size = base + 3, face = "bold"),
      plot.margin        = margin(5, 5, 5, 5)
    )
}

species_expr <- function(x) {
  lookup <- c(medicago = "italic('M. truncatula')", arabidopsis = "italic('A. thaliana')")
  unname(lookup[x])
}

save_fig <- function(plot, name, width_mm, height_mm, dir = fig_main_dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 300, device = ragg::agg_png,
         bg = "white")
  message("Wrote ", file.path(dir, name), ".{pdf,png}")
  invisible(plot)
}
