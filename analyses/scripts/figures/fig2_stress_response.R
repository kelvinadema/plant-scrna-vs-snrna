.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
out_dir <- file.path(fig_main_dir, "F2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(ggforce)
library(dplyr)
library(tidytext)

s1 <- "01_stress_response"

deg_counts <- read_tsv_q(section_table(s1, "fig_deg_counts.tsv")) |>
  filter(gate == "lfc_threshold_test") |>
  mutate(signed_n = if_else(direction == "up", n_genes, -n_genes)) |> 
  select(c(species, modality, contrast, direction, n_genes, genes_tested, fraction_of_tested))

deg_totals <- deg_counts |> 
  group_by(species, modality) |>
  summarise(n = sum(n_genes), tested = first(genes_tested), .groups = "drop") |>
  mutate(x_lab = sprintf("%s\n%s (%.1f%%)", lab_modality[modality],
                         format(n, big.mark = ","), 100 * n / tested))

deg_counts <- deg_counts |>
  left_join(deg_totals |> select(species, modality, x_lab), by = c("species", "modality")) |>
  mutate(species = factor(species, levels = c("medicago", "arabidopsis")),
         x_lab = factor(x_lab, levels = deg_totals |> arrange(desc(modality == "protoplast")) |>
                          pull(x_lab)))

y_max <- max(deg_counts$n_genes) * 1.15

pA <- ggplot(deg_counts, aes(x = x_lab, y = n_genes, fill = direction)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7), colour = "black") +
  facet_grid(~species, scales = "free_x", space = "free_x",
             labeller = labeller(species = as_labeller(species_expr, label_parsed))) +
  scale_fill_manual(values = c(up = "#009E73",
                               down = "#CC79A7")) +
  theme_bw(base_size = base_size) +
  theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = ink_primary, lineheight = 0.95)) +
  scale_y_continuous(labels = function(v) format(abs(v), big.mark = ","),
                     limits = c(0, y_max), expand = expansion(0)) + 
  xlab("") +
  ylab('Number of DEGs')

membership <- read_tsv_q(section_table(s1, "fig_overlap_membership.tsv")) |>
  filter(species == "medicago")

venn_counts <- membership |>
  distinct(gene, modality, direction) |>
  group_by(direction) |>
  summarise(
    protoplast_only = n_distinct(setdiff(gene[modality == "protoplast"], gene[modality == "nuclei"])),
    shared          = n_distinct(intersect(gene[modality == "protoplast"], gene[modality == "nuclei"])),
    nuclei_only     = n_distinct(setdiff(gene[modality == "nuclei"], gene[modality == "protoplast"])),
    .groups = "drop") |>
  mutate(panel = factor(if_else(direction == "up", "Induced", "Repressed"),
                        levels = c("Induced", "Repressed")))

r <- 0.8
d <- 0.50
circles <- expand_grid(panel = venn_counts$panel,
                       tibble(modality = c("protoplast", "nuclei"), x0 = c(-d, d)))

venn_labels <- venn_counts |>
  pivot_longer(c(protoplast_only, shared, nuclei_only), names_to = "region", values_to = "n") |>
  mutate(x = c(protoplast_only = -d - 0.45, shared = 0, nuclei_only = d + 0.45)[region],
         face = if_else(region == "shared", "bold", "plain"))

set_labels <- expand_grid(panel = venn_counts$panel,
                          tibble(modality = c("protoplast", "nuclei"), x = c(-0.12, 0.12),
                                 hjust = c(1, 0))) |>
  mutate(label = lab_modality[modality])

pB <- ggplot() + geom_circle(data = circles, aes(x0 = x0, y0 = 0, r = r, fill = modality), color = 'black', alpha = 0.25, linewidth = 0.3) +
  geom_text(
    data = venn_labels,
    aes(x, 0, label = format(n, big.mark=",")),
    size = if_else(venn_labels$region == "shared",(base_size + 1)/.pt, base_size/.pt), fontface = "bold") +
  geom_text(data = set_labels, aes(x, r + 0.1, label = label, colour = modality, hjust = hjust),
            size = (base_size + 1.5) / .pt, fontface = "bold") +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_manual(values = col_modality, guide = "none") +
  scale_colour_manual(values = col_modality, guide = "none") +
  coord_fixed(clip = "off") +
  theme_void(base_size = base_size + 2) +
  theme(strip.text    = element_text(face = "bold", size = base_size + 2, hjust = 0.5,
                                     margin = margin(b = 4)),
        plot.tag      = element_text(size = base_size + 3, face = "bold"),
        plot.margin   = margin(3, 3, 3, 3, "mm"))

go_terms <- read_tsv_q(file.path(bulk_dir, "cross-species-protoplast-strict/go_concordance_terms.tsv")) |>
  filter(direction == "up", ontology == "BP", agreement == "both")

padj_cap <- 20
go_long <- go_terms |>
  transmute(
    GO,
    term,
    mean_fold = (fold_enrich_Mt + fold_enrich_At) / 2,
    medicago_fold = fold_enrich_Mt,
    medicago_padj = padj_Mt,
    medicago_n = n_in_list_Mt,
    arabidopsis_fold = fold_enrich_At,
    arabidopsis_padj = padj_At,
    arabidopsis_n = n_in_list_At
  ) |>
  pivot_longer(
    starts_with(c("medicago_", "arabidopsis_")),
    names_to = c("species", ".value"),
    names_sep = "_"
  ) |>
  mutate(
    species = factor(species, levels = c("medicago", "arabidopsis")),
    mlog_p  = pmin(-log10(padj), padj_cap),
    term_lab = str_wrap(term, 42)
  )

go_plot <- go_long %>%
  mutate(
    theme = case_when(
      term %in% c(
        "regulation of defense response",
        "response to wounding",
        "response to chitin"
      ) ~ "Defence and immune response",
      
      term %in% c(
        "cellular response to lipid",
        "response to jasmonic acid",
        "cellular response to jasmonic acid stimulus",
        "jasmonic acid mediated signaling pathway",
        "response to fatty acid",
        "cellular response to fatty acid",
        "regulation of jasmonic acid mediated signaling pathway"
      ) ~ "Jasmonate and lipid signalling",
      
      term %in% c(
        "cellular response to endogenous stimulus",
        "cellular response to hormone stimulus",
        "hormone-mediated signaling pathway",
        "response to abscisic acid"
      ) ~ "Hormone signalling",
      
      term %in% c(
        "response to temperature stimulus",
        "response to osmotic stress",
        "response to salt stress",
        "response to heat"
      ) ~ "Stress response",
      
      term %in% c(
        "cellular response to oxygen-containing compound",
        "glutathione metabolic process",
        "response to hydrogen peroxide",
        "toxin metabolic process"
      ) ~ "Oxidative stress and detoxification",
      
      term %in% c(
        "secondary metabolic process",
        "cellular ketone metabolic process",
        "response to alcohol"
      ) ~ "Secondary metabolism",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(theme))

theme_levels <- c(
  "Defence and immune response",
  "Jasmonate and lipid signalling",
  "Hormone signalling",
  "Stress response",
  "Oxidative stress and detoxification",
  "Secondary metabolism"
)

go_plot <- go_plot %>%
  mutate(
    theme = factor(theme, levels = theme_levels, labels = str_wrap(theme_levels, 18)),
    term_f = reorder_within(term_lab, mean_fold, theme)
  )

pC <- ggplot(
  go_plot,
  aes(x = fold, y = term_f)
) +
  geom_point(
    aes(colour = mlog_p),
    size = 3.6
  ) +
  facet_grid(
    rows = vars(theme),
    cols = vars(species),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    labeller = labeller(
      species = as_labeller(species_expr, label_parsed)
    )
  ) +
  scale_y_reordered() +
  scale_colour_gradientn(
    colours = seq_blue[2:7],
    name = expression(-log[10]~italic(P)[adj]),
    breaks = c(2, 5, 10, 15, 20),
    labels = c("2", "5", "10", "15", "\u226520")
  ) +
  scale_x_continuous(
    name = "Fold enrichment",
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(y = NULL) +
  theme_bw(base_size = base_size) +
  theme(
    strip.placement = "outside",
    
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = base_size + 1,
      colour = ink_primary
    ),
    
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = base_size + 1,
      colour = ink_primary,
      margin = margin(r = 6, l = 6)
    ),
    
    axis.text.y = element_text(
      colour = ink_primary,
      size = base_size
    ),
    axis.text.x = element_text(
      colour = ink_primary,
      size = base_size
    ),
    
    panel.grid.major.y = element_line(
      colour = grid_col,
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    
    panel.spacing.y = unit(2.5, "mm"),
    panel.spacing.x = unit(2.5, "mm"),
    
    legend.position = "right")

pA <- pA + theme(axis.text.x = element_text(size = base_size - 1), plot.margin = margin(5, 5, 5, 5))
pB <- pB +theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
pC <- pC + theme(
  panel.spacing.y = unit(1, "mm"),
  panel.spacing.x = unit(1, "mm"),
  strip.text.y.left = element_text(size = base_size - 1), 
  axis.text.y = element_text(size = base_size - 1), 
  plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

pA <- pA + theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))
pC <- pC + theme(strip.text.y.left = element_text(angle = 0, size = base_size - 1))

ggsave(filename = 'p2A.png', pA, path = out_dir, width = 180, height = 100, units = "mm", dpi = 600)
ggsave(filename = 'p2B.png', pB, path = out_dir, width = 100, height = 50, units = "mm", dpi = 600)
ggsave(filename = 'p2C.png', pC, path = out_dir, width = 180, height = 180, units = "mm", dpi = 600)

left <- cowplot::plot_grid(pA, pB, ncol = 1, rel_heights = c(2, 1), labels = c("A", "B"), label_size = 10)
fig2 <- cowplot::plot_grid(left, pC, ncol = 2, rel_widths = c(1, 1.7), labels = c("", "C"), label_size = 10)
save_fig(fig2, "Fig2_stress_response", width_mm = 260, height_mm = 160, dir = out_dir)
