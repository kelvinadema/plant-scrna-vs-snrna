# plant-scRNA-vs-snRNA

Scripts and notebooks accompanying Chapter 3 of my PhD thesis, **Modality dependent trade-offs in reliability and information recovery for plant single-cell RNA sequencing**.

The study compares protoplast and nuclear transcriptomes from the *Medicago truncatula* root susceptible zone under mock, Nod factor and cytokinin treatments, together with published *Arabidopsis thaliana* datasets. The analyses cover preparation-associated responses, recovery and molecular yield, cell identity and treatment-response recovery.

## Contents

| Path | Contents |
| --- | --- |
| `analyses/notebooks/` | Eight analysis notebooks. |
| `analyses/scripts/analysis/` | Python versions of the notebooks. |
| `analyses/scripts/figures/` | R scripts for figures and supplementary tables. |
| `analyses/scripts/bulk/` | Bulk differential-expression and GO analyses. |
| [Figure map](analyses/figure_map.tsv) | Figures and supplementary tables mapped to scripts and expected input paths. |
| [Preprocessing](processing/README.md) | Pipeline version, cohort selection and model settings. |

## Data availability

The study data are not yet publicly available. The GEO submission is private pending release. This repository contains analysis source and environment definitions; sequencing files, count matrices, cell metadata, processed data, saved notebook outputs and result files are not included. Access details will be added when the data are released.

The published Arabidopsis datasets are GSE155304 and GSE123013 for the single-cell comparisons, and GSE123818 for the bulk comparison.

## Running the analyses

Use Python 3.12. From the repository root:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r analyses/requirements.txt
```

Execution requires the processed study inputs under `input/`, with the paths used in the scripts and [figure map](analyses/figure_map.tsv). Once those inputs are available:

```bash
cd analyses
python scripts/reproduce.py
```

The runner executes the eight analyses in order and writes results to `output/` at the project root. Paths in the figure map are also relative to the project root. To run selected analyses, use `python scripts/reproduce.py --section 04 08`. Section 01 copies precomputed bulk results into the output directory only when those files are absent, preserving any separately regenerated bulk results.

To use the notebooks, install `analyses/requirements-notebooks.txt` from the repository root, then launch `jupyter lab analyses/notebooks`. The notebooks and scripts implement the same analyses. A fresh clone can be read without data but cannot execute the analyses until the inputs are available.

## Figures and supplementary tables

After the Python analyses, run from `analyses/`:

```bash
conda env create -f environment-r.yml
conda run -n ch3-figures Rscript scripts/figures/render_all.R
```

Relative to the project root, figures are written to `output/figures/` and supplementary tables to `output/tables/supplementary/`. Figure 1 is manually assembled; other thesis figures may include manual layout adjustments.

The optional bulk workflow has a separate environment in `analyses/scripts/bulk/environment.yml`. Full Medicago GO analysis also requires the original gene annotation at `input/bulk/medicago_gene_annotation.tsv`.

## Workflow scope

These scripts start from processed inputs. Mapping, QC and reusable integration and annotation workflows are maintained in [scPlantTools](https://github.com/kelvinadema/scPlantTools). The [preprocessing notes](processing/README.md) describe the version and settings used for this study. Study-specific model fitting, manual filtering and generation of all intermediate inputs are not included here, so this is not a complete FASTQ-to-figure workflow.
