# Preprocessing

Mapping, QC and reusable integration and annotation workflows use [scPlantTools at commit 285f663a945efb0890a9519d32641da7dfa5ab35](https://github.com/kelvinadema/scPlantTools/tree/285f663a945efb0890a9519d32641da7dfa5ab35).

## Cohorts and mapping

The main Medicago comparison comprises `nMOCK`, `nNOD`, `nBAP`, `pMock`, `pNOD` and `pBAP`. The selected Arabidopsis models use GSE123013 protoplast replicates 1–3 and GSE155304 nuclear replicates 2, 4 and 5. Nuclear replicates 1 and 3 contribute to recovery/QC summaries but are excluded from those models. Additional preparations in the recovery analyses are selected explicitly in the analysis scripts.

The Medicago sample sheet incorporates corrected `nMOCK`/`nBAP` tube identities; original FASTQ directory names therefore do not always match treatment labels. Arabidopsis nuclear replicate 4 uses 10x v2 chemistry.

The following configuration paths are relative to the pinned scPlantTools repository:

| Stage | Configuration |
| --- | --- |
| Pipeline and cohort selection | `configs/pipeline_config.yaml` |
| Medicago libraries | `configs/sample_metadata.yaml` |
| Arabidopsis libraries | `configs/external/external_metadata.yaml` |
| Mapping and references | `configs/mapping_config.yaml`, `configs/index_config.yaml` |
| Chemistry | `configs/technology_config.yaml` |
| QC and filtering | `configs/qc_config.yaml`, `configs/filter_config.yaml` |
| Arabidopsis replicate-4 chemistry correction | `configs/chemistry_recheck_2026-09-21/pipeline.yaml`, `rep4_samples.yaml` in the same directory |
| Read-depth subsampling | `configs/C14/medicago_c14_f*.yaml`, `configs/C14/arabidopsis_c14_f*.yaml` |

The sample sheets also cover studies outside Chapter 3. Reference definitions include `Mt5-genome1.9.fasta`, `Mt5-annot_reformattedTID.gtf` and Arabidopsis genome/annotation files; the induced Medicago mapping profile additionally includes the *S. meliloti* reference. Reference files and sequencing inputs must be supplied separately, and cluster-specific paths must be adapted. Additional study-specific replicate-4 subsampling configurations are not included in the pinned pipeline.

## Integration and annotation

Processing proceeds through mapping and automatic QC, manual filtering where applied, integration using corrected consensus highly variable genes, marker identification, annotation and label refinement.

| Species | Selected SysVI model | Cycle weight | Learning rate |
| --- | --- | --- | --- |
| Medicago | `MtSZ_s_SV_E_l10_c50_HVGv2_CPU_l10_d2_c2` | 2 | 0.004221061640492432 |
| Arabidopsis | `At_tip_s_SV_E_l10_rep4v2_HVGv2_QCpass6_CPU_l10_d2_c50_lrdouble` | 50 | 0.005576094136353742 |

Both models use seed 42, 2,000 highly variable genes, 10 latent dimensions and two hidden layers. These are the selected fitting settings; base model recipes alone do not identify the final fits.

Marker and read-budget comparisons use 25 training and 25 test cells. `nNOD` has 47 eligible procambium cells and is excluded from this comparison because it does not meet the requirement for 50 disjoint training/test cells. It remains in the main six-library cohort.

## Downstream analyses

The Python scripts and notebooks begin with processed inputs organised by analysis under `input/`. Sections 02 and 03 share library summaries in `02_nuclear_recovery/`. Section 03 also uses per-cell QC and read-flow summaries. Subsequent sections use depth/yield summaries, refined labels and embeddings, cell-identity summaries, read-budget estimates and Nod-factor response inputs. Bulk RNA-seq and the Chapter 2 A17 Nod-factor reference are separate inputs.

Section 01 uses DESeq2's formal effect-size threshold test at adjusted P < 0.01 and |log2 fold change| > 1. Figure 7 uses corrected counts at AUROC 0.70; the Table S6 workflow also exports raw-count, alternative-target and second-seed estimates. The full input-to-figure mapping is in [figure_map.tsv](../analyses/figure_map.tsv).

Study-specific fitting configurations, manual barcode decisions, annotation-refinement records and all input-producing workflows are not distributed here. Reproduction from FASTQs therefore requires additional study materials beyond this repository and the pinned pipeline. Full Medicago bulk GO analysis also requires the original gene annotation; a gene/product annotation subset is not sufficient.
