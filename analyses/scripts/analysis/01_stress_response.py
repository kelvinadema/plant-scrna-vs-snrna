# %% [markdown]
# # Dissociation response — Figure 2, Tables S1–S2
# DEGs use DESeq2's formal test against |log2FC| > 1 at adjusted P < 0.01.
# Shared genes must change in the same direction. GO concordance uses induced BP terms.

# %%
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy import stats

if '__file__' in globals():
    ROOT = Path(__file__).resolve().parents[3]
else:
    ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
                if (p / 'input').is_dir() and (p / 'analyses/scripts').is_dir())
INPUT = ROOT / 'input/01_stress_response'
BULK = ROOT / 'output/bulk'
TABLES = ROOT / 'output/tables/01_stress_response'
TABLES.mkdir(parents=True, exist_ok=True)
# Seed missing results from supplied inputs; preserve optional R reruns.
for source in INPUT.rglob('*'):
    if source.is_file():
        destination = BULK / source.relative_to(INPUT)
        if not destination.exists():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)

def read(path):
    return pd.read_csv(path, sep='\t')

def save(frame, name):
    frame.to_csv(TABLES / name, sep='\t', index=False)

# %% [markdown]
# ## Differential expression

# %%
files = {
    ('medicago', 'protoplast'): 'isolation-effects-deg-nuclear-only/deseq2_protoplast_vs_tissue.tsv',
    ('medicago', 'nuclei'): 'isolation-effects-deg-nuclear-only/deseq2_nuclei_vs_tissue.tsv',
    ('arabidopsis', 'protoplast'): 'arabidopsis-isolation-effects-deg-nuclear-only/deseq2_protoplast_vs_tissue.tsv',
}
deg = pd.concat([read(BULK / path).assign(species=species, modality=modality,
                                         contrast=f'{modality}_vs_tissue')
                 for (species, modality), path in files.items()], ignore_index=True)
deg = deg.rename(columns={'baseMean': 'base_mean', 'log2FoldChange': 'log2fc',
                          'log2FoldChange_unshrunk': 'log2fc_unshrunk', 'lfcSE': 'lfc_se',
                          'padj_lfcThreshold': 'padj_lfc_threshold'})
deg['tested'] = deg.padj.notna()
deg['is_deg'] = deg.padj_lfc_threshold.lt(0.01)
deg['direction'] = np.where(~deg.is_deg, 'ns', np.where(deg.log2fc_unshrunk > 0, 'up', 'down'))
save(deg[['species', 'modality', 'contrast', 'gene', 'base_mean', 'log2fc', 'log2fc_unshrunk',
          'lfc_se', 'pvalue', 'padj', 'padj_lfc_threshold', 'tested', 'is_deg', 'direction']],
     'fig_lfc_per_gene.tsv')
rows = []
for (species, modality, contrast), group in deg.groupby(['species', 'modality', 'contrast']):
    for direction in ['up', 'down']:
        n = int(group.direction.eq(direction).sum())
        tested = int(group.tested.sum())
        rows.append(dict(species=species, modality=modality, contrast=contrast,
                         gate='lfc_threshold_test', direction=direction, n_genes=n,
                         genes_tested=tested, fraction_of_tested=n/tested))
counts = pd.DataFrame(rows)
save(counts, 'fig_deg_counts.tsv')
summary = counts.groupby(['species', 'modality'], as_index=False).agg(
    n_genes=('n_genes', 'sum'), genes_tested=('genes_tested', 'first'))
summary['percent_tested'] = 100 * summary.n_genes / summary.genes_tested
print(summary.to_string(index=False))

# %% [markdown]
# ## Overlap of protoplast and nuclear responses

# %%
mt = deg[deg.species.eq('medicago')]
directions = mt.pivot(index='gene', columns='modality', values='direction')
p, n = directions.protoplast, directions.nuclei
categories = pd.Series('not_called', index=directions.index)
for modality, own, other in [('protoplast', p, n), ('nuclei', n, p)]:
    mask = own.ne('ns') & other.eq('ns')
    categories.loc[mask] = modality + '_only_' + own[mask]
both = p.ne('ns') & n.ne('ns')
mask = both & p.eq(n)
categories.loc[mask] = 'shared_' + p[mask]
mask = both & p.ne(n)
categories.loc[mask] = 'opposite_protoplast_' + p[mask] + '_nuclei_' + n[mask]
overlap = categories.value_counts().rename_axis('category').rename('n_genes').reset_index()
overlap.insert(0, 'species', 'medicago')
save(overlap, 'fig_overlap_categories.tsv')
wide = directions.rename(columns=lambda c: f'direction_{c}').reset_index()
for column, prefix in [('log2fc', 'log2fc'), ('log2fc_unshrunk', 'log2fc_unshrunk')]:
    lfc = mt.pivot(index='gene', columns='modality', values=column).add_prefix(prefix + '_')
    wide = wide.merge(lfc, on='gene')
wide['species'] = 'medicago'
wide['category'] = wide.gene.map(categories)
save(wide, 'fig_lfc_mt_protoplast_vs_nuclei.tsv')
membership = mt[mt.is_deg][['species', 'gene', 'modality', 'direction']].copy()
membership['set'] = membership.modality + '_' + membership.direction
membership['category'] = membership.gene.map(categories)
save(membership, 'fig_overlap_membership.tsv')
print(overlap[overlap.category.str.startswith('shared')].to_string(index=False))

# %% [markdown]
# ## GO enrichment and cross-species concordance
# GO results use gene sets selected by the formal effect-size threshold test.

# %%
shared = read(BULK / 'isolation-effects-go-strict/go_enrichment_shared_unique.tsv')
shared['significant'] = shared.padj.lt(0.05)
save(shared, 'fig_go_shared_dotplot.tsv')
cross = read(BULK / 'cross-species-protoplast-strict/go_concordance_terms.tsv')
cross['gene_gate'] = 'lfc_threshold_test'
save(cross, 'fig_go_cross_species.tsv')
rows = []
for direction, group in cross.groupby('direction'):
    for ontology in ['all', 'BP']:
        d = group if ontology == 'all' else group[group.ontology.eq(ontology)]
        mt_sig, at_sig = d.padj_Mt.lt(0.05), d.padj_At.lt(0.05)
        a, b, c, e = [int(mask.sum()) for mask in
                      [mt_sig & at_sig, mt_sig & ~at_sig, ~mt_sig & at_sig, ~mt_sig & ~at_sig]]
        contingency = [[a, b], [c, e]]
        # Conditional estimate matches R fisher.test's odds ratio.
        odds = stats.contingency.odds_ratio(contingency, kind='conditional').statistic
        rows.append(dict(direction=direction, ontology=ontology, co_tested=len(d),
                         enriched_both=a, medicago_only=b, arabidopsis_only=c, neither=e,
                         odds_ratio=odds, fisher_p=stats.fisher_exact(contingency).pvalue,
                         spearman_rho_fold=stats.spearmanr(d.fold_enrich_Mt, d.fold_enrich_At).statistic))
go_summary = pd.DataFrame(rows)
save(go_summary, 'go_cross_species_summary.tsv')
print(go_summary.query("direction == 'up' and ontology == 'BP'").to_string(index=False))
