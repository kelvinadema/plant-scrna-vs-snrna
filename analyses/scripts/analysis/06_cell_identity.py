# %% [markdown]
# # Cell identity
# 25 training and 25 held-out cells per class. Primary depths: 1,500 UMIs (Medicago), 500 (Arabidopsis).


# %%
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'analyses/scripts').is_dir())

import numpy as np
import pandas as pd

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)
pd.set_option('display.max_rows', 200)

D = ROOT / 'input/06_cell_identity'
T = ROOT / 'output/tables/06_cell_identity'
T.mkdir(parents=True, exist_ok=True)

MT_NUC, MT_PROT = ['nMOCK', 'nNOD', 'nBAP'], ['pMock', 'pNOD', 'pBAP']
MT_TREAT = {'nMOCK': 'mock', 'pMock': 'mock', 'nNOD': 'nod_factor_3h', 'pNOD': 'nod_factor_3h',
            'nBAP': 'bap_3h', 'pBAP': 'bap_3h'}
SHORT = {'GSE123013_WT_rep1': 'Ryu rep1', 'GSE123013_WT_rep2': 'Ryu rep2', 'GSE123013_WT_rep3': 'Ryu rep3',
         'GSE155304_sNucRNA-seq_rep2': 'Farmer rep2', 'GSE155304_sNucRNA-seq_rep4': 'Farmer rep4',
         'GSE155304_sNucRNA-seq_rep5': 'Farmer rep5'}
SHORT.update({s: s for s in MT_NUC + MT_PROT})
LONG = {v: k for k, v in SHORT.items()}


def modality(sample):
    """nuclei / protoplast from the library name (Medicago n*/p*, Arabidopsis GSE155304 / GSE123013)."""
    s = LONG.get(sample, sample)
    if s in MT_NUC or s.startswith('GSE155304'):
        return 'nuclei'
    if s in MT_PROT or s.startswith('GSE123013'):
        return 'protoplast'
    raise KeyError(sample)


def replicate(sample):
    s = LONG.get(sample, sample)
    return s.rsplit('_', 1)[-1] if s.startswith('GSE') else ''


def maxdiff(a, b):
    """Maximum difference; unequal NaN positions fail."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        return np.inf
    m = ~np.isnan(a)
    return float(np.abs(a[m] - b[m]).max()) if m.any() else 0.0


# %% [markdown]
# ## 1. Held-out marker AUROC (C12)

# %%
perclass = pd.read_csv(D / 'markers/C12_perclass.tsv', sep='\t')
perclass = perclass[perclass.primary & ~perclass.pair_role.eq('cross_modality_other_treatment')].copy()
pairs12 = pd.read_csv(D / 'markers/C12_pairs.tsv', sep='\t')
pairs12 = pairs12[pairs12.primary & ~pairs12.pair_role.eq('cross_modality_other_treatment')].copy()
summary12 = pd.read_csv(D / 'markers/C12_summary.tsv', sep='\t')
summary12 = summary12[summary12.primary & ~summary12.pair_role.eq('cross_modality_other_treatment')].copy()

KEYS12 = ['species', 'depth_umis', 'n_train_per_class', 'n_test_per_class', 'n_classes', 'classes',
          'repeat', 'primary', 'source', 'target', 'pair_role', 'source_modality', 'target_modality',
          'source_treatment', 'target_treatment', 'validation', 'label_source']
re_pairs = perclass.groupby(KEYS12, dropna=False, observed=True).agg(
    heldout_gene_auc=('heldout_gene_auc', 'mean'), heldout_program_auc=('heldout_program_auc', 'mean')).reset_index()
chk = re_pairs.merge(pairs12[KEYS12 + ['heldout_gene_auc', 'heldout_program_auc']], on=KEYS12, suffixes=('_re', ''),
                     validate='one_to_one')
assert len(chk) == len(pairs12) == len(re_pairs)
d_pairs = max(maxdiff(chk.heldout_gene_auc_re, chk.heldout_gene_auc), maxdiff(chk.heldout_program_auc_re, chk.heldout_program_auc))
print(f'C12_pairs rebuilt from C12_perclass: {len(chk)} rows, max |diff| = {d_pairs:.2e}')
assert d_pairs < 1e-3, 'perclass -> pairs does not reproduce'

SUMKEYS = [k for k in KEYS12 if k != 'repeat']
re_summary = pairs12.groupby(SUMKEYS, dropna=False, observed=True).agg(
    heldout_gene_auc=('heldout_gene_auc', 'mean'), heldout_gene_auc_split_sd=('heldout_gene_auc', 'std'),
    heldout_gene_auc_split_min=('heldout_gene_auc', 'min'), heldout_gene_auc_split_max=('heldout_gene_auc', 'max'),
    heldout_program_auc=('heldout_program_auc', 'mean'), heldout_program_auc_split_sd=('heldout_program_auc', 'std'),
    heldout_program_auc_split_min=('heldout_program_auc', 'min'), heldout_program_auc_split_max=('heldout_program_auc', 'max'),
    n_splits=('repeat', 'nunique')).reset_index()
chk = re_summary.merge(summary12[SUMKEYS + ['heldout_gene_auc', 'heldout_program_auc']], on=SUMKEYS,
                       suffixes=('_re', ''), validate='one_to_one')
assert len(chk) == len(summary12)
d_sum = max(maxdiff(chk.heldout_gene_auc_re, chk.heldout_gene_auc), maxdiff(chk.heldout_program_auc_re, chk.heldout_program_auc))
print(f'C12_summary rebuilt from C12_pairs: {len(chk)} rows, max |diff| = {d_sum:.2e}')
assert d_sum < 1e-9


# %%
TEST_TYPE = {'own_library_heldout': 'own', 'within_modality_library': 'within_modality',
             'matched_modality_pair': 'cross_modality', 'cross_modality_other_treatment': 'cross_modality_other_treatment'}


def tidy_markers(df, per_split):
    base = pd.DataFrame(dict(
        species=df.species, depth_umis=df.depth_umis, n_train_per_class=df.n_train_per_class,
        n_test_per_class=df.n_test_per_class, classes=df.classes, primary=df.primary,
        source_sample=df.source, source_label=df.source.map(SHORT), source_modality=df.source_modality,
        source_treatment=df.source_treatment, source_replicate=df.source.map(replicate),
        target_sample=df.target, target_label=df.target.map(SHORT), target_modality=df.target_modality,
        target_treatment=df.target_treatment, target_replicate=df.target.map(replicate),
        test_type=df.pair_role.map(TEST_TYPE)))
    out = []
    for ms, col in (('single_gene', 'heldout_gene_auc'), ('programme', 'heldout_program_auc')):
        b = base.copy()
        b['marker_set'] = ms
        if per_split:
            b.insert(b.columns.get_loc('test_type') + 1, 'split', df['repeat'].values)
            b['auroc'] = df[col].values
        else:
            b['auroc_mean'] = df[col].values
            b['auroc_split_sd'] = df[f'{col}_split_sd'].values
            b['auroc_split_min'] = df[f'{col}_split_min'].values
            b['auroc_split_max'] = df[f'{col}_split_max'].values
            b['n_splits'] = df['n_splits'].values
        out.append(b)
    return pd.concat(out, ignore_index=True)


fig_markers_split = tidy_markers(pairs12, per_split=True)
fig_markers_summary = tidy_markers(re_summary, per_split=False)
fig_markers_split.to_csv(T / 'fig_markers_auroc_per_split.tsv', sep='\t', index=False)
fig_markers_summary.to_csv(T / 'fig_markers_auroc_summary.tsv', sep='\t', index=False)

pc = perclass[perclass.primary]
fig_markers_class = pd.DataFrame(dict(
    species=pc.species, source_sample=pc.source, source_modality=pc.source_modality, source_treatment=pc.source_treatment,
    target_sample=pc.target, target_modality=pc.target_modality, target_treatment=pc.target_treatment,
    test_type=pc.pair_role.map(TEST_TYPE), split=pc['repeat'], cell_type=pc.label,
    n_markers=pc.n_markers, n_markers_detected_in_target_test=pc.n_markers_detected_in_target_test,
    single_gene_auroc=pc.heldout_gene_auc, programme_auroc=pc.heldout_program_auc))
fig_markers_class.to_csv(T / 'fig_markers_auroc_per_class.tsv', sep='\t', index=False)
print(len(fig_markers_split), len(fig_markers_summary), len(fig_markers_class))

prim = re_summary[re_summary.primary].copy()
prim['test_type'] = prim.pair_role.map(TEST_TYPE)
view = prim[['species', 'source', 'target', 'test_type', 'heldout_gene_auc', 'heldout_program_auc', 'n_splits']].copy()
view['source'] = view.source.map(SHORT)
view['target'] = view.target.map(SHORT)
view.sort_values(['species', 'test_type', 'source', 'target'])


# %% [markdown]
# ## 2. Neighbour agreement in the V2 joint label space

# %%
for raw, out in [('labelspace_purity_by_type.tsv','fig_knn_purity_by_type.tsv'),
                 ('labelspace_pooled_by_direction.tsv','fig_knn_pooled_by_direction.tsv'),
                 ('labelspace_pooled_ratio.tsv','fig_knn_pooled_ratio.tsv')]:
    df = pd.concat([pd.read_csv(D/'knn'/sp/raw, sep='\t') for sp in ('medicago','arabidopsis')], ignore_index=True)
    df = df[df.condition.eq('base')].copy()
    df = df.rename(columns={'label':'cell_type', 'label_ok_both_sides':'shared_both_sides', 'null_purity_mean':'null_purity'})
    df['geometry'] = 'V2_species_joint_X_scVI'
    if 'source_sample' in df:
        df['source_label'] = df.source_sample.map(SHORT)
        df['target_label'] = df.target_sample.map(SHORT)
    df.to_csv(T/out, sep='\t', index=False)


# %% [markdown]
# ## 3. Centroid correlations and margins

# %%
def decorate(df):
    df = df[df.condition.eq('centred') & df.layer.eq('corrected') & df.exclusion.eq('none')].copy()
    df['primary_arm'] = df.condition.eq('centred') & df.layer.eq('corrected') & df.exclusion.eq('none')
    df['exclusion_arm'] = df.condition.eq('centred') & df.layer.eq('corrected') & df.exclusion.eq('c9_either_up')
    for suffix in ('row','col'):
        col = 'member_' + suffix
        if col in df:
            df['sample_'+suffix] = df[col]
            df['label_'+suffix] = df[col].map(SHORT)
            df['modality_'+suffix] = df[col].map(modality)
    return df
for raw, out in [('effects','per_type'), ('matrix','matrix')]:
    df = pd.concat([pd.read_csv(D/'centroids'/f'C16_{prefix}{raw}.tsv', sep='\t') for prefix in ('',)], ignore_index=True)
    df = decorate(df).rename(columns={'label':'cell_type'})
    if raw == 'matrix':
        df['diagonal'] = df.type_row.eq(df.type_col)
    df.to_csv(T/f'fig_centroid_{out}.tsv', sep='\t', index=False)
df = decorate(pd.read_csv(D/'centroids/C16_matrix_level.tsv', sep='\t'))
df.to_csv(T/'fig_centroid_pair_level.tsv', sep='\t', index=False)
assert (fig_markers_summary.n_train_per_class == 25).all()
assert (fig_markers_summary.n_test_per_class == 25).all()
assert 'nNOD' not in set(fig_markers_summary.source_sample)
fig_markers_summary['analysis_variant'] = 'V2_25_25'
fig_markers_summary.to_csv(T/'fig_markers_auroc_summary.tsv', sep='\t', index=False)
print('Marker aggregates verified; neighbour and centroid tables exported.')
