# %% [markdown]
# # Read yield and matched-depth detection
# Full-depth yield, matched reads per cell and matched corrected-UMI detection.


# %%
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'scripts').is_dir())
import json
import numpy as np
import pandas as pd

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)

DATA = ROOT / 'input/04_read_yield'
TABLES = ROOT / 'output/tables/04_read_yield'
TABLES.mkdir(parents=True, exist_ok=True)

core = pd.read_csv(DATA / 'C12_chapter_core_metrics.tsv', sep='\t')
elig = pd.read_csv(DATA / 'C12_detection_eligibility.tsv', sep='\t')
fate = pd.read_csv(DATA / 'C1/C1_readfate_perlib.tsv', sep='\t')
depth = pd.read_csv(DATA / 'C1/C1_depth_perlib.tsv', sep='\t')
det = pd.read_csv(DATA / 'C1/C1_detection_matched.tsv', sep='\t')
rare = pd.read_csv(DATA / 'C1/C1_rarefaction.tsv', sep='\t')
mdepth = pd.read_csv(DATA / 'C1/C1_matched_depth.tsv', sep='\t')
fc_perlib = pd.read_csv(DATA / 'C14_fixedcell_perlib.tsv', sep='\t')
fc_curve = pd.read_csv(DATA / 'C14_fixedcell_curve.tsv', sep='\t')
my_pairs_stored = pd.read_csv(DATA / 'C14_matched_yield/C14_matched_yield_pairs.tsv', sep='\t')
my_curve_stored = pd.read_csv(DATA / 'C14_matched_yield/C14_yield_curve.tsv', sep='\t')

MT_N, MT_P = ['nMOCK', 'nNOD', 'nBAP'], ['pMock', 'pNOD', 'pBAP']
AT_N = ['GSE155304_sNucRNA-seq_rep2', 'GSE155304_sNucRNA-seq_rep4', 'GSE155304_sNucRNA-seq_rep5']
AT_P = ['GSE123013_WT_rep1', 'GSE123013_WT_rep2', 'GSE123013_WT_rep3']
PRIMARY = MT_N + MT_P + AT_N + AT_P


def label(s):
    return s.replace('GSE155304_sNucRNA-seq_', 'Farmer ').replace('GSE123013_WT_', 'Ryu ')


def replicate(s):
    return int(s.rsplit('rep', 1)[1]) if 'rep' in s else pd.NA


def annotate(df, col='sample'):
    """Add label/replicate/primary flag columns next to `sample`."""
    df = df.copy()
    df['sample_label'] = df[col].map(label)
    df['replicate'] = pd.array([replicate(s) if s.startswith('GSE') else pd.NA for s in df[col]], dtype='Int64')
    df['primary_comparison'] = df[col].isin(PRIMARY)
    return df


core_i = core.set_index('sample')
print(core[['sample', 'species', 'modality', 'treatment', 'n_cells', 'frac_retained', 'matched_depth',
            'genes_at_matched', 'frac_eligible']].to_string(index=False))


# %% [markdown]
# ## 1. Full-depth yield

# %%
f = fate.set_index('sample')
assert np.allclose(core_i.frac_retained, f.loc[core_i.index, 'frac_retained'], rtol=0, atol=1e-15)
part_reads = f.frac_unmapped + f.frac_mapped_outside_permit_barcodes + f.frac_umi_duplicates + f.frac_molecules
part_mol = f.frac_removed_cellbender + f.frac_molecules_not_retained + f.frac_retained
print('max |reads partition - 1|      :', float((part_reads - 1).abs().max()))
print('max |molecule partition - mol| :', float((part_mol - f.frac_molecules).abs().max()))
assert (part_reads - 1).abs().max() < 1e-9 and (part_mol - f.frac_molecules).abs().max() < 1e-9

pct = core_i.frac_retained * 100
print((pct.loc[PRIMARY]).round(3).to_string())


# %% [markdown]
# ### Tables: per-library yield and read fate

# %%
yield_tab = core[['sample', 'species', 'modality', 'treatment', 'tissue', 'n_cells', 'total_reads', 'frac_retained',
                  'reads_per_retained_umi', 'sequencing_saturation', 'reads_per_retained_cell',
                  'median_total_corrected', 'cellbender_removed_frac_median']].copy()
yield_tab = yield_tab.merge(fate[['sample', 'frac_unmapped', 'frac_mapped_outside_permit_barcodes', 'frac_umi_duplicates',
                                  'frac_molecules', 'frac_removed_cellbender', 'frac_molecules_not_retained',
                                  'frac_retained_spliced_resolved', 'frac_retained_unspliced_resolved',
                                  'frac_retained_ambiguous']], on='sample', validate='one_to_one')
yield_tab = annotate(yield_tab)
yield_tab[yield_tab['sample'].isin(PRIMARY)].to_csv(TABLES / 'fig_yield_perlib.tsv', sep='\t', index=False)

FATES = [('unmapped', 'frac_unmapped'), ('mapped_outside_permit_barcodes', 'frac_mapped_outside_permit_barcodes'),
         ('umi_duplicates', 'frac_umi_duplicates'), ('removed_by_cellbender', 'frac_removed_cellbender'),
         ('molecules_not_in_retained_cells', 'frac_molecules_not_retained'), ('retained', 'frac_retained')]
stack = fate[fate['sample'].isin(PRIMARY)].melt(
    id_vars=['sample', 'species', 'modality', 'treatment'], value_vars=[c for _, c in FATES],
    var_name='source_column', value_name='fraction_of_input_reads')
stack['fate'] = stack.source_column.map({c: n for n, c in FATES})
stack['fate_order'] = stack.source_column.map({c: i for i, (_, c) in enumerate(FATES)})
stack = annotate(stack).sort_values(['species', 'modality', 'sample', 'fate_order'])
assert np.allclose(stack.groupby('sample').fraction_of_input_reads.sum(), 1)
stack[stack['sample'].isin(PRIMARY)].to_csv(TABLES / 'fig_readfate_stack.tsv', sep='\t', index=False)
print(len(yield_tab), 'rows fig_yield_perlib;', len(stack), 'rows fig_readfate_stack')


# %% [markdown]
# ## 2. Depth-matched yield

# %%
def c14_load(perlib, curve):
    flags = curve[curve.depth_mode == 'native'][['sample', 'fraction', 'seed_no', 'cell_set_intact', 'cellbender_anomaly']]
    p = perlib.merge(flags, on=['sample', 'fraction', 'seed_no'], how='left', validate='one_to_one')
    p = p[p.seed_no == 1].copy()
    p['yield'] = p.total_corrected_umis / p.output_pairs
    p['reads_per_retained_umi'] = p.output_pairs / p.total_corrected_umis
    p['valid'] = p.cell_set_intact.astype(bool) & ~p.cellbender_anomaly.astype(bool)
    return p


def yield_at(p, sample, reads):
    part = p[(p['sample'] == sample) & p.valid].sort_values('reads_per_fixed_cell')
    x, y = np.log(part.reads_per_fixed_cell.to_numpy()), np.log(part['yield'].to_numpy())
    if not x[0] - 1e-9 <= np.log(reads) <= x[-1] + 1e-9:
        raise ValueError(f'{sample}: {reads:.0f} reads/cell outside valid range')
    return float(np.exp(np.interp(np.log(reads), x, y)))


def short(s):
    return s.replace('GSE155304_sNucRNA-seq_', 'Farmer_').replace('GSE123013_WT_', 'Ryu_')


def pair_row(p, c1, species, nuc, pro, treatment):
    full = p[p.fraction == 1.0].set_index('sample')
    rn, rp = full.loc[nuc, 'reads_per_fixed_cell'], full.loc[pro, 'reads_per_fixed_cell']
    matched = min(rn, rp)
    yn, yp = yield_at(p, nuc, matched), yield_at(p, pro, matched)
    fn, fp = full.loc[nuc, 'yield'], full.loc[pro, 'yield']
    return dict(species=species, nuclei=short(nuc), protoplast=short(pro), treatment=treatment,
                nuclei_sample=nuc, protoplast_sample=pro,
                reads_per_cell_full_n=rn, reads_per_cell_full_p=rp, depth_ratio_full=rn / rp,
                saturation_n=c1.loc[nuc, 'sequencing_saturation'], saturation_p=c1.loc[pro, 'sequencing_saturation'],
                dup_frac_full_n=c1.loc[nuc, 'frac_umi_duplicates'], dup_frac_full_p=c1.loc[pro, 'frac_umi_duplicates'],
                yield_full_n=fn, yield_full_p=fp, yield_ratio_full=fp / fn,
                matched_reads_per_cell=matched,
                matched_at='protoplast full depth' if rp <= rn else 'nuclei full depth',
                yield_matched_n=yn, yield_matched_p=yp, yield_ratio_matched=yp / yn,
                share_of_full_gap_from_depth=1 - np.log(yp / yn) / np.log(fp / fn))


p14 = c14_load(fc_perlib, fc_curve)
full = p14[p14.fraction == 1.0].set_index('sample')
full_diff = (full['yield'] - f.loc[full.index, 'frac_retained']).abs()
print('full-depth C14 yield vs C1 frac_retained, max |diff|:', float(full_diff.max()))
assert (full_diff < 1e-4).all()

rows = [pair_row(p14, f, 'medicago', n, p, t) for n, p, t in zip(MT_N, MT_P, ['mock', 'nod_factor_3h', 'bap_3h'])]
rows += [pair_row(p14, f, 'arabidopsis', n, p, 'untreated') for n in AT_N for p in AT_P]
pairs = pd.DataFrame(rows)

# Check agreement with the supplied matched-yield results.
cmp = pairs.merge(my_pairs_stored, on=['species', 'nuclei', 'protoplast'], suffixes=('', '_stored'), validate='one_to_one')
num = [c for c in my_pairs_stored.columns if my_pairs_stored[c].dtype.kind == 'f']
max_rel = max(
    float(((cmp[c] - cmp[c + '_stored']).abs() / cmp[c + '_stored'].abs().clip(lower=1e-12)).max())
    for c in num
)
print(f'recomputed pairs vs stored C14_matched_yield_pairs.tsv: {len(cmp)} rows, max relative diff {max_rel:.2e}')
assert len(cmp) == 12 and max_rel < 1e-9

curve = p14[['species', 'sample', 'modality', 'treatment', 'fraction', 'reads_per_fixed_cell', 'total_corrected_umis',
             'output_pairs', 'yield', 'reads_per_retained_umi', 'cell_set_intact', 'cellbender_anomaly', 'valid']
            ].sort_values(['species', 'sample', 'fraction'])
cc = curve.merge(my_curve_stored, on=['species', 'sample', 'fraction'], suffixes=('', '_stored'), validate='one_to_one')
assert len(cc) == len(curve) == len(my_curve_stored)
assert np.allclose(cc['yield'], cc['yield_stored'], rtol=1e-12) and (cc.valid == cc.valid_stored).all()
print('invalid curve points:', curve.loc[~curve.valid, ['sample', 'fraction']].to_dict('records'))
print(pairs[['species', 'nuclei', 'protoplast', 'reads_per_cell_full_n', 'reads_per_cell_full_p', 'matched_reads_per_cell',
             'yield_ratio_full', 'yield_ratio_matched', 'share_of_full_gap_from_depth']].round(3).to_string(index=False))


# %%
pairs.to_csv(TABLES / 'fig_matched_yield_pairs.tsv', sep='\t', index=False)
curve_out = annotate(curve)
curve_out.to_csv(TABLES / 'fig_yield_curve.tsv', sep='\t', index=False)
print(len(pairs), 'rows fig_matched_yield_pairs;', len(curve_out), 'rows fig_yield_curve;')


# %% [markdown]
# ## 3. Gene detection at equal UMI depth

# %%
mt_det = det[det.definition.eq('corrected') & det.species.eq('medicago') & det['set'].eq('all')].set_index('sample')
assert mt_det.matched_depth.eq(1500).all()
assert np.allclose(mt_det.loc[MT_N + MT_P, 'genes_at_matched'], core_i.loc[MT_N + MT_P, 'genes_at_matched'])
assert np.allclose(mt_det.loc[MT_N + MT_P, 'frac_eligible'], core_i.loc[MT_N + MT_P, 'frac_eligible'])
# 1,500 UMIs is the largest grid depth with >=50% eligible cells in every Medicago library.
mt_rule = mdepth[mdepth.species.eq('medicago') & mdepth['set'].eq('all') & mdepth.definition.eq('corrected')]
print(mt_rule.to_string(index=False))
mt_grid = rare[rare.definition.eq('corrected') & rare['sample'].isin(MT_N + MT_P)].pivot(index='depth', columns='sample', values='frac_eligible')
mt_supported = int(mt_grid.index[mt_grid.ge(.5).all(axis=1)].max())
print('largest Medicago grid depth with >=50% eligible in all six libraries:', mt_supported)
assert mt_supported == 1500
print((mt_det.loc[MT_N + MT_P, ['n_cells', 'n_cells_eligible', 'frac_eligible', 'genes_at_matched', 'genes_at_matched_median']]).round(3).to_string())

# %%
at_grid = elig.pivot(index='depth', columns='sample', values='frac_eligible').reindex(columns=AT_N + AT_P)
at_supported = int(at_grid.index[at_grid.notna().all(axis=1) & at_grid.ge(.5).all(axis=1)].max())
print(at_grid.loc[[400, 500, 1000, 1500]].round(3).to_string())
assert at_supported == 500 and elig.loc[elig.selected_primary_depth, 'depth'].eq(500).all()
at_sel = elig[elig.depth.eq(at_supported)].set_index('sample')
assert (at_sel.loc[AT_N + AT_P, 'n_cells_eligible'] == at_sel.loc[AT_N + AT_P, 'n_cells']).all()
assert np.allclose(at_sel.loc[AT_N + AT_P, 'genes_detected'], core_i.loc[AT_N + AT_P, 'genes_at_matched'])
all_eligible = bool(at_sel.loc[AT_N + AT_P, 'frac_eligible'].eq(1).all())
print('every retained cell eligible at 500 in all six libraries:', all_eligible)
print(at_sel.loc[AT_N + AT_P, ['n_cells', 'n_cells_eligible', 'frac_eligible', 'genes_detected', 'genes_detected_median_cells']].round(2).to_string())

# %% [markdown]
# ### Tables: matched-depth detection and full rarefaction curves

# %%
at_det_c1 = det[det.definition.eq('corrected') & det.species.eq('arabidopsis') & det['set'].eq('all')].set_index('sample')
rows = []
for s in MT_N + MT_P + AT_N + AT_P:
    if s in mt_det.index:
        r = mt_det.loc[s]
        src = 'C1_detection_matched.tsv'
        rows.append(dict(sample=s, matched_depth=r.matched_depth, n_cells=r.n_cells, n_cells_eligible=r.n_cells_eligible,
                         frac_eligible=r.frac_eligible, genes_mean=r.genes_at_matched, genes_median=r.genes_at_matched_median,
                         genes_q25=r.genes_at_matched_q25, genes_q75=r.genes_at_matched_q75, entropy_mean=r.entropy_at_matched,
                         quantile_source=src))
    else:
        g = at_sel.loc[s]
        c1 = at_det_c1.loc[s] if s in at_det_c1.index else None
        if c1 is not None:
            assert np.isclose(c1.genes_at_matched, g.genes_detected, rtol=1e-6), s
        rows.append(dict(sample=s, matched_depth=at_supported, n_cells=g.n_cells, n_cells_eligible=g.n_cells_eligible,
                         frac_eligible=g.frac_eligible, genes_mean=g.genes_detected,
                         genes_median=c1.genes_at_matched_median if c1 is not None else g.genes_detected_median_cells,
                         genes_q25=c1.genes_at_matched_q25 if c1 is not None else np.nan,
                         genes_q75=c1.genes_at_matched_q75 if c1 is not None else np.nan,
                         entropy_mean=g.entropy,
                         quantile_source='C1_detection_matched.tsv' if c1 is not None else
                         'C12_detection_eligibility.tsv (genes_detected_median_cells; no IQR for rep4)'))
dm = pd.DataFrame(rows).merge(core[['sample', 'species', 'modality', 'treatment']], on='sample')
dm = annotate(dm)[['species', 'modality', 'sample', 'sample_label', 'treatment', 'replicate', 'matched_depth', 'n_cells',
                   'n_cells_eligible', 'frac_eligible', 'genes_mean', 'genes_median', 'genes_q25', 'genes_q75',
                   'entropy_mean', 'quantile_source']]
dm[dm['sample'].isin(PRIMARY)].to_csv(TABLES / 'fig_detection_matched.tsv', sep='\t', index=False)
print(dm.round(2).to_string(index=False))

rc = rare[rare.species.isin(['medicago', 'arabidopsis'])].copy()
rc = annotate(rc)[['species', 'modality', 'sample', 'sample_label', 'treatment', 'replicate', 'primary_comparison',
                   'definition', 'depth', 'genes_detected', 'genes_detected_sd_reps', 'genes_detected_median_cells',
                   'entropy', 'n_cells', 'n_cells_eligible', 'frac_eligible']]
rc[rc['sample'].isin(PRIMARY)].to_csv(TABLES / 'fig_rarefaction_curves.tsv', sep='\t', index=False)
print(len(dm), 'rows fig_detection_matched;', len(rc), 'rows fig_rarefaction_curves')

