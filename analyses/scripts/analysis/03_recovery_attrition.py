# %% [markdown]
# # Read flow and QC losses
# Barcode ranks, candidate QC fractions and losses across processing stages.


# %%
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'scripts').is_dir())
import numpy as np
import pandas as pd

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 40)

DATA = ROOT / 'input/03_recovery_attrition'
QC = DATA / 'readflow_qc'
LIBRARIES = ROOT / 'input/02_nuclear_recovery'
TABLES = ROOT / 'output/tables/03_recovery_attrition'
TABLES.mkdir(parents=True, exist_ok=True)

pl = pd.read_csv(LIBRARIES / 'per_library.tsv', sep='\t')
fl = pd.read_csv(QC / 'qc_flags_among_candidates.tsv', sep='\t')
sl = pd.read_csv(QC / 'qc_stage_losses.tsv', sep='\t')
interp = pd.read_csv(QC / 'nuclear_library_interpretation.tsv', sep='\t')
audit_val = pd.read_csv(QC / 'validation.tsv', sep='\t')
rel = pd.read_csv(LIBRARIES / 'C12_reliability_perlib.tsv', sep='\t')
percell = pd.read_csv(DATA / 'qc_percell.tsv.gz', sep='\t')
curves = pd.read_csv(DATA / 'barcode_rank_curves.tsv.gz', sep='\t')
export_val = pd.read_csv(DATA / 'export_qc_percell_validation.tsv', sep='\t')

assert len(pl) == 36 and pl['sample'].is_unique
assert audit_val.drop(columns='sample').all().all(), 'audit validation flags'
ok = ['n_count_gate_matches', 'threshold_matches', 'raw_umis_match', 'corrected_umis_match', 'quantiles_match',
      'rank_ratio_matches', 'stage_losses_match']
assert len(export_val) == 36 and export_val[ok].all().all(), 'per-cell export does not reproduce the audit'
print('audit + export validation: all 36 libraries pass')


def label(s):
    return s.replace('GSE155304_sNucRNA-seq_', 'Farmer ').replace('GSE123013_WT_', 'Ryu ')


def group(r):
    sp = 'Medicago' if r.species == 'medicago' else 'Arabidopsis'
    if r.modality == 'protoplast':
        return f'{sp} protoplast'
    return f'{sp} nuclei retained' if r.chapter_floor_status == 'pass' else f'{sp} nuclei below floor'


meta = pl[['sample', 'species', 'modality', 'technology', 'arm', 'chapter_floor_status', 'recorded_automatic_gates_status',
           'main_analysis_usable']].merge(rel[['sample', 'treatment']], on='sample', how='left', validate='one_to_one')
meta['group'] = meta.apply(group, axis=1)
meta['sample_label'] = meta['sample'].map(label)
meta['replicate'] = pd.array([int(s.rsplit('rep', 1)[1]) if s.startswith('GSE') else pd.NA for s in meta['sample']], dtype='Int64')
meta['plot_order'] = np.arange(len(meta)) + 1   # audit order: species, modality, sample
META_COLS = ['species', 'modality', 'sample', 'sample_label', 'treatment', 'replicate', 'arm', 'group',
             'chapter_floor_status', 'plot_order']
print(meta.groupby('group').size().to_string())


# %%
P = pl.set_index('sample')
failed_nuc = P[(P.modality == 'nuclei') & (P.chapter_floor_status == 'fail')]
retained_nuc = P[(P.modality == 'nuclei') & (P.chapter_floor_status == 'pass')]
print('failed nuclear libraries  :', list(failed_nuc.index))
print('retained nuclear libraries:', list(retained_nuc.index))


# %% [markdown]
# ## 1. Mapping rates do not separate failed from retained nuclear libraries

# %%
fm = failed_nuc.frac_mapping * 100
rm = retained_nuc.frac_mapping * 100
print(fm.round(2).sort_values().to_string())
print()
print(rm.round(2).sort_values().to_string())
mt_fail = fm[failed_nuc.species == 'medicago']
mt_facs = rm[retained_nuc.species == 'medicago']


# %% [markdown]
# ### Tables: per-library mapping and read fate, with read and UMI partitions in long format

# %%
READ_FATES = [('unmapped', None), ('barcode_loss', 'frac_read_barcode_loss'), ('duplicate_reads', 'frac_read_duplicates'),
              ('deduplicated_molecules', 'frac_read_to_raw_umi')]
UMI_FATES = ['outside_cb_model', 'removed_within_cb_model', 'corrected_below_count_gate', 'corrected_later_qc_removed', 'retained']
wide_cols = ['total_reads', 'frac_mapping', 'frac_read_barcode_loss', 'frac_read_duplicates', 'frac_read_to_raw_umi',
             'sequencing_saturation', 'frac_read_to_retained_umi', 'raw_umis', 'corrected_umis', 'retained_corrected_umis',
             'whole_matrix_removed_fraction'] + [f'raw_umi_fraction_{k}' for k in UMI_FATES] + [
             'barcodes_total', 'barcodes_cb_analyzed', 'barcodes_count_gate', 'cells_automatic', 'cells_retained',
             'count_gate_threshold', 'raw_rank_1000_to_10000_ratio', 'raw_top1000_umi_share']
wide = meta[META_COLS].merge(pl[['sample'] + wide_cols], on='sample')
wide.to_csv(TABLES / 'fig_mapping_readfate_perlib.tsv', sep='\t', index=False)

rows = []
for _, r in pl.iterrows():
    for i, (name, col) in enumerate(READ_FATES):
        v = 1 - r.frac_mapping if col is None else r[col]
        rows.append(dict(sample=r['sample'], partition='input_reads', fate=name, fate_order=i + 1, fraction=v))
    for i, name in enumerate(UMI_FATES):
        rows.append(dict(sample=r['sample'], partition='raw_umis', fate=name, fate_order=i + 1,
                         fraction=r[f'raw_umi_fraction_{name}']))
stack = meta[META_COLS].merge(pd.DataFrame(rows), on='sample')
chk = stack.groupby(['sample', 'partition']).fraction.sum()
print('partition sums: min', round(chk.min(), 12), 'max', round(chk.max(), 12))
assert np.allclose(chk, 1, atol=1e-9)
stack.to_csv(TABLES / 'fig_readfate_stack.tsv', sep='\t', index=False)
print(len(wide), 'rows fig_mapping_readfate_perlib;', len(stack), 'rows fig_readfate_stack')


# %% [markdown]
# ## 2. Failure mode A — barcode-rank curves without a knee

# %%
shape = meta[META_COLS].merge(pl[['sample', 'raw_rank_1000_to_10000_ratio', 'raw_umi_fraction_outside_cb_model',
                                  'barcodes_count_gate', 'cells_retained']], on='sample')
print(shape.sort_values(['group', 'raw_rank_1000_to_10000_ratio'])[['group', 'sample', 'raw_rank_1000_to_10000_ratio',
      'raw_umi_fraction_outside_cb_model', 'barcodes_count_gate', 'cells_retained']].round(3).to_string(index=False))
A_LIBS = ['1_S_zone', '2_tip', '3_nod', '4_BAP', 'N0h']       # count gate already below the 500-cell floor
facs = ['nMOCK', 'nNOD', 'nBAP']
ra = P.loc[A_LIBS, 'raw_rank_1000_to_10000_ratio']
rf = P.loc[facs, 'raw_rank_1000_to_10000_ratio']
prot = P[P.modality == 'protoplast'].raw_rank_1000_to_10000_ratio
oa = P.loc[A_LIBS, 'raw_umi_fraction_outside_cb_model'] * 100
of = P.loc[facs, 'raw_umi_fraction_outside_cb_model'] * 100
supported = ra.max() < rf.min() and ra.max() < prot.min() and oa.min() > of.max()


# %% [markdown]
# ### Table: barcode-rank curves

# %%
cur = meta[META_COLS].merge(curves.drop(columns=['species', 'modality', 'group']), on='sample', validate='one_to_many')
# check the curve reproduces the audit's scalar summary where rank 1,000 falls on a drawn point
for s, g in cur.groupby('sample'):
    assert g['rank'].is_monotonic_increasing and g.raw_umis.dropna().is_monotonic_decreasing
    assert (g.count_gate_threshold == P.loc[s, 'count_gate_threshold']).all()
cur.to_csv(TABLES / 'fig_barcode_rank_curves.tsv', sep='\t', index=False)
print(len(cur), 'rows fig_barcode_rank_curves;', cur['sample'].nunique(), 'libraries')


# %% [markdown]
# ## 3. Failure mode B — organellar contamination removes candidates after the count gate

# %%
cg = percell[percell.stage == 'count_gate_candidate']
ret = percell[percell.stage == 'retained']
med = cg.groupby('sample')[['mito_fraction', 'chloroplast_fraction', 'rrna_fraction']].median()
for col, audit_col in [('mito_fraction', 'count_gate_mitochondrial_fraction_corrected_counts_median'),
                       ('rrna_fraction', 'count_gate_ribosomal_rna_fraction_corrected_counts_median'),
                       ('chloroplast_fraction', 'count_gate_chloroplast_fraction_corrected_counts_median')]:
    assert np.allclose(med[col], P.loc[med.index, audit_col], rtol=0, atol=1e-12), col
B_QUOTED = ['N24h', 'N48h', 'N72h']
mb = P.loc[B_QUOTED, 'count_gate_mitochondrial_fraction_corrected_counts_median'] * 100
mr = P.loc[facs, 'count_gate_mitochondrial_fraction_corrected_counts_median'] * 100
mr_ret = P.loc[facs, 'retained_mitochondrial_fraction_corrected_counts_median'] * 100
s5 = P.loc['5_s_zone_exp', 'count_gate_mitochondrial_fraction_corrected_counts_median'] * 100
print('count-gate median mito % (quoted):', mb.round(3).to_dict())
print('retained FACS nuclei, count-gate median mito %:', mr.round(3).to_dict(), '| retained-cell medians:', mr_ret.round(3).to_dict())
print(f'5_s_zone_exp: {s5:.3f}% (all {int(P.loc["5_s_zone_exp", "barcodes_count_gate"])} candidates removed)')
all_mt_nuc_fail_high = (P.loc[failed_nuc.index, 'count_gate_mitochondrial_fraction_corrected_counts_median'] > 0.02)
print('failed nuclear libraries with count-gate median mito above the 2 % cap:', list(all_mt_nuc_fail_high[all_mt_nuc_fail_high].index))


caps = P.groupby('modality')[['mitochondrial_absolute_cap', 'chloroplast_absolute_cap', 'ribosomal_rna_absolute_cap']].agg(['min', 'max'])
print(caps.to_string())
nc, pc = P[P.modality == 'nuclei'].mitochondrial_absolute_cap, P[P.modality == 'protoplast'].mitochondrial_absolute_cap


# %% [markdown]
# ## 4. Failure mode C — DTT preparations exceed the rRNA cap

# %%
DTT = ['N_DTT_12', 'N_DTT_50']
rr = P.loc[DTT, 'count_gate_ribosomal_rna_fraction_corrected_counts_median'] * 100
cap = P.loc[DTT, 'ribosomal_rna_absolute_cap']
RRNA = 'flag_absolute_ribosomal_rna_fraction_corrected_counts'
dtt_rows = []
for s in DTT:
    c = cg[cg['sample'] == s]
    at_caps = c[c.removal_stage == 'absolute_caps']
    only = at_caps[(at_caps[RRNA] == 1) & (at_caps.flag_absolute_mitochondrial_fraction_corrected_counts == 0)
                   & (at_caps.flag_absolute_chloroplast_fraction_corrected_counts == 0)]
    n_flag = int((c[RRNA] == 1).sum())
    assert n_flag == int(fl[(fl['sample'] == s) & (fl.flag == RRNA)].n_flagged_among_count_gate.iloc[0])
    dtt_rows.append(dict(sample=s, count_gate_candidates=len(c), rrna_cap_flagged=n_flag,
                         removed_at_absolute_caps=len(at_caps), removed_at_caps_with_rrna_flag=int((at_caps[RRNA] == 1).sum()),
                         removed_at_caps_rrna_only=len(only), removed_mad_miqc=int((c.removal_stage == 'MAD_or_miQC').sum()),
                         retained=int(c.retained.sum())))
dtt = pd.DataFrame(dtt_rows).set_index('sample')
print(dtt.to_string())
print('rRNA median %:', rr.round(3).to_dict(), ' cap:', cap.to_dict())


# %% [markdown]
# ## 5. Failure-mode annotation per nuclear library

# %%
def mode(r):
    if r.chapter_floor_status == 'pass':
        return 'retained'
    if r.barcodes_count_gate < 500:
        return 'no_droplet_separation'
    if r.count_gate_mitochondrial_fraction_corrected_counts_median > r.mitochondrial_absolute_cap:
        return 'organellar'
    if r.count_gate_ribosomal_rna_fraction_corrected_counts_median > r.ribosomal_rna_absolute_cap:
        return 'rrna'
    return 'adaptive_attrition'


nuc = pl[pl.modality == 'nuclei'].copy()
nuc['failure_mode'] = nuc.apply(mode, axis=1)
fmodes = meta[META_COLS].merge(nuc[['sample', 'failure_mode', 'frac_mapping', 'raw_rank_1000_to_10000_ratio',
                                    'raw_umi_fraction_outside_cb_model', 'barcodes_count_gate', 'cells_retained',
                                    'fraction_count_gate_retained', 'n_removed_MAD_or_miQC', 'n_removed_absolute_caps',
                                    'n_removed_regression', 'n_removed_other', 'n_removed_manual',
                                    'count_gate_mitochondrial_fraction_corrected_counts_median',
                                    'count_gate_chloroplast_fraction_corrected_counts_median',
                                    'count_gate_ribosomal_rna_fraction_corrected_counts_median',
                                    'retained_mitochondrial_fraction_corrected_counts_median',
                                    'mitochondrial_absolute_cap', 'chloroplast_absolute_cap', 'ribosomal_rna_absolute_cap']],
                               on='sample').merge(interp.rename(columns={'interpretation': 'audit_note'}), on='sample', how='left')
fmodes = fmodes.rename(columns={'count_gate_mitochondrial_fraction_corrected_counts_median': 'count_gate_mito_median',
                                'count_gate_chloroplast_fraction_corrected_counts_median': 'count_gate_chloroplast_median',
                                'count_gate_ribosomal_rna_fraction_corrected_counts_median': 'count_gate_rrna_median',
                                'retained_mitochondrial_fraction_corrected_counts_median': 'retained_mito_median',
                                'mitochondrial_absolute_cap': 'mito_cap', 'chloroplast_absolute_cap': 'chloroplast_cap',
                                'ribosomal_rna_absolute_cap': 'rrna_cap'})
assert fmodes.audit_note.notna().all() and len(fmodes) == 19
fmodes.to_csv(TABLES / 'fig_failure_modes.tsv', sep='\t', index=False)
print(fmodes[['sample', 'failure_mode', 'barcodes_count_gate', 'cells_retained', 'count_gate_mito_median',
              'count_gate_rrna_median', 'audit_note']].to_string(index=False))


# %% [markdown]
# ## 6. QC-stage losses and per-cell distributions (supplementary figure tables)

# %%
# sequential losses after the count gate, as fractions of count-gate candidates
STAGES = [('MAD_or_miQC', 'n_removed_MAD_or_miQC'), ('absolute_caps', 'n_removed_absolute_caps'),
          ('regression', 'n_removed_regression'), ('other', 'n_removed_other'), ('manual', 'n_removed_manual'),
          ('retained', 'cells_retained')]
rows = []
for _, r in pl.iterrows():
    for i, (name, col) in enumerate(STAGES):
        rows.append(dict(sample=r['sample'], stage=name, stage_order=i + 1, n_cells=int(r[col]),
                         count_gate_candidates=int(r.barcodes_count_gate),
                         fraction_of_candidates=r[col] / r.barcodes_count_gate if r.barcodes_count_gate else np.nan))
losses = meta[META_COLS].merge(pd.DataFrame(rows), on='sample')
assert (losses.groupby('sample').n_cells.sum() == P.loc[losses.groupby('sample').n_cells.sum().index, 'barcodes_count_gate']).all()
chk = losses.merge(sl[['sample', 'stage', 'removed', 'applied_flags']], on=['sample', 'stage'], how='left')
assert (chk.dropna(subset=['removed']).eval('n_cells == removed')).all()
losses = chk.rename(columns={'removed': 'audit_removed'}).drop(columns='audit_removed')
# the same partition from the per-cell table
pc_counts = cg.groupby(['sample', 'removal_stage']).size()
for _, r in losses.iterrows():
    assert pc_counts.get((r['sample'], r.stage), 0) == r.n_cells, (r['sample'], r.stage)
losses.to_csv(TABLES / 'fig_qc_stage_losses.tsv', sep='\t', index=False)

flags = meta[META_COLS].merge(fl, on='sample')
flags['fraction_of_candidates'] = flags.n_flagged_among_count_gate / flags.count_gate_n.replace(0, np.nan)
flags['flag_family'] = flags.flag.str.extract(r'^flag_(bcrank|mad|miqc|absolute|regression)')[0]
def feature(flag):
    for f in ['mitochondrial', 'chloroplast', 'ribosomal_rna', 'cell_complexity', 'n_genes', 'total_counts', 'regression']:
        if f in flag:
            return f
    return flag


flags['feature'] = flags.flag.map(feature)
flags.to_csv(TABLES / 'fig_qc_flags_among_candidates.tsv', sep='\t', index=False)
print(len(losses), 'rows fig_qc_stage_losses;', len(flags), 'rows fig_qc_flags_among_candidates')


# %% [markdown]
# ### Per-cell distributions (violin / ECDF panels)

# %%
keep = ['barcode', 'stage', 'retained', 'removal_stage', 'raw_umis', 'corrected_umis', 'n_genes', 'removed_fraction',
        'mito_fraction', 'chloroplast_fraction', 'rrna_fraction', 'rspliced_fraction',
        'flag_absolute_mitochondrial_fraction_corrected_counts', 'flag_absolute_chloroplast_fraction_corrected_counts',
        'flag_absolute_ribosomal_rna_fraction_corrected_counts', 'flag_regression', 'flags_set']
mad = [c for c in percell.columns if c.startswith(('flag_mad_', 'flag_miqc_'))]
pc2 = percell[['sample'] + keep].copy()
pc2['flag_mad_or_miqc_any'] = percell[mad].fillna(0).astype(int).max(axis=1)   # same row order as percell
pco = meta[META_COLS].merge(pc2, on='sample', validate='one_to_many')
pco['mito_cap'] = pco['sample'].map(P.mitochondrial_absolute_cap)
pco['rrna_cap'] = pco['sample'].map(P.ribosomal_rna_absolute_cap)
pco['chloroplast_cap'] = pco['sample'].map(P.chloroplast_absolute_cap)
pco.to_csv(TABLES / 'fig_qc_percell.tsv.gz', sep='\t', index=False)

METRICS = ['mito_fraction', 'chloroplast_fraction', 'rrna_fraction', 'corrected_umis', 'n_genes', 'rspliced_fraction',
           'removed_fraction']
qs = pco.groupby(['sample', 'stage'])[METRICS].quantile([.1, .25, .5, .75, .9]).rename_axis(['sample', 'stage', 'q']).reset_index()
summ = qs.melt(id_vars=['sample', 'stage', 'q'], var_name='metric').pivot_table(index=['sample', 'stage', 'metric'], columns='q', values='value')
summ.columns = ['p10', 'p25', 'median', 'p75', 'p90']
summ = summ.join(pco.groupby(['sample', 'stage'])[METRICS].mean().melt(ignore_index=False, var_name='metric', value_name='mean')
                 .reset_index().set_index(['sample', 'stage', 'metric']))
summ = summ.reset_index().merge(pco.groupby(['sample', 'stage']).size().rename('n_cells').reset_index(), on=['sample', 'stage'])
summ = meta[META_COLS].merge(summ, on='sample')
# reproduces the audit per-library quantiles
for mcol, acol in [('mito_fraction', 'mitochondrial_fraction_corrected_counts'), ('rrna_fraction', 'ribosomal_rna_fraction_corrected_counts'),
                   ('chloroplast_fraction', 'chloroplast_fraction_corrected_counts'), ('corrected_umis', 'corrected_counts_total_counts')]:
    for st, pre in [('count_gate_candidate', 'count_gate'), ('retained', 'retained')]:
        sub = summ[(summ.metric == mcol) & (summ.stage == st)].set_index('sample')
        for q in ['p10', 'median', 'p90']:
            assert np.allclose(sub[q], P.loc[sub.index, f'{pre}_{acol}_{q}'], rtol=1e-12, atol=1e-12), (mcol, st, q)
summ.to_csv(TABLES / 'fig_qc_percell_summary.tsv', sep='\t', index=False)
print(len(pco), 'rows fig_qc_percell (gz);', len(summ), 'rows fig_qc_percell_summary')
print(summ[(summ.metric == 'mito_fraction') & summ['sample'].isin(B_QUOTED + facs + ['5_s_zone_exp'])]
      [['sample', 'stage', 'n_cells', 'p10', 'median', 'p90']].round(4).to_string(index=False))

