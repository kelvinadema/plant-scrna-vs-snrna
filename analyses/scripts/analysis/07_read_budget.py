# %% [markdown]
# # Sequencing requirements
# Corrected counts, fixed cells and measured AUROC 0.70 crossings; nNOD is not assessable.


# %%
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'scripts').is_dir())
import numpy as np
import pandas as pd
D = ROOT / 'input/07_read_budget'
T = ROOT / 'output/tables/07_read_budget'
T.mkdir(parents=True, exist_ok=True)
curve = pd.read_csv(D/'read_budget/C14_fixedcell_curve.tsv', sep='\t')
exchange = pd.read_csv(D/'read_budget/C14_exchange_rate.tsv', sep='\t')
curve = curve[curve.depth_mode.eq('native')].copy()
exchange = exchange[exchange.layer.str.startswith('corrected_counts') & exchange.target_auc.eq(0.70)].copy()
curve['analysis_variant'] = 'V2_25_25'
curve['layer'] = curve.depth_mode.map({'native':'corrected', 'native_raw':'raw', 'rarefied_500':'corrected_rarefied_500'})
assert curve.layer.notna().all(), set(curve.depth_mode)
curve = curve.rename(columns={'reads_per_fixed_cell':'reads_per_cell', 'heldout_gene_auc':'auroc',
                             'heldout_gene_auc_split_sd':'auroc_split_sd', 'heldout_program_auc':'programme_auroc',
                             'heldout_program_auc_split_sd':'programme_auroc_split_sd', 'name':'object_name',
                             'output_pairs':'output_read_pairs'})
curve['used_for_crossing'] = curve.supported & (curve.layer.eq('raw') | (curve.cell_set_intact & ~curve.cellbender_anomaly))
assert not curve.duplicated(['sample','depth_mode','fraction','seed_no']).any()
assert 'nNOD' not in set(curve['sample'])
curve.to_csv(T/'fig_auroc_vs_reads.tsv', sep='\t', index=False)


# %% [markdown]
# ## Recompute AUROC 0.70 crossings
# Interpolate linearly in log10(reads); preserve missing crossings and bounds.

# %%
def crossing(reads, auc, target):
    reads, auc = np.asarray(reads, float), np.asarray(auc, float)
    order = np.argsort(reads)
    reads, auc = reads[order], auc[order]
    valid = np.isfinite(auc)
    reads, auc = reads[valid], auc[valid]
    if len(auc) == 0:
        return np.nan, 'no curve'
    if auc.max() < target:
        return np.nan, f'not reached at {reads[-1]:.0f} reads/cell (max AUROC {auc.max():.3f})'
    i = int(np.argmax(auc >= target))
    if i == 0:
        return reads[0], f'reached at or below the lowest fraction ({reads[0]:.0f} reads/cell); upper bound on required reads'
    x0, x1 = np.log10(reads[i - 1]), np.log10(reads[i])
    y0, y1 = auc[i - 1], auc[i]
    return float(10 ** (x0 + (target - y0) * (x1 - x0) / (y1 - y0))), 'interpolated'

valid_curve = curve[curve.used_for_crossing]
for i, row in exchange.iterrows():
    required, statuses, alternate = {}, {}, {}
    for suffix, sample in [('n', row.nuclei_library), ('p', row.protoplast_library)]:
        points = valid_curve[valid_curve['sample'].eq(sample)]
        main = points[~points.technical_draw]
        required[suffix], statuses[suffix] = crossing(main.reads_per_cell, main.auroc, row.target_auc)
        second = points[points.technical_draw & points.fraction.eq(0.25)]
        if len(second):
            alt = pd.concat([main[~main.fraction.eq(0.25)], second])
            alternate[suffix] = crossing(alt.reads_per_cell, alt.auroc, row.target_auc)[0]
        assert np.isclose(required[suffix], row['reads_' + suffix], equal_nan=True)
        assert statuses[suffix] == row['status_' + suffix]
        exchange.loc[i, 'reads_' + suffix] = required[suffix]
    rn, rp = required['n'], required['p']
    ratio = rn / rp if np.isfinite(rn) and np.isfinite(rp) else np.nan
    censored_n = statuses['n'].startswith('reached at or below')
    censored_p = statuses['p'].startswith('reached at or below')
    if censored_n and censored_p:
        ratio = np.nan
    assert np.isclose(ratio, row.ratio_measured, equal_nan=True)
    exchange.loc[i, 'ratio_measured'] = ratio
    if alternate:
        an, ap = alternate.get('n', rn), alternate.get('p', rp)
        alt_ratio = an / ap if np.isfinite(an) and np.isfinite(ap) else np.nan
        assert np.isclose(alt_ratio, row.ratio_measured_seed2, equal_nan=True)
        exchange.loc[i, 'ratio_measured_seed2'] = alt_ratio
print('Recomputed and checked', len(exchange), 'threshold comparisons.')

assert len(exchange) == 11, 'Expected two Medicago and nine Arabidopsis pairs'
exchange['analysis_variant'] = 'V2_25_25'
exchange['layer_source'] = exchange.layer
exchange['layer'] = exchange.layer_source.map(lambda x: 'raw' if x.startswith('raw') else 'corrected')
exchange['assessable_both'] = exchange.reached_n & exchange.reached_p
exchange['ratio_nuclei_over_protoplast'] = exchange.ratio_measured
# A finite measured ratio requires both threshold crossings and positive measured read requirements.
finite = np.isfinite(exchange.ratio_measured)
assert exchange.loc[finite, 'assessable_both'].all()
assert (exchange.loc[finite, ['reads_n','reads_p']] > 0).all().all()
assert np.allclose(exchange.loc[finite,'ratio_measured'], exchange.loc[finite,'reads_n']/exchange.loc[finite,'reads_p'])
exchange.to_csv(T/'fig_exchange_ratio.tsv', sep='\t', index=False)
# Long crossing table keeps not-reached and lowest-depth upper-bound statuses visible.
long=[]
for suffix, modality in [('n','nuclei'),('p','protoplast')]:
    x = exchange[['species','pair','treatment','layer','primary','target','target_auc','analysis_variant']].copy()
    x['sample'] = exchange['nuclei_library' if suffix=='n' else 'protoplast_library']
    x['modality'] = modality
    x['reads_per_cell'] = exchange['reads_'+suffix]
    x['reached'] = exchange['reached_'+suffix]
    x['status'] = exchange['status_'+suffix]
    x['auc_full_depth'] = exchange['auc_full_depth_'+suffix]
    x['reads_full_depth'] = exchange['reads_full_depth_'+suffix]
    long.append(x)
pd.concat(long, ignore_index=True).to_csv(T/'fig_target_crossings.tsv', sep='\t', index=False)
for name in ('C14_support.tsv', 'C14_medicago_pair_support.tsv', 'C14_fixedcell_perlib.tsv', 'C14_subsample_manifest.tsv'):
    p = D/'read_budget'/name
    if p.exists():
        pd.read_csv(p, sep='\t').assign(analysis_variant='V2_25_25').to_csv(T/name.replace('C14_', 'support_'), sep='\t', index=False)
print('Exported',len(curve),'curve rows and',len(exchange),'exchange rows; missing crossings preserved.')
print(exchange.loc[exchange.primary, ['species','pair','target_auc','ratio_measured','status_n','status_p']].to_string(index=False))

# %% [markdown]
# ## Full current Table S6
# Preserve the supplied supplement's raw/corrected counts, both targets and
# second-seed estimates. Figure 7 above remains corrected counts at AUROC 0.70.

# %%
supplement = pd.read_csv(D/'read_budget/C14_exchange_rate.tsv', sep='\t')
assert len(supplement) == 44
supplement['analysis_variant'] = 'V2_25_25'
supplement['layer_source'] = supplement.layer
supplement['layer'] = supplement.layer_source.map(lambda x: 'raw' if x.startswith('raw') else 'corrected')
supplement.to_csv(T/'supplementary_exchange_ratio.tsv', sep='\t', index=False)
supplement_long = []
for suffix, modality in [('n', 'nuclei'), ('p', 'protoplast')]:
    x = supplement[['species','pair','treatment','layer','primary','target','target_auc','analysis_variant']].copy()
    x['sample'] = supplement['nuclei_library' if suffix == 'n' else 'protoplast_library']
    x['modality'] = modality
    x['reads_per_cell'] = supplement['reads_'+suffix]
    x['reached'] = supplement['reached_'+suffix]
    x['status'] = supplement['status_'+suffix]
    x['auc_full_depth'] = supplement['auc_full_depth_'+suffix]
    x['reads_full_depth'] = supplement['reads_full_depth_'+suffix]
    supplement_long.append(x)
pd.concat(supplement_long, ignore_index=True).to_csv(T/'supplementary_target_crossings.tsv', sep='\t', index=False)
