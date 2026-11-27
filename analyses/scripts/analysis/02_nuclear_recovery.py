# %% [markdown]
# # Nuclear recovery
# Chapter floor: 500 retained cells; 10x libraries only.


# %%
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'analyses/scripts').is_dir())
import numpy as np
import pandas as pd

DATA = ROOT / 'input/02_nuclear_recovery'
TABLES = ROOT / 'output/tables/02_nuclear_recovery'
TABLES.mkdir(parents=True, exist_ok=True)


lib = pd.read_csv(DATA / 'C12_reliability_perlib.tsv', sep='\t')
lib = lib[lib.comparison_scope.eq('chapter_reliability')].copy()
print(len(lib), 'libraries')
lib.groupby(['species', 'modality', 'technology']).size()


# %% [markdown]
# ## Tidy table: one row per library

# %%
def prep(row):
    arm = row.arm.lower()
    if 'facs' in arm:
        return 'facs_sorted'
    if 'unsorted' in arm:
        return 'unsorted'
    if 'dtt' in arm:
        return 'tip_unfixed_or_dtt'
    return 'standard'


tidy = pd.DataFrame({
    'sample': lib['sample'], 'species': lib.species, 'modality': lib.modality, 'study': lib.study,
    'arm': lib.arm, 'technology': lib.technology, 'prep': lib.apply(prep, axis=1),
    'treatment': lib.treatment.fillna('untreated'), 'tissue': lib.tissue_from_summary,
    'comparison_scope': lib.comparison_scope,
    'in_direct_comparison': lib.comparison_scope.eq('chapter_reliability'),
    'cells_barcode_rank': lib.cells_barcode_rank, 'cells_automatic': lib.cells_automatic,
    'cells_final_manual': lib.cells_final, 'retained_cells': lib.retained_cells,
    'retained_stage': lib.retained_stage, 'chapter_floor': lib.chapter_floor,
    'floor_status': lib.chapter_floor_status, 'retained': lib.chapter_floor_status.eq('pass'),
    'automatic_gate_status': lib.recorded_automatic_gates_status,
    'main_analysis_selected': lib.main_analysis_selected.astype(bool),
    'total_reads': lib.mapper_read_count, 'mapped_fraction': lib.mapper_mapped_fraction,
    'reads_per_retained_cell': lib.reads_per_retained_cell,
    'median_umis_retained': lib.median_umis_retained,
})
# The floor status must follow from the retained count.
assert (tidy.retained == tidy.retained_cells.ge(tidy.chapter_floor)).all()
tidy = tidy.sort_values(['species', 'modality', 'arm', 'sample']).reset_index(drop=True)
tidy.to_csv(TABLES / 'fig_library_recovery.tsv', sep='\t', index=False)

long = tidy.melt(id_vars=['sample', 'species', 'modality', 'study', 'arm', 'technology', 'prep', 'treatment',
                          'in_direct_comparison', 'retained', 'main_analysis_selected'],
                 value_vars=['cells_barcode_rank', 'cells_automatic', 'cells_final_manual', 'retained_cells'],
                 var_name='stage', value_name='n_cells')
long['stage'] = long.stage.str.replace('cells_', '', regex=False).replace({'retained_cells': 'retained'})
long.to_csv(TABLES / 'fig_library_stage_cells.tsv', sep='\t', index=False)
tidy[['sample', 'species', 'modality', 'arm', 'technology', 'cells_barcode_rank', 'cells_automatic',
      'retained_cells', 'floor_status', 'automatic_gate_status', 'main_analysis_selected']]


# %% [markdown]
# ### Cross-check against the read-flow audit

# %%
audit = pd.read_csv(DATA / 'per_library.tsv', sep='\t').set_index('sample')
x = tidy.set_index('sample').join(audit[['cells_retained']], how='inner')
print(len(x), '10x libraries in both tables; disagreements:', int((x.retained_cells != x.cells_retained).sum()))
assert (x.retained_cells == x.cells_retained).all()


# %%
mt_n = tidy[tidy.species.eq('medicago') & tidy.modality.eq('nuclei')]
n_fail = int((~mt_n.retained).sum())
src = 'tables/fig_library_recovery.tsv'
mt_n[['sample', 'prep', 'retained_cells', 'floor_status', 'automatic_gate_status']]


# %%
kept = mt_n[mt_n.retained]
lo, hi = int(kept.retained_cells.min()), int(kept.retained_cells.max())
only_facs = set(kept.prep) == {'facs_sorted'} and len(kept) == 3
kept[['sample', 'retained_cells', 'reads_per_retained_cell']]


# %%
far = tidy[tidy.study.eq('GSE155304')].sort_values('sample')
far_kept = far[far.retained]
vals = [int(v) for v in far_kept.retained_cells]
far[['sample', 'technology', 'cells_barcode_rank', 'cells_automatic', 'retained_cells', 'floor_status', 'automatic_gate_status']]


# %%
p10 = tidy[tidy.modality.eq('protoplast') & tidy.technology.isin(['10xv2', '10xv3'])]
n_pass = int((p10.retained & p10.automatic_gate_status.eq('pass')).sum())


# %% [markdown]
# ## Per-arm summary (and check against the producer's summary table)

# %%
arm = (tidy.groupby(['species', 'modality', 'arm', 'technology'], as_index=False)
       .agg(n_libraries=('sample', 'size'), n_retained=('retained', 'sum'),
            n_automatic_gate_pass=('automatic_gate_status', lambda s: int(s.eq('pass').sum())),
            n_main_analysis=('main_analysis_selected', 'sum'),
            retained_min=('retained_cells', 'min'), retained_max=('retained_cells', 'max')))
arm.to_csv(TABLES / 'recovery_by_arm.tsv', sep='\t', index=False)
ref = pd.read_csv(DATA / 'C12_reliability_summary.tsv', sep='\t').set_index('arm')
ref = ref.loc[arm.arm.unique()]
a = arm.groupby('arm').agg(n_libraries=('n_libraries', 'sum'), n_floor_pass=('n_retained', 'sum'),
                           retained_min=('retained_min', 'min'), retained_max=('retained_max', 'max'))
assert set(a.index) == set(ref.index)
assert (a.values.astype(float) == ref.loc[a.index, a.columns].values.astype(float)).all(), \
    'per-arm summary disagrees with C12_reliability_summary.tsv'

