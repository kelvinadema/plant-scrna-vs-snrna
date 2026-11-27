# %% [markdown]
# # Cell-type composition
# Supplied labels: six retained libraries per species; denominator includes unassigned cells.


# %%
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'analyses/scripts').is_dir())
import re
import numpy as np
import pandas as pd
import yaml

DATA = ROOT / 'input/05_composition'
LABELS = DATA / 'labels'
CONFIG = DATA / 'sample_metadata.json'
TABLES = ROOT / 'output/tables/05_composition'
TABLES.mkdir(parents=True, exist_ok=True)
PRIMARY = {
    'medicago': ['nMOCK', 'nNOD', 'nBAP', 'pMock', 'pNOD', 'pBAP'],
    'arabidopsis': ['GSE155304_sNucRNA-seq_rep2', 'GSE155304_sNucRNA-seq_rep4', 'GSE155304_sNucRNA-seq_rep5',
                    'GSE123013_WT_rep1', 'GSE123013_WT_rep2', 'GSE123013_WT_rep3'],
}


# %% [markdown]
# ## Recompute composition from the label files

# %%
labels = pd.concat([pd.read_csv(LABELS / f'C2_labels_{s}.tsv', sep='\t') for s in PRIMARY], ignore_index=True)
labels['species'] = labels.species_x
labels['cell_type'] = labels.label.map(lambda x: 'pericycle' if str(x).strip().lower() in {'pericycle', 'pericycle cell'}
                                       else str(x).strip())
meta = (labels.groupby('sample')
        .agg(species=('species', 'first'), modality=('modality', 'first'), treatment=('treatment', 'first'),
             study=('study_origin', 'first'), technology=('technology', 'first'), prep_variant=('prep_variant', 'first'),
             n_total=('cell', 'size'))
        .reset_index())
assert labels.groupby('sample')[['species', 'modality', 'treatment']].nunique().eq(1).all().all()


def replicate(sample):
    m = re.search(r'rep(\d+)$', sample)
    return int(m.group(1)) if m else 1


meta['replicate'] = meta['sample'].map(replicate)
meta['primary_comparison'] = meta['sample'].isin(sum(PRIMARY.values(), []))

# Library metadata must agree with the chapter config's sample sheet.
cfg = yaml.safe_load(CONFIG.read_text())['samples']
for r in meta.itertuples():
    c = cfg[r.sample]
    assert (c['species'], c['modality']) == (r.species, r.modality), r.sample

# %%
counts = labels.groupby(['sample', 'cell_type']).size().rename('n_cells').reset_index()
# Complete the grid per species so that a type not called in a library is an explicit zero.
grid = []
for species, m in meta.groupby('species'):
    types = sorted(counts[counts['sample'].isin(m['sample'])].cell_type.unique())
    idx = pd.MultiIndex.from_product([m['sample'], types], names=['sample', 'cell_type'])
    grid.append(counts.set_index(['sample', 'cell_type']).reindex(idx, fill_value=0).reset_index())
comp = pd.concat(grid, ignore_index=True).merge(meta, on='sample')
comp['proportion'] = comp.n_cells / comp.n_total
comp['called'] = comp.n_cells.gt(0)
comp['denominator'] = 'all_retained_cells_including_Unassigned'
assert np.allclose(comp.groupby('sample').proportion.sum(), 1)
comp = comp[['species', 'modality', 'sample', 'study', 'technology', 'prep_variant', 'treatment', 'replicate',
             'primary_comparison', 'cell_type', 'n_cells', 'n_total', 'proportion', 'called', 'denominator']]
comp = comp.sort_values(['species', 'modality', 'sample', 'cell_type']).reset_index(drop=True)
comp.to_csv(TABLES / 'fig_composition.tsv', sep='\t', index=False)
comp.pivot_table(index='cell_type', columns='sample', values='proportion').round(4)


# %% [markdown]
# ## Supplied UMAP embedding and library totals

# %%
import shutil
expected = {'medicago': 24161, 'arabidopsis': 18859}
assert labels.groupby('species').size().to_dict() == expected
assert meta.primary_comparison.all() and len(meta) == 12
assert not labels.duplicated(['sample', 'cell']).any()
for name in ('fig_umap_joint.tsv', 'fig_umap_palette.tsv'):
    shutil.copyfile(DATA / 'umap' / name, TABLES / name)
umap = pd.read_csv(TABLES / 'fig_umap_joint.tsv', sep='\t')
assert set(umap['sample']) == set(meta['sample']) and len(umap) == len(labels)
assert not umap.duplicated(['sample','cell']).any()
a = labels[['sample', 'cell', 'cell_type']].merge(
    umap[['sample', 'cell', 'cell_type']], on=['sample', 'cell'],
    suffixes=('_labels', '_umap'), validate='one_to_one')
assert len(a) == len(labels)
assert a.cell_type_labels.equals(a.cell_type_umap.replace({'pericycle cell':'pericycle'}))
meta.to_csv(TABLES / 'fig_library_totals.tsv', sep='\t', index=False)
comp.assign(percent=comp.proportion * 100).to_csv(TABLES / 'fig_composition.tsv', sep='\t', index=False)
print(meta[['species','sample','n_total']].to_string(index=False))
print('Validated V2 labels, all-cell composition denominators and matching UMAP coordinates.')
