# %% [markdown]
# # NF and cytokinin response
# 903 cells × 1,500 UMIs per library, ten technical draws. Remove genes induced by either preparation.


# %%
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import json, zlib, hashlib
from pathlib import Path
START = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
ROOT = next(p for p in (START, *START.parents) if (p / 'input').is_dir() and (p / 'analyses/scripts').is_dir())
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, fisher_exact

D = ROOT / 'input/08_nf_response'
T = ROOT / 'output/tables/08_nf_response'
T.mkdir(parents=True, exist_ok=True)
SPECIES = 'Medicago truncatula'
print('numpy', np.__version__, '| pandas', pd.__version__)


# %% [markdown]
# ## Metrics

# %%
def auc(sc, pos):
    n = int(pos.sum())
    m = len(pos) - n
    return float((rankdata(sc)[pos].sum() - n * (n + 1) / 2) / (n * m)) if n >= 3 and m >= 3 else np.nan

def hm(A, n1, n2):
    if not np.isfinite(A) or n1 < 3:
        return np.nan
    Q1 = A / (2 - A)
    Q2 = 2 * A * A / (1 + A)
    return float(np.sqrt((A * (1 - A) + (n1 - 1) * (Q1 - A * A) + (n2 - 1) * (Q2 - A * A)) / (n1 * n2)))

cpm = lambda v: v / v.sum() * 1e6              # run_C15_A17_reference.py
cpm0 = lambda v: v / max(v.sum(), 1) * 1e6      # layer / refree producers (guards an empty vector)
lfc = lambda a, b: np.log2((cpm(a) + 1) / (cpm(b) + 1))


# %% [markdown]
# ## Inputs: reference, gene sets, exclusion, pseudobulks

# %%
REF = D / 'reference_ch2_A17_3h'
SIG = {d: set(x.strip() for x in (REF / f'NF_response_A17_{d}.txt').read_text().split() if x.strip()) for d in ['up', 'down']}
lfc_ref = pd.read_csv(REF / 'NF_vs_mock_3h_A17.tsv', sep='\t').set_index('gene_id')
assert len(SIG['up']) == 311 and len(SIG['down']) == 42, {k: len(v) for k, v in SIG.items()}
meta = json.loads((REF / 'NF_response_A17_up.meta.json').read_text())
print(meta['model'], meta['contrast'], 'padj <', meta['padj_gate'], '|LFC| >=', meta['lfc_gate'], '| genes in LFC table', len(lfc_ref))

prep = pd.read_csv(D / 'C9_isolation_genesets.tsv', sep='\t')
st = {k: set(g.gene) for k, g in prep.groupby('set')}
EXC_SETS = {'either_up': st['protoplast_up'] | st['nuclei_up']}
print({k: len(v) for k, v in EXC_SETS.items()})

W = D / 'whole_library_pseudobulks'
LIBS = ['nMOCK', 'nNOD', 'pMock', 'pNOD']
wb = {s: pd.read_csv(W / f'{s}.tsv.gz', sep='\t', index_col='gene') for s in LIBS}
g = wb['nMOCK'].index
assert all(x.index.equals(g) for x in wb.values())
SCOPES = ['all_cells', 'without_root_hair']
MODS = {'nuclei': ('nNOD', 'nMOCK'), 'protoplast': ('pNOD', 'pMock')}
col = lambda s, d, scope: wb[s][('native__' if d == 0 else f'd{d}__') + scope].to_numpy(float)

R = lfc_ref.reindex(g)
bulk = R.log2FoldChange.to_numpy()
padj = R.padj.to_numpy()
up = np.asarray(g.isin(SIG['up']))
down = np.asarray(g.isin(SIG['down']))
EXC = {k: np.asarray(g.isin(v)) for k, v in EXC_SETS.items()}
print('genes', len(g), '| A17 up in gene index', up.sum(), '| down', down.sum())


# %% [markdown]
# ### Budget check: 903 cells × 1,500 UMIs

# %%
elig = pd.read_csv(D / 'C15_eligibility.tsv', sep='\t')
n_budget = int(elig[elig['sample'].isin(LIBS) & ~elig.label.eq('root hair cell')].groupby('sample').n_eligible.sum().min())
binding = elig[elig['sample'].isin(LIBS) & ~elig.label.eq('root hair cell')].groupby('sample').n_eligible.sum().idxmin()
DEPTH = 1500
man = {s: json.loads((W / f'{s}_manifest.json').read_text()) for s in LIBS}
for s in LIBS:
    assert man[s]['n_matched'] == n_budget and man[s]['depth'] == DEPTH
    for scope in SCOPES:
        for d in range(1, 11):
            assert col(s, d, scope).sum() == n_budget * DEPTH, (s, scope, d)
print('budget', n_budget, 'cells x', DEPTH, 'UMIs; binding library', binding)
print({s: (man[s]['job_id'], man[s]['script_sha256'][:8]) for s in LIBS})

comp = pd.concat([pd.read_csv(W / f'{s}_composition.tsv', sep='\t') for s in LIBS], ignore_index=True)
comp.insert(0, 'species', SPECIES)
comp.insert(2, 'modality', comp['sample'].str[0].map({'n': 'nuclei', 'p': 'protoplast'}))
comp.insert(3, 'treatment', comp['sample'].str[1:].str.upper().map({'MOCK': 'mock', 'NOD': 'NF'}))
comp.to_csv(T / 'fig_budget_composition.tsv', sep='\t', index=False)
rh = (comp[(comp.arm == 'matched_cells_UMIs') & (comp.scope == 'all_cells') & (comp.label == 'root hair cell')]
      .groupby('sample').fraction.mean())
print('root-hair share of the matched cells (all_cells scope, mean over draws):', rh.round(3).to_dict())


# %% [markdown]
# ## A17 recall: AUROC, direction, median LFC and absolute error (whole-library arm)

# %%
def a17_metrics(l, uni):
    kn = uni & np.isfinite(bulk)
    sig = uni & (up | down)
    return dict(
        auroc_up=auc(l[uni], up[uni]), auroc_down=auc(-l[uni], down[uni]),
        rho_transcriptome=spearmanr(l[kn], bulk[kn]).statistic,
        rho_sig=spearmanr(l[sig], bulk[sig]).statistic if sig.sum() > 10 else np.nan,
        dir_agree=float(np.mean(np.sign(l[sig]) == np.sign(bulk[sig]))) if sig.sum() else np.nan,
        med_up=float(np.median(l[uni & up])) if (uni & up).sum() else np.nan,
        med_down=float(np.median(l[uni & down])) if (uni & down).sum() > 2 else np.nan,
        abs_err=float(np.median(np.abs(l[sig] - bulk[sig]))) if sig.sum() else np.nan,
        n_up=int((uni & up).sum()), n_down=int((uni & down).sum()),
        n_universe=int(uni.sum()),
        n_transcriptome=int(kn.sum()))

ARMS = [('matched_cells_UMIs', list(range(1, 11)))]
rows = []
for scope in SCOPES:
    for arm, draws in ARMS:
        for exc, ex in [('either_up', EXC['either_up'])]:
            for d in draws:
                V = {s: col(s, d, scope) for s in LIBS}
                common = np.logical_and.reduce([(V[a] + V[b]) >= 10 for a, b in MODS.values()]) & ~ex
                for mod, (a, b) in MODS.items():
                    l = lfc(V[a], V[b])
                    for universe, uni in [('per_modality', ((V[a] + V[b]) >= 10) & ~ex)]:
                        rows.append(dict(species=SPECIES, arm=arm, scope=scope, exclusion=exc, universe=universe,
                                         modality=mod, sample_treated=a, sample_mock=b, treatment='NF', draw=d,
                                         **a17_metrics(l, uni)))
recall_pd = pd.DataFrame(rows)
KEYS = ['species', 'arm', 'scope', 'exclusion', 'universe', 'modality', 'sample_treated', 'sample_mock', 'treatment']
num = [c for c in recall_pd.columns if c not in KEYS + ['draw']]
recall_sm = recall_pd.groupby(KEYS, sort=False)[num].mean().reset_index()
sds = recall_pd.groupby(KEYS, sort=False).agg(auroc_up_draw_sd=('auroc_up', 'std'), rho_draw_sd=('rho_transcriptome', 'std'),
                                               dir_agree_draw_sd=('dir_agree', 'std'), n_draws=('draw', 'nunique')).reset_index()
recall_sm = recall_sm.merge(sds, on=KEYS, validate='one_to_one')
for a_, n1, se in [('auroc_up', 'n_up', 'auroc_up_gene_se'), ('auroc_down', 'n_down', 'auroc_down_gene_se')]:
    recall_sm[se] = [hm(A, n, N - n) for A, n, N in zip(recall_sm[a_], recall_sm[n1], recall_sm.n_universe)]
recall_pd.to_csv(T / 'fig_recall_auroc_direction_perdraw.tsv', sep='\t', index=False)
recall_sm.to_csv(T / 'fig_recall_auroc_direction_summary.tsv', sep='\t', index=False)

sel = lambda **k: recall_sm.query(' and '.join(f'{c} == {v!r}' for c, v in k.items()))
Q = dict(arm='matched_cells_UMIs', exclusion='either_up', universe='per_modality')
print(sel(**Q)[['scope', 'modality', 'auroc_up', 'dir_agree', 'med_up', 'abs_err', 'n_up', 'n_universe']].round(4).to_string(index=False))


# %% [markdown]
# ### Validation against the stored A17 table

# %%
val_rows = []
def validate(name, mine, stored, keys, cols, atol=1e-9):
    stored = stored.merge(mine[keys].drop_duplicates(), on=keys, how='inner', validate='one_to_one')
    m = mine.merge(stored, on=keys, suffixes=('', '__stored'), how='outer', indicator=True)
    assert (m._merge == 'both').all(), (name, m._merge.value_counts().to_dict())
    for c in cols:
        a, b = m[c].to_numpy(float), m[c + '__stored'].to_numpy(float)
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.where(both_nan, 0.0, np.abs(a - b))
        val_rows.append(dict(stored_table=name, column=c, n_rows=len(m), max_abs_diff=float(np.nanmax(diff)) if len(diff) else np.nan,
                             nan_mismatch=int((np.isnan(a) != np.isnan(b)).sum()),
                             agrees=bool((np.nan_to_num(diff, nan=np.inf) <= atol).all() and (np.isnan(a) == np.isnan(b)).all())))
    return m

wl_stored = pd.read_csv(D / 'a17_arm/A17_whole_library.tsv', sep='\t')
mine = recall_sm[recall_sm.universe == 'per_modality'].drop(columns=['species', 'universe', 'sample_treated', 'sample_mock', 'treatment'])
validate('a17_arm/A17_whole_library.tsv', mine, wl_stored, ['arm', 'scope', 'exclusion', 'modality'],
         [c for c in wl_stored.columns if c in mine.columns and c not in ['arm', 'scope', 'exclusion', 'modality']])
print(pd.DataFrame(val_rows).to_string(index=False))


# %% [markdown]
# ## Within-set Spearman rho against bulk log2FC

# %%
def boot_ci(x, y, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    k = len(x)
    s = [spearmanr(x[i], y[i]).statistic for i in (rng.integers(0, k, (n, k)))]
    return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))

SETS = [('induced_311_signature', up)]
ex = EXC['either_up']
rrows, rsum = [], []
for scope in SCOPES:
    for arm, draws in [('matched_cells_UMIs', list(range(1, 11)))]:
        for universe in ['per_modality']:
            for mod, (a, b) in MODS.items():
                for sname, S in SETS:
                    rs = []
                    last = None
                    for d in draws:
                        V = {s: col(s, d, scope) for s in LIBS}
                        base = ((V[a] + V[b]) >= 10) if universe == 'per_modality' else \
                            np.logical_and.reduce([(V[x] + V[y]) >= 10 for x, y in MODS.values()])
                        uni = base & ~ex & S & np.isfinite(bulk)
                        if uni.sum() < 10:
                            continue
                        l = lfc(V[a], V[b])
                        r = spearmanr(l[uni], bulk[uni]).statistic
                        rs.append(r)
                        last = (l[uni], bulk[uni])
                        rrows.append(dict(species=SPECIES, arm=arm, scope=scope, universe=universe, exclusion='either_up',
                                          modality=mod, sample_treated=a, sample_mock=b, treatment='NF', gene_set=sname,
                                          draw=d, rho=float(r), n_genes=int(uni.sum())))
                    if not rs:
                        continue
                    lo, hi = boot_ci(*last) if universe == 'per_modality' else (np.nan, np.nan)
                    rsum.append(dict(species=SPECIES, arm=arm, scope=scope, universe=universe, exclusion='either_up',
                                     modality=mod, sample_treated=a, sample_mock=b, treatment='NF', gene_set=sname,
                                     rho=float(np.mean(rs)), rho_draw_sd=float(np.std(rs, ddof=1)) if len(rs) > 1 else np.nan,
                                     rho_draw_sd_ddof0=float(np.std(rs)), n_genes_mean=float(np.mean([x['n_genes'] for x in rrows[-len(rs):]])),
                                     n_genes_last_draw=len(last[0]), boot_ci_low_last_draw=lo, boot_ci_high_last_draw=hi, n_draws=len(rs)))
rho_pd = pd.DataFrame(rrows)
rho_sm = pd.DataFrame(rsum)
rho_pd.to_csv(T / 'fig_within_set_rho_perdraw.tsv', sep='\t', index=False)
rho_sm.to_csv(T / 'fig_within_set_rho_summary.tsv', sep='\t', index=False)
print(rho_sm.query("arm == 'matched_cells_UMIs' and universe == 'per_modality'")[
    ['scope', 'modality', 'gene_set', 'rho', 'rho_draw_sd_ddof0', 'n_genes_mean', 'boot_ci_low_last_draw', 'boot_ci_high_last_draw']].round(3).to_string(index=False))


# %% [markdown]
# ### Validation against stored within-induced rho and common-universe metrics

# %%
rs_ = rho_sm.query("arm == 'matched_cells_UMIs' and universe == 'per_modality'").set_index(['scope','modality','gene_set']).rho
v2_check = pd.read_csv(D / 'NF_within_induced_and_common_universe_summary.tsv', sep='\t')
v2_keys = ['scope','modality','universe']
current = recall_sm.query("arm == 'matched_cells_UMIs' and exclusion == 'either_up'").merge(
    rho_sm.query("arm == 'matched_cells_UMIs' and gene_set == 'induced_311_signature'")[v2_keys+['rho']].rename(columns={'rho':'induced_rho'}),
    on=v2_keys, validate='one_to_one')
validate('V2_audit/NF_within_induced_and_common_universe_summary.tsv', current, v2_check,
         v2_keys, [c for c in v2_check if c not in v2_keys])


# %% [markdown]
# ## Per-gene table: single-cell NF log2FC (nuclei, protoplast) against A17 bulk log2FC

# %%
gene_rows = []
a17_set = np.where(up, 'up', np.where(down, 'down', 'none'))
for scope in SCOPES:
    Ls, U, Cu = {}, {}, []
    for d in range(1, 11):
        V = {s: col(s, d, scope) for s in LIBS}
        Cu.append(np.logical_and.reduce([(V[a] + V[b]) >= 10 for a, b in MODS.values()]) & ~EXC['either_up'])
        for mod, (a, b) in MODS.items():
            Ls.setdefault(mod, []).append(lfc(V[a], V[b]))
            U.setdefault(mod, []).append(((V[a] + V[b]) >= 10) & ~EXC['either_up'])
    Cu = np.vstack(Cu)
    for mod, (a, b) in MODS.items():
        Lm = np.vstack(Ls[mod])
        Um = np.vstack(U[mod])
        na, nb = col(a, 0, scope), col(b, 0, scope)
        cts = {s: np.vstack([col(s, d, scope) for d in range(1, 11)]).mean(0) for s in (a, b)}
        gene_rows.append(pd.DataFrame(dict(
            species=SPECIES, gene=g, a17_set=a17_set, bulk_log2fc=bulk, bulk_padj=padj,
            either_up_excluded=EXC['either_up'], scope=scope, modality=mod, sample_treated=a, sample_mock=b, treatment='NF',
            sc_log2fc_mean=Lm.mean(0), sc_log2fc_sd=Lm.std(0, ddof=1), sc_log2fc_native=lfc(na, nb),
            mean_counts_treated=cts[a], mean_counts_mock=cts[b],
            frac_draws_in_universe=Um.mean(0), in_universe_all_draws=Um.all(0), in_universe_any_draw=Um.any(0),
            in_common_universe_all_draws=Cu.all(0))))
genes_tab = pd.concat(gene_rows, ignore_index=True)
genes_tab.to_csv(T / 'fig_gene_lfc_single_cell_vs_bulk.tsv', sep='\t', index=False)
print(genes_tab.shape)
print(genes_tab[genes_tab.a17_set == 'up'].groupby(['scope', 'modality']).agg(
    in_uni_mean=('frac_draws_in_universe', 'sum'), median_mean_lfc=('sc_log2fc_mean', 'median')).round(2))


# %% [markdown]
# ## Spliced/unspliced decomposition

# %%
L_DIR = D / 'layers'
LIBS6 = ['nMOCK', 'nNOD', 'pMock', 'pNOD', 'nBAP', 'pBAP']
Z = {s: np.load(L_DIR / f'{s}.npz', allow_pickle=True) for s in LIBS6}
genesL = pd.Index(Z['nMOCK']['genes'])
for s in LIBS6:
    assert list(Z[s]['genes']) == list(genesL)
RL = lfc_ref.reindex(genesL)
bulkL = RL.log2FoldChange.to_numpy()
padjL = np.nan_to_num(RL.padj.to_numpy(), nan=1.0)
upL = np.asarray(genesL.isin(SIG['up']))
sigL = np.asarray(genesL.isin(SIG['up'] | SIG['down']))
exL = np.asarray(genesL.isin(EXC_SETS['either_up']))
seed_of = lambda *k: zlib.crc32('|'.join(map(str, k)).encode())

def l_metrics(A, B):
    uni = ((A + B) >= 10) & ~exL
    l = np.log2((cpm0(A) + 1) / (cpm0(B) + 1))
    kn = uni & np.isfinite(bulkL)
    su = uni & upL & np.isfinite(bulkL)
    sg = uni & sigL
    return dict(n_universe=int(uni.sum()), n_up=int((uni & upL).sum()),
        auroc_up=auc(l[uni], upL[uni]),
        dir_agree=float(np.mean(np.sign(l[sg]) == np.sign(bulkL[sg]))) if sg.sum() else np.nan,
        med_up=float(np.median(l[uni & upL])) if (uni & upL).sum() else np.nan,
        rho_tx=float(spearmanr(l[kn], bulkL[kn]).statistic) if kn.sum() > 10 else np.nan,
        rho_within_up=float(spearmanr(l[su], bulkL[su]).statistic) if su.sum() > 10 else np.nan,
        molecules_treated=int(A.sum()), molecules_mock=int(B.sum()))

def thin(v, Dm, rng):
    return rng.multivariate_hypergeometric(v.astype(np.int64), Dm) if v.sum() > Dm else v.astype(np.int64)

CON = {'NOD': {'nuclei': ('nNOD', 'nMOCK'), 'protoplast': ('pNOD', 'pMock')},
       'BAP': {'nuclei': ('nBAP', 'nMOCK'), 'protoplast': ('pBAP', 'pMock')}}
getL = lambda s, sc, d, c: Z[s][f'{sc}__d{d}__{c}']
lrows = []
FAM = {'resolved': ['tot', 'res_s', 'res_u', 'raw_rs', 'raw_ru']}
Dtab = []
for sc in ['without_root_hair']:
    for fam, comps in FAM.items():
        Dm = min(int(getL(s, sc, d, c).sum()) for s in LIBS6 for d in range(1, 11) for c in comps if c != 'tot')
        Dm = int(Dm * 0.98)
        Dtab.append(dict(scope=sc, family=fam, D=Dm))
        for cn, pairs in [('NOD', CON['NOD'])]:
            for mod, (a, b) in pairs.items():
                for comp_ in ['res_s', 'res_u']:
                    for d in range(1, 11):
                        rng = np.random.default_rng(seed_of(sc, fam, cn, mod, comp_, d))
                        A = thin(getL(a, sc, d, comp_), Dm, rng)
                        B = thin(getL(b, sc, d, comp_), Dm, rng)
                        lrows.append(dict(arm=f'B_matched_{fam}', scope=sc, contrast=cn, modality=mod, component=comp_, draw=d, D=Dm,
                                          **l_metrics(A, B)))
layer_pd = pd.DataFrame(lrows)
lkeys = ['arm', 'scope', 'contrast', 'modality', 'component']
lnum = [c for c in layer_pd.columns if c not in lkeys + ['draw']]
layer_sm = layer_pd.groupby(lkeys, sort=False)[lnum].mean().reset_index()
lsd = layer_pd.groupby(lkeys, sort=False)[['auroc_up', 'rho_within_up', 'rho_tx']].std().add_suffix('_draw_sd').reset_index()
layer_sm = layer_sm.merge(lsd, on=lkeys)
layer_sm['auroc_up_gene_se'] = [hm(a, n1, n - n1) for a, n1, n in zip(layer_sm.auroc_up, layer_sm.n_up, layer_sm.n_universe)]

# Validate against stored results before adding labels.
validate('layers/layer_perdraw.tsv', layer_pd, pd.read_csv(L_DIR / 'layer_perdraw.tsv', sep='\t'), lkeys + ['draw'],
         [c for c in layer_pd.columns if c not in lkeys + ['draw']])
validate('layers/layer_summary.tsv', layer_sm, pd.read_csv(L_DIR / 'layer_summary.tsv', sep='\t'), lkeys,
         [c for c in layer_sm.columns if c not in lkeys])
validate('layers/layer_armB_budgets.tsv', pd.DataFrame(Dtab), pd.read_csv(L_DIR / 'layer_armB_budgets.tsv', sep='\t'), ['scope', 'family'], ['D'])
print(pd.DataFrame(val_rows).tail(40).to_string(index=False))

COMP_NAMES = {'res_s': ('resolved', 'spliced'), 'res_u': ('resolved', 'unspliced')}
def decorate_layers(df):
    df = df.copy()
    df.insert(0, 'species', SPECIES)
    df.insert(4, 'treatment', df.contrast.map({'NOD': 'NF', 'BAP': 'BAP'}))
    df['layer_family'] = df.component.map(lambda c: COMP_NAMES[c][0])
    df['layer'] = df.component.map(lambda c: COMP_NAMES[c][1])
    df['sample_treated'] = [CON[c][m][0] for c, m in zip(df.contrast, df.modality)]
    df['sample_mock'] = [CON[c][m][1] for c, m in zip(df.contrast, df.modality)]
    return df
decorate_layers(layer_pd).to_csv(T / 'fig_layers_auroc_perdraw.tsv', sep='\t', index=False)
decorate_layers(layer_sm).to_csv(T / 'fig_layers_auroc_summary.tsv', sep='\t', index=False)
q = layer_sm.query("arm == 'B_matched_resolved' and scope == 'without_root_hair' and contrast == 'NOD'")
print(q[['modality', 'component', 'auroc_up', 'auroc_up_gene_se', 'dir_agree', 'med_up']].round(3).to_string(index=False))


# %% [markdown]
# ## Reference-free cross-modality overlap

# %%
F_DIR = D / 'reference_free'
ZR = {s: np.load(F_DIR / f'{s}.npz', allow_pickle=True) for s in LIBS6}
genesR = pd.Index(ZR['nMOCK']['genes'])
for s in LIBS6:
    assert list(ZR[s]['genes']) == list(genesR)
EXCR = {'either_up': np.asarray(genesR.isin(EXC_SETS['either_up']))}
lf0 = lambda a, b: np.log2((cpm0(a) + 1) / (cpm0(b) + 1))
def sp_(a, b):
    return float(spearmanr(a, b).statistic) if len(a) >= 10 else np.nan
RULES = ['2fold_reproducible', 'top100', 'top200', 'top500']
frows = []
for sc in ['without_root_hair']:
    for tr, pairs in CON.items():
        for exn, exr in [('either_up', EXCR['either_up'])]:
            for d in range(1, 11):
                qv = lambda s, k: ZR[s][f'{sc}__d{d}__Q{k}']
                full = lambda s: sum(qv(s, k) for k in range(1, 5))
                uni = np.logical_and.reduce([(full(t) + full(m)) >= 10 for t, m in pairs.values()]) & ~exr
                P = {}
                for mod, (t, m) in pairs.items():
                    P[mod] = dict(Q={k: lf0(qv(t, k), qv(m, k)) for k in range(1, 5)},
                                  H=lf0(qv(t, 3) + qv(t, 4), qv(m, 3) + qv(m, 4)),
                                  ok=uni & ((qv(t, 1) + qv(m, 1)) >= 5) & ((qv(t, 2) + qv(m, 2)) >= 5))
                N = int(uni.sum())
                for direction, sg in [('induced', 1.0)]:
                    for rule in RULES:
                        calls = {}
                        for mod, p in P.items():
                            x1, x2 = sg * p['Q'][1], sg * p['Q'][2]
                            if rule == '2fold_reproducible':
                                c = p['ok'] & (x1 > 1) & (x2 > 1)
                            else:
                                k = int(rule[3:])
                                score = np.where(p['ok'], np.minimum(x1, x2), -np.inf)
                                idx = np.argsort(-score)[:k]
                                c = np.zeros(len(genesR), bool)
                                c[idx] = True
                                c &= p['ok'] & (np.minimum(x1, x2) > 0)
                            calls[mod] = c
                        cn, cp = calls['nuclei'], calls['protoplast']
                        both = cn & cp
                        either = cn | cp
                        exp = cn.sum() * cp.sum() / max(N, 1)
                        tab = [[int(both.sum()), int((cn & ~cp).sum())], [int((~cn & cp & uni).sum()), int((uni & ~cn & ~cp).sum())]]
                        r = dict(scope=sc, treatment=tr, exclusion=exn, direction=direction, rule=rule, draw=d, n_universe=N,
                                 n_called_nuclei=int(cn.sum()), n_called_protoplast=int(cp.sum()), n_called_both=int(both.sum()),
                                 n_both_expected=exp, enrichment=both.sum() / exp if exp > 0 else np.nan,
                                 jaccard=both.sum() / either.sum() if either.sum() else np.nan,
                                 frac_nuclei_calls_also_protoplast=both.sum() / cn.sum() if cn.sum() else np.nan,
                                 frac_protoplast_calls_also_nuclei=both.sum() / cp.sum() if cp.sum() else np.nan,
                                 odds_ratio=fisher_exact(tab)[0] if min(cn.sum(), cp.sum()) > 0 else np.nan)
                        for src, tgt in [('nuclei', 'protoplast'), ('protoplast', 'nuclei')]:
                            c = calls[src]
                            bg = uni & ~c
                            r[f'{src}_calls_same_dir_in_{tgt}_heldout'] = float(np.mean(sg * P[tgt]['H'][c] > 0)) if c.sum() else np.nan
                            r[f'background_same_dir_in_{tgt}_heldout__vs_{src}'] = float(np.mean(sg * P[tgt]['H'][bg] > 0))
                            r[f'{src}_calls_same_dir_in_own_heldout'] = float(np.mean(sg * P[src]['H'][c] > 0)) if c.sum() else np.nan
                        for setn, S in [('both', both), ('either', either)]:
                            r[f'rank_{setn}_n'] = int(S.sum())
                            r[f'rank_{setn}_rho_cross_heldout'] = sp_(P['nuclei']['H'][S], P['protoplast']['H'][S])
                            r[f'rank_{setn}_rho_cross_quarter'] = np.nanmean([sp_(P['nuclei']['Q'][3][S], P['protoplast']['Q'][4][S]),
                                                                              sp_(P['nuclei']['Q'][4][S], P['protoplast']['Q'][3][S])])
                            r[f'rank_{setn}_rho_within_nuclei_quarter'] = sp_(P['nuclei']['Q'][3][S], P['nuclei']['Q'][4][S])
                            r[f'rank_{setn}_rho_within_protoplast_quarter'] = sp_(P['protoplast']['Q'][3][S], P['protoplast']['Q'][4][S])
                            r[f'rank_{setn}_median_heldout_nuclei'] = float(np.median(sg * P['nuclei']['H'][S])) if S.sum() else np.nan
                            r[f'rank_{setn}_median_heldout_protoplast'] = float(np.median(sg * P['protoplast']['H'][S])) if S.sum() else np.nan
                        frows.append(r)
ref_pd = pd.DataFrame(frows)
fk = ['scope', 'treatment', 'exclusion', 'direction', 'rule']
ref_sm = ref_pd.groupby(fk, sort=False).mean(numeric_only=True).drop(columns='draw')
rsd = ref_pd.groupby(fk, sort=False)[['enrichment', 'jaccard', 'rank_both_rho_cross_heldout', 'rank_both_rho_cross_quarter',
                                      'rank_either_rho_cross_quarter']].std().add_suffix('_draw_sd')
ref_sm = ref_sm.join(rsd).reset_index()
validate('reference_free/refree_perdraw.tsv', ref_pd, pd.read_csv(F_DIR / 'refree_perdraw.tsv', sep='\t'), fk + ['draw'],
         [c for c in ref_pd.columns if c not in fk + ['draw']])
validate('reference_free/refree_summary.tsv', ref_sm, pd.read_csv(F_DIR / 'refree_summary.tsv', sep='\t'), fk,
         [c for c in ref_sm.columns if c not in fk])

def decorate_refree(df):
    df = df.rename(columns={'treatment': 'contrast'}).copy()
    df.insert(0, 'species', SPECIES)
    df.insert(3, 'treatment', df.contrast.map({'NOD': 'NF', 'BAP': 'BAP'}))
    df.insert(4, 'modality_pair', 'nuclei_vs_protoplast')
    return df
decorate_refree(ref_pd).to_csv(T / 'fig_refree_overlap_perdraw.tsv', sep='\t', index=False)
decorate_refree(ref_sm).to_csv(T / 'fig_refree_overlap_summary.tsv', sep='\t', index=False)
q = ref_sm.query("scope == 'without_root_hair' and exclusion == 'either_up' and direction == 'induced'")
print(q[['treatment', 'rule', 'n_called_nuclei', 'n_called_protoplast', 'n_called_both', 'enrichment',
         'frac_nuclei_calls_also_protoplast', 'frac_protoplast_calls_also_nuclei', 'rank_both_rho_cross_quarter',
         'rank_both_rho_within_nuclei_quarter', 'rank_both_rho_within_protoplast_quarter']].round(3).to_string(index=False))


# %% [markdown]
# ## Validation summary

# %%
val = pd.DataFrame(val_rows)
val.to_csv(T / 'validation_vs_stored.tsv', sep='\t', index=False)
print(val.groupby('stored_table').agg(columns=('column', 'size'), all_agree=('agrees', 'all'), max_abs_diff=('max_abs_diff', 'max')).to_string())
bad = val[~val.agrees]
print('disagreeing columns:', 'none' if bad.empty else bad.to_string(index=False))
assert bad.empty, 'V2 pseudobulk recomputation disagrees with stored producer tables'
