#!/usr/bin/env python3
import math
import os
from pathlib import Path
import pandas as pd

base_env = os.environ.get('SC_PCQTL_GIMAP_MIXING_ROOT', '')
if not base_env:
    raise SystemExit('Set SC_PCQTL_GIMAP_MIXING_ROOT.')
BASE = Path(base_env).expanduser().resolve()
(BASE / 'coloc_susie').mkdir(parents=True, exist_ok=True)
manifest = pd.read_csv(BASE / 'data/tensorqtl_job_manifest.tsv', sep='\t')
scenario_summary = pd.read_csv(BASE / 'data/gimap_pseudobulk_scenario_summary.tsv', sep='\t')

def find_first(patterns):
    for pat in patterns:
        hits = sorted(BASE.glob(pat))
        if hits:
            return hits[0]
    return None

def read_table(path):
    if path is None or not path.exists():
        return None
    if path.suffix == '.parquet':
        return pd.read_parquet(path)
    return pd.read_csv(path, sep='\t')

def classify(pid):
    if '__e__' in str(pid):
        return 'single_gene'
    return 'cluster_PC'

def short_name(pid):
    pid = str(pid)
    if '__e__' in pid:
        return pid.split('__e__', 1)[1]
    return pid.rsplit('__', 1)[-1]

phenotype_rows = []
sumstat_frames = []
for _, row in manifest.iterrows():
    sc = row['scenario']; pset = row['phenotype_set']
    cis_path = find_first([f'tensorqtl/cis/{sc}/{pset}*.cis_qtl.txt.gz', f'tensorqtl/cis/{sc}/{pset}*.cis_qtl.txt'])
    cis = read_table(cis_path)
    if cis is not None and len(cis):
        for _, r in cis.iterrows():
            pid = r.get('phenotype_id', r.get('phenotype', None))
            phenotype_rows.append({
                'scenario': sc, 'phenotype_set': pset, 'phenotype_id': pid,
                'readout_type': classify(pid), 'readout_name': short_name(pid),
                'cis_status': 'complete',
                'lead_variant_cis': r.get('variant_id', r.get('variant', None)),
                'pval_nominal_cis': r.get('pval_nominal', math.nan),
                'pval_perm': r.get('pval_perm', r.get('pval_beta', math.nan)),
                'qval': r.get('qval', math.nan),
                'slope_cis': r.get('slope', math.nan),
            })
    else:
        phenotype_rows.append({'scenario': sc, 'phenotype_set': pset, 'phenotype_id': None, 'readout_type': None, 'readout_name': None, 'cis_status': 'missing'})

    nom_path = find_first([f'tensorqtl/cis_nominal/{sc}/{pset}*.cis_qtl_pairs.*.parquet', f'tensorqtl/cis_nominal/{sc}/{pset}*.cis_qtl_pairs.*.txt.gz', f'tensorqtl/cis_nominal/{sc}/{pset}*.cis_qtl_pairs.*.txt'])
    nom = read_table(nom_path)
    if nom is not None and len(nom):
        if pset in {'pc_only', 'single_gene_only'}:
            variant_col = 'variant_id' if 'variant_id' in nom.columns else 'variant'
            required = {'phenotype_id', variant_col, 'af', 'slope', 'slope_se'}
            missing = sorted(required.difference(nom.columns))
            if missing:
                raise ValueError(f'{nom_path} lacks full-summary columns: {missing}')
            full = nom[['phenotype_id', variant_col, 'af', 'slope', 'slope_se']].copy()
            full = full.rename(columns={variant_col: 'variant_id'})
            full.insert(0, 'readout', full['phenotype_id'].map(short_name))
            full.insert(0, 'qtl_type', 'pcQTL' if pset == 'pc_only' else 'eQTL')
            full.insert(0, 'scenario', sc)
            full = full.drop(columns='phenotype_id')
            sumstat_frames.append(full)

phenotype_df = pd.DataFrame(phenotype_rows)

# Add BH q-value fallback from TensorQTL permutation p-values when native qval is unavailable.
def bh_qvalues(vals):
    import numpy as np
    vals = pd.to_numeric(pd.Series(vals), errors='coerce')
    q = pd.Series([math.nan] * len(vals), index=vals.index, dtype='float64')
    ok = vals.notna()
    if ok.sum() == 0:
        return q
    p = vals[ok].astype(float)
    order = p.sort_values().index
    ranked = p.loc[order]
    m = len(ranked)
    raw = ranked.values * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(raw[::-1])[::-1]
    q.loc[order] = np.minimum(adj, 1.0)
    return q
if len(phenotype_df) and 'pval_perm' in phenotype_df.columns:
    phenotype_df['qval_native'] = phenotype_df.get('qval', math.nan)
    phenotype_df['qval_bh_from_perm'] = math.nan
    for (sc, ps), idx in phenotype_df.groupby(['scenario', 'phenotype_set']).groups.items():
        phenotype_df.loc[idx, 'qval_bh_from_perm'] = bh_qvalues(phenotype_df.loc[idx, 'pval_perm']).values
    phenotype_df['qval'] = pd.to_numeric(phenotype_df.get('qval_native', math.nan), errors='coerce')
    phenotype_df['qval'] = phenotype_df['qval'].where(phenotype_df['qval'].notna(), phenotype_df['qval_bh_from_perm'])
phenotype_df.to_csv(BASE / 'data/gimap_pseudobulk_tensorqtl_phenotype_level_results.tsv', sep='\t', index=False)
if not sumstat_frames:
    raise RuntimeError('No pc_only or single_gene_only cis-nominal outputs were found.')
sumstats = pd.concat(sumstat_frames, ignore_index=True)
n_map = scenario_summary.set_index('scenario')['n_donors']
sumstats['qtl_n'] = sumstats['scenario'].map(n_map)
if sumstats['qtl_n'].isna().any():
    missing_scenarios = sorted(sumstats.loc[sumstats['qtl_n'].isna(), 'scenario'].unique())
    raise ValueError(f'Missing donor count for scenarios: {missing_scenarios}')
sumstats.to_csv(BASE / 'coloc_susie/mix_all_sumstats.tsv', sep='\t', index=False)
print('Wrote TensorQTL phenotype summaries and full cis-nominal coloc input')
