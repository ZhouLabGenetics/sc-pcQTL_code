# Module 03 — pcQTL summaries and manuscript follow-up

Final downstream counting and targeted follow-up of cluster-PC QTLs from module
02.

- **Environment variables:** `SC_PCQTL_WORKFLOW_ROOT`,
  `SC_PCQTL_PCQTL_SUMMARY_ROOT`, `SC_PCQTL_PCQTL_EGENE_ROOT`,
  `SC_PCQTL_UPSTREAM_ROOT`, `ONEK1K_EQTL_ROOT`, `SMR_BIN`.
- **Upstream inputs:** pcQTL outputs (module 02) and explicitly staged inputs
  for the targeted GIMAP SMR analysis.

## Sub-pipelines

| Subdirectory | Purpose |
|---|---|
| `pcqtl_summary/` | `01_collect_results` followed by `02_identify_qtls`; applies cis-variant BH separately within each successfully tested cluster-PC phenotype. `03_extract_primary_egenes` applies the same within-phenotype BH rule to single-gene eQTLs, and `04_summarize_pcqtl_egene_overlap` produces the eGene-overlap summary. |
| `smr_followup/` | Harmonize raw FinnGen R12 endpoint 3019198 to the OneK1K hg19/LD allele frame, then run the two manuscript-targeted GIMAP SMR/HEIDI outcomes (`prepare_finngen3019198_ma.py`; `run_targeted_gimap_smr.py`). |

## Running

Run `pcqtl_summary/run_all.sh` after the pcQTL outputs from module 02
are complete and set `ONEK1K_EQTL_ROOT` for the matched eGene comparison. Set
`SC_PCQTL_RUN_EGENE_OVERLAP=0` only when regenerating count-only intermediate
outputs. Build the FinnGen `.ma` with
`smr_followup/prepare_finngen3019198_ma.py`, then run
`smr_followup/run_targeted_gimap_smr.py` with the final GIMAP PC3, FinnGen R12
endpoint 3019198 QC JSON, eQTL BESD, and OneK1K chromosome 7 LD inputs.

## Outputs

The pcQTL phenotype-level table, significant subset, per-cell-type counts,
pcQTL/eGene overlap summaries deposited as Supplementary Data, ACAT values used
by Supplementary Table S4, and targeted SMR summaries assembled into
Supplementary Figure S10.
