# Module 07 — Strict graph signal grouping

Groups colocalized QTL/GWAS signals into strict graph connected components at
several PPH4 thresholds, producing the connected-component class counts that
underlie the headline pcQTL-versus-eQTL colocalization comparison.

- **Environment variables:** `SC_PCQTL_FORMAL_COLOC_ROOT`.
- **Upstream inputs:** the `coloc.susie` outputs from module 06.

The scripts load the reviewed module-06 configuration/helpers from this code
checkout while reading and writing analysis data only under
`SC_PCQTL_FORMAL_COLOC_ROOT`. No copied scripts are required in the external
work directory.

## Scripts (run order)

| Script | Purpose |
|---|---|
| `01_compute_qtl_qtl_signal_edges.R` | Build QTL–QTL shared-signal edges from credible sets. |
| `02_build_strict_graph_groups.R` | Connected-component grouping at each PPH4 threshold. |

## Running

Run `01` → `02` after module 06 colocalization outputs are available, with
`SC_PCQTL_FORMAL_COLOC_ROOT` pointing at the colocalization tree.

## Outputs

Strict signal groups and connected-component class counts feeding the main
colocalization figure and the strict-graph supplementary tables (assembled in
module 11).
