# Module 08 — PIP-weighted gene effects and regulatory annotation

Characterizes strict colocalized signal groups using PIP-weighted nominal gene
effects and regulatory annotations.

- **Environment variables:** `SC_PCQTL_FORMAL_COLOC_ROOT`, `SC_PCQTL_UPSTREAM_ROOT`, `SC_PCQTL_CLUSTER_ENRICHMENT_ROOT`.
- **Upstream inputs:** strict signal groups (module 07), `coloc.susie` outputs (module 06), and QTL credible-set nominal-effect and regulatory-annotation summaries.

The scripts load configuration from module 06 and use
`SC_PCQTL_FORMAL_COLOC_ROOT` as the external analysis-data root.

## Scripts

| Script | Purpose |
|---|---|
| `01_strict_gene_effects.R` | Max absolute and cross-gene-CV PIP-weighted nominal effects by strict class and the BH-adjusted pairwise Wilcoxon results used by the final figure builder. |
| `02_strict_regulatory_annotation.R` | PIP-weighted regulatory-annotation distribution by strict class. |

## Running

Run with `SC_PCQTL_FORMAL_COLOC_ROOT` pointing at the colocalization tree, after
modules 06–07. Final figure composition is handled by module 11.

## Outputs

The pcQTL-specific loading, effect, and regulatory-annotation tables used to
build the two gene-effect panels in module 11.
