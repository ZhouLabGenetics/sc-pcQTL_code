# Module 01 — Pairwise hurdle co-expression screen and local cluster calling

For each cell type, screens retained within-chromosome gene pairs with the
`fasthurdle` co-expression model and calls local gene clusters with a sliding window. This is
the entry stage of the workflow: its clusters define the pcQTL phenotypes built
in module 02 and the units annotated in module 04.

Primary analysis begins with the central `../../config/celltype_eligibility.tsv`
contract. Cell types with fewer than 10,000 cells are excluded before the
pairwise screen, rather than removed from clusters or QTL results downstream.

- **Environment variables:** `SC_PCQTL_CELL_TYPE`, `SC_PCQTL_CELLTYPE_MANIFEST`, `COQTL_COUNT_FILE`, `COQTL_DATA_ROOT`, `COQTL_GENE_INFO_FILE`, `SC_PCQTL_HURDLE_WORK_DIR`.
- **Requires:** the custom R package `fasthurdle (>= 1.1.1)` (see `../../docs/dependencies.md`).
- **Upstream inputs:** OneK1K per-cell-type count tables with library-size/covariate columns and a gene coordinate table.

## Scripts (run order)

| Script | Purpose |
|---|---|
| `config.example.R` → copy to `config.R` | Per-cell-type paths and parameters. |
| `load_config.R` | Loads `config.R` / environment variables. |
| `step1_filter_sparse_genes.R` | Drops genes below the 1% nonzero-cell cutoff. |
| `calculate_chunks.R`, `compute_library_size_cache.R` | Chunk planning (keeps jobs < 6 h) and library-size cache. |
| `step2_calculate_sc_associations_single_chr.R` / `step2_calculate_sc_associations_chunked.R` | `fasthurdle` pairwise associations per chromosome (single or chunked). |
| `step2b_merge_chunks.R` | Merges chunked association outputs per chromosome. |
| `step3_identify_clusters_single_chr.R` | Sliding-window cluster calling per chromosome. |
| `step4_merge_clusters.R` | Merges clusters across chromosomes. |

## Running

Copy `config.example.R` to `config.R`, set `SC_PCQTL_CELL_TYPE`,
`COQTL_COUNT_FILE`, and `COQTL_GENE_INFO_FILE`, then submit the Slurm wrappers in
order: `submit_step1.sh` → `submit_step2_chunked.sh` (or `submit_step2.sh`) →
`submit_step2b_merge.sh` → `submit_step3.sh`, then run `step4_merge_clusters.R`.
`submit_all_chunked.sh` chains the whole pipeline with Slurm dependencies. The
R scripts can also be run directly after staging inputs. The wrapper skips
ineligible cell types, and direct module-02 execution rejects them explicitly.

## Outputs

Per-cell-type cluster definitions and assigned-gene tables, consumed by
module 02 (cluster-PC QTL) and module 04 (annotation enrichment).
