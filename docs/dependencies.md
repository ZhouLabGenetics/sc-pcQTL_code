# Software Dependencies and Environment

This release is **source code only**; it does not bundle interpreters, packages,
containers, or data. The pipeline was developed and run on a Slurm-based Linux
cluster. Most analysis steps are R; shell scripts orchestrate Slurm jobs and
containerized tools; a few helpers are Python.

Record your own environment with `sessionInfo()` (R) and `pip freeze` (Python)
when you reproduce results. Exact study versions are given below where they were
captured in the manuscript reproduction logs.

## Runtimes

| Runtime | Study version | Notes |
|---|---|---|
| R | 4.5.1 | Primary analysis and plotting language. |
| Python | 3 (cluster environment) | Helper scripts only; see below. |
| bash | POSIX / GNU | Job submission and tool orchestration. |
| Slurm | cluster scheduler | `submit_*.sh`, `*.sbatch`, `*.slurm` wrappers assume Slurm; the underlying R/Python can also be run directly. |

## R packages

Install the release's direct non-base dependencies from CRAN:

```r
install.packages(c(
  "data.table", "dplyr", "stringr", "jsonlite", "optparse",
  "ggplot2", "ggrepel", "scales", "gridExtra", "cowplot",
  "susieR", "coloc"
))
```

Base R packages also used directly: `parallel`, `grid`, `stats`, `tools`,
`utils`, `grDevices`, `graphics`.

| Purpose | Packages |
|---|---|
| Data wrangling / IO | data.table, dplyr, stringr, jsonlite |
| CLI argument parsing | optparse |
| Statistical genetics | **fasthurdle** (see below), susieR, coloc |
| Plotting | ggplot2, ggrepel, scales, gridExtra, cowplot |

Recorded study versions (from the manuscript figure-reproduction logs): R 4.5.1,
data.table 1.18.2.1, ggplot2 4.0.3, cowplot 1.2.0.

### Required custom package: `fasthurdle`

The hurdle co-expression screen (`modules/01`), hurdle simulations
(`modules/05`), and joint-score sensitivity (`modules/12`) require
the R package **`fasthurdle` (>= 1.1.1)**, a fast
C++/Rcpp hurdle-regression implementation authored by
Masahiro Kanai. It is **not** part of this repository and is licensed
separately (GPL-2). It is publicly available at
<https://github.com/mkanai/fasthurdle>; install it before running those
modules, e.g.:

```r
# requires a C++ toolchain
remotes::install_github("mkanai/fasthurdle")
```

`fasthurdle` itself depends on `fastglm (>= 0.0.4)`, `Rcpp (>= 1.0.14)`, and
links to `Rcpp`, `RcppArmadillo`, and `roptim` (a C++ toolchain is required to
build it). Both module-01 association workers enforce the minimum version at
runtime.

The real-data global-permutation null and its matched within-donor
diagnostic used `fasthurdle` 1.1.1 with `fastglm` 0.0.4. The production
pairwise screen used `fasthurdle` 1.2.0. The within-donor diagnostic enforces
the two published package versions and records the versions loaded for each
regenerated run. Module 12 requires exactly `fasthurdle` 1.2.0 because its
score-test cache and batch interfaces use that build.

## Python packages

The repository's own Python scripts use **NumPy**, **pandas**, and **SciPy** plus the standard library
(`argparse`, `csv`, `gzip`, `json`, `pathlib`, `re`, `subprocess`, `tarfile`,
`urllib`, ...). Install with:

```bash
pip install -r requirements.txt
```

`tensorQTL` and `LDSC` are used as **external tools** with their own
environments (not pinned in `requirements.txt`); see the next section.

## External command-line tools and containers

Configure executable/image paths in `config/example.env`.

| Tool | Used by | Purpose | Source |
|---|---|---|---|
| SAIGE-QTL 0.3.4 (Singularity image) | module 02 | Cell-type cluster-PC cis mapping | https://github.com/weizhou0/qtl |
| Singularity / Apptainer | module 02 | Container runtime for SAIGE-QTL | https://apptainer.org |
| PLINK (1.9 and/or 2) | modules 02, 03, 06, 10 | Genotype handling and LD reference | https://www.cog-genomics.org/plink/ |
| tensorQTL | module 10 (pseudobulk mixing) | Donor-level pseudobulk cis-QTL | https://github.com/broadinstitute/tensorqtl |
| SMR / HEIDI 1.4.0 | module 03 (targeted GIMAP follow-up) | Descriptive gene-level SMR and HEIDI | https://yanglab.westlake.edu.cn/software/smr/ |
| LDSC (ldsc.py, munge_sumstats.py) | module 09 | Stratified LD score regression (S-LDSC) | https://github.com/bulik/ldsc |
| liftOver | modules 06, 10 | hg19 → hg38 coordinate mapping | https://genome.ucsc.edu/ |

## Reference resources (not software)

S-LDSC needs baselineLD v2.2 annotations and 1000 Genomes Phase 3 EUR
weights/frequencies; colocalization needs FinnGen R12 official SuSiE
fine-mapping; locus and annotation steps need GENCODE gene models and the
genomic annotation resources listed in `input_inventory.md`. These are obtained
from their original providers and are not redistributed here.
