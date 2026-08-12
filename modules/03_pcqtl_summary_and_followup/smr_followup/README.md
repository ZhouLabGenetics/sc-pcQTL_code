# Targeted GIMAP SMR follow-up

This directory contains only the descriptive SMR/HEIDI analysis used for the
representative CD8-naive-T GIMAP locus. It does not run a genome-wide driver-gene
screen and does not contribute to strict SuSiE signal-group definitions.

`run_targeted_gimap_smr.py` runs two outcomes against the same seven-gene CD8
naive T-cell eQTL BESD:

1. the GIMAP cluster PC3 pcQTL summary;
2. FinnGen R12 lymphocyte-count endpoint 3019198.

Both runs use the OneK1K chromosome 7 genotype LD reference and the executable
distributed in the SMR 1.4.0 Linux release archive. The executable's own banner
reports version 1.3.2 (build 27 June 2025); its SHA-256 and banner metadata are
recorded separately in `run_provenance.json`. Runs use `--peqtl-smr 0.05`,
`--peqtl-heidi 0.01`,
`--heidi-min-m 3`, `--disable-freq-ck`, and the SMR HEIDI test. The HEIDI
candidate-eQTL threshold is set explicitly to 0.01 for every gene and both
outcomes so that descriptive HEIDI follow-up is not dropped solely because the
default candidate threshold leaves fewer than three variants. Figure labels use
the seven-gene Bonferroni reference (`0.05 / 7`) and HEIDI consistency reference
`p_HEIDI > 0.01`.

## Required staged inputs

- a harmonized SMR `.ma` file for GIMAP PC3;
- a harmonized SMR `.ma` file made specifically from FinnGen R12 endpoint
  3019198, not a substituted lymphocyte-count meta-GWAS;
- a GIMAP-only CD8-naive-T eQTL BESD prefix (`.besd`, `.epi`, `.esi`);
- the OneK1K chromosome 7 PLINK LD prefix (`.bed`, `.bim`, `.fam`).

`prepare_finngen3019198_ma.py` constructs the FinnGen input reproducibly. It
reads the raw GRCh38 FinnGen laboratory-value summary, reverse-maps chromosome 7
through the OneK1K hg19-to-hg38 position map, aligns alleles to the
OneK1K PLINK `.bim`, flips beta and allele frequency when needed, removes
palindromic SNPs, and records hashes and filtering counts in a QC JSON. The
endpoint sample size is 183,481, as reported in the FinnGen R12 lab-value
manifest.

```bash
python prepare_finngen3019198_ma.py \
  --finngen-gwas <finngen_R12_3019198.gz> \
  --posmap <chr7.hg19_to_hg38.posmap> \
  --ld-bim <ONEK1K_CHR7_PLINK_PREFIX.bim> \
  --output-ma <FINNGEN_R12_3019198.ma> \
  --qc-json <FINNGEN_R12_3019198.ma.qc.json>
```

The runner validates the `.ma` columns and refuses to run unless the supplied
QC JSON identifies endpoint 3019198, GRCh38-to-hg19 conversion, sample size
183,481, and the exact `.ma` SHA-256 hash.

```bash
python run_targeted_gimap_smr.py \
  --cluster-pc3-ma <GIMAP_PC3.ma> \
  --finngen-gwas-ma <FINNGEN_R12_3019198.ma> \
  --finngen-ma-qc <FINNGEN_R12_3019198.ma.qc.json> \
  --eqtl-besd <CD8_NAIVE_GIMAP_BESD_PREFIX> \
  --ld-prefix <ONEK1K_CHR7_PLINK_PREFIX> \
  --smr-bin <SMR_1.4.0_BINARY> \
  --out-dir <OUTPUT_DIRECTORY>
```

Use `SC_PCQTL_GIMAP_SMR_ROOT` from `config/example.env` as
`<OUTPUT_DIRECTORY>` so module 11 consumes these summaries.

Outputs are `gimap_cluster_pcqtl_smr_summary.tsv`,
`gimap_gwas_smr_summary.tsv`, a combined summary, per-gene SMR output/logs, and
`run_provenance.json`. The retained workflow rerun yielded 2,734
allele-compatible variants after SMR's intersection
of FinnGen, GIMAP eQTL BESD, and OneK1K LD inputs. Under the uniform main HEIDI
setting, every targeted gene returns a HEIDI result for both outcomes. HEIDI is
reported as an auxiliary consistency diagnostic and is not used to assign a
sole mediator.
