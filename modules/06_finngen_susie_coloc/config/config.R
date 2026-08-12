# Formal cis-pcQTL/eQTL-GWAS colocalization configuration.

.config_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
.module_dir <- if (nzchar(.config_file)) dirname(dirname(.config_file)) else getwd()
.release_root <- normalizePath(file.path(.module_dir, "..", ".."), mustWork = TRUE)
source(file.path(.release_root, "config", "celltype_eligibility.R"))
.celltype_eligibility <- load_celltype_eligibility(file.path(.release_root, "config", "celltype_eligibility.tsv"))

ROOT_DIR <- Sys.getenv("SC_PCQTL_FORMAL_COLOC_ROOT", unset = "")
if (!nzchar(ROOT_DIR)) {
  stop("Set SC_PCQTL_FORMAL_COLOC_ROOT to an external formal-colocalization work directory.")
}
ROOT_DIR <- normalizePath(ROOT_DIR, mustWork = FALSE)
WORKFLOW_ROOT <- Sys.getenv("SC_PCQTL_WORKFLOW_ROOT", unset = dirname(ROOT_DIR))

PCQTL_ROOT <- Sys.getenv(
  "SC_PCQTL_UPSTREAM_ROOT",
  unset = file.path(WORKFLOW_ROOT, "03_analysis_celltypes/01_upstream_main_pipeline_add_cov")
)
EQTL_ROOT  <- Sys.getenv("ONEK1K_EQTL_ROOT", unset = "")
FINNGEN_ROOT <- Sys.getenv("FINNGEN_ROOT", unset = file.path(ROOT_DIR, "resources", "finngen"))

ONEK1K_RAW_GENOTYPE_PREFIX <- Sys.getenv("ONEK1K_RAW_GENOTYPE_PREFIX", unset = "")
ONEK1K_MAF005_GENOTYPE_PREFIX <- Sys.getenv("ONEK1K_MAF005_GENOTYPE_PREFIX", unset = "")
ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN <- Sys.getenv(
  "ONEK1K_QTL_LIFTOVER_POSMAP_PATTERN",
  unset = file.path(ROOT_DIR, "resources/onek1k_liftover/chr%s.hg19_to_hg38.posmap")
)
ONEK1K_QTL_COORDINATE_BUILD <- "hg38_lifted_from_OneK1K_hg19_like_positions"

PLINK <- Sys.getenv("PLINK", "plink")

FINNGEN_R12_MANIFEST_URL <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"

# Official FinnGen SuSiE fine-mapping files are the main GWAS fine-mapping
# source. GWAS LD is therefore FinnGen's in-sample LD embedded in the official
# fine-mapping, not OneK1K LD.
FINNGEN_OFFICIAL_SUSIE_DIR <- Sys.getenv(
  "FINNGEN_OFFICIAL_SUSIE_DIR",
  unset = file.path(FINNGEN_ROOT, "finemapping", "r12_susie")
)

ALL_CELLTYPE_EQTL_MAP <- setNames(
  paste0("cis_", .celltype_eligibility$eqtl_celltype, ".tar.gz"),
  .celltype_eligibility$celltype
)
CELLTYPE_EQTL_MAP <- ALL_CELLTYPE_EQTL_MAP[.celltype_eligibility[include_primary == TRUE, celltype]]

# Reviewable endpoint universe: broad immune/inflammation/blood/infection scope.
# This list is applied to FinnGen phenotype/category text before any QTL/GWAS
# signal inspection to avoid winner's curse in endpoint selection.
FINNGEN_PRIORITY_REGEX <- paste(c(
  "immune", "autoimmune", "immun", "inflamm", "infection", "infectious",
  "blood", "haemat", "hemat", "anaemia", "anemia", "leuk", "lymph",
  "myelo", "neutro", "eosin", "thrombo", "platelet", "coag",
  "allerg", "asthma", "atopic", "dermatitis", "psoriasis", "arthritis",
  "rheumat", "lupus", "scleroderma", "vasculitis", "crohn",
  "ulcerative colitis", "colitis", "celiac", "coeliac", "thyroid",
  "diabetes", "multiple sclerosis", "sarcoid", "sepsis", "tuberc",
  "viral", "bacterial", "fungal", "parasit"
), collapse = "|")

MAIN_MIN_CASES <- 1000L

SUSIE_L <- 10L
SUSIE_COVERAGE <- 0.95
CS_MIN_ABS_CORR <- 0.5
COLOC_PPH4_CUTOFF <- 0.75
COLOC_PPH4_SENSITIVITY_CUTOFF <- 0.8
COLOC_SHARED_ALPHA_MASS_MIN <- 0.9

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
