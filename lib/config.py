"""
Single source of truth for paths, constants, and thresholds.
"""

from pathlib import Path

# ==============================================================================
# Project paths
# ==============================================================================

PROJECT_DIR = Path('/data/parks34/projects/4germicb')
DATA_DIR = PROJECT_DIR / 'data'
RESULTS_DIR = PROJECT_DIR / 'results'
SCRIPTS_DIR = PROJECT_DIR / 'scripts'

# Upstream WES pipeline
WES_DIR = Path('/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED')

# AlphaGenome results
ALPHAGENOME_DIR = RESULTS_DIR / 'alphagenome'

# eQTL validation results
EQTL_VALIDATION_DIR = RESULTS_DIR / 'eqtl_validation'

# Prioritized results
PRIORITIZED_DIR = RESULTS_DIR / 'prioritized'

# ==============================================================================
# eQTL reference paths
# ==============================================================================

DICE_HG38_DIR = DATA_DIR / 'eqtl_references' / 'dice' / 'hg38'
ONEK1K_PATH = DATA_DIR / 'eqtl_references' / 'oneK1K' / 'esnp_table.tsv.gz'
CIMA_EQTL_DIR = Path('/data/Jiang_Lab/Data/Seongyong/CIMA/xQTL')

# ==============================================================================
# Thresholds
# ==============================================================================

PVAL_SUGGESTIVE = 0.05
PVAL_STRINGENT = 0.001
FDR_THRESHOLD = 0.05

ALPHAGENOME_HIGH_THRESHOLD = 0.5
ALPHAGENOME_MEDIUM_THRESHOLD = 0.2

# ==============================================================================
# AlphaGenome API settings
# ==============================================================================

CHECKPOINT_INTERVAL = 100
MAX_RETRIES = 5
INITIAL_BACKOFF = 1.0  # seconds

# ==============================================================================
# Directories to skip when scanning WES cohort dirs
# ==============================================================================

SKIP_DIRS = {'summary_all', 'logis_batch', 'ridge', 'Transcriptome', 'DNAFormer', 'geneformer'}

# ==============================================================================
# Cohorts
# ==============================================================================

COHORTS = [
    'Allen_Melanoma',
    'Braun_RCC',
    'Cho_NSCLC',
    'Gide_Melanoma',
    'Hugo_Melanoma',
    'IMvigor210_UC',
    'Jung_Melanoma',
    'Kim_GC',
    'Lauss_Melanoma',
    'Liu_Melanoma',
    'Mariathasan_UC',
    'Miao_Melanoma',
    'Miao_RCC',
    'Riaz_Melanoma',
    'Snyder_UC',
    'VanAllen_Melanoma',
]
