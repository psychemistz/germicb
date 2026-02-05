# AlphaGenome + eQTL Database Validation Plan for Immunotherapy Germline Variants

## Executive Summary

This document outlines a validated approach for prioritizing germline variants associated with immunotherapy outcomes using AlphaGenome predictions and large-scale immune cell eQTL databases (CIMA, DICE, GTEx).

**Key Insight from Cytokine Atlas Project (2026-02-05):** AlphaGenome is most effective for predicting regulatory effects of variants with **unknown** effects, NOT for prioritizing variants that already have measured eQTL data. For immunotherapy-associated germline variants from WES, this is an **ideal use case**.

---

## Lessons Learned from Cytokine Atlas AlphaGenome Analysis

### What Works vs. What Doesn't

| Approach | Use Case | Effectiveness |
|----------|----------|---------------|
| AlphaGenome → predict unknown variant effects | **Immunotherapy WES variants** | ✅ IDEAL |
| AlphaGenome → prioritize known eQTLs | CIMA eQTL prioritization | ❌ Worse than naive |
| Direct eQTL database matching | Validate regulatory predictions | ✅ HIGHLY EFFECTIVE |

### Key Findings

1. **AlphaGenome API returns only 3 RNA-seq tracks** (not chromatin):
   - GTEx Lymphocytes (EBV-transformed)
   - GTEx Spleen
   - GTEx Whole Blood

2. **Direct eQTL validation achieves high concordance**:
   - CIMA vs DICE: 83.5% concordance (r=0.685)
   - CIMA vs GTEx: 68.9% concordance (r=0.370)

3. **DICE provides best immune cell validation** (cell-type-specific > bulk tissue)

4. **Genome build matters**: DICE is hg19, must liftover to hg38

---

## Data Resources Available

### Immunotherapy WES Cohorts (14 cohorts)

From `/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED/`:

| Cancer Type | Cohorts | Key Studies |
|-------------|---------|-------------|
| Melanoma | 6 | Snyder2014, VanAllen2015, Hugo2016, Roh2017, Riaz2017, Liu2019 |
| RCC | 3 | Miao2018, McDermott2018 |
| NSCLC | 2 | Rizvi2015, Ravi2023 |
| Other | 5 | HCC, Gastric, Urothelial, GBM, PanCancer |

**Data Structure per Cohort:**
```
cohort/
├── 4merge_vcf/          # Consensus variants (GATK ∩ DeepVariant ∩ FreeBayes)
│   └── filt_var.txt     # ~20K filtered variants with CADD scores
├── 6logis_batch/        # Therapy response association (Firth logistic)
│   └── output/fin_tbl.txt  # z-scores, coefficients, p-values
├── 7OS/, 8PFS/          # Survival analysis
└── 9eQTL/               # Local eQTL analysis
```

### Immune Cell eQTL Databases

| Database | Cell Types | N eQTLs | Genome Build | Location |
|----------|------------|---------|--------------|----------|
| **CIMA** | 69 immune subtypes | 48,627 (cytokine genes) | hg38 | `/data/Jiang_Lab/Data/Seongyong/CIMA/xQTL/` |
| **DICE** | 12 immune cell types | 3,706,208 | hg19 → **hg38 lifted** | `results/alphagenome/dice_data/hg38/` |
| **GTEx v10** | Whole blood, spleen, lymphocytes | 2,985,690 | hg38 | `results/alphagenome/gtex_data/` |

---

## Proposed Workflow

### Phase 1: Extract Therapy-Associated Variants

```
For each cohort:
    1. Load association results (6logis_batch/output/fin_tbl.txt)
    2. Filter: p < 0.05 (suggestive) or p < 0.001 (stringent)
    3. Extract variant info: chr, pos, ref, alt, gene, beta, pvalue

Meta-analysis:
    1. Combine across cohorts using inverse-variance weighting
    2. Identify variants significant in ≥2 cohorts
    3. Flag variants with consistent effect direction
```

**Expected yield**: ~500-2,000 therapy-associated variants per cohort

### Phase 2: AlphaGenome Regulatory Prediction

```
For therapy-associated variants (unknown regulatory effects):
    1. Format for AlphaGenome API (chr:pos REF>ALT, hg38)
    2. Query AlphaGenome API
    3. Extract immune-relevant track predictions:
       - Lymphocyte RNA-seq
       - Spleen RNA-seq
       - Whole blood RNA-seq
    4. Compute impact score = mean(|track_differences|)
```

**This is the CORRECT use of AlphaGenome** - predicting effects for variants without prior eQTL data.

### Phase 3: Validate Against eQTL Databases

```
For each variant with AlphaGenome predictions:
    1. Match to DICE eQTLs (by hg38 position)
       - If match: variant IS a known immune cell eQTL
       - Record: dice_beta, dice_pval, dice_celltype

    2. Match to GTEx eQTLs (by position)
       - Record: gtex_slope, gtex_pval

    3. Match to CIMA eQTLs (if cytokine/secreted protein gene)
       - Record: cima_beta, cima_pval, cima_celltype
```

### Phase 4: Prioritization Scoring

```python
# Composite prioritization score
priority_score = (
    therapy_association_score +      # -log10(p) from logistic/survival
    alphagenome_impact_score +       # Predicted regulatory effect
    eqtl_validation_score            # Bonus if validated in DICE/GTEx
)

# Tier assignment
Tier 1: Therapy-associated + AlphaGenome-predicted + eQTL-validated
Tier 2: Therapy-associated + AlphaGenome-predicted (no eQTL match)
Tier 3: Therapy-associated + eQTL-validated (no AlphaGenome prediction)
Tier 4: Therapy-associated only
```

---

## Expected Outcomes

### Validation Rates (Based on Cytokine Atlas Results)

| Comparison | Expected Matches | Expected Concordance |
|------------|------------------|----------------------|
| AlphaGenome-predicted → DICE | ~30-40% of variants | 70-80% direction match |
| AlphaGenome-predicted → GTEx | ~20-30% of variants | 60-70% direction match |
| Therapy-associated → any eQTL | ~10-20% of variants | Variable |

### Interpretation

- **High AlphaGenome impact + eQTL validated**: Strong regulatory evidence
- **High AlphaGenome impact + no eQTL**: May be context-specific (tumor/therapy)
- **Low AlphaGenome impact + eQTL validated**: May have regulatory effect AlphaGenome missed
- **No AlphaGenome + no eQTL**: Lowest priority, may be coding or structural

---

## Implementation Scripts

### Script 1: Extract Therapy-Associated Variants

```python
# scripts/01_extract_therapy_variants.py
"""
Extract therapy-associated variants from WES cohorts.
Input: 6logis_batch/output/fin_tbl.txt per cohort
Output: therapy_variants_all.csv
"""

import pandas as pd
from pathlib import Path

WES_DIR = Path('/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED')
COHORTS = [d.name for d in WES_DIR.iterdir() if d.is_dir() and not d.name.startswith(('summary', 'logis', 'ridge'))]

def load_cohort_associations(cohort_dir, pval_threshold=0.05):
    """Load association results for a cohort."""
    fin_tbl = cohort_dir / '6logis_batch/output/fin_tbl.txt'
    if not fin_tbl.exists():
        return None

    df = pd.read_csv(fin_tbl, sep='\t')
    df['cohort'] = cohort_dir.name

    # Filter by p-value
    df = df[df['pvalue'] < pval_threshold]
    return df

# Combine all cohorts
all_variants = []
for cohort in COHORTS:
    df = load_cohort_associations(WES_DIR / cohort)
    if df is not None:
        all_variants.append(df)

combined = pd.concat(all_variants, ignore_index=True)
combined.to_csv('therapy_variants_all.csv', index=False)
```

### Script 2: Query AlphaGenome

```python
# scripts/02_query_alphagenome.py
"""
Query AlphaGenome API for therapy-associated variants.
This is the CORRECT use case - predicting effects for unknown variants.
"""

# Use existing infrastructure from:
# /vf/users/parks34/projects/2secactpy/scripts/08_alphagenome_stage3_predict.py

# Key parameters:
# - Input: therapy_variants_all.csv
# - Output: alphagenome_predictions.h5ad
# - Checkpoint: Every 100 variants for resume
# - Rate limiting: Handle RESOURCE_EXHAUSTED errors
```

### Script 3: Validate Against eQTL Databases

```python
# scripts/03_validate_eqtl.py
"""
Validate AlphaGenome predictions against DICE, GTEx, CIMA.
Based on validated approach from cytokine atlas project.
"""

# Load lifted DICE data (hg38)
DICE_HG38_DIR = Path('/vf/users/parks34/projects/2secactpy/results/alphagenome/dice_data/hg38')

# Load GTEx v10
GTEX_PATH = '/data/parks34/projects/2secactpy/results/alphagenome/gtex_data/GTEx_Analysis_v10_eQTL_updated/Whole_Blood.v10.eQTLs.signif_pairs.parquet'

# Match by position (chr:pos)
def match_to_eqtl(variants_df, eqtl_df):
    """Match variants to eQTL database by position."""
    variants_df['pos_key'] = variants_df['CHROM'] + ':' + variants_df['POS'].astype(str)
    eqtl_df['pos_key'] = eqtl_df['chrom'] + ':' + eqtl_df['pos'].astype(str)

    matched = variants_df.merge(eqtl_df, on='pos_key', how='left')
    return matched
```

### Script 4: Prioritize and Report

```python
# scripts/04_prioritize_variants.py
"""
Compute priority scores and generate final report.
"""

def compute_priority_score(row):
    """Compute composite priority score."""
    score = 0

    # Therapy association (higher = more significant)
    if row['therapy_pval'] < 0.001:
        score += 3
    elif row['therapy_pval'] < 0.01:
        score += 2
    elif row['therapy_pval'] < 0.05:
        score += 1

    # AlphaGenome prediction
    if row['alphagenome_impact'] > 0.001:
        score += 2
    elif row['alphagenome_impact'] > 0:
        score += 1

    # eQTL validation
    if row['dice_matched']:
        score += 2
    if row['gtex_matched']:
        score += 1

    return score
```

---

## Comparison with Previous (Incorrect) Approach

| Aspect | Previous (Cytokine Atlas) | Proposed (Immunotherapy WES) |
|--------|---------------------------|------------------------------|
| Input variants | Known eQTLs (measured effects) | Therapy-associated (unknown effects) |
| AlphaGenome role | Prioritize known effects | **Predict unknown effects** |
| Validation | Circular (AlphaGenome trained on GTEx) | **Independent** (therapy outcomes) |
| Expected value | Low (75% vs 90% naive) | **High** (fills knowledge gap) |

---

## Files to Create

```
/data/parks34/projects/4germicb/
├── RESEARCH.md                         # Existing planning document
├── ALPHAGENOME_EQTL_VALIDATION_PLAN.md # This document
├── scripts/
│   ├── 01_extract_therapy_variants.py
│   ├── 02_query_alphagenome.py
│   ├── 03_validate_eqtl.py
│   └── 04_prioritize_variants.py
├── data/
│   ├── therapy_variants_all.csv
│   └── cohort_metadata.csv
└── results/
    ├── alphagenome_predictions.h5ad
    ├── eqtl_validation.csv
    └── prioritized_variants.csv
```

---

## Key Differences from Cytokine Atlas Approach

### Why This Will Work Better

1. **Correct use case**: Predicting effects for variants WITHOUT prior eQTL data
2. **Independent validation**: Therapy outcomes are independent of eQTL databases
3. **Clear value-add**: AlphaGenome fills the gap between association and mechanism
4. **Tiered evidence**: Variants validated in eQTLs have stronger regulatory evidence

### Expected Success Metrics

| Metric | Target |
|--------|--------|
| Therapy variants with AlphaGenome predictions | >80% |
| AlphaGenome-predicted variants matching DICE | 30-40% |
| Direction concordance (AlphaGenome vs DICE) | >70% |
| Tier 1 variants (therapy + AlphaGenome + eQTL) | 5-10% of input |

---

## Conclusion

This approach leverages AlphaGenome's true strength: **predicting regulatory effects for variants with unknown mechanisms**. By combining:

1. **Therapy association** (clinical relevance)
2. **AlphaGenome predictions** (regulatory potential)
3. **eQTL validation** (biological evidence)

We can identify high-confidence regulatory variants driving immunotherapy outcomes across independent cohorts.

---

## References

- AlphaGenome: DeepMind, January 2026
- DICE database: Schmiedel et al., Cell 2018
- GTEx v10: GTEx Consortium 2022
- CIMA: Single-cell immune atlas eQTLs
- Cytokine Atlas AlphaGenome Analysis: `/vf/users/parks34/projects/2secactpy/docs/pipelines/alphagenome_plan.md`
