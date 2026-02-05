# Integrating foundational genomics models with underpowered immunotherapy GWAS

> **UPDATE 2026-02-05**: See `ALPHAGENOME_EQTL_VALIDATION_PLAN.md` for validated implementation approach based on lessons from Cytokine Atlas AlphaGenome analysis. Key insight: AlphaGenome is most effective for variants with **unknown** regulatory effects (ideal for therapy-associated WES variants), NOT for prioritizing variants with known eQTL effects.

**AlphaGenome and related deep learning models can meaningfully prioritize your suggestive germline variants (p~10⁻⁵), but statistical power limitations with N=1,640 require combining meta-analysis, Bayesian fine-mapping with functional priors, and polygenic approaches to maximize discovery.** The four variants you've identified (rs3783947/TSHR, rs737321/PARP12, rs909723/HLA-F-AS1, rs3806268/NLRP3) represent novel candidates—none appear in published immunotherapy association studies—making functional prioritization particularly valuable. Recent advances, including the GeRI study's validated polygenic risk score for autoimmune diseases predicting ICI discontinuation (HR=1.24), demonstrate that germline genetics can predict immunotherapy outcomes even from moderately-sized cohorts.

---

## AlphaGenome represents a major advance for regulatory variant interpretation

Google DeepMind released AlphaGenome in January 2026, establishing it as the most comprehensive sequence-based model for variant effect prediction. Built on a transformer architecture accepting **1 million base pairs** of DNA sequence with single-nucleotide resolution, AlphaGenome predicts gene expression, splicing, chromatin accessibility, histone modifications, transcription factor binding, and 3D chromatin contacts simultaneously. It achieves state-of-the-art performance on **25 of 26 variant effect prediction benchmarks** and outperforms SpliceAI on 6 of 7 splicing tasks.

For your GWAS prioritization, AlphaGenome's key capability is scoring individual germline variants by comparing predictions between reference and alternate sequences. The model is accessible via a free API for non-commercial research (suitable for thousands of variants, though not millions). Each prediction requires approximately 1 second on NVIDIA H100 hardware. AlphaGenome's primary limitations include challenges modeling ultra-distant regulatory interactions (>100kb), inability to predict effects in cell types absent from training data, and API rate limits during high-demand periods.

**Enformer** (200kb context, 128bp resolution) remains valuable as a complementary tool with pre-computed scores available for all common 1000 Genomes variants, enabling rapid genome-wide screening. **Borzoi** (524kb context, 32bp resolution) specifically predicts RNA-seq coverage including splicing patterns. **SpliceAI** provides focused splicing predictions with delta scores interpretable as splice-altering probabilities—variants with delta scores >0.5 validate at 75% the rate of canonical splice site disruptions.

For integrative scoring, **CADD v1.7** combines 60+ annotations including protein language model scores (ESM-1v), with pre-computed PHRED-scaled scores available for all 9 billion possible SNVs. **Sei** uniquely classifies variants into 40 interpretable regulatory "sequence classes" (stem cell enhancers, brain enhancers, B-cell enhancers, etc.), enabling tissue-specific interpretation relevant to immune cell biology.

---

## Statistical approaches must address both small samples and cancer heterogeneity

Your cohort structure—23 studies with N ranging from 13-216 patients, heterogeneous cancer types, and case-control ratio of 494:1,146—creates specific analytical challenges requiring targeted methodological choices.

**Meta-analysis strategy**: Fixed-effects meta-analysis with METAL maximizes power for discovery when assuming shared genetic effects across cancer types (appropriate for immune-intrinsic germline variants). However, given cancer-type heterogeneity, random-effects analysis with GWAMA should serve as sensitivity analysis, reporting Cochran's Q and I² statistics for top signals. For your 23 cohorts, MR-MEGA's trans-ethnic meta-regression approach could model cancer type as a source of heterogeneity, though this requires at least 9 studies with sufficient sample sizes.

**Power reality**: Your effective sample size (N_eff = 4 × 494 × 1,146 / 1,640 ≈ 1,380) enables detection of odds ratios ≥**1.7-2.0** for variants with MAF 20-30% at genome-wide significance (5×10⁻⁸). At your observed suggestive threshold (p~10⁻⁵), detectable effect sizes are approximately 15% smaller. This power limitation makes functional prioritization of sub-threshold signals essential rather than optional.

**Firth's penalized regression** is critical for your individual cohort analyses, particularly for cohorts with N<100 or when analyzing lower-frequency variants (MAF 1-5%). Firth's method adds Jeffreys invariant prior to the likelihood, producing finite estimates even with quasi-complete separation (when a genotype class perfectly predicts response). Implementation options include:

- **REGENIE**: Uses approximate Firth correction (60x faster than exact) when p<0.05, optimal for your setting
- **SAIGE**: Combines saddlepoint approximation with Firth correction for case-control imbalance
- **R package `logistf`**: Exact Firth regression for detailed single-cohort analyses

For individual cohort analysis, REGENIE's two-step approach (ridge regression on genotype blocks followed by leave-one-chromosome-out prediction) handles relatedness and population structure while enabling Firth correction for rare variants and small samples.

---

## Only one genome-wide significant germline variant has been replicated for ICI outcomes

The field's most robust finding emerged from Groha et al.'s 2022 Nature Medicine GWAS: **rs16906115 near IL7** achieved combined P = 3.6×10⁻¹¹ (HR = 2.1 for all-grade immune-related adverse events) and replicated across three independent cohorts totaling 1,751 patients. Mechanistically, this variant creates a cryptic exon affecting lymphocyte homeostasis—carriers exhibited increased lymphocyte stability after ICI initiation, correlating with both irAEs and improved survival, suggesting shared genetic architecture between toxicity and response.

**HLA associations** show consistent effects: HLA-I heterozygosity at the A/B/C loci improves overall survival (HR = 1.40, P = 0.036 in 1,535 patients at MSK). Specific alleles show directional effects—**HLA-B44 supertype** associates with extended melanoma survival while **HLA-B62 supertype** and **HLA-A*03** associate with poor outcomes. HLA-B27 supertype correlates with grade 3 hepatitis and pneumonitis (67% vs 2.7% in non-carriers, q=0.007).

**Your candidate variants have not been published in ICI studies.** Searches for rs3783947 (TSHR), rs737321 (PARP12), rs909723 (HLA-F-AS1), and rs3806268 (NLRP3) in immunotherapy contexts yielded no specific associations—these represent potentially novel candidates. However, their biological pathways are mechanistically relevant:

- **NLRP3**: PD-L1/NLRP3 inflammasome signaling drives adaptive resistance to anti-PD-1; NLRP3 inhibition enhances efficacy
- **TSHR**: Thyroid irAEs occur in ~10% of ICI patients; genetic predisposition to thyroid autoimmunity shapes checkpoint blockade response
- **HLA-F-AS1**: Located in the HLA region with established importance, though non-classical HLA molecules are less characterized
- **PARP12**: PARP inhibitors induce immunogenic cell death and cGAS/STING pathway activation

---

## Functional prioritization workflow for sub-genome-wide-significant variants

A systematic variant prioritization pipeline integrating deep learning predictions with functional genomics data can substantially improve the signal-to-noise ratio in underpowered studies.

**Step 1—Score all suggestive variants (p<5×10⁻⁵) with multiple models:**

| Model | Output | Threshold for prioritization |
|-------|--------|------------------------------|
| AlphaGenome API | Multi-modal regulatory scores | Top 1% across predicted tracks |
| SpliceAI | Delta scores (0-1) | >0.2 (sensitive) or >0.5 (specific) |
| CADD v1.7 | PHRED score | >20 (top 1% genome-wide) |
| RegulomeDB v2 | Category + probability | Categories 1-2 (strong evidence) |
| Sei | Sequence class assignment | Immune-relevant classes (E4/B-cell, E5/monocyte) |

**Step 2—Colocalization with immune cell eQTLs:** Use COLOC or eCAVIAR to test whether GWAS signals share causal variants with eQTLs from immune-relevant tissues. Critical resources include:

- **DICE database**: eQTLs across 15 immune cell types (CD4+/CD8+ T cells, Tregs, NK cells, monocytes)—41% of immune eGenes show cell-type specificity
- **GTEx v8**: Whole blood, spleen, EBV-transformed lymphocytes
- **Blueprint Epigenome**: Primary immune cell types with ~38-41% colocalization rates for immune disease loci

Colocalization posterior probability H4 ≥0.75 indicates strong evidence for shared causal variant.

**Step 3—Functionally-informed fine-mapping with PolyFun:** This critical step uses stratified LD-score regression across 75+ functional annotations to compute per-SNP prior causal probabilities, then performs Bayesian fine-mapping with SuSiE or FINEMAP. Studies show **32% improvement** in identifying causal variants (PIP>0.95) compared to non-functional approaches. Deep learning scores can be incorporated as additional annotations in the S-LDSC framework.

**Step 4—Pathway enrichment testing:** Use MAGMA or E-MAGMA (eQTL-informed) for gene-based analysis aggregating variant associations to gene level. Test enrichment in immune pathways (Reactome immune system, KEGG cytokine signaling, MSigDB immunologic signatures). Significant enrichment in biologically relevant pathways strengthens confidence in sub-threshold signals.

---

## Polygenic risk scores offer an alternative discovery paradigm

The GeRI study (1,302 NSCLC patients, 2025) validated that **polygenic risk scores for autoimmune diseases (PRS_AD)** predict ICI discontinuation due to irAEs: HR = 1.24 per standard deviation (P = 0.004), with top-quintile patients showing 4.8% vs 2% discontinuation by 3 months. For combination PD-1/CTLA-4 therapy, the difference was 30.8% vs 14.3% discontinuation.

This approach circumvents single-variant power limitations by aggregating thousands of small-effect variants. For your cohort:

**PRS construction strategy:**
1. **Use PRS-CS with pseudovalidation**—performs only 3% worse than cross-validation when no external validation sample exists
2. **Calculate autoimmune disease PRSs** from published GWAS (rheumatoid arthritis, type 1 diabetes, hypothyroidism, psoriasis) as proxies for immune responsiveness
3. **Test association with response/OS/PFS** adjusting for standard covariates and cancer type

**Method selection based on benchmarking:**
- **LDpred2**: Highest overall accuracy (65.1% AUC average), best for autoimmune traits
- **PRS-CS**: Most robust without validation data
- **PUMAS-ensemble**: Combining multiple methods consistently improves prediction

**Cross-cancer considerations:** If genetic effects differ by cancer type, multi-trait methods (MTAG, mtPGS) can borrow information across cancer-type-stratified GWAS while modeling heterogeneity explicitly.

---

## Validation strategies without additional cohorts

Given the challenge of obtaining independent replication cohorts, several internal validation approaches provide robust evidence:

**Leave-one-cohort-out meta-analysis**: For each of your 23 cohorts, compute meta-analysis excluding that cohort and test whether the excluded cohort shows consistent effect direction and nominal significance. The R package **MetaSubtract** enables this analytically without re-running full meta-analyses. Consistent signals across leave-one-out analyses indicate robustness rather than single-cohort artifacts.

**MAMBA (Model-based Assessment of Meta-analysis replicaBility)**: Calculates replication probability based on consistency across datasets without requiring external validation—directly applicable to your multi-cohort structure.

**K-fold cross-validation for PRS**: Split patients into 10 folds, train PRS weights on 9 folds, evaluate prediction in held-out fold. Bootstrap bias-corrected confidence intervals provide unbiased performance estimates.

**Permutation testing for top signals**: Generate empirical null distributions by permuting phenotypes 10,000+ times. This controls type I error without parametric assumptions, particularly valuable for sub-threshold variants where asymptotic p-values may be unreliable.

---

## Recommended analytical pipeline integrating all approaches

**Phase 1—Quality control and single-cohort analysis:**
```
For each cohort:
  1. Sample QC: call rate >98%, sex check, heterozygosity outliers, relatedness (KING)
  2. Variant QC: MAF >1%, genotyping rate >95%, HWE P>10⁻⁶
  3. Population stratification: PCA projection onto 1000 Genomes, ancestry assignment
  4. Association: REGENIE with Firth correction, covariates = age, sex, PC1-10, cancer type
```

**Phase 2—Meta-analysis:**
```
1. Primary: Fixed-effects (METAL) inverse-variance weighted
2. Sensitivity: Random-effects (GWAMA) with heterogeneity metrics
3. Stratified: Cancer-type groups if sample sizes permit (melanoma/NSCLC vs. others)
4. Leave-one-out: MetaSubtract for internal validation
```

**Phase 3—Functional prioritization of suggestive signals (p<5×10⁻⁵):**
```
1. Annotate with FUMA (automated eQTL mapping, chromatin states, MAGMA)
2. Score with AlphaGenome API, SpliceAI, CADD
3. Colocalize with DICE/GTEx immune cell eQTLs (COLOC)
4. Fine-map with PolyFun + SuSiE using functional priors
5. Rank by composite score: PIP × (eQTL_coloc + DL_score + pathway_membership)
```

**Phase 4—Polygenic analysis:**
```
1. Calculate autoimmune disease PRSs using PRS-CS
2. Test association with response/survival in pooled and stratified analyses
3. Validate with 10-fold CV and leave-one-cohort-out
4. Compare predictive performance to clinical variables alone
```

**Phase 5—Biological interpretation:**
```
1. Prioritize variants with: PIP >0.1 AND (eQTL coloc OR CADD >20 OR SpliceAI >0.2)
2. Map to genes via positional + eQTL + chromatin interaction evidence
3. Test pathway enrichment (MAGMA gene-set analysis)
4. Cross-reference with drug-gene interactions (DGIdb)
```

---

## Essential software and resources

| Analysis | Tool | Access |
|----------|------|--------|
| GWAS association | REGENIE v3.3+ | github.com/rgcgithub/regenie |
| Meta-analysis | METAL, GWAMA | csg.sph.umich.edu/abecasis/metal |
| Fine-mapping | PolyFun + SuSiE | github.com/omerwe/polyfun |
| Colocalization | COLOC, eCAVIAR | CRAN, github.com/cew54/coloc |
| eQTL integration | S-PrediXcan | github.com/hakyimlab/MetaXcan |
| Functional annotation | FUMA | fuma.ctglab.nl |
| Deep learning scores | AlphaGenome API | deepmind.google/science/alphagenome |
| Splicing prediction | SpliceAI | github.com/Illumina/SpliceAI |
| Integrative scoring | CADD v1.7 | cadd.gs.washington.edu |
| PRS construction | PRS-CS, LDpred2 | github.com/getian107/PRScs |
| Immune eQTLs | DICE | dice-database.org |
| Validation | MetaSubtract, MAMBA | CRAN |

---

## Conclusion

Your 23 immunotherapy cohorts with N=1,640 represent a valuable resource that, while underpowered for conventional genome-wide discovery, can yield meaningful insights through the analytical framework outlined here. The integration of AlphaGenome and complementary deep learning models with Bayesian fine-mapping using immune-cell-specific functional priors provides a principled approach to prioritizing your suggestive variants. Your four candidate SNPs (TSHR, PARP12, HLA-F-AS1, NLRP3) occupy biologically plausible pathways but lack prior validation—making them genuinely novel candidates worth pursuing through this prioritization framework.

The most actionable near-term analysis may be constructing autoimmune disease PRSs to test whether polygenic immune predisposition predicts response in your cohort, following the GeRI study's validated approach. Simultaneously, functionally prioritizing your existing suggestive signals using the composite scoring framework will identify the strongest candidates for mechanistic follow-up. Cross-cancer meta-analysis maximizes power for immune-intrinsic variants, while cancer-type-stratified analyses and heterogeneity assessment will reveal any tumor-context-dependent effects.

---

## Addendum: Lessons from Cytokine Atlas AlphaGenome Analysis (2026-02-05)

### Critical Finding: Correct vs Incorrect Use of AlphaGenome

Analysis of AlphaGenome usage in the CIMA Cytokine Atlas eQTL project revealed important limitations:

| Approach | DICE Concordance | Verdict |
|----------|------------------|---------|
| Top by \|eQTL beta\| | **89.9%** | Best |
| Top by eQTL p-value | **89.2%** | Best |
| Random sample | 82.5% | Baseline |
| **AlphaGenome-filtered** | **75.1%** | **Worst** |

**Root cause**: AlphaGenome was used to prioritize variants that **already had measured eQTL effects**. For such variants, the eQTL statistics (beta, p-value) are more informative than AlphaGenome predictions.

### When AlphaGenome IS Effective

AlphaGenome excels at predicting regulatory effects for variants with **unknown** effects:

| Use Case | AlphaGenome Value |
|----------|-------------------|
| GWAS variants without eQTL data | ✅ High |
| Therapy-associated WES variants | ✅ High (this project) |
| Fine-mapping causal variants | ✅ High |
| Prioritizing known eQTLs | ❌ Low (use eQTL statistics instead) |

### AlphaGenome API Limitations

- Returns **only 3 RNA-seq tracks** (GTEx lymphocytes, spleen, whole blood)
- **No chromatin tracks** (ATAC, H3K27ac, etc.) for immune filtering
- ~50% RESOURCE_EXHAUSTED rate limiting errors
- GTEx validation may be circular (AlphaGenome trained on GTEx)

### Recommended Approach for This Project

1. **Extract therapy-associated variants** from WES cohorts (unknown regulatory effects)
2. **Query AlphaGenome** to predict regulatory impact
3. **Validate against DICE/GTEx** to identify variants that ARE immune cell eQTLs
4. **Prioritize by composite score**: therapy_association + alphagenome_impact + eqtl_validation

This workflow uses AlphaGenome for its intended purpose: **predicting effects of unknown variants**, then validating predictions against independent eQTL databases.

### Implementation

See `ALPHAGENOME_EQTL_VALIDATION_PLAN.md` for detailed implementation plan and scripts.

### Key eQTL Validation Resources

| Database | Cell Types | Genome Build | Expected Concordance |
|----------|------------|--------------|----------------------|
| DICE | 12 immune types | hg38 (lifted) | ~75-80% |
| GTEx v10 | Bulk blood | hg38 | ~65-70% |
| CIMA | 69 immune subtypes | hg38 | ~80-85% |

Direct eQTL matching provides stronger validation than AlphaGenome filtering when variants already have measured effects.
