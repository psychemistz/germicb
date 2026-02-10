# WES Immunotherapy Germline Variant Pipeline

Upstream preprocessing and association analysis pipeline for whole-exome sequencing
data from immunotherapy-treated cancer cohorts. Produces the therapy-associated variant
tables consumed by the 4germicb AlphaGenome prioritization pipeline.

**Location:** `/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED/`

## Overview

The pipeline processes raw WES FASTQ files through 9 stages: read preprocessing,
consensus variant calling (GATK + DeepVariant + FreeBayes), functional annotation,
and multiple association analyses (logistic regression, survival, eQTL, sQTL). Only
variants called by all 3 callers are retained (consensus intersection).

**Reference genome:** GRCh38/hg38 (`Homo_sapiens_assembly38.fasta`)

```
Stage 01: fastp → bwa → samtools → GATK MarkDuplicates → BQSR
Stage 1:  GATK HaplotypeCaller → GenomicsDBImport → GenotypeGVCFs → VQSR
Stage 2:  DeepVariant (WES mode)
Stage 3:  FreeBayes
Stage 4:  bcftools isec (consensus) → snpEff → SnpSift (annotation)
Stage 5:  plink2 (genotype matrix)
Stage 6a: logis_batch (Firth logistic regression)
Stage 6b: Wilcoxon rank-sum test
Stage 7:  Cox PH (overall survival)
Stage 8:  Cox PH (progression-free survival)
Stage 9a: MatrixeQTL (eQTL)
Stage 9b: MatrixeQTL (sQTL)
```

## Data Sources

**WES data:** `/data/Jiang_Lab/Data/Therapy/WES/{cohort}/`
**RNA-seq data:** `/data/Jiang_Lab/Data/Therapy/RNASeq/{cohort}/`
**Exome intervals:** `/data/parks34/projects/WES_ImmunoPredict/bedfiles/calling_regions.bed`

### 16 Cohorts

| Cohort | Cancer Type | Therapy | Publication |
|--------|-------------|---------|-------------|
| CCRCC_ICB_Miao2018 | Clear cell RCC | ICB | Miao et al. 2018 |
| GBM_PD1_Zhao2019 | Glioblastoma | Anti-PD1 | Zhao et al. 2019 |
| Hepatocellular_Atezolizumab_Finn2020 | HCC | Atezolizumab | Finn et al. 2020 |
| Melanoma_CTLA4_Snyder2014 | Melanoma | Anti-CTLA4 | Snyder et al. 2014 |
| Melanoma_ICB_Roh2017 | Melanoma | ICB (mixed) | Roh et al. 2017 |
| Melanoma_Ipilimumab_VanAllen2015 | Melanoma | Ipilimumab | Van Allen et al. 2015 |
| Melanoma_Nivolumab_Riaz2017 | Melanoma | Nivolumab | Riaz et al. 2017 |
| Melanoma_PD1_Hugo2016 | Melanoma | Anti-PD1 | Hugo et al. 2016 |
| Melanoma_PD1_Liu2019 | Melanoma | Anti-PD1 | Liu et al. 2019 |
| mGC_Pembrolizumab_Kim2018 | Gastric cancer | Pembrolizumab | Kim et al. 2018 |
| mRCC_Atezolizumab_McDermott2018 | Metastatic RCC | Atezolizumab | McDermott et al. 2018 |
| NSCLC_ICB_Ravi2023 | NSCLC | ICB | Ravi et al. 2023 |
| NSCLC_Pembrolizumab_Rizvi2015 | NSCLC | Pembrolizumab | Rizvi et al. 2015 |
| PanCancer_PD1_Cristescu2018 | Pan-cancer | Anti-PD1 | Cristescu et al. 2018 |
| PanCancer_Pembrolizumab_Yang2021 | Pan-cancer | Pembrolizumab | Yang et al. 2021 |
| Urothelial_Atezolizumab_Mariathasan2018 | Urothelial | Atezolizumab | Mariathasan et al. 2018 |

## Stage 01: Read Preprocessing

**Tools:** fastp 0.23.2, bwa 0.7.17, samtools 1.17, GATK 4.4.0.0

### 1. Read trimming (fastp)

```
fastp -i R1.fastq.gz -I R2.fastq.gz -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz \
  -j sample_fastp.json -h sample_fastp.html -w 2
```

### 2. Alignment (bwa mem)

```
bwa mem -v 2 -M -t 32 -p -R "@RG\tID:sample\tSM:sample\tPL:illumina\tLB:lib1" \
  Homo_sapiens_assembly38.fasta trimmed_R1.fq.gz trimmed_R2.fq.gz \
  | samtools sort -@ 12 -m 1707M -O BAM -o sample.bam
```

- Reference: `/fdb/GATK_resource_bundle/hg38-v0/Homo_sapiens_assembly38.fasta`
- 32 threads for alignment, 12 threads for sorting

### 3. Mark duplicates (GATK MarkDuplicatesSpark)

```
gatk --java-options "-Xms60G -Xmx60G -XX:ParallelGCThreads=2" \
  MarkDuplicatesSpark -I sample.bam -O sample_markdup.bam
```

### 4. Base quality score recalibration (BQSR)

```
gatk BaseRecalibrator \
  -I sample_markdup.bam \
  -R Homo_sapiens_assembly38.fasta \
  --known-sites dbsnp_146.hg38.vcf.gz \
  --known-sites Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -ip 50 -O recal_data.table

gatk ApplyBQSR \
  -I sample_markdup.bam -bqsr recal_data.table \
  -O sample_markdup_bqsr.bam
```

**Known sites:**
- dbSNP 146 (hg38)
- Known indels (Homo sapiens assembly 38)
- Mills & 1000G gold standard indels (hg38)

**Output:** `01preprocess/sample_markdup_bqsr.bam` (analysis-ready BAM)

**SLURM:** 34 CPUs, 65 GB RAM, 14h, 250 GB scratch

## Stage 1: GATK HaplotypeCaller Pipeline

**Tool:** GATK 4.4.0.0

### Per-sample GVCF generation

```
gatk HaplotypeCaller \
  -I sample_markdup_bqsr.bam \
  -R Homo_sapiens_assembly38.fasta \
  -O sample.g.vcf.gz \
  -ERC GVCF \
  -L calling_regions.bed -ip 50
```

### Joint genotyping (per chromosome)

```
gatk GenomicsDBImport \
  --batch-size 50 \
  --genomicsdb-shared-posixfs-optimizations true \
  --max-num-intervals-to-import-in-parallel 3 \
  -V sample1.g.vcf.gz -V sample2.g.vcf.gz ... \
  --genomicsdb-workspace-path chr1_gdb/ -L chr1

gatk GenotypeGVCFs \
  -R Homo_sapiens_assembly38.fasta \
  -V gendb://chr1_gdb/ -O chr1.vcf.gz
```

Per-chromosome VCFs merged with `picard GatherVcfs` into `merged.vcf.gz`.

### VQSR

**SNP recalibration:**

```
gatk VariantRecalibrator \
  -R Homo_sapiens_assembly38.fasta -V merged.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  --tranches-file output_SNP1.tranches \
  --rscript-file output_SNP1.plots.R \
  -mode SNP -O merged_SNP1.recal
```

**INDEL recalibration:**

```
gatk VariantRecalibrator \
  ... \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  -mode INDEL -O merged_indel1.recal
```

**ApplyVQSR:** Truth sensitivity threshold 99.9%

**Output:** `1GATK/vcf/indel.SNP.recalibrated_99.9.vcf.gz`

## Stage 2: DeepVariant

**Tool:** DeepVariant 1.5.0

```
run_deepvariant \
  --model_type=WES \
  --ref=Homo_sapiens_assembly38.fasta \
  --reads=sample_markdup_bqsr.bam \
  --output_vcf=sample.final.vcf.gz \
  --output_gvcf=sample.final.g.vcf.gz \
  --regions=calling_regions.bed \
  --num_shards=$(nproc) \
  --dry_run=false
```

- SLURM array job: 1 job per sample, 16 CPUs, 64 GB, 24h
- Joint calling: GLnexus (`--config DeepVariant`)

**Output:** `2DeepVariant/dv_all_variants_w_bed.vcf.gz`

## Stage 3: FreeBayes

**Tool:** FreeBayes (v1.3+)

```
freebayes \
  --fasta-reference Homo_sapiens_assembly38.fasta \
  --min-alternate-count 2 \
  --min-alternate-qsum 40 \
  --pvar 0.0001 \
  --use-mapping-quality \
  --site-selection-max-iterations 3 \
  --genotyping-max-iterations 25 \
  --gvcf \
  sample_markdup_bqsr.bam | bgzip -c > sample.vcf.gz
```

- SLURM array job: 1 job per sample, 16 CPUs, 64 GB, 24h
- Merged: `bcftools concat` into `3FreeBayes/merged.vcf.gz`

**Output:** `3FreeBayes/merged.vcf.gz`

## Stage 4: Consensus Variant Calling + Annotation

### Consensus intersection (require all 3 callers)

```
bcftools isec -n=3 \
  1GATK/vcf/indel.SNP.recalibrated_99.9.vcf.gz \
  2DeepVariant/dv_all_variants_w_bed.vcf.gz \
  3FreeBayes/merged.vcf.gz \
  -o output_snps.vcf.gz -o output_indels.vcf.gz
```

Only variants called by GATK AND DeepVariant AND FreeBayes are retained.

### Functional annotation

```
# snpEff: gene-level annotation
java -Xmx4096m -jar snpEff.jar ann -v -noLog -noStats -noLof GRCh38.105 \
  consensus.vcf > consensus_eff.vcf

# SnpSift: dbSNP rsIDs
java -Xmx4096m -jar SnpSift.jar annotate -v -id \
  /fdb/dbSNP/organisms/human_9606_b150_GRCh38p7/00-All.vcf.gz \
  consensus_eff.vcf > consensus_eff_dbsnp.vcf

# SnpSift: dbNSFP functional predictions
java -Xmx4096m -jar SnpSift.jar dbnsfp -v -m \
  -db /fdb/dbNSFP/dbNSFP3.3a.txt.gz \
  consensus_eff_dbsnp.vcf > consensus_eff_dbsnp_dbnsfp.vcf
```

**Annotation databases:**
- snpEff: GRCh38.105
- dbSNP: v150 (b150, GRCh38p7)
- dbNSFP: v3.3a

**Output:** `4merge_vcf/output_snps.vcf.gz`, `4merge_vcf/output_indels.vcf.gz`

## Stage 5: PLINK Genotype Matrix

**Tool:** plink2

```
# Filter to biallelic SNPs
bcftools view -m2 -M2 -v snps ../4merge_vcf/output_snps.vcf.gz > snps.vcf

# Convert to PLINK binary format
plink2 --vcf snps.vcf --make-bed --out ./input_data/snps
```

Phenotype and covariate files prepared from clinical data:
- `clin_prep/pheno.txt` (treatment response: 0/1)
- `clin_prep/cov_updated.txt` (Age, Gender)

**Output:** `5plink/input_data/snps.bed`, `.bim`, `.fam`

## Stage 6: Association Analysis

The choice of association method depends on the outcome type available for each cohort:
- **Stage 6a (logis_batch):** Cohorts with discrete binary response (responder/non-responder)
- **Stage 6b (rank-sum) / ridge:** Cohorts with continuous outcomes (e.g., RECIST measurements, tumor burden) where binary classification is not appropriate

### Stage 6a: Logistic Regression (Treatment Response)

**Tool:** logis_batch (custom C++, Peng Jiang 2014, data2intelligence/data_significance)

**Model:** `logit(Response) ~ Age + Gender + SNP_i` with Firth bias correction
**Used for:** Cohorts with discrete binary outcome (responder vs non-responder)

```
logis_batch \
  -B background_matrix \
  -X genotype_matrix \
  -Y phenotype_matrix \
  -out output \
  -maxiter 1000 \
  -tol 1e-5 \
  -delta 1.0 \
  -correction 1 \
  -cntthres 5 \
  -filterX 1 \
  -intercept 1 \
  -background 1
```

**Key parameters:**
- Newton-Raphson with Firth-penalized likelihood
- Convergence: tolerance 1e-5, max 1000 iterations
- MAF filter: `rowSums(gt==2) >= max(5, ceil(N/10))`
- Firth correction reduces bias in small-sample/rare-event logistic regression
- Implementation: Fisher information matrix via X'WX, Cholesky decomposition

**Output:** `6logis_batch/output/fin_tbl.txt` (SNP, gene, beta, SE, zscore, pvalue)

### Stage 6b: Wilcoxon Rank-Sum Test

**Tool:** R (Wilcoxon rank-sum / Mann-Whitney U)

**Used for:** Cohorts with continuous outcomes where logistic regression on a binary
endpoint is not applicable (GBM, mGC, mRCC, PanCancer).

**MAF filter:** `rowSums(gt==2) >= ceil(N/4)` (more aggressive than logis_batch;
appropriate since rank-sum has less power and benefits from higher MAF)

**Output:** `6wilcox/output/fin_tbl.txt`

## Stage 7: Overall Survival (Cox PH)

**Tool:** R survival package (`coxph`)

**Model:** `Surv(OS, OS.Event) ~ Age + SNP_i`

```r
fit <- coxph(Surv(OS, OS.Event) ~ Age + variant, data = surv_data)
```

Per-variant Cox proportional hazards model. Gender is NOT included as a covariate
(unlike Stage 6a logistic regression — see Issues).

**Output:** `7OS/output/` (per-variant hazard ratios and p-values)

## Stage 8: Progression-Free Survival (Cox PH)

**Tool:** R survival package (`coxph`)

**Model:** `Surv(PFS, PFS.Event) ~ Age + SNP_i`

Same structure as Stage 7. Available for cohorts with PFS data (Melanoma, HCC, mRCC,
NSCLC).

**Output:** `8PFS/output/`

## Stage 9a: eQTL Analysis

**Tool:** MatrixeQTL (R package)

```r
results <- Matrix_eQTL_main(
  snps = SlicedData$new(genotype_matrix),
  gene = SlicedData$new(expression_matrix),
  cvrt = SlicedData$new(covariate_matrix),  # Age, Gender, PCs
  snpspos = snp_annot,
  genepos = gene_annot
)
```

**Model:** `Expression ~ Age + Gender + SNP`
**Correction:** Bonferroni (FDR column in output)

**Output:** `9eQTL/output/fin_eqtl_tbl.txt` (SNP, gene, beta, statistic, pvalue, FDR)

**Note:** The `snp_annot.txt` file used by MatrixeQTL contains **hg19 coordinates**
despite the variant calling being done in hg38. This requires liftover when matching
against hg38 reference databases. See Issues #7.

## Stage 9b: sQTL Analysis

**Tool:** MatrixeQTL (adapted for splicing phenotypes)

**Model:** `SpliceScore ~ Age + Gender + SNP`
**Correction:** Bonferroni

Same framework as eQTL but with splice junction ratios as the phenotype.

**Output:** `9sQTL/output/fin_sqtl_tbl.txt`

## Output Summary

The primary output consumed by the 4germicb pipeline is:

```
{cohort}/6logis_batch/output/fin_tbl.txt
```

Columns: `CHROM POS REF ALT ID GENE coef zscore pvalue`

This file contains association results for all tested variants (no FDR correction
applied — see Issues #1). The 4germicb `01_extract_therapy_variants.py` script
aggregates these across all cohorts and filters at p < 0.05.

## Issues and Limitations

### 1. No multiple testing correction in logistic regression (Stage 6a)

Raw p-values are used in `fin_tbl.txt`. No Bonferroni or FDR correction is applied
at this stage. Users apply their own thresholds downstream (p < 0.05 suggestive,
p < 0.001 stringent in the 4germicb pipeline).

### 2. No multiple testing correction in survival analyses (Stages 7-8)

Same issue as Stage 6a. OS and PFS analyses report raw p-values.

### 3. Different MAF filters across association methods

- **logis_batch (Stage 6a):** `max(5, ceil(N/10))` — moderately stringent
- **Wilcoxon rank-sum (Stage 6b):** `ceil(N/4)` — more aggressive filter

These methods are applied to different cohorts based on outcome type (discrete vs
continuous), not as alternative analyses on the same data. The stricter rank-sum
filter is appropriate because the non-parametric test has less statistical power
and benefits from higher MAF variants. For a cohort of N=40, logis_batch requires
5 carriers while rank-sum requires 10.

### 4. Inconsistent variant_id formats across 4germicb scripts

- `01_extract_therapy_variants.py` creates `chrom_pos` format
- `02_query_alphagenome.py` creates `chrom:pos_ref>alt` format

The `lib/variants.py` module now standardizes to `chrom:pos_ref>alt`.

### 5. Gender not included in survival models (Stages 7-8)

OS and PFS models adjust for Age only (`Surv ~ Age + SNP`), while logistic
regression and eQTL analyses adjust for both Age and Gender. This is an
inconsistency that could confound sex-linked variant associations in survival
outcomes.

### 6. Inf/-Inf test statistics in eQTL output

`fin_eqtl_tbl.txt` contains `-Inf` statistics from quasi-complete separation in
MatrixeQTL linear regression. These arise when a genotype perfectly predicts
expression (typically due to very low MAF). The 4germicb `03_validate_eqtl.py`
handles this by filtering `np.isfinite(beta)`.

### 7. hg19 coordinates in cohort eQTL annotations

The `snp_annot.txt` files used by MatrixeQTL (Stage 9a) contain hg19 coordinates
despite the variant calling being done against the hg38 reference. This is likely
an oversight from pipeline migration. The 4germicb `03_validate_eqtl.py` uses
`pyliftover` (hg19 -> hg38) to correct coordinates before matching against DICE/OneK1K.

### 8. Ridge regression not used in 4germicb

Some cohorts have ridge regression results (`ridge/` directories) with L2-regularized
association analysis (lambda=10,000, 1000 permutations). These results are not
currently consumed by the 4germicb pipeline, which uses only `6logis_batch` output.
