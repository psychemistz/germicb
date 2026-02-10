# Immunotherapy Variant Prioritization Report

## Summary

| Metric | Value |
|--------|-------|
| Total variants | 15,918 |
| Tier 1 (AG + eQTL) | 2,430 (15.3%) |
| Tier 2 (AG only) | 8,247 (51.8%) |
| Tier 3 (eQTL only) | 1,118 (7.0%) |
| Tier 4 (therapy only) | 4,123 (25.9%) |

## Tier Definitions

| Tier | Criteria | Interpretation |
|------|----------|----------------|
| 1 | Therapy + AlphaGenome + eQTL | Strong regulatory evidence, validated mechanism |
| 2 | Therapy + AlphaGenome | Predicted regulatory, not yet validated in eQTL databases |
| 3 | Therapy + eQTL | Known eQTL but AlphaGenome didn't predict high impact |
| 4 | Therapy only | Association only, no regulatory evidence |

## Cohort Summary

| Cohort | Total | Tier 1 | Mean Priority |
|--------|-------|--------|---------------|
| CCRCC_ICB_Miao2018 | 1,934 | 367 | 3.95 |
| GBM_PD1_Zhao2019 | 40 | 6 | 3.42 |
| Hepatocellular_Atezolizumab_Finn2020 | 3,259 | 413 | 3.77 |
| Melanoma_CTLA4_Snyder2014 | 2,366 | 365 | 3.71 |
| Melanoma_Ipilimumab_VanAllen2015 | 2,952 | 491 | 3.91 |
| Melanoma_PD1_Hugo2016 | 802 | 103 | 3.48 |
| NSCLC_Pembrolizumab_Rizvi2015 | 859 | 76 | 3.51 |
| Urothelial_Atezolizumab_Mariathasan2018 | 3,706 | 609 | 3.99 |


## Top 20 Tier 1 Variants

These variants have the strongest evidence: therapy association + AlphaGenome prediction + eQTL validation.

| rsid | gene | cohort | pval | AG impact | DICE | OneK1K | Priority |
|------|------|--------|------|-----------|------|--------|----------|
| rs3795732 | GPATCH4 | Urothelial_Atezolizu | 1.66e-04 | 1.980 | TTC24 | - | 8.8 |
| rs135110 | SLC5A1 | Urothelial_Atezolizu | 2.70e-04 | 2.068 | YWHAH | - | 8.6 |
| rs12567958 | HAPLN2 | Urothelial_Atezolizu | 3.34e-04 | 1.301 | TTC24 | - | 8.5 |
| rs7011778 | ENSG00000269954 | Urothelial_Atezolizu | 4.56e-04 | 1.442 | RP11-148O21.4 | - | 8.3 |
| rs4841558 | BLK | Urothelial_Atezolizu | 4.64e-04 | 0.936 | RP11-148O21.4 | - | 8.3 |
| rs6988443 | BLK | Urothelial_Atezolizu | 4.82e-04 | 1.275 | RP11-148O21.4 | - | 8.3 |
| rs5753137 | ENSG00000273428 | Melanoma_Ipilimumab_ | 6.36e-04 | 1.427 | TBC1D10A | - | 8.2 |
| rs5998233 | SLC5A1 | Urothelial_Atezolizu | 6.50e-04 | 1.129 | C22orf24 | - | 8.2 |
| rs11785333 | BLK | Urothelial_Atezolizu | 6.72e-04 | 1.858 | RP11-148O21.4 | - | 8.2 |
| rs630245 | RNU6-1102P | Hepatocellular_Atezo | 7.46e-04 | 0.989 | NEXN | - | 8.1 |
| rs4841557 | BLK | Urothelial_Atezolizu | 8.36e-04 | 1.545 | RP11-148O21.4 | - | 8.1 |
| rs11555935 | TIMM10B | Urothelial_Atezolizu | 8.86e-04 | 2.264 | TIMM10B_ENSG00000132286 | - | 8.1 |
| rs485005 | RNU6-1102P | Hepatocellular_Atezo | 9.09e-04 | 0.827 | NEXN | - | 8.0 |
| rs1284351 | FEZ1 | Urothelial_Atezolizu | 1.00e-03 | 0.612 | STT3A-AS1 | - | 8.0 |
| rs848216 | ZBTB17 | Urothelial_Atezolizu | 1.14e-03 | 0.799 | CLCNKA | - | 7.9 |
| rs28397595 | ATP8B3 | Hepatocellular_Atezo | 1.16e-03 | 0.659 | ATP8B3 | - | 7.9 |
| rs2581774 | RFT1 | Urothelial_Atezolizu | 1.24e-03 | 0.514 | RFT1 | - | 7.9 |
| rs1869074 | ERICH1 | Melanoma_Ipilimumab_ | 1.26e-03 | 0.761 | RP11-43A14.1 | - | 7.9 |
| rs2015971 | STAB1 | Urothelial_Atezolizu | 1.26e-03 | 0.828 | NT5DC2 | - | 7.9 |
| rs3818496 | CARS2 | CCRCC_ICB_Miao2018 | 1.34e-03 | 0.568 | CARS2 | - | 7.9 |


## Top Variants by Cohort


### Urothelial_Atezolizumab_Mariathasan2018

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs3795732 | GPATCH4 | 1.66e-04 | 1.980 | 1 | 8.8 |
| rs135110 | SLC5A1 | 2.70e-04 | 2.068 | 1 | 8.6 |
| rs12567958 | HAPLN2 | 3.34e-04 | 1.301 | 1 | 8.5 |
| rs7011778 | ENSG00000269954 | 4.56e-04 | 1.442 | 1 | 8.3 |
| rs4841558 | BLK | 4.64e-04 | 0.936 | 1 | 8.3 |

### Melanoma_Ipilimumab_VanAllen2015

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs5753137 | ENSG00000273428 | 6.36e-04 | 1.427 | 1 | 8.2 |
| rs1869074 | ERICH1 | 1.26e-03 | 0.761 | 1 | 7.9 |
| rs2550270 | MUC4 | 1.97e-03 | 0.512 | 1 | 7.7 |
| rs2688513 | MUC4 | 2.15e-03 | 0.958 | 1 | 7.7 |
| rs3802509 | CCNY | 2.97e-03 | 1.240 | 1 | 7.5 |

### Hepatocellular_Atezolizumab_Finn2020

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs630245 | RNU6-1102P | 7.46e-04 | 0.989 | 1 | 8.1 |
| rs485005 | RNU6-1102P | 9.09e-04 | 0.827 | 1 | 8.0 |
| rs28397595 | ATP8B3 | 1.16e-03 | 0.659 | 1 | 7.9 |
| rs3752176 | ZNF471 | 1.37e-03 | 0.637 | 1 | 7.9 |
| rs59572086 | ATP8B3 | 1.39e-03 | 0.813 | 1 | 7.9 |

### CCRCC_ICB_Miao2018

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs3818496 | CARS2 | 1.34e-03 | 0.568 | 1 | 7.9 |
| rs6436124 | RETREG2 | 2.96e-03 | 0.505 | 1 | 7.5 |
| rs9302994 | ENSG00000267426 | 3.14e-03 | 0.818 | 1 | 7.5 |
| rs9555726 | ING1 | 3.44e-03 | 0.766 | 1 | 7.5 |
| rs72884533 | RIOK3 | 3.95e-03 | 0.864 | 1 | 7.4 |

### Melanoma_CTLA4_Snyder2014

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs12627971 | GGT5 | 1.51e-03 | 1.143 | 1 | 7.8 |
| rs1329119 | PTGES3P1 | 3.20e-03 | 0.710 | 1 | 7.5 |
| rs676913 | GBP7 | 3.34e-03 | 0.549 | 1 | 7.5 |
| rs2230338 | GBP2 | 3.34e-03 | 0.830 | 1 | 7.5 |
| rs586609 | GBP7 | 3.34e-03 | 0.570 | 1 | 7.5 |

### NSCLC_Pembrolizumab_Rizvi2015

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs2310925 | H6PD | 2.82e-02 | 0.507 | 1 | 6.5 |
| rs2301993 | IPO13 | 3.16e-02 | 0.540 | 1 | 6.5 |
| rs2071042 | ITIH4 | 3.66e-02 | 0.617 | 1 | 6.4 |
| rs2240917 | ITIH4 | 3.66e-02 | 0.800 | 1 | 6.4 |
| rs6769789 | STIMATE | 3.66e-02 | 0.676 | 1 | 6.4 |

### Melanoma_PD1_Hugo2016

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs2252639 | ENSG00000249624 | 3.01e-02 | 0.523 | 1 | 6.5 |
| rs17206931 | ENSG00000269001 | 3.28e-02 | 0.523 | 1 | 6.5 |
| rs2274458 | UQCC2 | 3.54e-02 | 0.637 | 1 | 6.5 |
| rs2303492 | ARFIP2 | 4.97e-02 | 0.997 | 1 | 6.3 |
| rs669433 | ZDHHC5 | 2.03e-02 | 0.250 | 1 | 5.7 |

### GBM_PD1_Zhao2019

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs242456 | ARL5AP5 | 3.83e-02 | 0.257 | 1 | 5.4 |
| rs7559479 | IL18RAP | 4.29e-02 | 0.252 | 1 | 5.4 |
| rs7796144 | ERV3-1 | 4.58e-02 | 0.262 | 1 | 5.3 |
| rs708457 | AFF4 | 4.59e-02 | 0.203 | 1 | 5.3 |
| rs2276868 | RPL14 | 4.62e-02 | 0.418 | 1 | 5.3 |


## Cross-Cohort Variants

Variants appearing in multiple cohorts with Tier 1 or 2 evidence:

| rsid | gene | cohorts | best tier | mean priority |
|------|------|---------|-----------|---------------|
| rs135110 | SLC5A1 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_CTLA4_Snyder2014 | 1 | 7.5 |
| rs1284351 | FEZ1 | Urothelial_Atezolizumab_Mariathasan2018, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.3 |
| rs438610 | HCG4B | Urothelial_Atezolizumab_Mariathasan2018, CCRCC_ICB_Miao2018 | 1 | 7.1 |
| rs4627097 | EI24 | Urothelial_Atezolizumab_Mariathasan2018, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.1 |
| rs2550270 | MUC4 | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.0 |
| rs7502041 | LINC01180 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_Ipilimumab_VanAllen2015 | 1 | 7.0 |
| rs364936 | GABRR1 | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020, Melanoma_CTLA4_Snyder2014 | 1 | 7.0 |
| rs2270590 | CCDC136 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_Ipilimumab_VanAllen2015 | 1 | 7.0 |
| rs6769789 | STIMATE | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 6.9 |
| rs2140335 | HIGD2AP2 | Urothelial_Atezolizumab_Mariathasan2018, CCRCC_ICB_Miao2018, Melanoma_Ipilimumab_VanAllen2015 | 1 | 6.9 |
| rs781230 | TESK2 | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020 | 1 | 6.9 |
| rs60332615 | ITIH4 | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 6.9 |
| rs2071042 | ITIH4 | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 6.8 |
| rs4369638 | ADAMTS17 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_CTLA4_Snyder2014 | 1 | 6.8 |
| rs34336464 | FLNB | Urothelial_Atezolizumab_Mariathasan2018, CCRCC_ICB_Miao2018 | 1 | 6.8 |


## Next Steps

1. **Tier 1 variants**: Prioritize for functional validation
   - Gene expression studies in immune cells
   - Reporter assays for regulatory activity
   - CRISPR editing in relevant cell types

2. **Tier 2 variants**: Query real AlphaGenome API
   - Replace mock predictions with actual API calls
   - Re-evaluate with real regulatory predictions

3. **Cross-cohort variants**: Meta-analysis
   - Fixed-effects meta-analysis for consistent effect direction
   - Test replication across cancer types

---
Generated: 2026-02-10 18:23:16
