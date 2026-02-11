# Immunotherapy Variant Prioritization Report

## Summary

| Metric | Value |
|--------|-------|
| Total variants | 15,918 |
| Tier 1 (AG + eQTL) | 3,681 (23.1%) |
| Tier 2 (AG only) | 6,996 (44.0%) |
| Tier 3 (eQTL only) | 1,720 (10.8%) |
| Tier 4 (therapy only) | 3,521 (22.1%) |

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
| CCRCC_ICB_Miao2018 | 1,934 | 537 | 4.24 |
| GBM_PD1_Zhao2019 | 40 | 9 | 3.62 |
| Hepatocellular_Atezolizumab_Finn2020 | 3,259 | 655 | 3.98 |
| Melanoma_CTLA4_Snyder2014 | 2,366 | 531 | 3.96 |
| Melanoma_Ipilimumab_VanAllen2015 | 2,952 | 754 | 4.16 |
| Melanoma_PD1_Hugo2016 | 802 | 149 | 3.68 |
| NSCLC_Pembrolizumab_Rizvi2015 | 859 | 153 | 3.71 |
| Urothelial_Atezolizumab_Mariathasan2018 | 3,706 | 893 | 4.23 |


## Top 20 Tier 1 Variants

These variants have the strongest evidence: therapy association + AlphaGenome prediction + eQTL validation.

| rsid | gene | cohort | pval | AG impact | DICE | OneK1K | Priority |
|------|------|--------|------|-----------|------|--------|----------|
| rs3795732 | GPATCH4 | Urothelial_Atezolizu | 1.66e-04 | 1.980 | TTC24 | APOA1BP | 9.8 |
| rs135110 | SLC5A1 | Urothelial_Atezolizu | 2.70e-04 | 2.068 | YWHAH | YWHAH | 9.6 |
| rs12567958 | HAPLN2 | Urothelial_Atezolizu | 3.34e-04 | 1.301 | TTC24 | APOA1BP | 9.5 |
| rs7011778 | ENSG00000269954 | Urothelial_Atezolizu | 4.56e-04 | 1.442 | RP11-148O21.4 | BLK | 9.3 |
| rs4841558 | BLK | Urothelial_Atezolizu | 4.64e-04 | 0.936 | RP11-148O21.4 | BLK | 9.3 |
| rs6988443 | BLK | Urothelial_Atezolizu | 4.82e-04 | 1.275 | RP11-148O21.4 | BLK | 9.3 |
| rs5753137 | ENSG00000273428 | Melanoma_Ipilimumab_ | 6.36e-04 | 1.427 | TBC1D10A | RP4-539M6.22 | 9.2 |
| rs5998233 | SLC5A1 | Urothelial_Atezolizu | 6.50e-04 | 1.129 | C22orf24 | YWHAH | 9.2 |
| rs11785333 | BLK | Urothelial_Atezolizu | 6.72e-04 | 1.858 | RP11-148O21.4 | BLK | 9.2 |
| rs630245 | RNU6-1102P | Hepatocellular_Atezo | 7.46e-04 | 0.989 | NEXN | NEXN | 9.1 |
| rs4841557 | BLK | Urothelial_Atezolizu | 8.36e-04 | 1.545 | RP11-148O21.4 | BLK | 9.1 |
| rs485005 | RNU6-1102P | Hepatocellular_Atezo | 9.09e-04 | 0.827 | NEXN | NEXN | 9.0 |
| rs848216 | ZBTB17 | Urothelial_Atezolizu | 1.14e-03 | 0.799 | CLCNKA | RP11-169K16. | 8.9 |
| rs2581774 | RFT1 | Urothelial_Atezolizu | 1.24e-03 | 0.514 | RFT1 | RFT1 | 8.9 |
| rs2015971 | STAB1 | Urothelial_Atezolizu | 1.26e-03 | 0.828 | NT5DC2 | PPM1M | 8.9 |
| rs3818496 | CARS2 | CCRCC_ICB_Miao2018 | 1.34e-03 | 0.568 | CARS2 | CARS2 | 8.9 |
| rs12405992 | HAPLN2 | Urothelial_Atezolizu | 1.35e-04 | 0.407 | TTC24 | TTC24 | 8.9 |
| rs11214598 | ANKK1 | Hepatocellular_Atezo | 1.39e-03 | 0.686 | ANKK1 | TTC12 | 8.9 |
| rs2550270 | MUC4 | Melanoma_Ipilimumab_ | 1.97e-03 | 0.512 | MUC4 | MUC20 | 8.7 |
| rs7860949 | ENDOG | Urothelial_Atezolizu | 2.07e-03 | 0.933 | PKN3 | ENDOG | 8.7 |


## Top Variants by Cohort


### Urothelial_Atezolizumab_Mariathasan2018

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs3795732 | GPATCH4 | 1.66e-04 | 1.980 | 1 | 9.8 |
| rs135110 | SLC5A1 | 2.70e-04 | 2.068 | 1 | 9.6 |
| rs12567958 | HAPLN2 | 3.34e-04 | 1.301 | 1 | 9.5 |
| rs7011778 | ENSG00000269954 | 4.56e-04 | 1.442 | 1 | 9.3 |
| rs4841558 | BLK | 4.64e-04 | 0.936 | 1 | 9.3 |

### Melanoma_Ipilimumab_VanAllen2015

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs5753137 | ENSG00000273428 | 6.36e-04 | 1.427 | 1 | 9.2 |
| rs2550270 | MUC4 | 1.97e-03 | 0.512 | 1 | 8.7 |
| rs2688513 | MUC4 | 2.15e-03 | 0.958 | 1 | 8.7 |
| rs3802509 | CCNY | 2.97e-03 | 1.240 | 1 | 8.5 |
| rs4759794 | ENSG00000256299 | 3.95e-03 | 0.991 | 1 | 8.4 |

### Hepatocellular_Atezolizumab_Finn2020

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs630245 | RNU6-1102P | 7.46e-04 | 0.989 | 1 | 9.1 |
| rs485005 | RNU6-1102P | 9.09e-04 | 0.827 | 1 | 9.0 |
| rs11214598 | ANKK1 | 1.39e-03 | 0.686 | 1 | 8.9 |
| rs4938016 | ANKK1 | 2.30e-03 | 0.908 | 1 | 8.6 |
| rs4348317 | CENPQ | 2.63e-03 | 0.943 | 1 | 8.6 |

### CCRCC_ICB_Miao2018

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs3818496 | CARS2 | 1.34e-03 | 0.568 | 1 | 8.9 |
| rs6436124 | RETREG2 | 2.96e-03 | 0.505 | 1 | 8.5 |
| rs9555726 | ING1 | 3.44e-03 | 0.766 | 1 | 8.5 |
| rs72884533 | RIOK3 | 3.95e-03 | 0.864 | 1 | 8.4 |
| rs885937 | HLA-V | 5.28e-03 | 0.571 | 1 | 8.3 |

### Melanoma_CTLA4_Snyder2014

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs1329119 | PTGES3P1 | 3.20e-03 | 0.710 | 1 | 8.5 |
| rs676913 | GBP7 | 3.34e-03 | 0.549 | 1 | 8.5 |
| rs586609 | GBP7 | 3.34e-03 | 0.570 | 1 | 8.5 |
| rs2230338 | GBP2 | 3.34e-03 | 0.830 | 1 | 8.5 |
| rs1803632 | GBP2 | 3.34e-03 | 1.218 | 1 | 8.5 |

### NSCLC_Pembrolizumab_Rizvi2015

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs2310925 | H6PD | 2.82e-02 | 0.507 | 1 | 7.5 |
| rs60332615 | ITIH4 | 3.66e-02 | 0.771 | 1 | 7.4 |
| rs2071042 | ITIH4 | 3.66e-02 | 0.617 | 1 | 7.4 |
| rs6769789 | STIMATE | 3.66e-02 | 0.676 | 1 | 7.4 |
| rs2240917 | ITIH4 | 3.66e-02 | 0.800 | 1 | 7.4 |

### Melanoma_PD1_Hugo2016

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs2252639 | ENSG00000249624 | 3.01e-02 | 0.523 | 1 | 7.5 |
| rs17206931 | ENSG00000269001 | 3.28e-02 | 0.523 | 1 | 7.5 |
| rs2274458 | UQCC2 | 3.54e-02 | 0.637 | 1 | 7.5 |
| rs669433 | ZDHHC5 | 2.03e-02 | 0.250 | 1 | 6.7 |
| rs6978892 | CICP11 | 2.54e-02 | 0.341 | 1 | 6.6 |

### GBM_PD1_Zhao2019

| rsid | gene | pval | AG impact | Tier | Priority |
|------|------|------|-----------|------|----------|
| rs7559479 | IL18RAP | 4.29e-02 | 0.252 | 1 | 6.4 |
| rs2276868 | RPL14 | 4.62e-02 | 0.418 | 1 | 6.3 |
| rs1736982 | HLA-V | 4.76e-02 | 0.300 | 1 | 6.3 |
| rs242456 | ARL5AP5 | 3.83e-02 | 0.257 | 1 | 5.4 |
| rs3817465 | IL18RAP | 4.29e-02 | 0.158 | 3 | 5.4 |


## Cross-Cohort Variants

Variants appearing in multiple cohorts with Tier 1 or 2 evidence:

| rsid | gene | cohorts | best tier | mean priority |
|------|------|---------|-----------|---------------|
| rs135110 | SLC5A1 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_CTLA4_Snyder2014 | 1 | 8.5 |
| rs2550270 | MUC4 | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020 | 1 | 8.0 |
| rs6769789 | STIMATE | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 7.9 |
| rs781230 | TESK2 | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.9 |
| rs60332615 | ITIH4 | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 7.9 |
| rs2071042 | ITIH4 | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 7.8 |
| rs2297995 | L2HGDH | Melanoma_Ipilimumab_VanAllen2015, Melanoma_CTLA4_Snyder2014 | 1 | 7.8 |
| rs17851970 | DIS3L | Melanoma_Ipilimumab_VanAllen2015, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.8 |
| rs2240917 | ITIH4 | Urothelial_Atezolizumab_Mariathasan2018, NSCLC_Pembrolizumab_Rizvi2015 | 1 | 7.7 |
| rs3750207 | ZC3H3 | Urothelial_Atezolizumab_Mariathasan2018, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.7 |
| rs2274458 | UQCC2 | Urothelial_Atezolizumab_Mariathasan2018, Melanoma_PD1_Hugo2016 | 1 | 7.6 |
| rs378971 | HCG9 | CCRCC_ICB_Miao2018, Urothelial_Atezolizumab_Mariathasan2018 | 1 | 7.3 |
| rs1284351 | FEZ1 | Urothelial_Atezolizumab_Mariathasan2018, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.3 |
| rs35162188 | ENSG00000267672 | Urothelial_Atezolizumab_Mariathasan2018, Hepatocellular_Atezolizumab_Finn2020 | 1 | 7.2 |
| rs1866561 | ABHD1 | Melanoma_Ipilimumab_VanAllen2015, Urothelial_Atezolizumab_Mariathasan2018 | 1 | 7.2 |


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
Generated: 2026-02-10 19:40:14
