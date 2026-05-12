# Option 1: 4-Way Analysis (Manuscript + Child) — Comparison Report

**Date:** 2026-05-12
**Setup:** 14 internal samples from authentic manuscript count matrix + 3 Child samples quantified with our v3 pipeline.

## Setup Details

- **14 manuscript samples:** Use `manuscript_count_matrix_19samples.txt` columns directly (raw counts, manuscript original)
- **3 Child samples:** Use our v3 quantification (`per_sample_v3/Child{1,2,3}.counts.txt`)
- **Combined matrix:** 62,153 genes × 17 samples

Library sizes:
- Manuscript samples: 20-33M
- Child v3 samples: 25-32M (similar range — no major depth issue)

## 4 Analyses

| Label | Samples | Controls | Design | Purpose |
|---|---|---|---|---|
| **A** | 14 (manuscript) | AA-PRO, AA-KEW (2) | `~ group` | Manuscript replica (gold standard) |
| **B** | 17 (+ Child) | 5 mixed (naive) | `~ group` | What happens if we just add public |
| **C** | 17 (+ Child) | 5 mixed (batch-aware) | `~ cohort + group` | Proper way to add public |
| **D** | 15 (no internal ctrl) | 3 Child only | `~ group` | Public-only baseline |

## TABLE 1: DEG counts

| Analysis | G-AA vs Ctrl total | G-AA mRNA | G-AA lncRNA | U-AA mRNA | U-AA lncRNA | G-AA vs U-AA |
|---|---|---|---|---|---|---|
| **A** (manuscript replica) | **4158** | **2081** | **1168** | 1335 | 1007 | 11 |
| B (17 naive) | 4591 | 1619 | 1645 | 983 | 1494 | 21 |
| C (17 batch-aware) | 3923 | 2023 | 1081 | 1310 | 915 | 17 |
| D (3 Child only) | 5008 | 2045 | 1647 | 1968 | 1672 | 13 |

**Manuscript paper reports:** 2078 mRNA + 1167 lncRNA (G-AA vs Ctrl). A reproduces this within ±3 genes.

## TABLE 2: Gene-level overlap with manuscript (A)

| Contrast | A vs B | **A vs C** | A vs D |
|---|---|---|---|
| G-AA vs Ctrl | 54.7% (2274/4158) | **89.7% (3730/4158)** ✅ | 45.4% (1887/4158) |
| U-AA vs Ctrl | 53.0% (1595/3011) | **90.8% (2733/3011)** ✅ | 49.9% (1502/3011) |
| G-AA vs U-AA | 54.5% (6/11) | 72.7% (8/11) | 45.5% (5/11) |

## TABLE 3: Manuscript-highlighted genes status

### ⭐ 6 Strong DE-lncRNA biomarkers (manuscript Figure 2)

| Gene | A: 14ms | B: 17 naive | C: 17 batch-aware | D: 3 Child only |
|---|---|---|---|---|
| HCG11 | ✓ LFC 9.48 | **✗ 0.85** (n.s.) | ✓ 9.50 | **✗ 0.10** |
| HCP5 | ✓ 9.07 | **✗ 0.09** | ✓ 9.09 | **✗ -0.65** |
| SNHG32 | ✓ 9.65 | **✗ 0.74** | ✓ 9.69 | **✗ 0.00** |
| PSMB8-AS1 | ✓ 11.42 | **✗ -0.13** | ✓ 11.45 | **✗ -0.88** |
| FAM30A | ✓ 5.46 | **✗ 0.40** | ✓ 5.50 | **✗ -0.33** |
| MIR22HG | ✓ 9.89 | **✗ -0.12** | ✓ 9.90 | **✗ -0.86** |

**6/6 biomarkers LOST in naive B and D analyses.** Only batch-aware C preserves all 6.

### 5 Borderline DE-lncRNAs

| Gene | A | B | C | D |
|---|---|---|---|---|
| ATP1A1-AS1 | ✓ 2.34 | ✓ 2.35 | ✓ 2.36 | ✓ 2.34 |
| USP3-AS1 | ✓ 2.60 | ✓ 3.27 | ✓ 2.62 | ✓ 3.94 |
| TAGAP-AS1 | ✓ 4.36 | ✓ 3.56 | ✓ 4.39 | ✓ 3.19 |
| LINC01036 | ✓ 3.25 | ✓ 3.87 | ✓ 3.28 | ✓ 4.46 |
| MALAT1 | ✓ -6.56 | ✓ -6.30 | ✓ -6.54 | ✓ -6.13 |

All borderline lncRNAs remain sig across all 4 analyses (less sensitive to control composition).

### G-AA-specific genes

| Gene | A | B | C | D |
|---|---|---|---|---|
| TEN1-CDK3 | ✓ -4.57 | ✓ -3.50 | ✓ -4.54 | **✗ -1.95** |
| FANCA | ✓ -2.24 | ✓ -1.68 | ✓ -2.22 | **✗ -1.17** |

D analysis loses G-AA-specific signatures.

## Biological Interpretation

### Why Strong DE-lncRNAs disappear in B (naive)

The 6 manuscript biomarkers (HCG11, HCP5, SNHG32, PSMB8-AS1, FAM30A, MIR22HG) are:
- HIGH in BMF patient samples
- LOW in internal aspirate controls (whole bone marrow, mostly erythroid/granulocytic)
- HIGH in Child MNC samples (lymphoid-enriched MNC isolate)

When you naively pool Child MNC with internal aspirate as "Control":
- Control mean becomes HIGH (because Child contributes high values)
- Patient mean stays HIGH
- → No difference → not significant

This is a **textbook example of biology-confounded-with-batch**.

### Why C (batch-aware) preserves them

The `~ cohort + group` design separates the cohort effect (MNC vs aspirate) from the group effect (patient vs control). The group estimate is computed WITHIN cohort, so:
- Within internal aspirate samples: patient HIGH vs Control LOW → preserved
- Cohort term absorbs the MNC-baseline elevation

→ **Batch-aware is the ONLY valid way to add public Child controls.**

## Recommendations

### For Manuscript Revision Response

1. **Primary analysis: A (14 ms, 2 ctrl) — manuscript as published**
2. **Sensitivity analysis: C (17, batch-aware) — adds 3 public Child controls**
   - Demonstrates robustness (90% DEG overlap)
   - Addresses reviewer 2's concern about n=2 controls
3. **Show B as cautionary** — naive pooling loses 6 key biomarkers (dramatic visualization)
4. **Do NOT use D as substitute** — different cohort baseline

### For Future Sample Addition

When adding more samples to extend the dataset:
- Use `PIPELINE/06_add_new_samples.sh` (v3 quantification)
- Always use batch-aware design if mixing with manuscript samples
- Test naive vs batch-aware as sensitivity check

## Output files

`/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/option1_analyses/`:
- `combined_14ms_3child_v3.txt` — 17-sample combined count matrix
- `table1_deg_counts.tsv` — DEG counts per analysis
- `table2_overlap_with_A.tsv` — gene-level overlaps
- `table3_key_genes_GvC.tsv` — 19 highlighted gene status across 4 analyses
- `all_dds_objects.rds` — saved DESeq2 objects for re-use

Reproduce: `Rscript hypothesis_tests/option1_combined_analysis.R`
