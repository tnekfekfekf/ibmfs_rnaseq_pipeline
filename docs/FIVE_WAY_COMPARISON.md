# 5-Way DEG Comparison — Final Reference Table

**Date:** 2026-05-12
**Question:** How do 5 different analytical approaches compare for DEG identification?

## 5 Analyses

| # | Setup | Counts | Design | n |
|---|---|---|---|---|
| **1** | Manuscript paper (published) | — | — | 14 |
| **2** | Manuscript raw counts → DESeq2 | `manuscript_count_matrix_19samples.txt` (14 cols) | `~ group` | 14 |
| **3** | v3 pipeline 14 samples | `fc_manuscript_v3.txt` | `~ group` | 14 |
| **4** | **Hybrid: MS counts + Child v3** | MS 14 + Child v3 3 | `~ cohort + group` | 17 |
| **5** | v3 pipeline 17 samples | `fc_v3_17samples.txt` | `~ cohort + group` | 17 |

## TABLE 1: DEG counts (padj<0.05, |LFC|>1)

| Contrast | Metric | (1) Paper | (2) MS counts | (3) v3 14 | (4) Hybrid | (5) v3 17 |
|---|---|---|---|---|---|---|
| G-AA vs Ctrl | mRNA | **2078** | 2078 ✅ | 1425 | 2023 | 1453 |
| G-AA vs Ctrl | lncRNA | **1167** | 1167 ✅ | 262 | 1081 | 265 |
| U-AA vs Ctrl | mRNA | **1315** | 1315 ✅ | 906 | 1310 | 908 |
| U-AA vs Ctrl | lncRNA | **992** | 992 ✅ | 203 | 915 | 199 |
| G-AA vs U-AA | mRNA | **4** | 4 ✅ | 4 | 5 | 2 |
| G-AA vs U-AA | lncRNA | **4** | 4 ✅ | 3 | 7 | 2 |

## TABLE 2: Total + Up/Down breakdown

| Contrast | Analysis | Total | mRNA | lncRNA | Up | Down |
|---|---|---|---|---|---|---|
| G-AA vs Ctrl | (2) MS DESeq | 4158 | 2078 | 1167 | 3350 | 808 |
| G-AA vs Ctrl | (3) v3 14 | 2014 | 1425 | 262 | 1428 | 586 |
| G-AA vs Ctrl | (4) Hybrid | 3923 | 2023 | 1081 | 3173 | 750 |
| G-AA vs Ctrl | (5) v3 17 | 2038 | 1453 | 265 | 1478 | 560 |
| U-AA vs Ctrl | (2) MS DESeq | 2968 | 1315 | 992 | 2822 | 146 |
| U-AA vs Ctrl | (3) v3 14 | 1277 | 906 | 203 | 1191 | 86 |
| U-AA vs Ctrl | (4) Hybrid | 2831 | 1310 | 915 | 2711 | 120 |
| U-AA vs Ctrl | (5) v3 17 | 1293 | 908 | 199 | 1214 | 79 |
| G-AA vs U-AA | (2) MS DESeq | 11 | 4 | 4 | 1 | 10 |
| G-AA vs U-AA | (3) v3 14 | 13 | 4 | 3 | 3 | 10 |
| G-AA vs U-AA | (4) Hybrid | 17 | 5 | 7 | 0 | 17 |
| G-AA vs U-AA | (5) v3 17 | 9 | 2 | 2 | 2 | 7 |

## TABLE 3: 13 manuscript-highlighted genes (G-AA vs Ctrl)

| Gene | (1) Paper | (2) MS | (3) v3 14 | (4) Hybrid | (5) v3 17 |
|---|---|---|---|---|---|
| **HCG11** | 9.49 sig | 9.49 sig | 9.52 sig | 9.50 sig | 9.52 sig |
| **HCP5** | 9.07 sig | 9.07 sig | 9.16 sig | 9.09 sig | 9.16 sig |
| **SNHG32** | 9.65 sig | 9.65 sig | 9.75 sig | 9.69 sig | 9.76 sig |
| **PSMB8-AS1** | 11.43 sig | 11.43 sig | 12.65 sig | 11.45 sig | 12.65 sig |
| **FAM30A** | 5.46 sig | 5.46 sig | 5.60 sig | 5.50 sig | 5.61 sig |
| **MIR22HG** | 9.89 sig | 9.89 sig | 9.95 sig | 9.90 sig | 9.94 sig |
| **ATP1A1-AS1** ⭐ | 2.34 sig | 2.34 sig | **0.22 n.s.** ❌ | 2.36 sig | **0.22 n.s.** ❌ |
| **USP3-AS1** ⭐ | 2.60 sig | 2.60 sig | **0.16 n.s.** ❌ | 2.62 sig | **0.16 n.s.** ❌ |
| **TAGAP-AS1** | 4.36 sig | 4.36 sig | **0.58 n.s.** ❌ | 4.39 sig | **0.58 n.s.** ❌ |
| **LINC01036** | 3.25 sig | 3.25 sig | **1.76 n.s.** ❌ | 3.28 sig | **1.77 n.s.** ❌ |
| **MALAT1** | -6.56 sig | -6.56 sig | **-0.63 n.s.** ❌ | -6.54 sig | **-0.63 n.s.** ❌ |
| **TEN1-CDK3** | -4.57 sig | -4.57 sig | -4.62 sig | -4.54 sig | -4.61 sig |
| **FANCA** | -2.24 sig | -2.24 sig | -2.45 sig | -2.22 sig | -2.45 sig |

⭐ = RT-qPCR validated

## Key observations

### (1) ↔ (2): Perfect match — manuscript pipeline verified
Manuscript paper values 100% reproduce when running DESeq2 on the authentic count matrix.

### (3) vs (2): v3 vs MS quantification — count level difference
The v3 pipeline (different featureCounts options) gives ~69% mRNA recall, ~22% lncRNA recall. 5 borderline lncRNAs lose significance.

### (4) ↔ (2): Hybrid (MS + Child) preserves manuscript
Adding 3 Child controls with `~ cohort + group` design + manuscript counts for the 14 = 97% mRNA recall, 93% lncRNA recall. **All 6 strong biomarkers + 5 borderline lncRNAs preserved.**

### (3) ↔ (5): v3 framework consistent
Adding Child to v3 framework doesn't help recover the missing borderline lncRNAs because the v3 count matrix itself doesn't contain those signals.

## Conclusion: which analysis to use?

| Goal | Best analysis |
|---|---|
| **Reproduce manuscript exactly** | (2) MS counts + DESeq2 |
| **Address reviewer 2 with public controls** | (4) Hybrid + batch-aware |
| **Future samples added to manuscript baseline** | (4)-style: MS counts as anchor, new samples v3 |
| **Complete v3-only pipeline (less faithful but consistent)** | (5) v3 17 samples |

## Files

- `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/five_way_comparison/`
  - `table1_deg_counts.tsv`
  - `table2_total_updown.tsv`
  - `table3_key_genes.tsv`
- `hypothesis_tests/five_way_comparison.R` — reproduce
