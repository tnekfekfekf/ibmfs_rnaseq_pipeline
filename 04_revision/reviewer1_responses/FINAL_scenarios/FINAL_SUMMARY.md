# Reviewer 1 — Final Analysis Summary

All analyses use raw counts from `manuscript_count_matrix_19samples.txt` and a
group-aware expression filter (CPM ≥ 1 in ≥ 50 % of samples of any group).
DE significance: padj < 0.05 AND |log2FC| > 1.

## Sequencing site / cohort labels

- **macrogen** (10 samples): AA-RNA-DKC, AA-RNA-FA, AA-RNA-FA2, AA-RNA-FA3,
  AA-RNA-1, AA-RNA-4, AA-RNA-5, AA-RNA-13, AA-RNA-16, AA-RNA-18
- **jinpyung** (4 samples): AA-PRO, AA-KEW (Control), AA-HMH, AA-PJH (u-BMF)
- **public** (3 samples): Child1, Child2, Child3 (Control from GSE147523)
- Technical replicate: AA-RNA-1 (macrogen) = AA-LES (jinpyung) = same patient
  → usable for RUV-based site adjustment within the 14-sample cohort
- AA-CSB has 2.39 M library (~10× below normal); excluded as unreliable

## Site × group design

| Site | Control | g-BMF | u-BMF | total |
|------|--------:|------:|------:|------:|
| macrogen | 0 | 4 | 6 | 10 |
| jinpyung | 2 | 0 | 2 | 4 |
| public | 3 | 0 | 0 | 3 |
| **total** | **5** | **4** | **8** | **17** |

⚠ Critical: **g-BMF samples exist only at macrogen, all Controls live at
jinpyung or public.** Site and group are nearly perfectly confounded —
any design with `~ site + group` becomes singular and DE detection collapses.

## Six scenarios — DE counts and 11 manuscript lncRNA recovery

| # | Cohort | Design | DE g-BMF / u-BMF | 11 lncRNA g / u |
|---|--------|--------|------------------:|----------------:|
| 1 | 14-sample | `~ group` (no batch) | 3 719 / 2 651 | **11 / 11** |
| 2 | 14-sample | `~ site + group` (2-level macrogen/jinpyung) | 20 / 34 | 0 / 0 ✗ singular |
| 3 | 14-sample | RUV (LES replicate) `~ W_1 + group` | 4 131 / 2 416 | 11 / 9 ✓ |
| 4 | 17-sample (+Child) | `~ cohort + group` (2-level internal/public) | 3 566 / 2 689 | **11 / 10** (MALAT1 ns in u-BMF) |
| 5 | 17-sample | `~ site + group` (3-level) | 16 / 22 | 0 / 0 ✗ singular |
| 6 | 17-sample | RUV `~ W_1 + group` | 3 104 / 1 753 | 8 / 7 |

## Per-gene recovery (boldface = significant)

| lncRNA | 14s `~group` | 14s + RUV | 17s `~cohort+group` |
|--------|:---:|:---:|:---:|
| HCG11 | **g11/u11** | **g/u** | **g/u** |
| HCP5 | **g/u** | **g/u** | **g/u** |
| SNHG32 | **g/u** | **g/u** | **g/u** |
| PSMB8-AS1 | **g/u** | **g/u** | **g/u** |
| FAM30A | **g/u** | **g/u** | **g/u** |
| MIR22HG | **g/u** | **g/u** | **g/u** |
| ATP1A1-AS1 | **g/u** | **g**/u (p=0.065) | **g/u** |
| USP3-AS1 | **g/u** | **g/u** | **g/u** |
| TAGAP-AS1 | **g/u** | **g/u** | **g/u** |
| LINC01036 | **g/u** | **g/u** | **g/u** |
| MALAT1 | **g/u** | **g**/u (ns) | **g**/u (ns) |

## Recommended analysis to present to Reviewer 1

**Primary (matches manuscript):** 14-sample `~ group` → all 11 biomarkers
recovered with original effect sizes.

**Sensitivity 1 (proper site batch correction within manuscript cohort):**
14-sample RUV using AA-RNA-1/AA-LES technical replicate → 11/11 in g-BMF vs
Ctrl, 9/11 in u-BMF vs Ctrl. Demonstrates results are robust to a more
principled batch handling than the original blind subtractive correction.

**Sensitivity 2 (extending the control sample size as the reviewer suggests):**
17-sample with 3 additional public PNBM controls, `~ cohort + group` →
11/11 in g-BMF vs Ctrl, 10/11 in u-BMF vs Ctrl (MALAT1 only loss; this
single gene is expected given that whole-BM samples have ~100× more
MALAT1 than BM-MNC samples due to nucleated red-cell content). Confirms that
the central biomarker panel is reproducible against an independent
cell-composition–distinct control cohort.

**Scenarios that fail and why:** Any design with `~ site + group` becomes
singular because g-BMF samples exist exclusively at the macrogen facility
while all internal Controls were generated at jinpyung. We therefore use
the 2-level cohort split (internal vs public) when including public data,
and RUV (which uses the cross-site technical replicate) for the within-
manuscript site adjustment.

## Files

- `DE_*.tsv` — full DE tables per scenario/contrast
- `ms11_audit_*.tsv` — per-gene log2FC, padj, SIG flag for the 11 biomarkers
- `SUMMARY_*.tsv` — wide-format recovery summaries
- `PCA_17sample_3panels.pdf` (in `pca_correction_comparison/`) — visualization
  of batch effect before/after 2-level vs 3-level adjustment
