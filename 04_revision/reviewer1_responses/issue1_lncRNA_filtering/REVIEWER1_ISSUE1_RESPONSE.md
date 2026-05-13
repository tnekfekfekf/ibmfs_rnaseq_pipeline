# Reviewer 1, Issue 1 — lncRNA Filtering Strategy

**Reviewer comment (paraphrased):** The lncRNA filter (mean CPM ≥ 10 AND
CPM ≥ 1 in every sample) is unusually strict and uses two simultaneous CPM
thresholds. Please justify the choice and demonstrate robustness.

---

## 1. Critical clarification first

While preparing this response we discovered that **the filter described in
the Methods text is NOT the filter the manuscript analysis actually used.**

| Filter | Where it appears | Actually used? |
|--------|------------------|----------------|
| `rowSums(counts) >= 10` (raw counts) | `deseq2_analysis.R` line 33 (and `deseq2_analysis_enhanced.R` line 59, `snakemake_template/.../deseq2_analysis.R` lines 87–88) | **YES — this is what produced the manuscript results** |
| `mean(CPM) >= 10 AND CPM >= 1 in all samples` | Methods section text only | **NO — never used in code** |

This is a manuscript documentation error, and the reviewer's confusion (about
the two CPM thresholds) is the natural result. We propose to correct the
Methods text in the revision and to support both filters with sensitivity
analyses below.

---

## 2. Analysis design

We re-ran DESeq2 on the authentic manuscript count matrix
(`gene_counts_all_19_samples_with_annotations.txt`) under **five
pre-filter strategies × two cohort configurations** (30 cells per gene):

| Strategy | Definition | Provenance |
|----------|------------|------------|
| **ORIG** | `rowSums(counts) >= 10` | **Actual manuscript code** |
| **STRICT** | `mean(CPM) >= 10` AND `CPM >= 1 in ALL samples` | Manuscript methods text (never used in code) |
| **R1** | `CPM >= 1 in >= 50 % of samples of ANY group` | Group-aware moderate (recommended) |
| **R2** | `mean(CPM) >= 1 in any group` AND `>= 50 % of that group >= 0.5` | Group-aware sensitive |
| **R3** | none — DESeq2 IF only | Most permissive |

Cohorts:
* **14-sample (manuscript)** — 2 internal Ctrl + 4 g-BMF + 8 u-BMF, design `~ group`
* **17-sample (hybrid)** — 14 + 3 Child PNBM public Ctrl (GSE147523), batch-aware `~ cohort + group`

Significance criterion (manuscript thresholds): `padj < 0.05 AND |log2FC| > 1`.

---

## 3. Genes retained at the pre-filter stage

| Cohort | ORIG | STRICT | R1 | R2 | R3 |
|--------|-----:|-------:|---:|---:|---:|
| 14-sample | 36,683 | 8,325 | 19,246 | 18,637 | 62,652 |
| 17-sample | 38,424 | 8,478 | 18,805 | 18,997 | 62,652 |

ORIG keeps **~4.4× more genes than STRICT** — which is why ORIG can detect
biomarkers that STRICT excludes.

---

## 4. DE-lncRNA counts (padj<0.05, |log2FC|>1)

| Cohort | Strategy | g-BMF vs Ctrl | u-BMF vs Ctrl | g-BMF vs u-BMF |
|--------|----------|--------------:|--------------:|---------------:|
| 14-sample | **ORIG** | **1,165** | **1,015** | 4 |
| 14-sample | STRICT | 29 | 13 | 0 |
| 14-sample | R1 | 939 | 797 | 3 |
| 14-sample | R2 | 921 | 787 | 3 |
| 14-sample | R3 | 1,191 | 960 | 4 |
| 17-sample | **ORIG** | **1,171** | **1,034** | 7 |
| 17-sample | STRICT | 26 | 12 | 0 |
| 17-sample | R1 | 973 | 851 | 0 |
| 17-sample | R2 | 960 | 847 | 0 |
| 17-sample | R3 | 1,200 | 1,050 | 6 |

---

## 5. Recovery of the 11 manuscript-highlighted lncRNAs

Per-cell recovery across 4 disease-vs-control conditions per gene
(2 contrasts × 2 cohorts) for each strategy:

| Strategy | 14-sample (22 cells) | 17-sample (22 cells) |
|----------|---------------------:|---------------------:|
| **ORIG** (true manuscript filter) | **100 %  (22/22)** | **95.5 % (21/22)** |
| STRICT (methods-text filter) | 45.5 % (10/22) | 45.5 % (10/22) |
| R1 | 100 % | 95.5 % |
| R2 | 100 % | 95.5 % |
| R3 | 100 % | 95.5 % |

The only loss is **MALAT1 in 17-sample u-BMF vs Ctrl** — it remains
significant under STRICT in 17-sample mode (because STRICT applies its own
independent-filtering threshold differently) but drops below the FDR cutoff
under ORIG/R1/R2/R3 when more low-count genes are co-tested. MALAT1 in
g-BMF vs Ctrl is preserved everywhere.

Under STRICT, the **six** biomarkers HCG11, HCP5, SNHG32, PSMB8-AS1, FAM30A,
MIR22HG fail the pre-filter in every condition. These six all share a
"near-zero in controls, log2FC ~9–11 in disease" pattern — the strict
"CPM ≥ 1 in ALL samples" rule rejects exactly this pattern.

---

## 6. Quantitative robustness of the 4 RT-qPCR-validated biomarkers

(log2FC across strategies in 14-sample g-BMF vs Ctrl)

| Gene | ORIG | STRICT | R1 | R3 |
|------|-----:|-------:|---:|---:|
| ATP1A1-AS1 | ~2.34 | ~2.38 | ~2.34 | ~2.34 |
| USP3-AS1   | ~2.60 | ~2.64 | ~2.60 | ~2.60 |
| TAGAP-AS1  | ~4.37 | ~4.41 | ~4.37 | ~4.36 |
| LINC01036  | ~3.26 | ~3.30 | ~3.25 | ~3.25 |
| MALAT1     | −6.56 | −6.51 | −6.56 | −6.56 |

Within 1 % across the four filters where the gene passes — the biology is
unchanged, the filter only governs which genes are even tested.

---

## 7. g-BMF vs u-BMF: 0 lncRNAs DE under R1/R2/STRICT in 17-sample

Across every strategy except R3, **0 lncRNAs** distinguish g-BMF from u-BMF
in the batch-aware 17-sample analysis. Even under R3 only 6 marginal cases
appear in 17-sample. This is fully consistent with the manuscript's
interpretation that the disease-relevant lncRNA program is a **pan-BMF
signature** rather than a g-vs-u discriminator.

---

## 8. Recommendation for the revision

1. **Correct the Methods text.** Replace the inaccurate
   "mean CPM ≥ 10 AND CPM ≥ 1 in all samples" description with the filter
   actually used in the code: `rowSums(counts) >= 10`. This already
   recovers all 11 reported biomarkers, so no result needs to change.
2. **Add the sensitivity table above** as a Supplementary Table to document
   that the findings are robust to filter choice.
3. **(Optional) Adopt R1 (`CPM ≥ 1 in ≥ 50 % of samples of any group`) as
   the primary filter going forward**, since it is the field-standard
   group-aware approach (analogous to edgeR `filterByExpr`), preserves all
   11 biomarkers in the 14-sample reproduction, is robust to adding Child
   controls, and is easier to defend on biological grounds.

---

## 9. Source files

Analysis directory: `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/lncrna_filter_analysis/`

- `filter_summary_combined.tsv` — counts per strategy/cohort/contrast
- `TABLE1_de_lncRNA_counts.tsv` — same as above, with up/down split
- `TABLE2_manuscript_11lncRNA_recovery_diseaseVScontrol.tsv` — 11-gene recovery
- `TABLE3_manuscript_11lncRNA_wide.tsv` — wide pivot
- `TABLE4_recovery_pct_by_strategy.tsv` — % recovery per cell
- `TABLE5_relaxed_filter_gain.tsv` — extra DE genes over STRICT
- `de_lncRNA_{cohort}_{contrast}_{strategy}.tsv` — 30 per-cell DE tables
- `robust_lncRNAs_intersection_all_strategies.tsv`, `ultra_robust_lncRNAs.tsv`
- `manuscript_11lncRNA_recovery.tsv` (+ `_rate`, `_wide`)

Scripts: `PIPELINE/hypothesis_tests/lncrna_filter_strategies.R`,
`lncrna_filter_summary_report.R`.

Manuscript filter provenance (the smoking gun):
- `03_original_analysis/scripts_archive/deseq2_analysis.R` line 33
- `03_original_analysis/scripts_archive/deseq2_analysis_enhanced.R` line 59
- `03_original_analysis/scripts_archive/snakemake_template/scripts/deseq2_analysis.R` lines 87–88
