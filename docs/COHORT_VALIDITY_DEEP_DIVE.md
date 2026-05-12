# Public vs Internal Control — Deep Dive on Cohort Equivalence

**Date:** 2026-05-12
**Question:** "Public controls (MNC isolated) vs internal controls (aspirate) + different depth — is `~cohort + group` correction really sufficient?"

**Honest answer:** Partial. The design works for the sensitivity analysis, but interpretation requires nuance.

## 1. Depth (good news)

| | Avg library size |
|---|---|
| Internal aspirate (AA-PRO, AA-KEW) | 29.1 M assigned reads |
| Public MNC (Child1/2/3) | 29.1 M assigned reads |

Initial concern about 3x depth difference was misleading — Child BAMs are larger pre-quantification but after `-s 2 -t exon` filtering, both cohorts yield similar effective library sizes.

## 2. Cell composition differences (HUGE)

| Marker | Internal CPM | Public CPM | Ratio | Interpretation |
|---|---|---|---|---|
| **TRBC2** (T cell) | 1.24 | 81.44 | **65x ↑ public** | MNC enriched for T cells |
| TRBC1 (T cell) | 2.45 | 55.36 | 22x ↑ public | |
| CD3D, TRAC | low | high | 3-7x ↑ public | T cell signature |
| LYZ (Monocyte) | 1161 | 8983 | 7.7x ↑ public | MNC enriched for monocytes |
| S100A8/9 | 316/395 | 1624/1717 | 4-5x ↑ public | |
| **ALAS2** (Erythroid) | **10347** | 402 | **26x ↑ internal** | Aspirate has erythroid cells |
| **SLC4A1** (Erythroid) | **11190** | 872 | **13x ↑ internal** | |
| ITGA2B (Megakaryocyte) | 94 | 17 | 5.7x ↑ internal | Aspirate has megakaryocytes |

**Internal (aspirate) = bone marrow ALL cells: erythroid + myeloid + megakaryocyte + lymphoid**

**Public (MNC) = Ficoll-isolated mononuclear cells: enriched for lymphoid + monocytes, depleted erythroid/megakaryocytes**

These are fundamentally different sample types, not just technical batches.

## 3. Direct DEG comparison: public vs internal Control

Treated cohort as DEG variable:
- **3,718 genes** differentially expressed (padj<0.05, |LFC|>1) — **15.3%** of tested genes
- 3,211 with |LFC|>2 (13.2%)
- 1,579 with |LFC|>5 (6.5%)
- Top hits include HLA-DRA, RPS18, RPS25, HLA-DPA1 with LFC 15-17

This is far beyond "batch noise" — it's substantial biological difference.

## 4. Why does ~cohort + group design still yield A ≈ C?

The 92% gene-level overlap between Analysis A (14 samples) and Analysis C (17 samples, batch-aware) doesn't mean public ≡ internal. It reflects analysis structure:

**Analysis C composition:**
- G-AA: 4 internal aspirate
- U-AA: 8 internal aspirate
- Control: 2 internal aspirate + 3 public MNC

**What ~cohort + group does:**
- cohort term absorbs all systematic public-vs-internal differences (cell composition, etc.)
- group estimates are computed AFTER cohort-adjustment
- Effectively asks: "within the same cohort, how does G-AA differ from Control?"

**Key insight:** Since G-AA, U-AA, and 2 of 5 Control samples are all internal, the **within-cohort biology** is anchored by internal samples (14/17 = 82%). The 3 public samples augment Control group statistically but don't dictate biology.

→ **Adding public controls preserves the internal-anchored analysis structure.**

## 5. The flip side — Analysis D evidence

When public controls are used **alone** (no internal):
- Analysis D: 3 Child only as controls (15 samples total)
- DEG count: 1140 mRNA, 331 lncRNA G-AA vs Ctrl
- **Gene overlap with manuscript A: only 10% (205/2014)** ⚠️

→ This is the technical proof that public ≠ internal. Different baseline transcriptome.

## 6. Honest interpretation matrix

| Claim | Status |
|---|---|
| Public Child controls are biologically equivalent to internal controls | ❌ FALSE — cell composition very different |
| The cohort effect is just technical batch | ❌ FALSE — it's also cell-composition |
| `~cohort + group` removes all confounding | ⚠️ MOSTLY — depends on whether biology correlates with composition |
| Adding public controls (with batch-aware design) preserves manuscript conclusions | ✅ TRUE — but driven by internal-sample dominance |
| C analysis proves robustness to control choice | ✅ TRUE for this dataset structure (internal-dominant) |
| Manuscript could have used public-only controls | ❌ FALSE (D analysis only 10% overlap) |

## 7. Improved reviewer-response framing

```
Sample composition differences:
- Internal samples (n=14): bone marrow aspirate, contains all hematopoietic lineages
- Public samples (n=3, GSE147523): mononuclear cells isolated by density gradient,
  enriched for lymphoid and monocytic lineages, depleted for erythroid/megakaryocyte

Quantitative cohort effect:
- 15.3% of genes show |LFC|>1 difference between cohorts at padj<0.05
- Cell-type markers (TRBC2: 65x↑public, ALAS2: 26x↑internal) confirm composition difference

Statistical handling:
- We did NOT pool public + internal as equivalent controls
- Instead, used `~ cohort + group` design where cohort absorbs sample-type variance
- Public samples augment statistical power of Control group estimation
- Internal samples (12/17 of non-public) anchor within-cohort biological estimates
- This yielded 92% DEG overlap with original 2-control analysis (median |ΔlogFC| = 0.003)

Interpretation:
- We do not claim public controls are biologically equivalent to internal
- We claim the original biological conclusions are robust to sample-size augmentation
  under appropriate statistical handling of cell-composition differences
- Public-only analysis (sensitivity D, n=3 public controls only) yielded only 10% 
  DEG overlap with original, confirming the cohorts are genuinely different
- This demonstrates the original choice of internal IDA/HUS controls was appropriate
  given absence of true healthy pediatric bone marrow aspirate datasets
```

## 8. Bottom line

✅ **Analysis C 결과는 valid** — manuscript의 결론이 robust함을 입증
⚠️ **단 해석은 정확해야** — public이 internal과 같다는 게 아니라, 내부 anchor 분석에서 추가 control이 결과를 깨지 않는다는 것
❌ **C 결과를 "public controls work as substitutes for internal"** 라고 주장하면 잘못

## Files

- `docs/public_vs_internal_ctrl_DEG.tsv` — full DEG list between cohorts
- `docs/cohort_effect_markers.tsv` — cohort coefficient for cell-type markers in ~cohort+group fit
- `hypothesis_tests/verify_cohort_validity.R` — reproducible analysis
