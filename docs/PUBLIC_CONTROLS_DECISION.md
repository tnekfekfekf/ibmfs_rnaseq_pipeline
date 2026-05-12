# Public Healthy Controls (Child1/2/3) 추가 — 결론

**Date:** 2026-05-12
**Question (사용자):** Public Child 데이터셋의 healthy control 3명을 넣어도 매뉴스크립트 결과와 비슷한가? 아니면 batch effect가 너무 커서 못 넣는가?

**핵심 답:**

| 결과 | 답 |
|---|---|
| **Batch effect 존재?** | ✅ 있음 (PC2의 분산 17.6%를 cohort가 설명) |
| **단순 추가 (~group)로 넣으면?** | ❌ 결과 매우 바뀜 (DEG 60-80% 감소) |
| **Batch-aware (~cohort+group)로 넣으면?** | ✅ 매뉴스크립트와 거의 동일 (92% overlap) |

---

## 1. 4가지 분석 설계

| | n_samples | Controls | Design | 용도 |
|---|---|---|---|---|
| **A** | 14 | AA-PRO, AA-KEW (2 internal) | `~ group` | Manuscript replica |
| **B** | 17 | + Child1/2/3 (5 mixed) | `~ group` (naive) | "그냥 넣으면 어떻게 되나" |
| **C** | 17 | + Child1/2/3 (5 mixed) | `~ cohort + group` (batch-aware) | "제대로 넣었을 때" |
| **D** | 15 | Child1/2/3 만 (3 public) | `~ group` | "Public만 쓰면 어떻게 되나" |

---

## 2. DEG 결과 — 한 눈에

| | A: 2 internal | B: 5 ctrl naive | C: 5 ctrl batch-aware | D: 3 public only |
|---|---|---|---|---|
| **G-AA vs Ctrl mRNA** | **1425** | 535 ⬇⬇ | **1453** ≈A | 1140 |
| G-AA vs Ctrl lncRNA | 262 | 151 ⬇ | 265 ≈A | 331 ⬆ |
| **U-AA vs Ctrl mRNA** | **906** | 158 ⬇⬇⬇ | 908 ≈A | 1186 ⬆ |
| U-AA vs Ctrl lncRNA | 203 | 64 ⬇⬇ | 199 ≈A | 329 ⬆ |
| **G-AA vs U-AA total** | **13** | 14 | 9 | 5 |

### Manuscript (A)와 gene-level overlap

| 비교 | DEG 일치 (G-AA vs Ctrl) | 의미 |
|---|---|---|
| **A vs C (batch-aware)** | **1861 / 2014 = 92%** ✅ | **결과 robust** |
| A vs B (naive) | 360 / 2014 = 18% ⚠️ | 매우 다름 |
| A vs D (public only) | 205 / 2014 = 10% ⚠️ | 완전히 다른 결론 |

---

## 3. Batch effect — 정량적 증거

### PCA on top 2000 variable genes (VST normalized)

**Before batch correction:**

| PC | Variance | Cohort p-value | Group p-value | 해석 |
|---|---|---|---|---|
| PC1 | 58.3% | 0.245 (n.s.) | 0.501 (n.s.) | Sample-level noise |
| **PC2** | **17.6%** | **2×10⁻⁶ ⚠️** | 7×10⁻⁴ | **Cohort batch effect 강함** |
| PC3 | 6.3% | 0.71 | 0.59 | - |

**After ComBat-seq correction:**

| PC | Variance | Cohort p-value | Group p-value | 해석 |
|---|---|---|---|---|
| **PC1** | **80.1%** | 0.027 (mild) | **8×10⁻⁴** ✅ | **Biology dominates** |
| PC2 | 3.9% | 0.82 (n.s.) | 0.88 (n.s.) | Cohort 효과 사라짐 |

**시각화:** `batch_effect/PCA_before_correction.png` vs `PCA_after_combatseq.png`

---

## 4. 결론

### Q1: Batch effect가 너무 커서 못 넣는가?

❌ **아니오.** Batch effect는 있지만 (PC2의 17.6%), **ComBat-seq 또는 design에 cohort term 포함하면 잘 보정됨**.

### Q2: 그냥 넣으면 매뉴스크립트와 비슷한가?

❌ **아니오.** `~ group` naive design은:
- G-AA mRNA: 1425 → 535 (62% 감소)
- U-AA mRNA: 906 → 158 (83% 감소)
- Manuscript와 18%만 overlap → 사실상 다른 결과

→ **이렇게 넣으면 안 됨**

### Q3: 제대로 넣으면 매뉴스크립트와 비슷한가?

✅ **그렇다.** `~ cohort + group` batch-aware design:
- G-AA mRNA: 1425 → 1453 (≈동일)
- Manuscript와 **92% overlap** (1861/2014)
- G-AA vs U-AA contrast: 13 vs 9 (안정)

→ **이게 정답**

### Q4: 그러면 답변 어떻게?

**Reviewer에게 답할 수 있는 stance:**

> "We added 3 public PNBM controls (Child1/2/3) to expand sample size. Cohort batch effect was identified (17.6% of PC2 variance, ANOVA p=2×10⁻⁶). Using a cohort-aware design (`~ cohort + group`), our 17-sample analysis recovered **92% of the original DEG set**, validating the robustness of our biological conclusions. Naive pooling without batch correction reduced DEG counts by 60-83% and disagreed substantially with the original analysis (only 18% overlap), highlighting the importance of proper batch handling. The G-AA vs U-AA contrast (primary clinical comparison) remained stable across all sensitivity analyses (5-14 DEGs)."

---

## 5. 추천 분석 방향

✅ **Revision에서 사용할 것:**
1. **Primary analysis 그대로 유지** (A: 2 internal controls, manuscript 결과)
2. **Sensitivity analysis로 C 추가** (5 controls + `~ cohort + group`):
   - 결과 매우 robust 보임 (92% overlap)
   - Reviewer 2의 "n=2 너무 적음" 우려 해소

❌ **하지 말 것:**
- B (naive 5 control): 잘못된 통계, 결과 misleading
- D (public-only): cohort baseline 자체 다름, valid comparison 아님

---

## 6. 출력 파일

| 파일 | 내용 |
|---|---|
| `/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/control_comparison/DEG_summary_4analyses.tsv` | 4가지 분석 DEG 수치 |
| `dds_all.rds`, `deg_tables.rds` | DESeq2 객체 + DEG 테이블 |
| `batch_effect/PCA_before_correction.pdf` | Cohort 효과 시각화 |
| `batch_effect/PCA_after_combatseq.pdf` | ComBat-seq 보정 후 |
| `batch_effect/sample_distance_heatmap.pdf` | Sample-sample distance |
| `batch_effect/PC_anova_{before,after}.tsv` | ANOVA p-values |
| `batch_effect/counts_combatseq.rds` | Batch-corrected counts |

## 7. 재실행

스크립트:
- `hypothesis_tests/analyze_with_controls.R` — 4가지 DESeq2 분석
- `hypothesis_tests/batch_effect_analysis.R` — PCA + ComBat-seq + ANOVA
