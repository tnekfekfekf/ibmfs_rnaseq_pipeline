# IBMFS RNA-seq Pipeline

Pediatric Bone Marrow Failure RNA-seq analysis pipeline (Scientific Reports manuscript revision).

**Code backup rule:** Every change is auto-committed to this repo (see [CONVENTIONS.md](CONVENTIONS.md)).

## 🎯 Current strategy: Option 1 (manuscript matrix preserved)

We use the **authentic manuscript count matrix** as baseline and add new samples with consistent quantification on top.

### Authoritative count matrix

**Location:** `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt`

- 62,653 genes × 19 samples (14 manuscript + Child1/2/3 + AA-LES + AA-CSB)
- Raw integer counts
- **Reproduces 100% of manuscript published DEG values** (verified 2026-05-12)
- 4158/4158 sig DEG match
- All highlighted genes (HCG11, HCP5, SNHG32, ATP1A1-AS1, USP3-AS1, etc.) PERFECT match

See `MANUSCRIPT_COUNTS/README.md` (in data folder) for full details.

## Pipeline scripts

| Script | Purpose |
|---|---|
| `PIPELINE/05_use_manuscript_matrix.R` | **Primary workflow** — DESeq2 on authentic matrix (reproduces manuscript) |
| `PIPELINE/06_add_new_samples.sh` | Add NEW samples (quantify + append to matrix) |
| `PIPELINE/01_run_fc_per_sample.sh` | featureCounts (per-sample, used for new samples) |
| `PIPELINE/02_merge_counts.R` | Merge per-sample counts |
| `PIPELINE/03_run_deseq2_manuscript.R` | DESeq2 manuscript-style (eAT5.R replica) |
| `PIPELINE/04_validate_pipeline.R` | CPM-level validation vs manuscript |

## Quick reference

### Reproduce manuscript DEG values
```bash
Rscript PIPELINE/05_use_manuscript_matrix.R
# Output: /Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/deseq2_replica/DE_*.tsv
```

### Add new samples
```bash
bash PIPELINE/06_add_new_samples.sh /path/to/bams sample01 sample02
# Output: combined matrix with new sample columns
```

### Validated featureCounts command (for new samples)
```bash
featureCounts -T 8 -p -s 2 -t exon -g gene_id \
  -a gencode.v44.annotation.no_rRNA.gtf  BAM
```

## Folder structure

| Folder | Contents |
|---|---|
| `PIPELINE/` | Validated pipeline scripts (6 files) |
| `docs/` | Analysis reports + decisions log |
| `03_original/` | Manuscript-era scripts (Aug-Oct 2025) |
| `04_revision/scripts/` | Revision pipeline (HISAT2 + featureCounts + 8 R analyses) |
| `hypothesis_tests/` | featureCounts option tests (v1-v6), batch effect analysis |

## Key documents

- **`docs/PIPELINE_README.md`** — Full pipeline documentation
- **`docs/HEAD_TO_HEAD_COMPARISON.md`** — Manuscript original vs our reproduction
- **`docs/A_VS_C_COMPARISON.md`** — 2 controls vs 5 controls (with public Child added)
- **`docs/PUBLIC_CONTROLS_DECISION.md`** — Batch-aware analysis decision
- **`docs/COHORT_VALIDITY_DEEP_DIVE.md`** — Public MNC vs internal aspirate biology
- **`docs/DESEQ2_VERSION_TEST.md`** — DESeq2 v1.48 = v1.50 (version drift ruled out)
- **`docs/MANUSCRIPT_GENES_COMPARISON.md`** — 19 highlighted genes reproduction

## Validation summary

| Test | Result |
|---|---|
| DESeq2 1.48 vs 1.50 reproducibility | IDENTICAL ✅ |
| Manuscript matrix → DESeq2 → published DEG | 100% match ✅ |
| Our v3 reconstruction (independent quantify) | ~95% similarity, 70% recall at padj<0.05 |
| Adding public Child controls (batch-aware) | 92% DEG overlap with original |

## Reference data location

Large files on ExtremeSSD external drive:
- BAMs: `/Volumes/ExtremeSSD/ibmfs/01_raw_data/bams_macrogen_hisat2/`
- GTF: `/Volumes/ExtremeSSD/download_ssd/gencode.v44.annotation.gtf`
- **Manuscript count matrix: `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt`** ⭐

Not in git (see `.gitignore`).
