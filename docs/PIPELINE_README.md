# IBMFS RNA-seq Pipeline — Full Documentation

**Last updated:** 2026-05-12
**Strategy: Option 1 — Use authentic manuscript count matrix as baseline**

---

## Background

After extensive investigation we determined:
1. The manuscript's count matrix was a **patchwork** — different samples quantified with different featureCounts options at different times (Aug 16, Sept 12, etc.)
2. A single featureCounts command can NOT exactly reproduce the manuscript matrix (~95% similarity ceiling)
3. The authentic manuscript count matrix was recovered on 2026-05-12 and stored at `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt`
4. DESeq2 on this matrix reproduces manuscript DEG values 100%

**Strategy:** Use this authentic matrix as the source of truth. Add new samples via consistent v3 pipeline (`-T 8 -p -s 2 -t exon -g gene_id`).

---

## Files

### `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/`

```
manuscript_count_matrix_19samples.txt  ⭐ Authentic source (62k genes × 19 samples)
all_samples_expression_table_CPM_*.txt  Normalized CPM (multiple versions)
all_samples_expression_table_final.txt  log-CPM
aa_final_cpm_anchor.xlsx                Excel with all DEG/GSEA results
AA RNA supplementary 260108.xlsx        Supplementary tables (sTable 1-11)
README.md                               Detailed file documentation
```

### `PIPELINE/` (this repo)

```
01_run_fc_per_sample.sh   featureCounts per-sample (v3 options)
02_merge_counts.R          Merge per-sample outputs
03_run_deseq2_manuscript.R DESeq2 manuscript-style (eAT5.R replica, from quantification)
04_validate_pipeline.R     CPM-level validation vs manuscript CPM
05_use_manuscript_matrix.R ⭐ Use authentic matrix directly (Option 1 primary workflow)
06_add_new_samples.sh      Add new samples → append to manuscript matrix
```

---

## Critical decisions documented

### Library / sequencing facts
- **Library type:** TruSeq Stranded Total RNA + Ribo-Zero (REVERSE-stranded → featureCounts `-s 2`)
- **Aligner:** HISAT2 v2.1.0 with `--dta --rna-strandness RF` (Macrogen)
- **Reference:** GRCh38, Gencode v44 (some Aug scripts used Ensembl GRCh38.110)
- **rRNA removal:** Custom GTF `gencode.v44.annotation.no_rRNA.gtf`

### featureCounts options (validated)
```bash
featureCounts -T 8 -p -s 2 -t exon -g gene_id \
  -a gencode.v44.annotation.no_rRNA.gtf \
  -o output.txt  BAM
```

| Option | Value | Why |
|---|---|---|
| `-T` | 8 | threads |
| `-p` | (set) | paired-end |
| `-s` | **2** | reverse-stranded (TruSeq Stranded) — required |
| `-t` | exon | mature transcript focus |
| `-g` | gene_id | gene-level counts |
| `-a` | no_rRNA.gtf | exclude rRNA |

**Tested and rejected:**
- `-M --primary` (v4): made worse
- `-B -C` (v2): mostly same
- `-O` overlap-allowed (v6): +35% reads, too much
- `-Q 10` MAPQ filter (v7): no effect
- `-t gene` (revision): too many intronic reads
- Ensembl GRCh38.110 GTF (v5): identical to Gencode v44

### DESeq2 design (manuscript replica)
```r
# eAT5.R Oct 20 2025
counts <- counts[rowSums(counts) > 0, ]   # minimal filter
dds <- DESeqDataSetFromMatrix(counts, meta, ~ group)   # ~ group only
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", g1, g2))
# Post-filter padj != NA, then subset by gene_type for mRNA/lncRNA reporting
```

**Versions tested:**
- DESeq2 v1.48 (Bioc 3.21, manuscript era) = v1.50.2 (current) — IDENTICAL output

---

## Workflow

### Reproduce manuscript results
```bash
Rscript PIPELINE/05_use_manuscript_matrix.R
# Reproduces: 2078 mRNA DEG (G-AA vs Ctrl), 1167 lncRNA, etc.
# 100% match with manuscript published DEG files
```

### Add new samples
```bash
# 1. Place BAM files in some directory
mkdir /path/to/new_bams
cp new_sample.bam /path/to/new_bams/new_sample_sorted.bam

# 2. Run pipeline
bash PIPELINE/06_add_new_samples.sh /path/to/new_bams new_sample01 new_sample02

# 3. Output:
# /Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/combined_matrix_with_new.txt
# (manuscript 19 samples + new samples)

# 4. Re-analyze
Rscript -e "
suppressPackageStartupMessages({library(DESeq2)})
# ... DESeq2 with combined matrix and appropriate design
"
```

### Caveats for Option 1

1. **Mixed quantification:** Existing 19 samples are patchwork; new samples will be v3-consistent. Expect 5-7% systematic offset between old/new samples.
2. **Batch-aware design recommended** if adding many new samples (`~ cohort + group` with cohort as old/new indicator).
3. **Strong signals robust:** LFC > 5 reproduces 100%; borderline (LFC 2-6) may differ.

---

## Lessons learned (don't repeat)

1. **Read user's existing code FIRST** — manuscript Methods text isn't the source of truth
2. **Library type → strand option**: TruSeq Stranded = `-s 2`, check BAM `--rna-strandness` header
3. **Reproduction priority:** code files → input data → manuscript text
4. **When count matrix is missing:** CPM-level comparison validates pipeline
5. **Per-sample sequential featureCounts** for progress visibility
6. **Local SSD for BAM/GTF** to avoid external drive I/O bottleneck

---

## Validation history

| Date | Test | Result |
|---|---|---|
| 2026-05-11 | Pipeline reproducibility (v3 vs manuscript CPM) | 95-100% Pearson |
| 2026-05-12 | Manuscript count matrix recovered | 100% DEG match |
| 2026-05-12 | DESeq2 1.48 vs 1.50 | Identical output |
| 2026-05-12 | Adding Child controls (batch-aware) | 92% DEG overlap with original |

---

## References

- Manuscript: Scientific Reports submission (Jan 2026)
- DEG result files: `/Volumes/ExtremeSSD/ibmfs/03_original_analysis/deg_results/`
- Authentic count matrix: `/Volumes/ExtremeSSD/ibmfs/MANUSCRIPT_COUNTS/manuscript_count_matrix_19samples.txt`
- Manuscript master pipeline (Aug 5): `03_original/master_pipeline_Aug5_BTui.sh`
- Manuscript final DESeq2 (Oct 20): `03_original/manuscript_final_eAT5.R`
