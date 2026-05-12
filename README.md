# IBMFS RNA-seq Pipeline

Pediatric Bone Marrow Failure RNA-seq analysis pipeline for Scientific Reports manuscript revision.

**Code backup rule:** Every script/analysis change is auto-committed to this repo. See [CONVENTIONS.md](CONVENTIONS.md).

## Quick start

```bash
cd /Volumes/ExtremeSSD/ibmfs/PIPELINE   # actual data location
bash PIPELINE/01_run_fc_per_sample.sh   # featureCounts per-sample
Rscript PIPELINE/02_merge_counts.R       # merge into matrix
Rscript PIPELINE/03_run_deseq2_manuscript.R   # eAT5.R-style DEG
Rscript PIPELINE/04_validate_pipeline.R       # CPM comparison vs manuscript
```

See `docs/PIPELINE_README.md` for full documentation.

## Folder structure

| Folder | Contents |
|---|---|
| `PIPELINE/` | **Validated** featureCounts + DESeq2 scripts |
| `docs/` | Pipeline docs + control-comparison report |
| `03_original/` | Manuscript-era scripts (Aug-Oct 2025) |
| `04_revision/scripts/` | Revision pipeline (HISAT2 + featureCounts + 8 R analyses) |
| `hypothesis_tests/` | Tested featureCounts options (v1-v5), batch effect analysis |
| `CONVENTIONS.md` | Working rules (auto-backup, commit conventions) |

## Validated featureCounts command

```bash
featureCounts -T 8 -p -s 2 -t exon -g gene_id \
  -a gencode.v44.annotation.no_rRNA.gtf \
  -o output.txt  BAM
```

**Validated** (2026-05-12) against manuscript CPM:
- Pearson r = 0.9998 (4 samples — perfect)
- Pearson r = 0.93-0.99 (10 samples — 95% match)

Other tested options (all worse or no improvement):
- `-M --primary` (v4) → CPM correlation **worse**
- `-B -C` (v2) → slightly worse
- `-t gene` (revision) → much worse
- Ensembl GRCh38.110 GTF (v5) → same as Gencode v44 (identical assignment)

## Reference data (NOT in git)

Large files on ExtremeSSD:
- `/Volumes/ExtremeSSD/ibmfs/01_raw_data/bams_macrogen_hisat2/` — 14 manuscript BAMs
- `/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/aligned/Child{1,2,3}.bam` — 3 public controls
- `/Volumes/ExtremeSSD/download_ssd/gencode.v44.annotation.gtf` — reference GTF

## Key results

See `docs/CONTROL_COMPARISON_REPORT.md`:
- 4-way DEG analysis: 2 controls vs 5 controls (naive vs batch-aware) vs 3 public-only
- Conclusion: adding public controls with `~ cohort + group` design recovers 92% of manuscript DEGs (robust)
- Adding controls naively (no batch correction) drops DEG count 60-80%
