# IBMFS RNA-seq Pipeline

Pediatric Bone Marrow Failure RNA-seq analysis pipeline for Scientific Reports manuscript revision.

## Quick start

```bash
cd /Volumes/ExtremeSSD/ibmfs/PIPELINE   # actual data location
bash 01_run_fc_per_sample.sh            # featureCounts per-sample
Rscript 02_merge_counts.R               # merge into matrix
Rscript 03_run_deseq2_manuscript.R      # eAT5.R-style DEG
Rscript 04_validate_pipeline.R          # CPM comparison vs manuscript
```

See `docs/PIPELINE_README.md` for full pipeline documentation.

## Folder structure

| Folder | Contents |
|---|---|
| `PIPELINE/` | **Final validated pipeline** — featureCounts + DESeq2 scripts |
| `docs/` | Pipeline docs + control-comparison report |
| `03_original/` | Manuscript-era scripts (Aug-Oct 2025) |
| `04_revision/scripts/` | Revision pipeline (HISAT2 + featureCounts + R analyses) |
| `hypothesis_tests/` | featureCounts option tests (v1-v5), batch effect analysis |

## Validated featureCounts command

```bash
featureCounts -T 8 -p -s 2 -t exon -g gene_id \
  -a gencode.v44.annotation.no_rRNA.gtf \
  -o output.txt  BAM
```

Validated against manuscript CPM: Pearson r = 0.9998 (4 samples perfect), 0.93-0.99 (10 samples 95% match).

## Reference data location

Large data files (BAMs, GTFs, count matrices) are on **ExtremeSSD** external drive at:
- `/Volumes/ExtremeSSD/ibmfs/01_raw_data/bams_macrogen_hisat2/` — 14 manuscript BAMs
- `/Volumes/ExtremeSSD/ibmfs/04_revision_analysis/aligned/Child{1,2,3}.bam` — 3 public controls
- `/Volumes/ExtremeSSD/download_ssd/gencode.v44.annotation.gtf` — reference GTF

These are not in git (see `.gitignore`).

## Key analysis results

See `docs/CONTROL_COMPARISON_REPORT.md` for the manuscript replication + Child controls added analysis (4 designs compared, batch effect quantified).

