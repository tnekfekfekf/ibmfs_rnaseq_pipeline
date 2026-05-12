# Working Conventions

## Rule: Auto-backup to Git after every significant change

**Established:** 2026-05-12
**Repo:** https://github.com/tnekfekfekf/ibmfs_rnaseq_pipeline

### Claude's working protocol

After ANY of the following:
- New analysis script created
- Pipeline option tested/validated
- Significant config change
- New report/documentation written
- Major bug fix

→ **Auto-execute:**

```bash
cd /Users/jaeeunyoo/Documents/ibmfs_rnaseq_pipeline
# Copy/move new files into appropriate folders
git add -A
git commit -m "<descriptive message>"
git push
```

### What goes in git

✅ **Always commit:**
- Shell scripts (`.sh`)
- R scripts (`.R`)
- Python scripts (`.py`)
- Documentation (`.md`)
- Configuration (`.yaml`, `.toml`, etc.)
- Small reference TSV/CSV (<1 MB)

❌ **Never commit (in `.gitignore`):**
- BAM, FASTQ, GTF files (large data)
- Count matrices, RDS, RData (intermediate results)
- Plot PNGs (regenerable)
- `.DS_Store`, `._*`, swap files

### Folder structure to maintain

| Folder | Purpose |
|---|---|
| `PIPELINE/` | **Validated, current pipeline** (don't break this) |
| `hypothesis_tests/` | Experimental options (v1-vN tests) |
| `03_original/` | Manuscript-era reference scripts (read-only history) |
| `04_revision/` | Revision pipeline scripts |
| `docs/` | Reports, READMEs, decisions log |

### Commit message conventions

- `add:` new files (e.g. `add: v6 hypothesis test with -O overlap`)
- `fix:` bug fix
- `update:` modify existing (e.g. `update: PIPELINE_README with v5 results`)
- `validate:` empirical verification (e.g. `validate: v3 vs manuscript 95% match`)
- `doc:` documentation only

### When to remind user

If user asks me to:
- Run a long analysis → backup BEFORE starting (so config is saved)
- Modify a validated script → commit current state first
- Add new analysis → commit + push at completion

### Reproducibility checklist (before claiming "done")

1. ✓ Script runs from clean state
2. ✓ Inputs/outputs documented in script header
3. ✓ Validated against expected output (manuscript or known result)
4. ✓ Committed + pushed to git
5. ✓ README updated if user-facing
