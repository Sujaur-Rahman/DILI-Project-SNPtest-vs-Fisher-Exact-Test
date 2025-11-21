# DILI WES vs 1000G EUR: Fisher Exact Test vs SNPTEST

This repository contains a workflow to identify significant SNPs in
137 **DILI / DILI-ALH** whole-exome sequencing (WES) samples by
comparing them to **503 European (EUR)** samples from the 1000 Genomes Project.

Two statistical approaches are used and later compared:

1. **Fisher exact test** using manually derived positive/negative allele counts
2. **SNPTEST** (planned) using genotype likelihood–based association testing

At this stage, the repository documents the **Fisher exact test pipeline**
for **chr1** in both case and control cohorts, plus the structure for adding the
SNPTEST pipeline later.

---

## Project goals

- Filter and harmonize WES VCF data from:
  - 137 DILI + DILI-ALH cases
  - 503 1000G EUR controls
- For each SNP, compute:
  - **Positive allele count** (ALT allele, case/control)
  - **Negative allele count** (non-ALT, case/control)
- Build 2×2 contingency tables and perform **Fisher exact tests**
- Run **SNPTEST** on the same data and compare:
  - Overlap of significant SNPs
  - Effect size direction and p-values

---

## Data description

**Cases**  
- 137 DILI and DILI-ALH WES samples  
- Per-sample VCFs (one per individual), initially containing all chromosomes  
- In the pipeline we restrict to **chr1** to save time for method comparison.

**Controls (1000G EUR)**  
- 503 European samples from 1000 Genomes  
- VCFs merged and filtered to match the same SNP-level inclusion criteria as cases  
- Also restricted to **chr1** for the comparison step.

> **Note**  
> Raw VCFs and intermediate data files are **not** stored in this repository
> because of size and privacy. Only scripts and documentation are included.

---

## Dependencies

Core tools used in the Fisher pipeline:

- `bash`
- `awk`
- `grep`
- [`bcftools`](http://www.htslib.org/doc/bcftools.html)
- `bgzip` / `tabix` (from htslib)
- R or Python (for running Fisher exact tests on the final contingency tables;
  the exact script can be added later)

Planned for SNPTEST pipeline:

- [`SNPTEST`](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)

---

## Repository layout

```text
scripts/
  fisher/
    01_filter_DP_GQ.sh          # Per-sample DP/GQ filter
    02_bgzip.sh                 # Compress VCFs
    03_snp_filter.sh            # Keep biallelic SNPs only
    04_merge.sh                 # Merge samples into cohort VCF
    05_snp_filter_merged_cohort.sh
    05a_unzip.sh
    06_rm_metadata.sh           # Remove '##' metadata lines
    07_snp_count.sh             # Count samples carrying the SNP
    08_no_snp_count.sh          # Compute NO_SNP_COUNT per variant
    09_allele_count.sh          # Count AC=1 and AC=2 per variant
    10_conti_pos_neg.sh         # Derive POSITIVE / NEGATIVE allele counts
  snptest/
    (to be added)

docs/
  workflow_overview.md          # High-level overview
  fisher_pipeline_chr1.md       # Detailed Fisher pipeline notes

results/
  fisher_chr1/                  # Fisher test outputs (not tracked)
  snptest_chr1/                 # SNPTEST outputs (not tracked)

data/
  cases/                        # Case VCFs (not tracked)
  controls/                     # 1000G VCFs (not tracked)
