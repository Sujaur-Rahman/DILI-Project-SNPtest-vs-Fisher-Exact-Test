# DILI / DILI-ALH WES vs 1000G EUR Controls  
**Comparison of Fisher Exact Test and SNPTEST for Significant SNP Discovery (chr1)**

## Overview

This repository contains a workflow to identify significant SNPs in **137 DILI and DILI-ALH WES samples** by comparing them to **503 European (EUR) control samples** from the 1000 Genomes Project.

Two statistical approaches are used:

1. **Fisher exact test**  
   - Manual derivation of allele counts (positive vs negative) from case and control VCFs  
   - Construction of 2×2 contingency tables per SNP  
2. **SNPTEST**  
   - Generation of `.gen` and `.sample` files  
   - SNP-level association testing using SNPTEST  

To speed up method comparison, the initial analyses are performed on **chromosome 1 (chr1)** only.

The main goals are:

- To obtain significant SNPs in DILI / DILI-ALH vs 1000G EUR controls  
- To compare the overlap and concordance of significant SNPs between **Fisher exact test** and **SNPTEST**  

---

## Data

- **Cases:** 137 WES samples (DILI + DILI-ALH)
- **Controls:** 503 EUR samples from 1000 Genomes Project
- **Chromosome:** chr1 (subset of the full WES data used for comparison)

> Note: Raw VCF files and individual-level data are **not** stored in the repository due to size and privacy constraints. Only scripts and documentation are included.

---

## Fisher Exact Test Pipeline (chr1)

For the Fisher approach, allele frequencies are counted manually from VCF data to obtain:

- **POSITIVE** allele counts (ALT alleles)
- **NEGATIVE** allele counts (non-ALT alleles)

This process is performed separately for **cases** and **controls**, and then merged to build 2×2 contingency tables.

### Step 0 — Restrict to chr1

Before running the scripts below, input VCFs are restricted to **chr1** (e.g., using `bcftools view -r 1` or equivalent).

---

### Script 1 — Filter by DP and GQ (`1.filter_DP_GQ.sh`)

```bash
#!/bin/bash

# Script:1.filter_DP_GQ.sh

for vcf in *.final.vcf
do
  awk -F'\t' '
    /^#/ {print; next}  # Print all header lines
    {
      split($10, arr, ":");  # Assuming the first sample is in column 10
      dp = arr[3];           # Adjust index if DP is not the 3rd subfield
      gq = arr[4];           # Adjust index if GQ is not the 4th subfield
      if (dp >= 20 && gq >= 30) print  # Keep lines where DP >= 20 and GQ >= 30
    }
  ' "$vcf" > "${vcf%.vcf}.filtered.vcf"
done
