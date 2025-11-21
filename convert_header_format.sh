#cd /media/iscsi13T/home/mohammad/new_analysis_wes_4datasets/new_snptest/convert_missing_gt_to/work

# Create clean VCFv4.2 copies without touching originals
for F in 1KG.EUR.chr1.common.vcf.gz DILI.cases.chr1.shared.00.vcf.gz; do
  BASENAME=$(basename "$F" .vcf.gz)
  OUT="${BASENAME}.v42.vcf.gz"

  # Build a fresh header: enforce ##fileformat=VCFv4.2 on line 1,
  # keep all other header lines, and then append the body
  {
    echo "##fileformat=VCFv4.2"
    bcftools view -h "$F" | grep -v '^##fileformat='
    bcftools view -H "$F"
  } | sed 's/\r$//' | bgzip -c > "$OUT"

  tabix -f -p vcf "$OUT"
done

# Quick check (should show ##fileformat=VCFv4.2)
#bgzip -cd 1KG.EUR.chr1.shared.v42.vcf.gz | head -3
#bgzip -cd DILI.cases.chr1.shared.v42.vcf.gz | head -3
