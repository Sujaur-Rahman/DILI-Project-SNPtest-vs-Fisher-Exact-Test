#!/bin/bash

#Script:1.filter_DP_GQ.sh#

for vcf in *.final.vcf
do
  awk -F'\t' '
    /^#/ {print; next}  # Print all header lines
    {
      split($10, arr, ":");  # Assuming the first sample is in column 10
      dp = arr[3];  # Adjust index if DP is not the 3rd subfield
      gq = arr[4];  # Adjust index if GQ is not the 4th subfield
      if (dp >= 20 && gq >= 30) print  # Print lines where DP >= 20 and GQ >= 30
    }
  ' $vcf > "${vcf%.vcf}.filtered.vcf"
done 
  
#Script:2.bgzip.sh#

for vcf in *.filtered.vcf
do
  bgzip -c $vcf > "${vcf}.gz"
  rm $vcf  # Optional: remove the original file after compressing
done

#Script:3.snp_filter.sh#

for input_file in *.final.filtered.vcf.gz;
do
    output_file="${input_file%.vcf.gz}.snps.max2.vcf.gz"
    bcftools view -v snps -m2 -M2 $input_file -Oz -o $output_file
    bcftools index $output_file
    echo "Filtering complete for $input_file. Output saved to $output_file"
done

#Script:4.merge.sh#

bcftools merge *.filtered.snps.max2.vcf.gz -o eur_merged_cohort.vcf.gz -O z

#Script:5.snp_filter_merged_cohort.sh#

bcftools view -v snps -m2 -M2 eur_merged_cohort.vcf.gz -Oz -o filtered_eur_merged_cohort.vcf.gz
bcftools index filtered_eur_merged_cohort.vcf.gz
echo "Filtering complete"

#Script:5.1.unzip.sh#

gunzip filtered_eur_merged_cohort.vcf.gz


#Script:6.rm_metadata.sh#

# Input VCF file
input_file="filtered_eur_merged_cohort.vcf"

# Output file
output_file="no_metadata_eur_merged_cohort.vcf"

# Use grep to filter out lines starting with ##
grep -v '^##' "$input_file" > "$output_file"

echo "Lines starting with ## removed. Output saved to $output_file"


#Script:7.SNP_COUNT.sh#

# Your original VCF file
vcf_file="no_metadata_eur_merged_cohort.vcf"

# Temporary file to hold intermediate results
tmp_file=$(mktemp)

# Output file with the added column for SNP count
output_file="2.modified_${vcf_file}"

# Copy the header to the output file and modify the #CHROM line to add new column name
grep '^#' $vcf_file | while read -r line; do
    if [[ $line == *CHROM* ]]; then
        echo -e "${line}\tSNP_COUNT"
    else
        echo "$line"
    fi
done > $output_file

# Process the VCF file, skipping header lines
grep -v '^#' $vcf_file | while read -r line; do
    # Extract all sample fields into an array
    IFS=$'\t' read -r -a fields <<< "$line"
    samples=("${fields[@]:9}") # Extract samples starting from the 10th field
    snp_count=0
    for sample in "${samples[@]}"; do
        if [[ $sample != "./.:.:.:.:." ]]; then # Check if SNP is present
            ((snp_count++))
        fi
    done
    # Add the SNP count as a new column
    echo -e "$line\t$snp_count"
done > $tmp_file
# Append the modified lines to the output file
cat $tmp_file >> $output_file

# Remove the temporary file
rm $tmp_file

echo "Modified file is saved as $output_file"

#Script:8.NO_SNP_COUNT.sh#

# Your modified VCF file with SNP_COUNT, AC=1_COUNT, and AC=2_COUNT columns
vcf_file="2.modified_no_metadata_eur_merged_cohort.vcf"

# Final output file with an additional column for the count of patients without each SNP
enhanced_output_file="NO_SNP_${vcf_file}"

# Copy the header to the enhanced output file and modify the #CHROM line to add new column name
grep '^#' $vcf_file | while read -r line; do
    if [[ $line == *CHROM* ]]; then
        echo -e "${line}\tNO_SNP_COUNT" # Directly modify the #CHROM line to include the new column
    else
        echo "$line" # Echo other header lines unmodified
    fi
done > $enhanced_output_file

# Process the VCF file, skipping header lines
grep -v '^#' $vcf_file | while read -r line; do
    IFS=$'\t' read -r -a fields <<< "$line"
    # The SNP_COUNT is now directly extracted from the 107th column, assuming the first column is index 1
    snp_count="${fields[146]}" # Adjust the index by subtracting 1 because arrays are 0-indexed
    no_snp_count=$((137 - snp_count))
    # Add the NO_SNP_COUNT as a new column and output the modified line
    echo -e "$line\t$no_snp_count"
done >> $enhanced_output_file

echo "Enhanced file with NO_SNP_COUNT is saved as $enhanced_output_file"


#Script:9.Allele_Count.sh#

vcf_file="NO_SNP_2.modified_no_metadata_eur_merged_cohort.vcf"
final_output_file="Allele_Count_${vcf_file}"

# Extract the header and modify the #CHROM line to add new columns for AC=1_COUNT and AC=2_COUNT
grep '^#' $vcf_file | while read -r line; do
    if [[ $line == *CHROM* ]]; then
        echo -e "${line}\tAC=1_COUNT\tAC=2_COUNT" # Directly modify the #CHROM line to include new columns
    else
        echo "$line" # Echo other header lines unmodified
    fi
done > $final_output_file

# Process the file with awk to calculate and add AC=1_COUNT and AC=2_COUNT columns
awk -F'\t' '
    !/^#/ { # Skip header lines
        ac1_count=0
        ac2_count=0
        for (i=10; i<=NF; i++) { # Assuming genotype info starts at column 10
            split($i, arr, ":")
            genotype=arr[1]
            if (genotype == "0/1") ac1_count++
            if (genotype == "1/1") ac2_count++
        }
        print $0"\t"ac1_count"\t"ac2_count
    }
' $vcf_file >> $final_output_file

echo "Final modified file is saved as $final_output_file"

#Script:10.conti_pos_neg.sh#

# Define the input VCF file and the output file
vcf_file="Allele_Count_NO_SNP_2.modified_no_metadata_eur_merged_cohort.vcf"
output_file="2.final_with_pos_neg.vcf"

# Process the header to include the new columns for Positive and Negative counts
grep '^#' "$vcf_file" | while read -r line; do
    if [[ $line == *CHROM* ]]; then
        echo -e "${line}\tPOSITIVE\tNEGATIVE"
    else
        echo "$line"
    fi
done > "$output_file"

# Read through the VCF file, skipping header lines, and process each data line
awk -F'\t' 'BEGIN {OFS="\t"}
!/^#/ {
  # Extract the necessary counts
  snp_count=$147  # Number of patients with the SNP
  no_snp_count=$148  # Number of patients without the SNP
  ac1_count=$149  # Number of patients with AC=1
  ac2_count=$150  # Number of patients with AC=2

  # Calculate Positive and Negative
  positive=ac2_count * 2 + ac1_count
  negative=no_snp_count * 2 + ac1_count

  # Print the original line with the Positive and Negative counts appended
  print $0, positive, negative
}' "$vcf_file" >> "$output_file"

echo "The file with added POSITIVE and NEGATIVE columns has been saved as $output_file"
































