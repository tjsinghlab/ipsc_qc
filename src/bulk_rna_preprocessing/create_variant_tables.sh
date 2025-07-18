#!/bin/bash

input_dir="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling-outputs/RNA-variant-calling-gvcf-May2025-2"
output_dir="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/variant_tables_may2025"

mkdir -p "$output_dir"

for vcf in "$input_dir"/*.variant_filtered.vcf.gz; do
  base=$(basename "$vcf" .variant_filtered.vcf.gz)
  out_file="$output_dir/${base}_variants.txt"
  bcftools query -f '%CHROM\t%POS\t%DP\t[%AD{0}\t%AD{1}\n]' "$vcf" > "$out_file"
done