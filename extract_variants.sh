#!/bin/bash

# Usage: ./extract_variant.sh chr1:723052 sample.cram
# Output: prints VCF line to stdout if found

POS=$1
CRAM=$2

samtools mpileup -r chr1:723052-723053 ./crams/converted/Batch10-1710/Batch10-1710.converted.cram | cut -f5

bcftools mpileup -f "$REF" -r "$POS" -a FORMAT/DP -Ou "$CRAM" 2>/dev/null | \
bcftools call -mv -Ou 2>/dev/null | \
bcftools norm -f "$REF" -Oz 2>/dev/null | \
bcftools view -H -