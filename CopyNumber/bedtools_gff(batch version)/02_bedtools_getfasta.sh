#!/bin/bash
# set your path
fna_dir="/path/to/your/FNA_data"
bed_dir="/path/to/your/BED_data"
output_dir="./FASTA"

# create output folder（if not exist）
mkdir -p "$output_dir"

# traverse .fna files
for fna in "$fna_dir"/*.fna; do
    base=$(basename "$fna" .fna)
    bed="$bed_dir/${base}.bed"
    fasta="$output_dir/${base}.fasta"

    if [[ -f "$bed" ]]; then
        bedtools getfasta -fi "$fna" -bed "$bed" -fo "$fasta" -s -name
        echo "✅ Extracted: $base → $fasta"
    else
        echo "⚠️ No BED file for $base — skipped."
    fi
done
