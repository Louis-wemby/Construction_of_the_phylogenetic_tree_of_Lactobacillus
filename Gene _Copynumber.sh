#!/bin/bash
## Set up necessary environment variables
BIN=/home/your_hostname/miniconda3/envs/your_environment/bin
mkdir -p Genome_copynumber
cd Genome_copynumber
PROJECT_HOME=$PWD
mkdir -p Genmeo_data
mkdir -p Reference_genome
DATA=$PWD/Genome_data
RESULT=$PWD/results
mkdir -p $RESULT

## Set up the reference genome catalog
REFERENCE_GENOME=$PROJECT_HOME/Reference_genome/your_Reference_genome.fna
REFERENCE_GFF=$PROJECT_HOME/Reference_genome/your_Reference_genome.gff

## Alignment
mkdir -p $RESULT/aligned
for sample in $DATA/*.fna
do
    sample_name=$(basename $sample .fna)
    SAMPLE_READS=$DATA/${sample_name}.fna
    ALIGNMENT_OUTPUT=$RESULT/aligned/${sample_name}_aligned.bam
    bwa mem $REFERENCE_GENOME $SAMPLE_READS > $ALIGNMENT_OUTPUT.sam
    samtools view -Sb $ALIGNMENT_OUTPUT.sam > $ALIGNMENT_OUTPUT.bam
    samtools sort $ALIGNMENT_OUTPUT.bam -o $ALIGNMENT_OUTPUT.sorted.bam
    samtools index $ALIGNMENT_OUTPUT.sorted.bam
    bedtools coverage -a $REFERENCE_GFF -b $ALIGNMENT_OUTPUT.sorted.bam > $RESULT/copy_number/${sample_name}_coverage.txt
done

## Calculate the copy number of each gene in each sample
for coverage_file in $RESULT/copy_number/*_coverage.txt
do
    sample_name=$(basename "$coverage_file" _coverage.txt)
    avg_coverage=$(awk '$3 == "gene" {sum+=$11; count++} END {if (count > 0) print sum/count; else print "ERROR: No genes found"}' "$coverage_file")
    awk -v sample="$sample_name" -v avg_cov="$avg_coverage" '
    $3 == "gene" {
       split($9, arr, ";");
       gene_name = "unknown";
       for (i in arr) {
            if (arr[i] ~ /^Name=/) {
                gene_name = substr(arr[i], 6);
            }
       }
        coverage = $(NF-2);
        copy_number = coverage / avg_cov;
        print sample, gene_name, copy_number;
    }' "$coverage_file" > "$RESULT/copy_number/${sample_name}_copy_number.tsv"

    echo "Processed $sample_name"
done

## Merge the copy number results
INPUT_DIR="results/copy_number"
OUTPUT_FILE="results/merged_copy_number.tsv"
awk '1 {print $2}' $INPUT_DIR/*.tsv | sort | uniq > gene_names.txt
echo -e "Sample\t$(paste -sd '\t' gene_names.txt)" > $OUTPUT_FILE

for file in $INPUT_DIR/*.tsv; do
    sample_name=$(basename "$file" copy_number.tsv)
    declare -A copy_number_map
    awk '{copy[$2] = $4} END {for (gene in copy) print gene, copy[gene]}' "$file" > temp_copy.txt
    printf "%s" "$sample_name" > temp_row.txt
    while read -r gene; do
        value=$(grep -w "$gene" temp_copy.txt | awk '{print $2}')
        printf "\t%s" "${value:-0}" >> temp_row.txt
    done < gene_names.txt
    echo "" >> temp_row.txt

    cat temp_row.txt >> $OUTPUT_FILE
done

rm gene_names.txt temp_copy.txt temp_row.txt

echo "✅ mission complete！result file：$OUTPUT_FILE"
