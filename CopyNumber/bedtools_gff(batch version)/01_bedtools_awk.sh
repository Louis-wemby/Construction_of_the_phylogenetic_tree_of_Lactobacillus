#!/bin/bash
cd /path/to/your/Ref_Data

for gff in *.gff; do
    bed="${gff%.gff}.bed"
    awk -v OFS='\t' '!/^#/ && $3=="gene" {
        split($9, attr, ";");
        gene_id="unknown";
        for (i in attr) {
            if (attr[i] ~ /^ID=/) {
                split(attr[i], id, "=");
                gene_id=id[2];
                break;
            }
        }
        print $1, $4-1, $5, gene_id, ".", $7
    }' "$gff" > "$bed"
    echo "Converted $gff → $bed"
done
