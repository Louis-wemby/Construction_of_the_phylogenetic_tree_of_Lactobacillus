import os
from Bio import SeqIO
from collections import defaultdict
import pandas as pd

# ==== Parameter settings ====
gff_dir = "/your/full/path/to/annotation_gff"
fasta_dir = "/your/full/path/to/genome_fasta"
output_dir = "/your/full/path/to/output_results"

os.makedirs(output_dir, exist_ok=True)

# ==== GFF parsing function ====
def parse_gff_for_gene_names(gff_file):
    gene_id_to_name = {}
    cds_product_map = {}
    rrna_trna_product_map = {}

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]
            attr_dict = dict(
                [field.split("=", 1) for field in attributes.split(";") if "=" in field]
            )

            # Process gene, record gene ID → gene name (only gene field is allowed)
            if feature_type == "gene":
                gene_id = attr_dict.get("ID")
                gene_name = attr_dict.get("gene")  # Only the gene field is allowed
                gene_id_to_name[gene_id] = gene_name if gene_name else None

            # Process CDS and record Parent → product
            elif feature_type == "CDS":
                parent_id = attr_dict.get("Parent")
                product = attr_dict.get("product")
                if parent_id and product:
                    cds_product_map[parent_id] = product + " (product)"

            # Process rRNA/tRNA and record its Parent → product (here Parent refers to the gene level)
            elif feature_type in {"rRNA", "tRNA"}:
                gene_level_parent = attr_dict.get("Parent")
                product = attr_dict.get("product")
                if gene_level_parent and product:
                    rrna_trna_product_map[gene_level_parent] = product + " (product)"

    # Merge gene names, and if missing, add them from product
    for gene_id in list(gene_id_to_name.keys()):
        name = gene_id_to_name[gene_id]
        if not name:
            if gene_id in cds_product_map:
                gene_id_to_name[gene_id] = cds_product_map[gene_id]
            elif gene_id in rrna_trna_product_map:
                gene_id_to_name[gene_id] = rrna_trna_product_map[gene_id]
            else:
                del gene_id_to_name[gene_id]  # Delete if no valid names

    return gene_id_to_name

# The gene copy number for each sample was collected
all_sample_data = defaultdict(dict)

# ==== Iterate over all GFF/FASTA pairs ====
for filename in os.listdir(fasta_dir):
    if not filename.endswith(".fasta"):
        continue
    sample_id = filename.replace(".fasta", "")
    fasta_path = os.path.join(fasta_dir, f"{sample_id}.fasta")
    gff_path = os.path.join(gff_dir, f"{sample_id}.gff")

    if not os.path.exists(gff_path):
        print(f"❌ Missing GFF for {sample_id}, skipping.")
        continue

    print(f"✅ Processing {sample_id}...")

    id_to_name = parse_gff_for_gene_names(gff_path)
    sequence_counts = defaultdict(int)
    gene_sequences = {}

    # Parse FASTA file and calculate the copy number of each gene
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        sequence_counts[seq] += 1
        gene_sequences[record.id] = seq

    # Reverse mapping seq → gene_id (only one is kept)
    seen_seqs = {}
    for gene_id, seq in gene_sequences.items():
        if seq not in seen_seqs:
            seen_seqs[seq] = gene_id

    # Record the copy number of the sample into the dictionary
    for seq, count in sequence_counts.items():
        gene_id = seen_seqs.get(seq, "unknown").split("::")[0]
        gene_name = id_to_name.get(gene_id, gene_id)
        all_sample_data[gene_name][sample_id] = count

    # Output a single sample copy number file
    output_file = os.path.join(output_dir, f"{sample_id}_copy_number.tsv")
    with open(output_file, "w") as f:
        f.write("GeneID\tCopyNumber\n")
        for seq, count in sequence_counts.items():
            gene_id = seen_seqs.get(seq, "unknown").split("::")[0]
            gene_name = id_to_name.get(gene_id, gene_id)
            f.write(f"{gene_name}\t{count}\n")

# ==== Merge all sample data ====
# Fill the copy number data of all samples into the DataFrame
genes = list(all_sample_data.keys())
samples = sorted({sample for gene_data in all_sample_data.values() for sample in gene_data})

# Generate a DataFrame containing all samples and genes
df = pd.DataFrame(index=genes, columns=samples)

# Fill the DataFrame, if a sample is missing a gene, fill it with 0
for gene in genes:
    for sample in samples:
        df.loc[gene, sample] = all_sample_data[gene].get(sample, 0)

# Output merge result
merged_output = os.path.join(output_dir, "All_Gene_CopyNumbers.tsv")
df.to_csv(merged_output, sep="\t")

print(f"✅ Mission complete! Output saved to: {merged_output}")
