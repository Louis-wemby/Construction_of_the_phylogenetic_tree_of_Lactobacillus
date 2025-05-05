import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_gff_attributes(attribute_str):
    """Parse the attribute field of column 9 of GFF and return a dictionary"""
    attributes = {}
    for attr in attribute_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    return attributes


def extract_16s_with_copies(fna_path, gff_path, output_dir):
    """Extract 16S rRNA and keep only the first sequence"""
    try:
        # Read all sequences and store as a dictionary {seq_id: seq_record}
        seq_dict = {seq.id: seq for seq in SeqIO.parse(fna_path, "fasta")}
        print(f"✅ Read {len(seq_dict)} sequences: {list(seq_dict.keys())}")
        records = []

        with open(gff_path, 'r') as gff:
            rrna_found = False
            for line in gff:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 9:
                    continue

                # Check if it is rRNA type
                if cols[2].lower() == "rrna":
                    rrna_found = True
                    # Parsing attribute fields
                    attributes = parse_gff_attributes(cols[8])
                    product = attributes.get('product', '').lower()
                    strand = cols[6]
                    # DEBUG: Print detailed information about rRNA records
                    print(f"📝 Discovery of rRNA records: seq_id={cols[0]}, type={cols[2]}, strand={strand}, attributes={attributes}")

                    # Matching 16S rRNA
                    if '16s' in product or 'rrs' in product.lower():
                        seq_id = cols[0]
                        if seq_id not in seq_dict:
                            print(f"⚠️ Sequence ID does not match: {seq_id}")
                            continue

                        # Extract sequences
                        start = int(cols[3]) - 1  # GFF is 1-based, Python is 0-based
                        end = int(cols[4])
                        seq = seq_dict[seq_id].seq[start:end]
                        if strand == '-':
                            seq = seq.reverse_complement()  # Negative strand reverse complement

                        # Generate a unique identifier
                        locus_tag = attributes.get('locus_tag', 'unknown')
                        feature_id = f"{seq_id}_{locus_tag}_{start + 1}-{end}({strand})"
                        records.append(SeqRecord(seq, id=feature_id, description=""))
                        print(f"✅ Extract 16S rRNA: {feature_id}")

            if not rrna_found:
                print(f"⚠️ No rRNA records found in the GFF file: {os.path.basename(gff_path)}")

        # Only the first 16S rRNA sequence is retained
        if records:
            base_name = os.path.splitext(os.path.basename(fna_path))[0]
            selected_record = records[0]  # Keep the first sequence
            # Use base_name as sequence ID, suitable for MEGA
            selected_record.id = base_name
            selected_record.description = ""
            output_path = os.path.join(output_dir, f"{base_name}.fasta")  # Remove _16s suffix
            SeqIO.write([selected_record], output_path, "fasta")
            print(f"✅ 1 16S rRNA sequence ({len(records)} sequences found) has been retained to {output_path}")
            return selected_record
        else:
            print(f"⚠️ 16S rRNA not found: {os.path.basename(fna_path)}")
            return None

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(fna_path)}: {str(e)}")
        return None


def batch_extract_16s(fna_dir, gff_dir, output_dir="output"):
    """Batch process files in two directories and merge them into a single FASTA file"""
    os.makedirs(output_dir, exist_ok=True)
    paired_files = 0
    extracted_files = 0
    failed_files = 0

    # Get the base of the file name (without extension)
    fna_bases = {os.path.splitext(f)[0]: f for f in os.listdir(fna_dir) if f.endswith(('.fna', '.fasta'))}
    gff_bases = {os.path.splitext(f)[0]: f for f in os.listdir(gff_dir) if f.endswith('.gff')}

    total_fna_files = len(fna_bases)
    print(f"📊 {total_fna_files} .fna/.fasta files found")

    # Process all pairing files
    for base_name in fna_bases:
        if base_name in gff_bases:
            paired_files += 1
            print(f"📂 Processing Pairing Files: {base_name}")
            record = extract_16s_with_copies(
                fna_path=os.path.join(fna_dir, fna_bases[base_name]),
                gff_path=os.path.join(gff_dir, gff_bases[base_name]),
                output_dir=output_dir
            )
            if record:
                extracted_files += 1
                print(f"✅ Extracted Sequence: {base_name}")
            else:
                failed_files += 1
                print(f"⚠️ Sequence not extracted: {base_name}")
        else:
            print(f"⛔ Missing paired GFF file: {base_name}.gff")

    # Check and merge single strain FASTA files
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith(".fasta")]
    print(f"📊 {len(fasta_files)} single strain .fasta files found")
    merged_records = []
    seen_base_names = set()

    for fasta_file in fasta_files:
        fasta_path = os.path.join(output_dir, fasta_file)
        base_name = os.path.splitext(fasta_file)[0]  # Use the file name directly (without suffix)
        try:
            record = next(SeqIO.parse(fasta_path, "fasta"))
            # Handling duplicate base_name
            if base_name in seen_base_names:
                new_base_name = f"{base_name}_{len(seen_base_names)}"
                print(f"⚠️ Duplicate base_name found: {base_name}, renamed to {new_base_name}")
                record.id = new_base_name
            else:
                record.id = base_name
            record.description = ""
            merged_records.append(record)
            seen_base_names.add(base_name)
            print(f"✅ Add sequence from {fasta_file} to the merge list: {record.id}")
        except Exception as e:
            print(f"❌ Error reading {fasta_file}: {str(e)}")

    # Merge all sequences into one FASTA file
    if merged_records:
        merged_output_path = os.path.join(output_dir, "merged.fasta")  # 移除 _16s
        SeqIO.write(merged_records, merged_output_path, "fasta")
        print(f"✅ Merged {len(merged_records)} 16S rRNA sequences to {merged_output_path}")
    else:
        print(f"⚠️ No 16S rRNA sequences were found, unable to generate merge file")

    # Print Messages
    print(f"📊 Processing summary：")
    print(f"  - Total .fna/.fasta files: {total_fna_files}")
    print(f"  - Number of successfully paired files: {paired_files}")
    print(f"  - Number of files successfully extracted: {extracted_files}")
    print(f"  - Number of files that failed to be extracted: {failed_files}")
    print(f"  - Number of single strain .fasta files generated: {len(fasta_files)}")
    print(f"  - The number of sequences successfully merged: {len(merged_records)}")

    # Check for missing sequences
    if len(fasta_files) > len(merged_records):
        print(f"⚠️ Warning: {len(fasta_files)} single strain .fasta files were generated, but only {len(merged_records)} sequences were merged")
        fasta_base_names = {os.path.splitext(f)[0] for f in fasta_files}
        merged_base_names = {record.id for record in merged_records}
        missing_base_names = fasta_base_names - merged_base_names
        print(f"📉 Missing sequence (base_name): {missing_base_names}")


# Example
batch_extract_16s(
    fna_dir="./Sequence_data",
    gff_dir="./Reference_data",
    output_dir="./16s_results"
)
