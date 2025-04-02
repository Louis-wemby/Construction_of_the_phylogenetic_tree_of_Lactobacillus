import hashlib
import os
from Bio import SeqIO
from collections import defaultdict


def parse_prokka_tsv(tsv_path):
    locus_to_gene = {}
    with open(tsv_path, 'r', encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')
        locus_idx = headers.index('locus_tag')
        gene_idx = headers.index('gene')
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < max(locus_idx, gene_idx) + 1:
                continue  # 跳过不完整行
            locus_tag = cols[locus_idx].strip()
            gene = cols[gene_idx].strip()
            if gene:  # 只保留有geneID的条目
                locus_to_gene[locus_tag] = gene
    return locus_to_gene

def process_genome(ffn_path, tsv_path, output_csv="combined_gene_copies.csv"):
    genome_name = os.path.splitext(os.path.basename(ffn_path))[0]  # 从文件名提取菌种名
    gene_data = []

    # 统计单个基因组内的拷贝数
    hash_counts = defaultdict(int)
    seq_hash_map = {}  # 记录基因ID与哈希的映射
    for record in SeqIO.parse(ffn_path, "fasta"):
        gene_id = record.id
        seq = str(record.seq).upper().replace(' ', '')
        seq_hash = hashlib.md5(seq.encode()).hexdigest()
        seq_hash_map[gene_id] = seq_hash
        hash_counts[seq_hash] += 1

    # 获取Prokka基因标识符到geneID的映射
    prokka_to_geneid = parse_prokka_tsv(tsv_path)

    # 按基因ID整理数据并添加菌种信息
    for prokka_id, seq_hash in seq_hash_map.items():
        copy_number = hash_counts[seq_hash]
        gene_id = prokka_to_geneid.get(prokka_id)
        if not gene_id:
            continue
        gene_data.append({
            "Genome": genome_name,
            "GeneID": gene_id,
            "CopyNumber": copy_number
        })

    # 追加写入总表
    with open(output_csv, 'a') as f:
        if os.stat(output_csv).st_size == 0:  # 文件为空时写入表头
            f.write("Genome,GeneID,CopyNumber\n")
        for entry in gene_data:
            line = f"{entry['Genome']},{entry['GeneID']},{entry['CopyNumber']}\n"
            f.write(line)


if __name__ == '__main__':
    ffn_files = "./ffn_file"
    tsv_files = "./tsv_file"
    for ffn_file in os.listdir(ffn_files):
        if ffn_file.endswith(".ffn"):
            file_name_with_extension = os.path.basename(ffn_file)  # 获取文件名（包含扩展名）
            file_name, extension = os.path.splitext(file_name_with_extension)  # 分割文件名和扩展名
            ffn_path = os.path.join(ffn_files, ffn_file)
            process_genome(ffn_path, f"{tsv_files}/{file_name}.tsv")
