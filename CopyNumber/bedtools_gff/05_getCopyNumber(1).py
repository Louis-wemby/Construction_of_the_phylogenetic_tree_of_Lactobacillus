from Bio import SeqIO
from collections import defaultdict
import re

# 从GFF文件中构建基因ID到通用名称的映射字典
def parse_gff_for_gene_names(gff_file):
    id_to_name = {}
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            attr_dict = {}
            for attr in re.split('[;=]', attributes):
                attr = attr.strip()
                if not attr:
                    continue
                if attr in ["ID", "gene", "Name", "locus_tag"]:
                    pairs = re.findall(r'([^=]+)=([^;]+)', attributes)
                    for k, v in pairs:
                        attr_dict[k] = v
                    break
            gene_id = attr_dict.get("ID")
            gene_name = (
                attr_dict.get(";gene") or  # 优先geneID
                attr_dict.get(";Name") or  # 之后是Name
                attr_dict.get(";locus_tag") or  # 再之后是locus_tag
                gene_id  # 如果都没有就使用默认gene_id
            )
            if gene_id and gene_name:
                id_to_name[gene_id] = gene_name
    return id_to_name

# 从GFF文件中构建映射
id_to_name = parse_gff_for_gene_names("D:/Louis/大创＋生竞-乳酸杆菌进化树/Data/Lactiplantibacillus plantarum/Genome data/Ref_Data/DSM20174.gff")  # 替换为你的 GFF 路径
print(id_to_name)
# 统计序列出现次数
sequence_counts = defaultdict(int)
gene_sequences = {}
for record in SeqIO.parse("D:/Louis/大创＋生竞-乳酸杆菌进化树/Data/Lactiplantibacillus plantarum/Genome data/DSM20174.fasta", "fasta"):
    seq = str(record.seq).upper()
    sequence_counts[seq] += 1
    gene_sequences[record.id] = seq

# 反向建立 seq -> gene_id 的唯一映射，只保留一个 gene_id 输出
seen_seqs = {}
for gene_id, seq in gene_sequences.items():
    if seq not in seen_seqs:
        seen_seqs[seq] = gene_id

# 输出每个基因的拷贝数
with open("gene_copy_numbers.tsv", "w") as f:
    f.write("GeneID\tCopyNumber\n")
    for seq, count in sequence_counts.items():
        gene_id = seen_seqs.get(seq, "unknown").split("::")[0]
        gene_name = id_to_name.get(gene_id, gene_id)
        f.write(f"{gene_name}\t{count}\n")
