import pandas as pd

input_csv = "combined_gene_copies.csv"
output_csv = "combined_gene_copies_wide.csv"
df = pd.read_csv(input_csv)
# 长变宽
pivot_df = df.pivot_table(
    index='Genome',
    columns='GeneID',
    values='CopyNumber',
    fill_value=0,
).reset_index()  # 将索引转换为普通列
pivot_df.to_csv(output_csv, index=False)






