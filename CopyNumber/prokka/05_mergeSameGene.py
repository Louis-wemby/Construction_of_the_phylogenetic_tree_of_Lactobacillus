import pandas as pd

input_csv = "combined_gene_copies_wide.csv"
output_csv = "combined_gene_copies_wide_merge.csv"
df = pd.read_csv(input_csv)
# 合并相同基因
cols_to_remove = []  # 存储需要删除的列
for col in df.columns:  # 遍历所有列
    if '_' in col:
        parts = col.split('_')
        # 确保是a_b格式且只包含一个下划线
        if len(parts) == 2:
            prefix = parts[0]
            # 如果前缀列不存在则创建
            if prefix not in df.columns:
                df[prefix] = 0
            # 数值相加（处理NaN值）
            df[prefix] = df[prefix].fillna(0) + df[col].fillna(0)
            cols_to_remove.append(col)
df.drop(columns=cols_to_remove, inplace=True)  # 删除原始带下划线的列
df.to_csv(output_csv, index=False)
