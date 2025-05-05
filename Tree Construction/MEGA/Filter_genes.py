import pandas as pd

# Read copy number matrix
df = pd.read_csv("./All_Gene_CopyNumbers.tsv", sep='\t', index_col=0)

# Calculate the frequency of 0 in each gene (row)
zero_freq = (df == 0).mean(axis=1)

# Filtering condition: retain genes with 0 occurrence frequency ≤ N%
filtered_df = df[zero_freq <= N%] # Set as your own threshold

# Save the results
filtered_df.to_csv("./filtered_copy_number_matrix.tsv", sep='\t', index=True)

print(f"Number of genes before filtering: {len(df)}")
print(f"Number of genes after filtering: {len(filtered_df)}")
print("Results saved to filtered_gene_matrix.tsv")

# Check filtered data
print("The first 5 rows of data：")
print(df.head())
