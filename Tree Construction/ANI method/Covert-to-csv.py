import pandas as pd

ani_matrix = pd.read_csv("./ani_output/ANIb_percentage_identity.tab", sep="\t", index_col=0)
ani_matrix.to_csv("./ani_output/ani_matrix.csv")

print(ani_matrix.head())
