import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from itertools import combinations
from sklearn.model_selection import train_test_split

# Step 1. Input data
X = pd.read_csv("./filtered_copy_number_matrix.tsv", sep='\t', index_col=0)
if not X.index[0].isdigit():
    X = X.T

ani_df = pd.read_csv("./ani_output/ani_matrix.csv", index_col=0)

# Step 2. Construct strain pair features (L1 distance) and labels (ANI)
pair_data = []
labels = []

for i, j in combinations(X.index, 2):
    diff = np.abs(X.loc[i] - X.loc[j])  # L1 distance
    pair_data.append(diff.values)
    labels.append(ani_df.loc[i, j])  # The distance can be inverted using 100 - ani_df.loc[i, j]

X_pairs = np.array(pair_data)
y_pairs = np.array(labels)

# Step 3. Train the Random Forest Regression Model
X_train, X_test, y_train, y_test = train_test_split(X_pairs, y_pairs, random_state=42)
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
print("R² score:", r2_score(y_test, y_pred))

# Step 4. Get feature importance
importances = model.feature_importances_
genes = X.columns
top_n = 30  # Adjustable

top_idx = np.argsort(importances)[::-1][:top_n]
top_genes = genes[top_idx]
top_importances = importances[top_idx]

# Convert to Series with index
top_df = pd.Series(top_importances, index=top_genes, name="Importance")

# 5. Calculate the coexistence ratio (proportion of non-zero strains)
nonzero_ratio = (X[top_df.index] > 0).sum(axis=0) / X.shape[0]

# 6. Secondary screening: only retaining top genes present in ≥ N% of strains
filtered_df = top_df[nonzero_ratio >= N%]  # Set as your own threshold

print(f"\nOriginal Top Gene Number: {len(top_df)}")
print(f"Number of genes with coexistence ratio ≥ N%: {len(filtered_df)}")
print(filtered_df)

# 7. Visualize the conserved important coexisting genes
plt.figure(figsize=(8, max(4, len(filtered_df) * 0.3)))
sns.barplot(x=filtered_df.values, y=filtered_df.index, palette="viridis")
plt.title("Top Genes with High Prevalence")
plt.xlabel("Feature Importance")
plt.ylabel("Gene")
plt.tight_layout()
plt.show()
