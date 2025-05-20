import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

# Step 1: Read ANI Matrix
ani_df = pd.read_csv("./ani_output/ani_matrix.csv", index_col=0)
# Symmetric processing
ani_symmetric = (ani_df + ani_df.T) / 2
distance_df = 1 - ani_symmetric

# Step 2: Ensure that the matrix is ​​symmetric and has no missing values
assert (distance_df.values == distance_df.values.T).all(), "ANI matrix is ​​not symmetrical"
assert not distance_df.isnull().any().any(), "ANI matrix contains missing values"

# Step 3: Construct the tree using UPGMA clustering
linkage_matrix = linkage(squareform(distance_df), method='average')

# Step 4: Convert to Newick format
def linkage_to_newick(linkage_matrix, labels):
    tree = to_tree(linkage_matrix, rd=False)
    def _to_newick(node):
        if node.is_leaf():
            return labels[node.id]
        left = _to_newick(node.left) if node.left else ""
        right = _to_newick(node.right) if node.right else ""
        if left and right:
            return f"({left},{right}):{node.dist}"
        return left or right
    return f"{_to_newick(tree)};"
newick_str = linkage_to_newick(linkage_matrix, distance_df.index.tolist())

# Step 5: Read the tree with ete3
t = Tree(newick_str)

# Step 6: Standardize terminal branch lengths
max_dist = max([t.get_distance(leaf) for leaf in t.iter_leaves()])
for leaf in t.iter_leaves():
    current_dist = t.get_distance(leaf)
    additional_dist = max_dist - current_dist
    leaf.dist = additional_dist

# Step 7: Set the circular tree style
ts = TreeStyle()
ts.mode = "c"  # circle
ts.show_leaf_name = False
ts.show_scale = False
ts.branch_vertical_margin = 15
ts.rotation = 90
ts.arc_start = 0
ts.arc_span = 360

# Step 8: Set red branches
for node in t.traverse():
    if node.is_root():
        continue
    nstyle = NodeStyle()
    nstyle["hz_line_color"] = "#FF0000"
    nstyle["vt_line_color"] = "#FF0000"
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_width"] = 2
    node.set_style(nstyle)

# Step 9: Add strain names
for node in t.iter_leaves():
    node.add_face(TextFace(node.name, fsize=6), column=0, position="aligned")

# Step 10: Render the treemap
t.render("UPGMA_tree_circle.pdf", tree_style=ts, w=800, units="mm")

print("Circular UPGMA tree PDF generated.")
