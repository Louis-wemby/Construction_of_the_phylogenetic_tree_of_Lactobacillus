from ete3 import Tree, TreeStyle, TextFace, NodeStyle
from skbio import DistanceMatrix
from skbio.tree import nj
import pandas as pd
import numpy as np
from io import StringIO

# Step 1: Read the ANI percentage matrix and convert into a distance matrix
ani_df = pd.read_csv(".ani_output/ani_matrix.csv", index_col=0).astype(float)

# Step 2: Convert to symmetric distance matrix if needed
ani_sym = (ani_df + ani_df.T) / 2
distance_df = 1 - ani_sym
np.fill_diagonal(distance_df.values, 0)

# Step 3: Construct the NJ tree
dm = DistanceMatrix(distance_df.values, ids=distance_df.index.tolist())
tree_skbio = nj(dm)

# Step 4: Export the tree to Newick format and save
newick_io = StringIO()
tree_skbio.write(file=newick_io)
newick_str = newick_io.getvalue()
with open("./NJ_tree.nwk", "w") as f:
    f.write(newick_str)

# Step 5: Read the tree with ete3
t = Tree(newick_str)

# Step 6: Standardize terminal branch length
max_dist = max([t.get_distance(leaf) for leaf in t.iter_leaves()])
for leaf in t.iter_leaves():
    current_dist = t.get_distance(leaf)
    additional_dist = max_dist - current_dist
    leaf.dist = additional_dist

# Step 7: Set the tree style to circular
ts = TreeStyle()
ts.mode = "c"  # Set to circular mode (radial treemap)
ts.show_leaf_name = False  # Manually add strain name
ts.show_scale = False
ts.branch_vertical_margin = 15  # Control the node spacing
ts.rotation = 90  # Rotate the treemap so the top is facing up (optional, to adjust the visual effect)
ts.arc_start = 0  # The original angle of the circular tree
ts.arc_span = 360  # Circular tree covers a full 360 degrees

# Step 8: Change the color of the branch to red
for node in t.traverse():
    if node.is_root():
        continue
    nstyle = NodeStyle()
    nstyle["hz_line_color"] = "#FF0000"  # The horizontal line color is red
    nstyle["vt_line_color"] = "#FF0000"  # The vertical line color is red
    nstyle["hz_line_width"] = 2  # Branch Thickness
    nstyle["vt_line_width"] = 2
    node.set_style(nstyle)

# Step 9: Add strain names and adjust position
for node in t.iter_leaves():
    node.add_face(TextFace(node.name, fsize=8), column=0, position="aligned")  # Font size 8, name aligned

# Step 10: Render the treemap
t.render("NJ_tree_circle.pdf", tree_style=ts, w=800, units="mm")

print("Circular NJ tree PDF generated.")
