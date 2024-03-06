import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree  # Assuming you have ete3 installed, you can install it using: pip install ete3

# Read the tab-separated file into a DataFrame
df = pd.read_csv('combined_full_tables.tsv', sep='\t')

# Assuming 'phylogenetic_tree_file' is the file containing the Newick format tree
# You might need to adjust the file path accordingly
phylogenetic_tree_file = 'pruned-tree.txt'

# Read the Newick format phylogenetic tree from the file
with open(phylogenetic_tree_file, 'r') as file:
    phylogenetic_tree_string = file.read()

# Modify the Newick string to remove bootstrap values
# Assuming the format is (Species1:0.01[100],Species2:0.02[95]) - remove the [100] and [95]
phylogenetic_tree_string_cleaned = ''.join([c.split('[')[0] for c in phylogenetic_tree_string.split(']')])

# Assign the cleaned phylogenetic tree string to the DataFrame
df['phylogenetic_tree'] = phylogenetic_tree_string_cleaned

# Convert the 'status' column to a binary column (1 for presence, 0 for absence)
df['presence'] = 1

# Pivot the DataFrame to create a presence/absence matrix
presence_matrix = df.pivot_table(index='Gene_ID', columns=['Species', 'Status'], values='presence', aggfunc='max', fill_value=0)

# Filter the presence_matrix for genes with status 'Complete'
complete_genes_matrix = presence_matrix.xs(key='Complete', level='Status', axis=1, drop_level=True)

# Parse the tree and get the tree structure
tree = Tree(phylogenetic_tree_string_cleaned, format=1)

# Create a dictionary to map species to the corresponding phylogenetic position
phylo_position_dict = {leaf.name: leaf.get_distance(tree) for leaf in tree.iter_leaves()}

# Create a custom function to order the columns based on the phylogenetic position
def order_by_phylogeny(x):
    return phylo_position_dict[x]

# Order the columns based on phylogenetic position
column_order = sorted(complete_genes_matrix.columns.get_level_values(0), key=order_by_phylogeny)
complete_genes_matrix = complete_genes_matrix[column_order]

# Cluster rows and columns with inverted color map
clustered_matrix = sns.clustermap(complete_genes_matrix, cmap="viridis_r", method='average', metric='euclidean', figsize=(35, 35), vmin=0, vmax=1)

# Increase font size and rotate tick labels for X-axis
plt.setp(clustered_matrix.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=15)

# Increase font size for Y-axis
plt.setp(clustered_matrix.ax_heatmap.yaxis.get_majorticklabels(), fontsize=15)

# Display the plot
plt.title('Gene Presence/Absence Matrix - Status: Complete')

# Save the plot as a PDF
#plt.savefig('your_plot_filename.pdf', format='pdf')

# Show the plot
plt.show()
