import pandas as pd
import argparse

parser= argparse.ArgumentParser(add_help=False)
parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "Giving one file with 3 tab separated columns of busco results, where first is gene id, second is frag/missing/complete and third is species name, parse this file to find genes complete, missing or fragmented in all species. Can create lists of these genes as well, need to uncomment lines 28,29,30 for it.") 
parser.add_argument("-f", help= "-f: 3 column input file as described above", required = "True")
args = parser.parse_args()

# Assuming you have a file with three columns: 'gene_name', 'status', 'species'
file_path = args.f

# Read the file into a DataFrame
df = pd.read_csv(file_path, names=['gene_name', 'status', 'species'], sep="\t")

# Pivot the DataFrame to have 'gene_name' as rows, 'species' as columns, and 'status' as values
pivot_df = df.pivot_table(index='gene_name', columns='species', values='status', aggfunc='first')

# Identify genes that are missing, fragmented, or complete in all three species
missing_genes = pivot_df.loc[(pivot_df == 'Missing').all(axis=1)].index
fragmented_genes = pivot_df.loc[(pivot_df == 'Fragmented').all(axis=1)].index
complete_genes = pivot_df.loc[(pivot_df == 'Complete').all(axis=1)].index

# Convert indices to lists. 
missing_genes_list = missing_genes.tolist()
fragmented_genes_list = fragmented_genes.tolist()
complete_genes_list = complete_genes.tolist()

#complete_genes.to_series().to_csv('complete_genes-inters.txt', header=False, index=False, sep='\n')
#missing_genes.to_series().to_csv('missing_genes-inters.txt', header=False, index=False, sep='\n')
#fragmented_genes.to_series().to_csv('fragmented_genes-inters.txt', header=False, index=False, sep='\n')

# Create a DataFrame with the counts
result_df = pd.DataFrame({
    'Genes Complete in All Species': [len(complete_genes)],
    'Genes Missing in All Species': [len(missing_genes)],
    'Genes Fragmented in All Species': [len(fragmented_genes)],
})

# Save the result DataFrame
result_df.to_csv("Buscos_intersections.txt", index=False, sep="\t")
