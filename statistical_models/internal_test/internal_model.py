import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Function to read genes from a file
def read_genes(file_path):
    with open(file_path, 'r') as file:
        genes = [line.strip() for line in file.readlines()]
    return genes

# Function to read SRR names and generate file paths
def generate_file_paths(base_path, srr_file):
    with open(srr_file, 'r') as f:
        srr_names = [line.strip() for line in f.readlines() if line.strip()]
    file_paths = [base_path.replace("SRR", srr_name) for srr_name in srr_names]
    return file_paths

# Function to calculate weighted scores
def calculate_weighted_scores(files, gene_list, top_n=10):
    all_data = pd.DataFrame()
    gene_file_count = {}

    for file in files:
        if os.path.exists(file):
            data = pd.read_csv(file)
            data['SourceFile'] = os.path.basename(file)  # Add a column to identify the source file
            all_data = pd.concat([all_data, data], ignore_index=True)
            
            # Update gene file count
            for gene in data['Gene Name'].unique():
                if gene in gene_list:
                    if gene in gene_file_count:
                        gene_file_count[gene].add(file)
                    else:
                        gene_file_count[gene] = {file}
        else:
            print(f"Warning: {file} does not exist and will be skipped.")
    
    # Ensure column names are correct
    all_data.columns = ['Variant', 'Classification', 'Position', 'Gene Name', 'SourceFile']
    
    # Filter out relevant columns and genes in the gene list
    gene_data = all_data[all_data['Gene Name'].isin(gene_list)][['Gene Name', 'Classification', 'Position', 'SourceFile']]
    
    # Add placeholder for 'variant' and 'pathogenic_score'
    gene_data['variant'] = all_data['Variant']
    gene_data['pathogenic_score'] = all_data['Position']  # This is just an example; update based on your data
    
    # Group by gene and calculate counts and summary statistics within each file
    gene_stats = gene_data.groupby('Gene Name').agg({
        'variant': 'count',
        'pathogenic_score': ['mean', 'std'],
        'SourceFile': lambda x: len(x.unique())  # Count the number of unique files each gene appears in
    })

    gene_stats.columns = ['variant_count', 'pathogenic_mean', 'pathogenic_std', 'file_count']
    
    # Correctly update the file_count
    for gene in gene_stats.index:
        if gene in gene_file_count:
            gene_stats.at[gene, 'file_count'] = len(gene_file_count[gene])

    # Calculate weighted scores
    gene_stats['normalized_score'] = (
        (gene_stats['pathogenic_mean'] - gene_stats['pathogenic_mean'].min()) /
        (gene_stats['pathogenic_mean'].max() - gene_stats['pathogenic_mean'].min())
    )
    gene_stats['weighted_score'] = (
        gene_stats['variant_count'] * gene_stats['normalized_score'] * gene_stats['file_count']
    )

    #print(gene_stats[['variant_count', 'normalized_score', 'file_count', 'weighted_score']].head())

    min_weighted_score = gene_stats['weighted_score'].min()
    max_weighted_score = gene_stats['weighted_score'].max()
    gene_stats['normalized_weighted_score'] = (
        (gene_stats['weighted_score'] - min_weighted_score) /
        (max_weighted_score - min_weighted_score)
    )
    
    overall_scores = gene_stats['normalized_weighted_score']  # Directly use normalized weighted scores
    
    top_genes = overall_scores.nlargest(top_n)
    
    return overall_scores, gene_stats, top_genes

# Example usage
base_path = "/path/to/ngs/results/SRR/ANNOTATION/Franklin_Results_Without_Benign_LikelyBenign.csv"
srr_file = "/path/to/srr_names.txt"
gene_list_file = "/path/to/gene_list.txt"

# Read the gene list
gene_list = read_genes(gene_list_file)

# Generate the list of file paths
files = generate_file_paths(base_path, srr_file)

# Calculate scores and visualize results
overall_scores, gene_stats, top_genes = calculate_weighted_scores(files, gene_list, top_n=10)

#print("Top Genes:")
#print(top_genes)

# Directory for plots
plot_dir = r"/path/to/output/plots"

# Helper function to save plots in multiple formats
def save_plot(plot_name):
    for ext in ['pdf', 'png', 'jpeg']:
        plt.savefig(os.path.join(plot_dir, f"{plot_name}.{ext}"))

# Plot: Top genes by normalized weighted score
plt.figure(figsize=(12, 8))
top_genes.sort_values(ascending=True).plot(kind='barh', color='skyblue', edgecolor='black')
plt.xlabel('Normalized Weighted Score', fontsize=14)
plt.ylabel('Gene', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
save_plot('top_genes')
plt.close()

# Create a pivot table for the heatmap
pivot_table = gene_stats.pivot_table(index='Gene Name', values='normalized_weighted_score', aggfunc='mean')

# Enhanced Heatmap Plot
plt.figure(figsize=(12, 8))
sns.heatmap(pivot_table, cmap='coolwarm', annot=True, fmt=".2f", cbar_kws={'label': 'Normalized Weighted Score'})
plt.xlabel('Gene', fontsize=14)
plt.ylabel('Gene Name', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
save_plot('enhanced_heatmap')
plt.close()

# Plot histogram of normalized weighted scores
plt.figure(figsize=(12, 8))
plt.hist(overall_scores.astype(float), bins=20, color='skyblue', edgecolor='black', alpha=0.7)
plt.xlabel('Normalized Weighted Score', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
save_plot('histogram')
plt.close()

# Plot box plot of normalized weighted scores
plt.figure(figsize=(12, 8))
sns.boxplot(x=gene_stats['normalized_weighted_score'], color='skyblue')
plt.xlabel('Normalized Weighted Score', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
save_plot('box_plot')
plt.close()

# Plot violin plot of normalized weighted scores by gene
plt.figure(figsize=(12, 8))
sns.violinplot(data=gene_stats[['normalized_weighted_score']], color='skyblue', inner='quartile')
plt.ylabel('Normalized Weighted Score', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
save_plot('violin_plot')
plt.close()

# Clustered heatmap of gene-gene correlations
numeric_data = gene_stats.select_dtypes(include=[np.number])

plt.figure(figsize=(12, 8))
sns.heatmap(numeric_data.corr(), cmap='coolwarm', annot=True, fmt=".2f", cbar_kws={'label': 'Correlation'})
plt.tight_layout()
save_plot('clustered_heatmap')
plt.close()

# Save results to CSV
output_csv_path = r"/path/to/output/Final_result_of_analysis.csv"
gene_scores_df = pd.DataFrame({'Gene Name': overall_scores.index, 'Score': overall_scores.values})
gene_scores_df.to_csv(output_csv_path, index=False)
print(f"Gene scores saved to {output_csv_path}")

# Save results to new CSV with Exist_Number column
output_csv_path_new = r"/path/to/output/Myself_Final_result_of_analysis.csv"
gene_scores_df_new = pd.DataFrame({
    'Gene Name': overall_scores.index, 
    'Score': overall_scores.values,
    'Exist_Number': gene_stats['file_count'].values
})
gene_scores_df_new.to_csv(output_csv_path_new, index=False)
print(f"Gene scores saved to {output_csv_path_new}")

# Save results to a new CSV with Gene Name, Variant Count, Normalized Score, File Count
output_csv_path_weight_score = r"/path/to/output/Weight_Score.csv"

# Include gene names explicitly in the DataFrame
weight_score_df = gene_stats[['variant_count', 'normalized_score', 'file_count', 'weighted_score']]
weight_score_df['Gene Name'] = overall_scores.index  # Add gene names as a new column

# Reorder columns to match the desired format
weight_score_df = weight_score_df[['Gene Name', 'variant_count', 'normalized_score', 'file_count', 'weighted_score']]

# Save the DataFrame to CSV
weight_score_df.to_csv(output_csv_path_weight_score, index=False)
print(f"Weight scores saved to {output_csv_path_weight_score}")

# Sort overall_scores in descending order
overall_scores = overall_scores.sort_values(ascending=False)

# Save results to CSV
output_csv_path = r"/path/to/output/descending_Final_result_of_analysis.csv"
gene_scores_df = pd.DataFrame({'Gene Name': overall_scores.index, 'Score': overall_scores.values})
gene_scores_df.to_csv(output_csv_path, index=False)
print(f"Gene scores saved to {output_csv_path}")
