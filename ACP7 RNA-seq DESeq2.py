import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from adjustText import adjust_text

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

class FeatureCountsProcessor:
    def __init__(self, file_path):

        self.file_path = file_path
        self.data = None
        self.counts = None
        self.metadata = None
        self.res = None  
        self.dds = None
    
    def load_data(self):

        try:
            self.data = pd.read_csv(
                self.file_path,
                sep="\t",  # or ',' for CSV; '\t' for TSV
                comment="#",  # Skip lines starting with '#'
                header=0,
            )
            print("Data successfully loaded.")
        except Exception as e:
            print(f"Error reading the file: {e}")
            self.data = None

    def clean_column_names(self, df):

        if df is None:
            print("DataFrame is empty. Cannot clean column names.")
            return

        try:
            # Apply the transformation to each column name
            df.columns = [col.split('_')[0].replace('./', '') if './' in col else col for col in df.columns]
            print("Column names cleaned.")
        except Exception as e:
            print(f"Error cleaning column names: {e}")

    def extract_counts_table(self):

        if self.data is None:
            print("Data not loaded. Please run load_data() first.")
            return

        try:
            # Ensure 'GeneId' column exists
            gene_id_column = self.data.columns[0]  # Typically 'Geneid'

            # Extract columns G to AD (Excel-style), i.e., index 6 to 30
            count_columns = self.data.columns[6:31]  # 6 to 30 inclusive
            # Build counts dataframe
            self.counts = self.data[[gene_id_column] + list(count_columns)].copy()
            
            # Remove 'CELE_' prefix from GeneId
            self.counts[gene_id_column] = self.counts[gene_id_column].str.replace('CELE_', '', regex=False)
            
            self.counts.set_index(gene_id_column, inplace=True)
            
            print("Counts table successfully extracted.")
            
            # Clean up the column names in counts table
            self.clean_column_names(self.counts)  # Clean column names in counts
        except Exception as e:
            print(f"Error extracting counts table: {e}")
            self.counts = None
            
    def filter_zero_rows(self, min_reads=10):

        if self.counts is None:
            print("Counts table not found. Run extract_counts_table() first.")
            return

        # Before filtering
        initial_shape = self.counts.shape

        # Filter genes with less than 'min_reads' across all samples
        self.counts = self.counts[self.counts.sum(axis=1) >= min_reads]

        # After filtering
        final_shape = self.counts.shape

        print(f"Filtered counts table: {initial_shape[0]} → {final_shape[0]} genes (removed {initial_shape[0] - final_shape[0]})")

    def create_metadata(self):

        if self.counts is None:
            print("Counts table not found. Run extract_counts_table() first.")
            return None

        try:
            sample_names = self.counts.columns.tolist()
            num_samples = len(sample_names)

            if num_samples != 25:
                print(f" Expected 25 samples, but found {num_samples}. Please verify the counts columns.")
                return None

            # Define condition labels based on sample index
            conditions = (
                ['Wildtype'] * 5 +
                ['ACP7'] * 5 +
                ['ACP7-7A'] * 5 
            )

            metadata = pd.DataFrame({
                'Sample': sample_names,
                'Condition': conditions
            })

            print("Metadata successfully created.")
            metadata.set_index("Sample", inplace=True)

            # Store metadata as an instance variable
            self.metadata = metadata

            return metadata
        except Exception as e:
            print(f"Error creating metadata: {e}")
            return None
    
    def run_deseq2_all(self, design_factors='Condition'):

        if self.counts is None:
            print("Counts table not found. Run extract_counts_table() first.")
            return None
    
        if self.metadata is None:
            print("Metadata not found. Run create_metadata() first.")
            return None
    
        try:
            # --- Step 1: Clean up data ---
            if 'WormScarlet_gene' in self.counts.index:
                self.counts = self.counts.drop(index='WormScarlet_gene')
                print("Removed 'WormScarlet_gene' from counts.")
    
            # --- Step 2: Filter for desired conditions ---
            filtered_metadata = self.metadata[
                self.metadata[design_factors].isin(['Wildtype', 'ACP7', 'ACP7-7A'])
            ]
            filtered_counts = self.counts[filtered_metadata.index].T
            print(f"Using samples from conditions: {filtered_metadata[design_factors].unique().tolist()}")
    
            # --- Step 3: Create and fit DESeq2 dataset ---
            self.dds = DeseqDataSet(
                counts=filtered_counts,
                metadata=filtered_metadata,
                design_factors=design_factors
            )
            print("DESeq2 DeseqDataSet object created successfully.")
    
            self.dds.deseq2()
            print("DESeq2 model fitted across all conditions.")
    
            return self.dds
    
        except Exception as e:
            print(f"Error during DESeq2 fitting: {e}")
            return None
    
    
    def run_pairwise_deseq2(self, conditions, design_factors='Condition'):

        if self.dds is None:
            print("DESeq2 dataset not found. Run run_deseq2_all() first.")
            return None
    
        try:
            if not isinstance(conditions, (list, tuple)) or len(conditions) != 2:
                print("Please provide exactly two conditions, e.g., ('ACP7', 'ACP7-7A').")
                return None
    
            condition_col = self.metadata[design_factors]
            if not all(c in condition_col.unique() for c in conditions):
                print(f"One or both specified conditions not found in metadata['{design_factors}']. "
                      f"Available: {condition_col.unique().tolist()}")
                return None
    
            contrast = (design_factors, conditions[0], conditions[1])
            print(f"Running contrast: {conditions[0]} vs {conditions[1]}")
    
            stat_res = DeseqStats(self.dds, contrast=contrast)
            stat_res.summary()
            self.res = stat_res.results_df
    
            # Save results
            cond_1_clean = conditions[0].replace(" ", "_")
            cond_2_clean = conditions[1].replace(" ", "_")
            filename = f'deseq2_results_{cond_1_clean}_vs_{cond_2_clean}.csv'
            self.res.to_csv(filename)
            print(f"Differential expression results saved to '{filename}'.")
    
            return self.res
    
        except Exception as e:
            print(f"Error during pairwise DESeq2 contrast: {e}")
            return None

            
def plot_sig_gene_heatmap(
    dds, 
    res, 
    condition_col='Condition', 
    conditions=None, 
    log2fc_threshold=1, 
    padj_threshold=0.05, 
    save_path="significant_genes_heatmap.png"
):
    # Ensure variance-stabilized counts are available
    if "vst_counts" not in dds.layers:
        dds.vst()

    # --- Filter significant genes ---
    sig_genes = res[
        (res["padj"].fillna(1) < padj_threshold) &
        (res["log2FoldChange"].abs() > log2fc_threshold)
    ].index

    if len(sig_genes) == 0:
        print("No significant genes found with the given thresholds.")
        return

    # --- Subset samples for the two conditions of interest ---
    if conditions is not None:
        if not isinstance(conditions, (list, tuple)) or len(conditions) != 2:
            raise ValueError("Please provide exactly two conditions (e.g., ['ACP7', 'Wildtype']).")

        # Identify samples matching the specified conditions
        condition_mask = dds.obs[condition_col].isin(conditions)
        selected_samples = dds.obs_names[condition_mask]
        print(f"Plotting heatmap for samples in conditions: {conditions}")
    else:
        # Default: include all samples
        selected_samples = dds.obs_names
        print("No specific conditions provided; using all samples.")

    # --- Extract VST data for significant genes and selected samples ---
    heat_data = pd.DataFrame(
        dds[:, sig_genes].layers["vst_counts"].T,
        index=sig_genes,
        columns=dds.obs_names
    )[selected_samples]

    # --- Generate clustered heatmap (z-score by gene) ---
    g = sns.clustermap(heat_data, z_score=0, cmap="RdYlBu_r", yticklabels=False)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.savefig(save_path, dpi=300)
    plt.close()

    print(f"Heatmap saved as {save_path}")


def plot_volcano_with_labels_and_save(results_df, log2fc_threshold=1, padj_threshold=0.05,
                                      top_n=5, save_path='volcano_plot.png'):

    # Set global font to Arial-like
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'DejaVu Sans']

    # Drop rows with missing padj or log2FC
    results_df = results_df.dropna(subset=['padj', 'log2FoldChange']).copy()
    
    # Replace zero padj values with smallest positive float64 value
    results_df.loc[results_df['padj'] == 0, 'padj'] = np.finfo(float).tiny
    
    # Calculate -log10(padj) for plotting
    results_df['-log10(padj)'] = -np.log10(results_df['padj'])
    results_df = results_df[np.isfinite(results_df['-log10(padj)'])]

    # Define significant genes using thresholds (padj < threshold)
    results_df['significant'] = (
        (results_df['log2FoldChange'].abs() >= log2fc_threshold) &
        (results_df['padj'] < padj_threshold)
    )

    # Get top N significant genes for labeling
    top_genes = results_df[results_df['significant']].nsmallest(top_n, 'padj')
    
    # Start plotting
    plt.figure(figsize=(6, 6))
    ax = sns.scatterplot(
        data=results_df,
        x='log2FoldChange',
        y='-log10(padj)',
        hue='significant',
        hue_order=[False, True],
        palette=['lightgrey', 'black'],
        legend=None
    )

    # Add threshold lines
    ax.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=1.5)

    # Add gene name labels to top N significant genes
    texts = []
    for i in range(len(top_genes)):
        texts.append(plt.text(
            x=top_genes.iloc[i]['log2FoldChange'],
            y=top_genes.iloc[i]['-log10(padj)'],
            s=top_genes.index[i],
            fontsize=14,
            weight='bold'
        ))
    
    # Prevent text overlap
    adjust_text(
        texts, 
        arrowprops=dict(arrowstyle='-', color='k')
    )

    # Aesthetics
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width=2)
    plt.xticks(size=14, weight='bold')
    plt.yticks(size=14, weight='bold')

    plt.xlabel("$log_{2}$ fold change", size=17)
    plt.ylabel("-$log_{10}$ FDR", size=17)

    # Save the plot
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    
    print(f"Plot saved to {save_path}")
    plt.close()

def save_significant_genes(res_df1, res_df2=None, log2fc_threshold=1, padj_threshold=0.05, output_csv='significant_genes.csv'):

    def process_df(res_df):
        df = res_df.copy()
        df.loc[df['padj'] == 0, 'padj'] = np.finfo(np.float64).tiny
        return df[(df['padj'] < padj_threshold) & (df['log2FoldChange'].abs() >= log2fc_threshold)]

    sig1 = process_df(res_df1)
    sig2 = process_df(res_df2) if res_df2 is not None else None

    # Intersection: significant in both
    if sig2 is not None:
        common_genes = sig1.index.intersection(sig2.index)
        intersection_df = sig1.loc[common_genes]
        print(f"Number of genes significant in both datasets (intersection): {len(intersection_df)}")
    else:
        intersection_df = sig1
        print(f"Number of significant genes: {len(intersection_df)}")

    # Union: significant in either
    if sig2 is not None:
        union_df = pd.concat([sig1, sig2]).loc[~pd.concat([sig1, sig2]).index.duplicated(keep='first')]
        print(f"Number of genes significant in either dataset (union): {len(union_df)}")
    else:
        union_df = sig1

    # Save intersection
    intersection_df[['log2FoldChange']].copy().rename_axis("gene").to_csv(output_csv)
    print(f"Saved intersection to '{output_csv}'")

    # Save union
    union_csv = output_csv.replace('.csv', '_union.csv')
    union_df[['log2FoldChange']].copy().rename_axis("gene").to_csv(union_csv)
    print(f"Saved union to '{union_csv}'")
    
def plot_pca_and_save(dds):
    # Apply variance stabilizing transformation (stores result in dds.layers['vst_counts'])
    dds.vst()
    # Use vst counts for PCA
    dds.X = dds.layers["vst_counts"]
    # Perform PCA
    sc.tl.pca(dds)
    # Extract PCA loadings
    loadings = pd.DataFrame(dds.varm['PCs'], index=dds.var_names)
    # Save top 500 genes by absolute loading for each of the first 3 PCs
    for i in range(3):
        pc_col = loadings.iloc[:, i]
        top_genes = pc_col.abs().sort_values(ascending=False).head(500)
        top_genes_df = pd.DataFrame({
            f'PC{i+1}_loading': pc_col.loc[top_genes.index]
        })
        top_genes_df.to_csv(f'top500_pc{i+1}_loadings.csv')
        print(f"Saved top 500 loadings for PC{i+1} to 'top500_pc{i+1}_loadings.csv'")
    
    # Print the percent variance explained by the first few PCs
    variance_ratio = dds.uns['pca']['variance_ratio']
    for i, var in enumerate(variance_ratio[:5]):  # Show first 5 PCs
        print(f"PC{i+1} explains {var * 100:.2f}% of the variance")
        
    # Save full loadings matrix too
    loadings.to_csv('pca_loadings.csv')
    
    # Set matplotlib parameters for better axis display
    plt.rcParams.update({'font.size': 12})
    
    # Get PCA coordinates
    pca_coords = dds.obsm['X_pca']
    conditions = dds.obs['Condition']
    
    # Print PCA coordinates for each sample
    print("\nPCA Coordinates for each sample:")
    print("=" * 50)
    for i, sample_id in enumerate(dds.obs_names):
        condition = conditions.iloc[i]
        pc1 = pca_coords[i, 0]
        pc2 = pca_coords[i, 1]
        pc3 = pca_coords[i, 2] if pca_coords.shape[1] > 2 else 0
        print(f"{sample_id:15} | {condition:10} | PC1: {pc1:7.2f} | PC2: {pc2:7.2f} | PC3: {pc3:7.2f}")
    print("=" * 50)
    
    # PC1 vs PC2 plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create scatter plot with condition colors
    for condition in conditions.unique():
        mask = conditions == condition
        ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1], 
                  label=condition, s=80, alpha=0.7)
    
    # Add variance explained to axis labels
    pc1_var = variance_ratio[0] * 100
    pc2_var = variance_ratio[1] * 100
    ax.set_xlabel(f'PC1 ({pc1_var:.1f}%)', fontsize=14)
    ax.set_ylabel(f'PC2 ({pc2_var:.1f}%)', fontsize=14)
    
    # Ensure axis ticks are visible and properly formatted
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Force axis ticks to be shown
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)
    
    plt.tight_layout()
    plt.savefig('pca_pc1_pc2.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('pca_pc1_pc2.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # PC2 vs PC3 plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create scatter plot for PC2 vs PC3
    for condition in conditions.unique():
        mask = conditions == condition
        ax.scatter(pca_coords[mask, 1], pca_coords[mask, 2], 
                  label=condition, s=80, alpha=0.7)
    
    # Add variance explained to axis labels
    pc3_var = variance_ratio[2] * 100
    ax.set_xlabel(f'PC2 ({pc2_var:.1f}%)', fontsize=14)
    ax.set_ylabel(f'PC3 ({pc3_var:.1f}%)', fontsize=14)
    
    # Ensure axis ticks are visible and properly formatted
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Force axis ticks to be shown
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)
    
    plt.tight_layout()
    plt.savefig('pca_pc2_pc3.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('pca_pc2_pc3.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("PCA plots saved as 'pca_pc1_pc2.pdf/png' and 'pca_pc2_pc3.pdf/png'.")
    print("PCA loadings saved to 'pca_loadings.csv'.")
    
def main():
    # --- Path to featureCounts file ---
    file_path = "{PATH}/all_samples.featureCounts"
    
    # --- Initialize and preprocess data ---
    processor = FeatureCountsProcessor(file_path)
    processor.load_data()
    processor.extract_counts_table()
    processor.filter_zero_rows()
    metadata = processor.create_metadata()
    
    # --- Fit DESeq2 model ONCE across all 3 conditions ---
    dds_all = processor.run_deseq2_all()
    
    # --- Run pairwise contrasts from the pre-fitted model ---
    res_1 = processor.run_pairwise_deseq2(('ACP7', 'Wildtype'))
    res_2 = processor.run_pairwise_deseq2(('ACP7', 'ACP7-7A'))
    res_3 = processor.run_pairwise_deseq2(('ACP7-7A', 'Wildtype'))
    
    # --- Plot volcano plots for each comparison ---
    plot_volcano_with_labels_and_save(res_1, save_path='oe_v_wt_volcano_plot.pdf')
    plot_volcano_with_labels_and_save(res_2, save_path='oe_v_dead_volcano_plot.pdf')
    plot_volcano_with_labels_and_save(res_3, save_path='dead_v_wt_volcano_plot.pdf')
    
    # --- Plot significant gene heatmaps for each comparison ---
    plot_sig_gene_heatmap(
        dds_all,
        res_1,
        condition_col='Condition',
        conditions=('ACP7', 'Wildtype'),
        save_path="oe_v_wt_sig_heatmap.pdf"
    )
    
    plot_sig_gene_heatmap(
        dds_all,
        res_2,
        condition_col='Condition',
        conditions=('ACP7', 'ACP7-7A'),
        save_path="oe_v_dead_sig_heatmap.pdf"
    )
    
    plot_sig_gene_heatmap(
        dds_all,
        res_3,
        condition_col='Condition',
        conditions=('ACP7-7A', 'Wildtype'),
        save_path="dead_v_wt_sig_heatmap.pdf"
    )
    
    # --- Save significant gene lists and overlaps ---
    save_significant_genes(res_1, output_csv='oe_v_wt_sig_genes.csv')
    save_significant_genes(res_2, output_csv='oe_v_dead_sig_genes.csv')
    save_significant_genes(res_3, output_csv='dead_v_wt_sig_genes.csv')
    
    # (optional) Combine overlaps
    save_significant_genes(res_1, res_2, output_csv='union_sig_genes.csv')
    
    # --- PCA plot using the full model (all samples) ---
    plot_pca_and_save(dds_all)

if __name__ == "__main__":
    main()