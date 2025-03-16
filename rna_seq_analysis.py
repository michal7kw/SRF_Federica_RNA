#!/usr/bin/env python3

# %% Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc
from statsmodels.stats.multitest import fdrcorrection
import os
from pathlib import Path

# %% Set parameters
counts_file = 'results/counts/all_samples_counts.txt'
output_dir = Path('results/analysis')
output_dir.mkdir(parents=True, exist_ok=True)

# %% Load and process data
def load_data(counts_file):
    """Load and process the count data"""
    # Read the counts file, skipping the first line (featureCounts info)
    counts = pd.read_csv(counts_file, sep='\t', skiprows=1)
    
    # Set gene names as index
    counts.set_index('Geneid', inplace=True)
    
    # Remove unnecessary columns from featureCounts
    cols_to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']
    counts = counts.drop(columns=[col for col in cols_to_drop if col in counts.columns])
    
    # Clean up sample names
    counts.columns = [os.path.basename(col).split('_Aligned')[0] for col in counts.columns]
    
    return counts

counts = load_data(counts_file)

# Calculate basic statistics
total_counts = counts.sum()
genes_detected = (counts > 0).sum()

# %% Normalize counts
def normalize_counts(counts, method='cpm'):
    """
    Normalize count data
    
    Parameters:
    -----------
    method : str
        Normalization method ('cpm' or 'tmm')
    """
    if method == 'cpm':
        # CPM normalization
        norm_counts = counts.copy()
        norm_counts = (norm_counts * 1e6) / norm_counts.sum()
        return norm_counts
    else:
        raise ValueError(f"Normalization method '{method}' not implemented")

norm_counts = normalize_counts(counts)

# %% Filter low expression
def filter_low_expression(counts, norm_counts=None, min_counts=10, min_samples=2):
    """
    Filter out lowly expressed genes
    
    Parameters:
    -----------
    min_counts : int
        Minimum counts required
    min_samples : int
        Minimum number of samples that must meet min_counts
    """
    # Keep genes that have at least min_counts in at least min_samples samples
    mask = (counts >= min_counts).sum(axis=1) >= min_samples
    filtered_counts = counts.loc[mask]
    
    if norm_counts is not None:
        filtered_norm_counts = norm_counts.loc[mask]
        return filtered_counts, filtered_norm_counts
    return filtered_counts

counts, norm_counts = filter_low_expression(counts, norm_counts)

# %% Plot sample statistics
def plot_sample_stats(total_counts, genes_detected, output_dir):
    """Generate basic sample statistics plots"""
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Total counts per sample
    sns.barplot(x=total_counts.index, y=total_counts.values, ax=ax1)
    ax1.set_title('Total Counts per Sample')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
    ax1.set_ylabel('Total Counts')
    
    # Genes detected per sample
    sns.barplot(x=genes_detected.index, y=genes_detected.values, ax=ax2)
    ax2.set_title('Genes Detected per Sample')
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
    ax2.set_ylabel('Number of Genes')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sample_stats.png')
    plt.close()

plot_sample_stats(total_counts, genes_detected, output_dir)

# %% Plot expression distribution
def plot_expression_dist(counts, output_dir):
    """Plot expression distribution"""
    plt.figure(figsize=(10, 6))
    
    # Plot density of log2 transformed counts
    for col in counts.columns:
        sns.kdeplot(data=np.log2(counts[col] + 1), label=col)
        
    plt.title('Expression Distribution')
    plt.xlabel('log2(counts + 1)')
    plt.ylabel('Density')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / 'expression_dist.png')
    plt.close()

plot_expression_dist(counts, output_dir)

# %% Plot PCA
def plot_pca(data, output_dir):
    """Perform and plot PCA"""
    # Log transform
    data = np.log2(data + 1)
    
    # Perform PCA
    pca = sc.tl.pca(data.T)
    variance_ratio = pca.explained_variance_ratio_
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(pca.X[:, 0], pca.X[:, 1])
    
    # Add sample labels
    for i, sample in enumerate(data.columns):
        plt.annotate(sample, (pca.X[i, 0], pca.X[i, 1]))
        
    plt.xlabel(f'PC1 ({variance_ratio[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({variance_ratio[1]:.1%} variance)')
    plt.title('PCA Plot')
    plt.tight_layout()
    plt.savefig(output_dir / 'pca_plot.png')
    plt.close()

plot_pca(norm_counts, output_dir)

# %% Differential expression analysis
def run_differential_expression(norm_counts, group1, group2, output_dir, method='t-test'):
    """
    Perform differential expression analysis
    
    Parameters:
    -----------
    group1 : list
        List of sample names for group 1
    group2 : list
        List of sample names for group 2
    method : str
        Statistical method to use
    """
    results = []
    for gene in norm_counts.index:
        g1_expr = norm_counts.loc[gene, group1]
        g2_expr = norm_counts.loc[gene, group2]
        
        if method == 't-test':
            stat, pval = stats.ttest_ind(g1_expr, g2_expr)
            fc = np.log2(g1_expr.mean() / g2_expr.mean())
            results.append({
                'gene': gene,
                'log2FC': fc,
                'pvalue': pval,
                'stat': stat
            })
            
    # Create results DataFrame
    de_results = pd.DataFrame(results)
    de_results.set_index('gene', inplace=True)
    
    # Multiple testing correction
    de_results['padj'] = fdrcorrection(de_results['pvalue'])[1]
    
    # Sort by adjusted p-value
    de_results.sort_values('padj', inplace=True)
    
    # Save results
    de_results.to_csv(output_dir / 'differential_expression_results.csv')
    
    return de_results

# Example groups - replace with actual sample groups
group1 = ['1', '2', '3']  # control samples
group2 = ['4', '5', '6']  # treatment samples

de_results = run_differential_expression(norm_counts, group1, group2, output_dir)

# %% Plot volcano plot
def plot_volcano(de_results, output_dir, fc_threshold=1, pval_threshold=0.05):
    """
    Create volcano plot from differential expression results
    
    Parameters:
    -----------
    fc_threshold : float
        Log2 fold change threshold
    pval_threshold : float
        Adjusted p-value threshold
    """
    plt.figure(figsize=(10, 8))
    
    # Plot all points
    plt.scatter(
        de_results['log2FC'],
        -np.log10(de_results['padj']),
        alpha=0.5
    )
    
    # Highlight significant points
    significant = (abs(de_results['log2FC']) > fc_threshold) & \
                 (de_results['padj'] < pval_threshold)
    
    plt.scatter(
        de_results.loc[significant, 'log2FC'],
        -np.log10(de_results.loc[significant, 'padj']),
        color='red',
        alpha=0.5
    )
    
    plt.axhline(-np.log10(pval_threshold), color='gray', linestyle='--')
    plt.axvline(-fc_threshold, color='gray', linestyle='--')
    plt.axvline(fc_threshold, color='gray', linestyle='--')
    
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 Adjusted p-value')
    plt.title('Volcano Plot')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'volcano_plot.png')
    plt.close()

plot_volcano(de_results, output_dir)
