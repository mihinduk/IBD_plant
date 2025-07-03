#!/usr/bin/env python3
"""
Host Response Signatures Analysis for IBD Virome Study
======================================================

This script identifies viral-host interaction signatures by analyzing 
correlations between viral abundance and host inflammatory markers,
incorporating graph-based approaches for multi-viral signatures.

Author: Claude Code for Kathie Mihindukulasuriya
Date: 2025-06-12
"""

import pandas as pd
import numpy as np
import argparse
import json
import os
import sys
from pathlib import Path
import warnings
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.cluster import KMeans, DBSCAN
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import cross_val_score, StratifiedKFold
import networkx as nx
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for HTCF
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class HostResponseAnalyzer:
    """
    Main class for host response signature analysis
    """
    
    def __init__(self, output_dir, min_abundance=10, min_prevalence=0.1, 
                 correlation_threshold=0.3, fdr_alpha=0.05):
        """
        Initialize the analyzer
        
        Parameters:
        -----------
        output_dir : str
            Directory to save results
        min_abundance : int
            Minimum read count for inclusion
        min_prevalence : float
            Minimum prevalence across samples for inclusion
        correlation_threshold : float
            Minimum correlation strength for significance
        fdr_alpha : float
            FDR alpha level for multiple testing correction
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_abundance = min_abundance
        self.min_prevalence = min_prevalence
        self.correlation_threshold = correlation_threshold
        self.fdr_alpha = fdr_alpha
        
        # Initialize containers
        self.metadata = None
        self.counts = None
        self.taxonomy = None
        self.viral_families = None
        self.viral_species = None
        self.host_markers = {}
        self.correlation_results = {}
        self.viral_networks = {}
        self.response_signatures = {}
        
    def load_data(self, metadata_path, counts_path, taxonomy_path):
        """
        Load and preprocess all data files
        """
        print("Loading data files...")
        
        # Load metadata
        print(f"  Loading metadata from {metadata_path}")
        self.metadata = pd.read_csv(metadata_path, sep='\t')
        print(f"    Loaded {len(self.metadata)} samples")
        
        # Load taxonomy FIRST to identify viral contigs
        print(f"  Loading taxonomy from {taxonomy_path}")
        self.taxonomy = pd.read_csv(taxonomy_path, sep='\t')
        
        # Check for contig column name (might be 'contig_id' or 'contig' or 'Contig' or first column)
        if 'contig_id' not in self.taxonomy.columns:
            if 'Contig' in self.taxonomy.columns:
                self.taxonomy['contig_id'] = self.taxonomy['Contig']
            elif 'contig' in self.taxonomy.columns:
                self.taxonomy['contig_id'] = self.taxonomy['contig']
            else:
                # Assume first column is contig ID
                self.taxonomy['contig_id'] = self.taxonomy.iloc[:, 0]
                print(f"    Using first column '{self.taxonomy.columns[0]}' as contig_id")
        
        print(f"    Loaded taxonomy for {len(self.taxonomy)} contigs")
        
        # Pre-filter for viral contigs
        kingdom_col = 'superkingdom' if 'superkingdom' in self.taxonomy.columns else 'kingdom'
        viral_contigs = set(self.taxonomy[self.taxonomy[kingdom_col] == 'Viruses']['contig_id'].astype(str))
        print(f"    Identified {len(viral_contigs)} viral contigs")
        
        # Load counts - filter to viral contigs only
        print(f"  Loading counts from {counts_path}")
        
        # Check if this is a large file by reading just the header
        if counts_path.endswith('.gz'):
            test_df = pd.read_csv(counts_path, sep='\t', compression='gzip', nrows=5)
        else:
            test_df = pd.read_csv(counts_path, sep='\t', nrows=5)
        
        # Check if this is the long format (Sample, Contig, Reads, etc.)
        if 'Sample' in test_df.columns and 'Contig' in test_df.columns:
            print("    Detected long format - reading in chunks to filter viral contigs...")
            
            # Read in chunks and filter
            chunk_size = 10_000_000  # 10M rows at a time
            viral_counts = []
            total_rows = 0
            viral_rows = 0
            
            if counts_path.endswith('.gz'):
                reader = pd.read_csv(counts_path, sep='\t', compression='gzip', chunksize=chunk_size)
            else:
                reader = pd.read_csv(counts_path, sep='\t', chunksize=chunk_size)
            
            for chunk in reader:
                total_rows += len(chunk)
                # Filter to viral contigs
                viral_chunk = chunk[chunk['Contig'].astype(str).isin(viral_contigs)]
                viral_rows += len(viral_chunk)
                if len(viral_chunk) > 0:
                    viral_counts.append(viral_chunk)
                
                if total_rows % 50_000_000 == 0:
                    print(f"    Processed {total_rows:,} rows, found {viral_rows:,} viral records...")
            
            print(f"    Read {total_rows:,} total count records")
            print(f"    Filtered to {viral_rows:,} viral contig records")
            
            # Combine all viral chunks
            counts_df = pd.concat(viral_counts, ignore_index=True)
            del viral_counts  # Free memory
            
            print("    Converting long format counts to wide format...")
            # Pivot to create contig x sample matrix using Reads column
            self.counts = counts_df.pivot_table(
                index='Contig', 
                columns='Sample', 
                values='Reads',  # Use raw read counts
                fill_value=0
            )
            # Reset column name (remove the 'Sample' name from columns)
            self.counts.columns.name = None
            # Convert column names to strings to match metadata
            self.counts.columns = self.counts.columns.astype(str)
            print(f"    Created matrix: {self.counts.shape[0]} viral contigs x {self.counts.shape[1]} samples")
        else:
            # Assume it's already in wide format
            self.counts = counts_df.set_index(counts_df.columns[0])
            # Filter to viral contigs
            viral_in_counts = [c for c in viral_contigs if c in self.counts.index]
            self.counts = self.counts.loc[viral_in_counts]
            print(f"    Loaded {self.counts.shape[0]} viral contigs x {self.counts.shape[1]} samples")
        
        # Filter taxonomy to just viral contigs
        self.taxonomy = self.taxonomy[self.taxonomy['contig_id'].astype(str).isin(viral_contigs)]
        
        # Clean and align data
        self._clean_and_align_data()
        
        # Extract host response markers
        self._extract_host_markers()
        
    def _clean_and_align_data(self):
        """
        Clean data and align sample names across datasets
        """
        print("Cleaning and aligning data...")
        
        # Standardize sample names - use database_ID to match Sample column in counts
        metadata_samples = set(self.metadata['database_ID'].astype(str))
        count_samples = set(self.counts.columns.astype(str))
        
        print(f"  Metadata samples: {sorted(list(metadata_samples))}")
        print(f"  Count samples: {sorted(list(count_samples))}")
        
        # Find common samples
        common_samples = metadata_samples.intersection(count_samples)
        print(f"  Found {len(common_samples)} common samples: {sorted(list(common_samples))}")
        
        # Filter to common samples
        self.metadata = self.metadata[self.metadata['database_ID'].astype(str).isin(common_samples)]
        self.counts = self.counts[list(common_samples)]
        
        # Align contig IDs between counts and taxonomy
        count_contigs = set(self.counts.index.astype(str))
        tax_contigs = set(self.taxonomy['contig_id'].astype(str))
        common_contigs = count_contigs.intersection(tax_contigs)
        print(f"  Found {len(common_contigs)} contigs with both counts and taxonomy")
        
        # Filter to common contigs
        self.counts = self.counts.loc[list(common_contigs)]
        self.taxonomy = self.taxonomy[self.taxonomy['contig_id'].astype(str).isin(common_contigs)]
        
        # Sort to ensure alignment
        self.metadata = self.metadata.sort_values('database_ID')
        self.counts = self.counts.reindex(sorted(self.counts.columns))
        self.taxonomy = self.taxonomy.sort_values('contig_id')
        
        # Set index to database_ID for easier access
        self.metadata = self.metadata.set_index('database_ID')
        
    def _extract_host_markers(self):
        """
        Extract host inflammatory and response markers from metadata
        """
        print("Extracting host response markers...")
        
        # Define host response markers of interest
        marker_columns = {
            'inflammation': ['logFCP', 'CRP', 'ESR'],
            'intestinal_barrier': ['FOBT', 'calprotectin'],
            'immune_response': ['IL6', 'TNF', 'IFN'],
            'metabolic': ['albumin', 'hemoglobin'],
            'clinical': ['diagnosis', 'disease_activity', 'medications']
        }
        
        # Extract available markers
        available_columns = self.metadata.columns.tolist()
        
        for category, markers in marker_columns.items():
            available_markers = [m for m in markers if m in available_columns]
            if available_markers:
                self.host_markers[category] = available_markers
                print(f"  {category}: {available_markers}")
        
        # Create binary diagnosis variables
        if 'diagnosis' in available_columns:
            diagnosis_dummies = pd.get_dummies(self.metadata['diagnosis'], prefix='dx')
            self.metadata = pd.concat([self.metadata, diagnosis_dummies], axis=1)
            self.host_markers['diagnosis_binary'] = diagnosis_dummies.columns.tolist()
        
    def process_viral_data(self):
        """
        Process viral data at family and species levels
        """
        print("Processing viral data...")
        
        # Filter to viral contigs - check for kingdom/superkingdom column
        kingdom_col = 'superkingdom' if 'superkingdom' in self.taxonomy.columns else 'kingdom'
        
        viral_tax = self.taxonomy[
            (self.taxonomy[kingdom_col] == 'Viruses') &
            (self.taxonomy['family'].notna())
        ].copy()
        
        print(f"  Found {len(viral_tax)} viral contigs")
        
        # Process at family level
        family_counts = []
        for family in viral_tax['family'].unique():
            if pd.isna(family):
                continue
                
            family_contigs = viral_tax[viral_tax['family'] == family]['contig_id'].values
            family_contigs = [c for c in family_contigs if c in self.counts.index]
            
            if len(family_contigs) > 0:
                family_sum = self.counts.loc[family_contigs].sum(axis=0)
                family_counts.append({
                    'taxon': family,
                    'level': 'family',
                    'n_contigs': len(family_contigs),
                    'counts': family_sum
                })
        
        # Process at species level
        species_counts = []
        for species in viral_tax['species'].dropna().unique():
            if pd.isna(species):
                continue
                
            species_contigs = viral_tax[viral_tax['species'] == species]['contig_id'].values
            species_contigs = [c for c in species_contigs if c in self.counts.index]
            
            if len(species_contigs) > 0:
                species_sum = self.counts.loc[species_contigs].sum(axis=0)
                species_counts.append({
                    'taxon': species,
                    'level': 'species',
                    'n_contigs': len(species_contigs),
                    'counts': species_sum
                })
        
        # Create abundance matrices
        if family_counts:
            self.viral_families = pd.DataFrame({
                fc['taxon']: fc['counts'] for fc in family_counts
            }).fillna(0)
            # Set index to match sample IDs
            self.viral_families.index = self.counts.columns
            
            # Filter by abundance and prevalence
            abundance_filter = (self.viral_families >= self.min_abundance).sum(axis=0) >= len(self.viral_families) * self.min_prevalence
            self.viral_families = self.viral_families.loc[:, abundance_filter]
            print(f"  Retained {self.viral_families.shape[1]} viral families")
        
        if species_counts:
            self.viral_species = pd.DataFrame({
                sc['taxon']: sc['counts'] for sc in species_counts
            }).fillna(0)
            # Set index to match sample IDs
            self.viral_species.index = self.counts.columns
            
            # Filter by abundance and prevalence
            abundance_filter = (self.viral_species >= self.min_abundance).sum(axis=0) >= len(self.viral_species) * self.min_prevalence
            self.viral_species = self.viral_species.loc[:, abundance_filter]
            print(f"  Retained {self.viral_species.shape[1]} viral species")
        
    def compute_viral_host_correlations(self):
        """
        Compute correlations between viral abundance and host markers
        """
        print("Computing viral-host correlations...")
        
        # Merge metadata with viral data
        sample_data = self.metadata
        
        # Analyze correlations for each taxonomic level
        for level, viral_data in [('family', self.viral_families), ('species', self.viral_species)]:
            if viral_data is None or viral_data.empty:
                continue
                
            print(f"  Analyzing {level}-level correlations...")
            
            # Combine viral and host data
            combined_data = sample_data.join(viral_data.T, how='inner')
            
            level_results = []
            
            # Iterate through viral taxa and host markers
            for viral_taxon in viral_data.columns:
                viral_abundance = combined_data[viral_taxon]
                
                for marker_category, markers in self.host_markers.items():
                    for marker in markers:
                        if marker not in combined_data.columns:
                            continue
                        
                        host_values = combined_data[marker]
                        
                        # Remove missing values
                        valid_idx = ~(pd.isna(viral_abundance) | pd.isna(host_values))
                        if valid_idx.sum() < 10:  # Need at least 10 samples
                            continue
                        
                        viral_clean = viral_abundance[valid_idx]
                        host_clean = host_values[valid_idx]
                        
                        # Compute correlation
                        if marker_category == 'diagnosis_binary':
                            # Point-biserial correlation for binary diagnosis
                            corr, p_value = stats.pointbiserialr(host_clean, viral_clean)
                        else:
                            # Spearman correlation for continuous variables
                            corr, p_value = stats.spearmanr(viral_clean, host_clean)
                        
                        level_results.append({
                            'viral_taxon': viral_taxon,
                            'host_marker': marker,
                            'marker_category': marker_category,
                            'correlation': corr,
                            'p_value': p_value,
                            'n_samples': valid_idx.sum(),
                            'viral_prevalence': (viral_clean > 0).mean(),
                            'host_mean': host_clean.mean(),
                            'host_std': host_clean.std()
                        })
            
            if level_results:
                # Convert to DataFrame and apply multiple testing correction
                level_df = pd.DataFrame(level_results)
                
                # FDR correction
                _, level_df['p_adjusted'], _, _ = multipletests(
                    level_df['p_value'], alpha=self.fdr_alpha, method='fdr_bh'
                )
                
                # Filter for significant correlations
                significant = level_df[
                    (level_df['p_adjusted'] < self.fdr_alpha) &
                    (np.abs(level_df['correlation']) >= self.correlation_threshold)
                ].copy()
                
                self.correlation_results[level] = {
                    'all_results': level_df,
                    'significant': significant
                }
                
                print(f"    Found {len(significant)} significant correlations")
                
                # Save results
                level_df.to_csv(
                    self.output_dir / f'viral_host_correlations_{level}.csv',
                    index=False
                )
                
                significant.to_csv(
                    self.output_dir / f'significant_correlations_{level}.csv',
                    index=False
                )
    
    def build_viral_interaction_networks(self):
        """
        Build networks of viral-viral and viral-host interactions
        """
        print("Building viral interaction networks...")
        
        for level, viral_data in [('family', self.viral_families), ('species', self.viral_species)]:
            if viral_data is None or viral_data.empty:
                continue
                
            print(f"  Building {level}-level network...")
            
            # Create viral-viral correlation network
            viral_corr_matrix = viral_data.T.corr(method='spearman')
            
            # Create network graph
            G = nx.Graph()
            
            # Add viral taxa as nodes
            for taxon in viral_data.columns:
                G.add_node(taxon, node_type='viral', level=level)
            
            # Add significant viral-viral edges
            for i in range(len(viral_data.columns)):
                for j in range(i + 1, len(viral_data.columns)):
                    taxon1 = viral_data.columns[i]
                    taxon2 = viral_data.columns[j]
                    corr = viral_corr_matrix.loc[taxon1, taxon2]
                    
                    if not np.isnan(corr) and abs(corr) >= self.correlation_threshold:
                        G.add_edge(taxon1, taxon2, weight=corr, edge_type='viral_viral')
            
            # Add host markers as nodes and edges
            if level in self.correlation_results:
                significant_corrs = self.correlation_results[level]['significant']
                
                for _, row in significant_corrs.iterrows():
                    host_marker = row['host_marker']
                    viral_taxon = row['viral_taxon']
                    correlation = row['correlation']
                    
                    # Add host marker node if not already present
                    if host_marker not in G.nodes():
                        G.add_node(host_marker, 
                                 node_type='host', 
                                 category=row['marker_category'])
                    
                    # Add viral-host edge
                    G.add_edge(viral_taxon, host_marker, 
                             weight=correlation, 
                             edge_type='viral_host')
            
            self.viral_networks[level] = G
            
            # Save network
            nx.write_gml(G, self.output_dir / f'viral_host_network_{level}.gml')
            
            print(f"    Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    def identify_response_signatures(self):
        """
        Identify multi-viral response signatures using machine learning
        """
        print("Identifying host response signatures...")
        
        # Merge all data
        sample_data = self.metadata
        
        for level, viral_data in [('family', self.viral_families), ('species', self.viral_species)]:
            if viral_data is None or viral_data.empty:
                continue
                
            print(f"  Analyzing {level}-level signatures...")
            
            combined_data = sample_data.join(viral_data.T, how='inner')
            
            # Define target variables for prediction
            target_variables = []
            
            # Continuous inflammatory markers
            for category, markers in self.host_markers.items():
                if category in ['inflammation', 'intestinal_barrier']:
                    for marker in markers:
                        if marker in combined_data.columns:
                            target_variables.append((marker, 'continuous'))
            
            # Binary diagnosis outcomes
            if 'diagnosis_binary' in self.host_markers:
                for dx_var in self.host_markers['diagnosis_binary']:
                    if dx_var in combined_data.columns:
                        target_variables.append((dx_var, 'binary'))
            
            level_signatures = {}
            
            for target_var, var_type in target_variables:
                print(f"    Predicting {target_var}...")
                
                # Prepare data
                X = combined_data[viral_data.columns].fillna(0)  # Viral features
                y = combined_data[target_var].dropna()
                
                # Align X and y
                common_idx = X.index.intersection(y.index)
                X = X.loc[common_idx]
                y = y.loc[common_idx]
                
                if len(y) < 20:  # Need minimum samples
                    continue
                
                # Scale features
                scaler = RobustScaler()
                X_scaled = pd.DataFrame(
                    scaler.fit_transform(X),
                    columns=X.columns,
                    index=X.index
                )
                
                # Train model
                if var_type == 'continuous':
                    model = RandomForestRegressor(n_estimators=100, random_state=42)
                    # Cross-validation for regression
                    cv_scores = cross_val_score(model, X_scaled, y, cv=5, scoring='r2')
                else:
                    model = RandomForestClassifier(n_estimators=100, random_state=42)
                    # Cross-validation for classification
                    if len(np.unique(y)) > 1:  # Need both classes
                        cv_scores = cross_val_score(model, X_scaled, y, cv=5, scoring='roc_auc')
                    else:
                        continue
                
                # Fit final model
                model.fit(X_scaled, y)
                
                # Get feature importance
                feature_importance = pd.DataFrame({
                    'viral_taxon': X.columns,
                    'importance': model.feature_importances_
                }).sort_values('importance', ascending=False)
                
                # Identify top features (signature)
                top_features = feature_importance.head(10)
                
                level_signatures[target_var] = {
                    'model_type': var_type,
                    'cv_score_mean': cv_scores.mean(),
                    'cv_score_std': cv_scores.std(),
                    'feature_importance': feature_importance,
                    'top_features': top_features,
                    'n_samples': len(y)
                }
                
                # Save feature importance
                feature_importance.to_csv(
                    self.output_dir / f'feature_importance_{level}_{target_var.replace("/", "_")}.csv',
                    index=False
                )
            
            self.response_signatures[level] = level_signatures
    
    def create_visualizations(self):
        """
        Create comprehensive visualizations
        """
        print("Creating visualizations...")
        
        # 1. Correlation heatmaps
        for level, results in self.correlation_results.items():
            if 'significant' not in results or results['significant'].empty:
                continue
                
            # Pivot significant correlations for heatmap
            sig_corrs = results['significant']
            heatmap_data = sig_corrs.pivot(
                index='viral_taxon', 
                columns='host_marker', 
                values='correlation'
            ).fillna(0)
            
            if heatmap_data.shape[0] > 1 and heatmap_data.shape[1] > 1:
                plt.figure(figsize=(12, 8))
                sns.heatmap(heatmap_data, 
                           center=0, cmap='RdBu_r',
                           annot=True, fmt='.2f',
                           cbar_kws={'label': 'Correlation'})
                plt.title(f'Significant Viral-Host Correlations ({level.title()} Level)')
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig(
                    self.output_dir / f'correlation_heatmap_{level}.png',
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
        
        # 2. Network visualizations
        for level, G in self.viral_networks.items():
            if G.number_of_nodes() == 0:
                continue
                
            plt.figure(figsize=(14, 10))
            
            # Layout
            pos = nx.spring_layout(G, k=1, iterations=50)
            
            # Separate nodes by type
            viral_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'viral']
            host_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'host']
            
            # Draw nodes
            nx.draw_networkx_nodes(G, pos, nodelist=viral_nodes, 
                                 node_color='lightblue', node_size=300, alpha=0.8)
            nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, 
                                 node_color='lightcoral', node_size=300, alpha=0.8)
            
            # Draw edges
            viral_viral_edges = [(u, v) for u, v, d in G.edges(data=True) 
                               if d.get('edge_type') == 'viral_viral']
            viral_host_edges = [(u, v) for u, v, d in G.edges(data=True) 
                              if d.get('edge_type') == 'viral_host']
            
            nx.draw_networkx_edges(G, pos, edgelist=viral_viral_edges, 
                                 edge_color='gray', alpha=0.5)
            nx.draw_networkx_edges(G, pos, edgelist=viral_host_edges, 
                                 edge_color='red', alpha=0.7)
            
            # Labels for small networks
            if G.number_of_nodes() <= 20:
                nx.draw_networkx_labels(G, pos, font_size=8)
            
            plt.title(f'Viral-Host Interaction Network ({level.title()} Level)')
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(
                self.output_dir / f'interaction_network_{level}.png',
                dpi=300, bbox_inches='tight'
            )
            plt.close()
        
        # 3. Feature importance plots
        for level, signatures in self.response_signatures.items():
            for target_var, sig_data in signatures.items():
                top_features = sig_data['top_features']
                
                plt.figure(figsize=(10, 6))
                plt.barh(range(len(top_features)), top_features['importance'])
                plt.yticks(range(len(top_features)), top_features['viral_taxon'])
                plt.xlabel('Feature Importance')
                plt.title(f'Top Viral Predictors for {target_var} ({level.title()} Level)')
                plt.gca().invert_yaxis()
                plt.tight_layout()
                plt.savefig(
                    self.output_dir / f'feature_importance_{level}_{target_var.replace("/", "_")}.png',
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
        
        print("  Visualization creation complete")
    
    def save_summary_report(self):
        """
        Save comprehensive summary report
        """
        print("Generating summary report...")
        
        summary = {
            'analysis_date': pd.Timestamp.now().isoformat(),
            'data_summary': {
                'n_samples': len(self.metadata),
                'n_contigs': self.counts.shape[0],
                'n_viral_families': self.viral_families.shape[1] if self.viral_families is not None else 0,
                'n_viral_species': self.viral_species.shape[1] if self.viral_species is not None else 0,
                'host_marker_categories': list(self.host_markers.keys())
            },
            'correlation_analysis': {
                level: {
                    'total_correlations': len(results['all_results']),
                    'significant_correlations': len(results['significant'])
                }
                for level, results in self.correlation_results.items()
            },
            'network_analysis': {
                level: {
                    'nodes': G.number_of_nodes(),
                    'edges': G.number_of_edges(),
                    'viral_nodes': len([n for n in G.nodes() if G.nodes[n].get('node_type') == 'viral']),
                    'host_nodes': len([n for n in G.nodes() if G.nodes[n].get('node_type') == 'host'])
                }
                for level, G in self.viral_networks.items()
            },
            'signature_analysis': {
                level: {
                    'n_signatures': len(signatures),
                    'targets': list(signatures.keys())
                }
                for level, signatures in self.response_signatures.items()
            },
            'parameters': {
                'min_abundance': self.min_abundance,
                'min_prevalence': self.min_prevalence,
                'correlation_threshold': self.correlation_threshold,
                'fdr_alpha': self.fdr_alpha
            }
        }
        
        with open(self.output_dir / 'host_response_analysis_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"  Summary report saved to {self.output_dir / 'host_response_analysis_summary.json'}")


def main():
    """
    Main execution function
    """
    parser = argparse.ArgumentParser(
        description='Host Response Signatures Analysis for IBD Virome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--metadata', required=True,
                        help='Path to metadata file')
    parser.add_argument('--counts', required=True,
                        help='Path to contig count table (TSV, optionally gzipped)')
    parser.add_argument('--taxonomy', required=True,
                        help='Path to taxonomy file')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for results')
    parser.add_argument('--min-abundance', type=int, default=10,
                        help='Minimum read count for inclusion')
    parser.add_argument('--min-prevalence', type=float, default=0.1,
                        help='Minimum prevalence across samples')
    parser.add_argument('--correlation-threshold', type=float, default=0.3,
                        help='Minimum correlation strength for significance')
    parser.add_argument('--fdr-alpha', type=float, default=0.05,
                        help='FDR alpha level for multiple testing correction')
    parser.add_argument('--skip-visualizations', action='store_true',
                        help='Skip visualization generation (faster)')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = HostResponseAnalyzer(
        output_dir=args.output_dir,
        min_abundance=args.min_abundance,
        min_prevalence=args.min_prevalence,
        correlation_threshold=args.correlation_threshold,
        fdr_alpha=args.fdr_alpha
    )
    
    try:
        # Run analysis pipeline
        analyzer.load_data(args.metadata, args.counts, args.taxonomy)
        analyzer.process_viral_data()
        analyzer.compute_viral_host_correlations()
        analyzer.build_viral_interaction_networks()
        analyzer.identify_response_signatures()
        
        if not args.skip_visualizations:
            analyzer.create_visualizations()
        
        analyzer.save_summary_report()
        
        print(f"\nAnalysis complete! Results saved to {args.output_dir}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()