#!/usr/bin/env python3
"""
Contig-Level Temporal-Ecological Fingerprinting Analysis
========================================================

This script extends the family-level analysis to individual contigs,
identifying contigs that share similar temporal patterns regardless 
of taxonomic classification.

Author: Claude Code for Kathie Mihindukulasuriya
Date: 2025-06-18
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
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import networkx as nx

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for HTCF
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class ContigTemporalAnalyzer:
    """
    Contig-level temporal pattern analysis
    """
    
    def __init__(self, output_dir, min_samples=2, min_total_reads=10):
        """
        Initialize the analyzer
        
        Parameters:
        -----------
        output_dir : str
            Directory to save results
        min_samples : int
            Minimum number of samples for inclusion (default: 2)
        min_total_reads : int
            Minimum total reads across all samples for a contig (default: 50)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_samples = min_samples
        self.min_total_reads = min_total_reads
        
        # Initialize containers
        self.metadata = None
        self.counts = None
        self.taxonomy = None
        self.contig_patterns = None
        self.pattern_clusters = None
        
    def load_data(self, metadata_path, counts_path, taxonomy_path):
        """
        Load and preprocess all data files
        """
        print("Loading data files...")
        
        # Load metadata
        print(f"  Loading metadata from {metadata_path}")
        self.metadata = pd.read_csv(metadata_path, sep='\t')
        print(f"    Loaded {len(self.metadata)} samples")
        
        # Load counts
        print(f"  Loading counts from {counts_path}")
        if counts_path.endswith('.gz'):
            counts_df = pd.read_csv(counts_path, sep='\t', compression='gzip')
        else:
            counts_df = pd.read_csv(counts_path, sep='\t')
        
        # Check if this is the long format (Sample, Contig, Reads, etc.)
        if 'Sample' in counts_df.columns and 'Contig' in counts_df.columns:
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
        else:
            # Assume it's already in wide format
            self.counts = counts_df.set_index(counts_df.columns[0])
        
        print(f"    Loaded {self.counts.shape[0]} contigs x {self.counts.shape[1]} samples")
        
        # Ensure columns are strings for consistent matching
        self.counts.columns = self.counts.columns.astype(str)
        
        # Load taxonomy
        print(f"  Loading taxonomy from {taxonomy_path}")
        self.taxonomy = pd.read_csv(taxonomy_path, sep='\t')
        print(f"    Loaded taxonomy for {len(self.taxonomy)} contigs")
        
        # Standardize taxonomy contig ID column name
        if 'Contig' in self.taxonomy.columns:
            self.taxonomy = self.taxonomy.rename(columns={'Contig': 'contig_id'})
        
    def align_data(self):
        """
        Align metadata, counts, and taxonomy data
        """
        print("Cleaning and aligning data...")
        
        # Ensure database_ID is string for consistent matching
        self.metadata['database_ID'] = self.metadata['database_ID'].astype(str)
        
        # Find common samples
        metadata_samples = set(self.metadata['database_ID'].astype(str))
        count_samples = set(self.counts.columns.astype(str))
        common_samples = metadata_samples.intersection(count_samples)
        
        print(f"  Found {len(common_samples)} common samples: {sorted(common_samples)}")
        
        if len(common_samples) == 0:
            raise ValueError("No matching samples found between metadata and counts!")
        
        # Filter to common samples
        self.metadata = self.metadata[self.metadata['database_ID'].astype(str).isin(common_samples)]
        self.counts = self.counts[sorted(common_samples)]
        
        # Align taxonomy with counts (find contigs with both counts and taxonomy)
        count_contigs = set(self.counts.index)
        tax_contigs = set(self.taxonomy['contig_id'])
        common_contigs = count_contigs.intersection(tax_contigs)
        
        print(f"  Found {len(common_contigs)} contigs with both counts and taxonomy")
        
        # Filter to common contigs
        self.counts = self.counts.loc[list(common_contigs)]
        self.taxonomy = self.taxonomy[self.taxonomy['contig_id'].isin(common_contigs)]
        
    def create_abundance_pattern(self, values):
        """
        Create abundance pattern string from numeric values
        """
        if len(values) == 0 or all(v == 0 for v in values):
            return 'O' * len(values)
        
        # Calculate thresholds based on non-zero values
        non_zero_values = [v for v in values if v > 0]
        if len(non_zero_values) == 0:
            return 'O' * len(values)
        
        # Use quartiles of non-zero values for thresholds
        q25 = np.percentile(non_zero_values, 25)
        q75 = np.percentile(non_zero_values, 75)
        
        pattern = []
        for val in values:
            if val == 0:
                pattern.append('O')  # Absent
            elif val <= q25:
                pattern.append('L')  # Low
            elif val <= q75:
                pattern.append('M')  # Medium
            else:
                pattern.append('H')  # High
        
        return '-'.join(pattern)
    
    def classify_pattern_type(self, presence_pattern, abundance_pattern):
        """
        Classify the type of temporal pattern
        """
        n_present = presence_pattern.count('X')
        n_total = len(presence_pattern)
        
        if n_present == n_total:
            return 'stable_present'
        elif n_present == 0:
            return 'stable_absent'
        elif presence_pattern.startswith('O') and presence_pattern.endswith('X'):
            return 'emerging'
        elif presence_pattern.startswith('X') and presence_pattern.endswith('O'):
            return 'declining'
        else:
            return 'episodic'
    
    def analyze_contig_patterns(self):
        """
        Analyze temporal patterns for individual contigs
        """
        print("Analyzing contig-level temporal patterns...")
        
        # Filter contigs by minimum total reads
        total_reads = self.counts.sum(axis=1)
        filtered_contigs = total_reads[total_reads >= self.min_total_reads].index
        
        print(f"  Analyzing {len(filtered_contigs)} contigs with ≥{self.min_total_reads} total reads")
        
        pattern_results = []
        
        for contig_id in filtered_contigs:
            contig_counts = self.counts.loc[contig_id]
            
            # Calculate pattern metrics
            n_samples_present = (contig_counts > 0).sum()
            
            # Skip if not present in minimum samples
            if n_samples_present < self.min_samples:
                continue
            
            total_samples = len(contig_counts)
            prevalence_pct = (n_samples_present / total_samples) * 100
            
            # Create presence pattern (X = present, O = absent)
            presence_pattern = ''.join(['X' if x > 0 else 'O' for x in contig_counts])
            
            # Create abundance pattern (H = high, M = medium, L = low, O = absent)
            abundance_pattern = self.create_abundance_pattern(contig_counts.values)
            
            # Classify pattern type
            pattern_type = self.classify_pattern_type(presence_pattern, abundance_pattern)
            
            # Calculate summary statistics
            total_reads = contig_counts.sum()
            mean_abundance = contig_counts.mean()
            max_abundance = contig_counts.max()
            
            # Get taxonomy info if available
            contig_tax = self.taxonomy[self.taxonomy['contig_id'] == contig_id]
            if len(contig_tax) > 0:
                family = contig_tax.iloc[0].get('family', 'Unknown')
                kingdom = contig_tax.iloc[0].get('superkingdom', 
                         contig_tax.iloc[0].get('kingdom', 'Unknown'))
                species = contig_tax.iloc[0].get('species', 'Unknown')
            else:
                family = 'Unknown'
                kingdom = 'Unknown'
                species = 'Unknown'
            
            # Use pipe delimiter for contig IDs to avoid comma conflicts
            contig_id_safe = str(contig_id).replace(',', '|')
            
            pattern_results.append({
                'contig_id': contig_id_safe,
                'n_samples_present': n_samples_present,
                'total_samples': total_samples,
                'prevalence_pct': round(prevalence_pct, 1),
                'presence_pattern': presence_pattern,
                'abundance_pattern': abundance_pattern,
                'pattern_type': pattern_type,
                'total_reads': total_reads,
                'mean_abundance': round(mean_abundance, 2),
                'max_abundance': max_abundance,
                'taxonomy_kingdom': kingdom,
                'taxonomy_family': family,
                'taxonomy_species': species
            })
        
        # Convert to DataFrame
        self.contig_patterns = pd.DataFrame(pattern_results)
        
        print(f"  Analyzed {len(self.contig_patterns)} contigs meeting criteria")
        print(f"  Pattern types found: {self.contig_patterns['pattern_type'].value_counts().to_dict()}")
        
        return self.contig_patterns
    
    def cluster_similar_patterns(self):
        """
        Cluster contigs with similar temporal patterns
        """
        print("Clustering contigs by similar patterns...")
        
        if self.contig_patterns is None or len(self.contig_patterns) == 0:
            print("  No patterns to cluster")
            return
        
        # Initialize cluster_id column for all contigs
        self.contig_patterns['cluster_id'] = 0  # 0 = no cluster (singleton)
        
        # Group by identical patterns
        pattern_groups = self.contig_patterns.groupby(['presence_pattern', 'abundance_pattern'])
        
        cluster_results = []
        cluster_id = 1
        
        for (presence_pattern, abundance_pattern), group in pattern_groups:
            if len(group) >= 2:  # Only clusters with multiple contigs
                # Assign cluster_id to contigs in this cluster
                contig_indices = group.index
                self.contig_patterns.loc[contig_indices, 'cluster_id'] = cluster_id
                
                cluster_results.append({
                    'cluster_id': cluster_id,
                    'presence_pattern': presence_pattern,
                    'abundance_pattern': abundance_pattern,
                    'pattern_type': group.iloc[0]['pattern_type'],
                    'n_contigs': len(group),
                    'contig_ids': '|'.join(group['contig_id'].tolist()),
                    'taxonomic_families': '|'.join(group['taxonomy_family'].unique()),
                    'taxonomic_kingdoms': '|'.join(group['taxonomy_kingdom'].unique()),
                    'total_reads_sum': group['total_reads'].sum(),
                    'mean_prevalence': round(group['prevalence_pct'].mean(), 1)
                })
                cluster_id += 1
        
        self.pattern_clusters = pd.DataFrame(cluster_results)
        
        print(f"  Found {len(self.pattern_clusters)} clusters with ≥2 contigs")
        print(f"  Assigned cluster IDs to {(self.contig_patterns['cluster_id'] > 0).sum()} contigs")
        
        return self.pattern_clusters
    
    def create_visualizations(self):
        """
        Create visualization of contig patterns
        """
        print("Creating visualizations...")
        
        if self.contig_patterns is None or len(self.contig_patterns) == 0:
            print("  No data to visualize")
            return
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Contig-Level Temporal Pattern Analysis', fontsize=16, fontweight='bold')
        
        # Plot 1: Pattern type distribution
        pattern_counts = self.contig_patterns['pattern_type'].value_counts()
        axes[0, 0].bar(pattern_counts.index, pattern_counts.values)
        axes[0, 0].set_title('Distribution of Pattern Types')
        axes[0, 0].set_xlabel('Pattern Type')
        axes[0, 0].set_ylabel('Number of Contigs')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # Plot 2: Prevalence distribution
        axes[0, 1].hist(self.contig_patterns['prevalence_pct'], bins=20, alpha=0.7)
        axes[0, 1].set_title('Distribution of Contig Prevalence')
        axes[0, 1].set_xlabel('Prevalence (%)')
        axes[0, 1].set_ylabel('Number of Contigs')
        
        # Plot 3: Total reads vs prevalence
        axes[1, 0].scatter(self.contig_patterns['total_reads'], 
                          self.contig_patterns['prevalence_pct'], 
                          alpha=0.6)
        axes[1, 0].set_xlabel('Total Reads')
        axes[1, 0].set_ylabel('Prevalence (%)')
        axes[1, 0].set_title('Abundance vs Prevalence')
        axes[1, 0].set_xscale('log')
        
        # Plot 4: Kingdom distribution
        if 'taxonomy_kingdom' in self.contig_patterns.columns:
            kingdom_counts = self.contig_patterns['taxonomy_kingdom'].value_counts().head(10)
            axes[1, 1].barh(kingdom_counts.index, kingdom_counts.values)
            axes[1, 1].set_title('Top Taxonomic Kingdoms')
            axes[1, 1].set_xlabel('Number of Contigs')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'contig_pattern_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("  Visualization saved to contig_pattern_analysis.png")
    
    def _create_family_taxonomy_mapping_from_contigs(self):
        """
        Create family-taxonomy mapping from contig-level analysis results
        """
        print("  Creating family-taxonomy mapping from contig results...")
        
        if self.contig_patterns is None or len(self.contig_patterns) == 0:
            return
        
        # Group contigs by family
        family_groups = self.contig_patterns.groupby('taxonomy_family')
        
        mapping_results = []
        
        for family, group in family_groups:
            if pd.isna(family) or family == 'Unknown':
                continue
            
            # Get all contigs for this family
            family_contigs = group['contig_id'].tolist()
            
            # Use the first contig to get representative taxonomy
            representative = group.iloc[0]
            
            # Build full taxonomy string
            taxonomy_parts = []
            if representative['taxonomy_kingdom'] not in ['Unknown', '']:
                taxonomy_parts.append(f"kingdom:{representative['taxonomy_kingdom']}")
            if family not in ['Unknown', '']:
                taxonomy_parts.append(f"family:{family}")
            if representative['taxonomy_species'] not in ['Unknown', '']:
                taxonomy_parts.append(f"species:{representative['taxonomy_species']}")
            
            full_taxonomy = "; ".join(taxonomy_parts) if taxonomy_parts else "Unknown"
            
            # Join contigs with pipe delimiter
            contigs_str = "|".join(sorted(family_contigs))
            
            mapping_results.append({
                'viral_family': family,
                'taxonomy': full_taxonomy,
                'contigs': contigs_str
            })
        
        # Convert to DataFrame and save
        if mapping_results:
            mapping_df = pd.DataFrame(mapping_results)
            mapping_df = mapping_df.sort_values('viral_family')
            
            mapping_path = self.output_dir / 'viral_family_taxonomy_mapping.csv'
            mapping_df.to_csv(mapping_path, index=False)
            
            print(f"    Saved family-taxonomy mapping to {mapping_path}")
            print(f"    Mapped {len(mapping_df)} families")
        else:
            print("    No valid families found for mapping")
    
    def save_results(self):
        """
        Save all analysis results
        """
        print("Saving results...")
        
        # Save contig patterns
        if self.contig_patterns is not None:
            contig_output_path = self.output_dir / 'contig_temporal_patterns.csv'
            self.contig_patterns.to_csv(contig_output_path, index=False)
            print(f"  Contig patterns saved to {contig_output_path}")
            
            # Create family-taxonomy mapping from contig results
            self._create_family_taxonomy_mapping_from_contigs()
        else:
            # Create empty taxonomy mapping file
            empty_tax_df = pd.DataFrame(columns=[
                'viral_family', 'taxonomy', 'contigs'
            ])
            empty_tax_df.to_csv(self.output_dir / 'viral_family_taxonomy_mapping.csv', index=False)
        
        # Save pattern clusters
        if self.pattern_clusters is not None:
            cluster_output_path = self.output_dir / 'contig_pattern_clusters.csv'
            self.pattern_clusters.to_csv(cluster_output_path, index=False)
            print(f"  Pattern clusters saved to {cluster_output_path}")
        
        # Create summary statistics
        summary_stats = {
            'analysis_type': 'contig_level_temporal_patterns',
            'parameters': {
                'min_samples': self.min_samples,
                'min_total_reads': self.min_total_reads
            },
            'data_summary': {
                'n_samples': len(self.metadata) if self.metadata is not None else 0,
                'n_contigs_analyzed': len(self.contig_patterns) if self.contig_patterns is not None else 0,
                'n_pattern_clusters': len(self.pattern_clusters) if self.pattern_clusters is not None else 0
            }
        }
        
        if self.contig_patterns is not None:
            summary_stats['pattern_type_counts'] = self.contig_patterns['pattern_type'].value_counts().to_dict()
            summary_stats['kingdom_distribution'] = self.contig_patterns['taxonomy_kingdom'].value_counts().to_dict()
        
        # Save summary
        summary_path = self.output_dir / 'contig_analysis_summary.json'
        with open(summary_path, 'w') as f:
            json.dump(summary_stats, f, indent=2)
        
        print(f"  Summary saved to {summary_path}")
        
        # Create data dictionary
        self._create_data_dictionary()
    
    def _create_data_dictionary(self):
        """
        Create comprehensive data dictionary for all output files
        """
        print("  Creating data dictionary...")
        
        data_dict = {
            "Data Dictionary": "IBD Virome Contig-Level Temporal Pattern Analysis",
            "Analysis_Date": "Generated automatically during analysis",
            "Files_Described": [
                "contig_temporal_patterns.csv",
                "contig_pattern_clusters.csv", 
                "viral_family_taxonomy_mapping.csv",
                "contig_analysis_summary.json"
            ],
            
            "File_Descriptions": {
                "contig_temporal_patterns.csv": {
                    "description": "Temporal patterns for individual contigs across all samples",
                    "columns": {
                        "contig_id": "Unique contig identifier (pipe-delimited if needed to avoid CSV conflicts)",
                        "n_samples_present": "Number of samples where this contig has >0 reads",
                        "total_samples": "Total number of samples in the analysis",
                        "prevalence_pct": "Percentage of samples where contig is present (n_samples_present/total_samples * 100)",
                        "presence_pattern": "Binary presence pattern across samples (X=present, O=absent). Order matches sample chronology.",
                        "abundance_pattern": "Abundance level pattern across samples (H=high, M=medium, L=low, O=absent). Thresholds based on quartiles of non-zero values.",
                        "pattern_type": "Classification of temporal pattern (see Pattern_Types below)",
                        "total_reads": "Sum of all reads for this contig across all samples",
                        "mean_abundance": "Average read count per sample (including zeros)",
                        "max_abundance": "Maximum read count observed in any single sample",
                        "taxonomy_kingdom": "Taxonomic kingdom assignment (e.g., Viruses, Bacteria)",
                        "taxonomy_family": "Taxonomic family assignment",
                        "taxonomy_species": "Taxonomic species assignment",
                        "cluster_id": "Cluster identifier (>0 if contig belongs to a cluster with ≥2 members, 0 if singleton)"
                    }
                },
                
                "contig_pattern_clusters.csv": {
                    "description": "Groups of contigs that share identical temporal patterns",
                    "columns": {
                        "cluster_id": "Unique identifier for each pattern cluster",
                        "presence_pattern": "Shared binary presence pattern (X=present, O=absent)",
                        "abundance_pattern": "Shared abundance pattern (H=high, M=medium, L=low, O=absent)",
                        "pattern_type": "Classification of temporal pattern (see Pattern_Types below)",
                        "n_contigs": "Number of contigs in this cluster",
                        "contig_ids": "Pipe-delimited list of all contig IDs in this cluster",
                        "taxonomic_families": "Pipe-delimited list of unique taxonomic families represented",
                        "taxonomic_kingdoms": "Pipe-delimited list of unique taxonomic kingdoms represented",
                        "total_reads_sum": "Combined total reads from all contigs in cluster",
                        "mean_prevalence": "Average prevalence percentage across contigs in cluster"
                    }
                },
                
                "viral_family_taxonomy_mapping.csv": {
                    "description": "Mapping between viral families and their constituent contigs with full taxonomy",
                    "columns": {
                        "viral_family": "Viral family name",
                        "taxonomy": "Complete taxonomic hierarchy (semicolon-separated: kingdom:value; phylum:value; etc.)",
                        "contigs": "Pipe-delimited list of all contigs belonging to this family"
                    }
                }
            },
            
            "Pattern_Types": {
                "stable_present": "Contig present in all analyzed samples (100% prevalence)",
                "episodic": "Contig present intermittently - irregular pattern of presence/absence",
                "emerging": "Contig absent early in time series, becomes present in later samples", 
                "declining": "Contig present early in time series, becomes absent in later samples",
                "stable_absent": "Contig absent in all samples (should not appear in results)",
                "single_occurrence": "Contig present in only one sample",
                "highly_variable": "Contig alternates presence/absence in every consecutive sample"
            },
            
            "Abundance_Categories": {
                "H": "High abundance - above 75th percentile of non-zero values for this contig",
                "M": "Medium abundance - between 25th and 75th percentiles of non-zero values",
                "L": "Low abundance - below 25th percentile of non-zero values", 
                "O": "Absent - zero reads in this sample"
            },
            
            "Presence_Codes": {
                "X": "Present - contig has >0 reads in this sample",
                "O": "Absent - contig has 0 reads in this sample"
            },
            
            "Analysis_Parameters": {
                "min_samples": f"Minimum samples required for contig inclusion: {self.min_samples}",
                "min_total_reads": f"Minimum total reads required per contig: {self.min_total_reads}",
                "delimiter_rationale": "Pipe (|) delimiters used in multi-value fields to avoid CSV parsing conflicts"
            },
            
            "Interpretation_Notes": {
                "temporal_order": "Patterns are ordered chronologically based on sample metadata",
                "biological_significance": "Stable patterns may indicate core microbiome members, episodic patterns may reflect environmental exposures or dietary influences",
                "clinical_relevance": "Emerging/declining patterns may correlate with disease progression or treatment response",
                "clustering_utility": "Pattern clusters identify contigs with coordinated temporal behavior, potentially indicating functional relationships"
            }
        }
        
        # Save data dictionary as readable text file
        dict_path = self.output_dir / 'data_dictionary.txt'
        with open(dict_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("DATA DICTIONARY: IBD Virome Contig-Level Temporal Pattern Analysis\n")
            f.write("="*80 + "\n\n")
            f.write(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Analysis Parameters: min_samples={self.min_samples}, min_total_reads={self.min_total_reads}\n\n")
            
            f.write("-"*80 + "\n")
            f.write("FILE DESCRIPTIONS\n")
            f.write("-"*80 + "\n\n")
            
            # contig_temporal_patterns.csv
            f.write("1. contig_temporal_patterns.csv\n")
            f.write("   Description: Temporal patterns for individual contigs across all samples\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - contig_id: Unique contig identifier (pipe-delimited if needed to avoid CSV conflicts)\n")
            f.write("   - n_samples_present: Number of samples where this contig has >0 reads\n")
            f.write("   - total_samples: Total number of samples in the analysis\n")
            f.write("   - prevalence_pct: Percentage of samples where contig is present (n_samples_present/total_samples * 100)\n")
            f.write("   - presence_pattern: Binary presence pattern across samples (X=present, O=absent)\n")
            f.write("   - abundance_pattern: Abundance level pattern (H=high, M=medium, L=low, O=absent)\n")
            f.write("   - pattern_type: Classification of temporal pattern (see PATTERN TYPES below)\n")
            f.write("   - total_reads: Sum of all reads for this contig across all samples\n")
            f.write("   - mean_abundance: Average read count per sample (including zeros)\n")
            f.write("   - max_abundance: Maximum read count observed in any single sample\n")
            f.write("   - taxonomy_kingdom: Taxonomic kingdom assignment (e.g., Viruses, Bacteria)\n")
            f.write("   - taxonomy_family: Taxonomic family assignment\n")
            f.write("   - taxonomy_species: Taxonomic species assignment\n")
            f.write("   - cluster_id: Cluster identifier (>0 if contig belongs to a cluster with ≥2 members, 0 if singleton)\n\n")
            
            # contig_pattern_clusters.csv
            f.write("2. contig_pattern_clusters.csv\n")
            f.write("   Description: Groups of contigs that share identical temporal patterns\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - cluster_id: Unique identifier for each pattern cluster\n")
            f.write("   - presence_pattern: Shared binary presence pattern (X=present, O=absent)\n")
            f.write("   - abundance_pattern: Shared abundance pattern (H=high, M=medium, L=low, O=absent)\n")
            f.write("   - pattern_type: Classification of temporal pattern\n")
            f.write("   - n_contigs: Number of contigs in this cluster\n")
            f.write("   - contig_ids: Pipe-delimited list of all contig IDs in this cluster\n")
            f.write("   - taxonomic_families: Pipe-delimited list of unique taxonomic families represented\n")
            f.write("   - taxonomic_kingdoms: Pipe-delimited list of unique taxonomic kingdoms represented\n")
            f.write("   - total_reads_sum: Combined total reads from all contigs in cluster\n")
            f.write("   - mean_prevalence: Average prevalence percentage across contigs in cluster\n\n")
            
            # viral_family_taxonomy_mapping.csv
            f.write("3. viral_family_taxonomy_mapping.csv\n")
            f.write("   Description: Mapping between viral families and their constituent contigs with full taxonomy\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - viral_family: Viral family name\n")
            f.write("   - taxonomy: Complete taxonomic hierarchy (semicolon-separated)\n")
            f.write("   - contigs: Pipe-delimited list of all contigs belonging to this family\n\n")
            
            f.write("-"*80 + "\n")
            f.write("PATTERN TYPE DEFINITIONS\n")
            f.write("-"*80 + "\n\n")
            f.write("stable_present    = Contig present in all analyzed samples (100% prevalence)\n")
            f.write("episodic         = Contig present intermittently - irregular pattern of presence/absence\n")
            f.write("emerging         = Contig absent early in time series, becomes present in later samples\n")
            f.write("declining        = Contig present early in time series, becomes absent in later samples\n")
            f.write("single_occurrence = Contig present in only one sample\n")
            f.write("highly_variable  = Contig alternates presence/absence in every consecutive sample\n")
            f.write("stable_absent    = Contig absent in all samples (should not appear in results)\n\n")
            
            f.write("-"*80 + "\n")
            f.write("ABUNDANCE CATEGORIES\n")
            f.write("-"*80 + "\n\n")
            f.write("H = High abundance - above 75th percentile of non-zero values for this contig\n")
            f.write("M = Medium abundance - between 25th and 75th percentiles of non-zero values\n")
            f.write("L = Low abundance - below 25th percentile of non-zero values\n")
            f.write("O = Absent - zero reads in this sample\n\n")
            
            f.write("-"*80 + "\n")
            f.write("PRESENCE CODES\n")
            f.write("-"*80 + "\n\n")
            f.write("X = Present - contig has >0 reads in this sample\n")
            f.write("O = Absent - contig has 0 reads in this sample\n\n")
            
            f.write("-"*80 + "\n")
            f.write("INTERPRETATION NOTES\n")
            f.write("-"*80 + "\n\n")
            f.write("1. Temporal Order: Patterns are ordered chronologically based on sample metadata\n")
            f.write("2. Biological Significance: Stable patterns may indicate core microbiome members,\n")
            f.write("   while episodic patterns may reflect environmental exposures or dietary influences\n")
            f.write("3. Clinical Relevance: Emerging/declining patterns may correlate with disease\n")
            f.write("   progression or treatment response\n")
            f.write("4. Clustering Utility: Pattern clusters identify contigs with coordinated temporal\n")
            f.write("   behavior, potentially indicating functional relationships\n")
            f.write("5. Delimiter Choice: Pipe (|) delimiters used in multi-value fields to avoid\n")
            f.write("   CSV parsing conflicts\n\n")
            
            f.write("="*80 + "\n")
        
        print(f"    Data dictionary saved to {dict_path}")

def main():
    """
    Main function to run contig-level temporal analysis
    """
    parser = argparse.ArgumentParser(description='Contig-level temporal pattern analysis')
    parser.add_argument('--metadata', required=True, help='Path to metadata file')
    parser.add_argument('--counts', required=True, help='Path to counts file')
    parser.add_argument('--taxonomy', required=True, help='Path to taxonomy file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--min-samples', type=int, default=2, help='Minimum samples for inclusion')
    parser.add_argument('--min-total-reads', type=int, default=10, help='Minimum total reads per contig (field standard for virome studies)')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ContigTemporalAnalyzer(
        output_dir=args.output_dir,
        min_samples=args.min_samples,
        min_total_reads=args.min_total_reads
    )
    
    try:
        # Run analysis pipeline
        analyzer.load_data(args.metadata, args.counts, args.taxonomy)
        analyzer.align_data()
        analyzer.analyze_contig_patterns()
        analyzer.cluster_similar_patterns()
        analyzer.create_visualizations()
        analyzer.save_results()
        
        print("\nContig-level analysis complete!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()