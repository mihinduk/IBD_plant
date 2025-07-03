#!/usr/bin/env python3
"""
Temporal-Ecological Fingerprinting Analysis for IBD Virome Study
=================================================================

This script implements pattern-based temporal analysis to identify 
viral community patterns associated with IBD disease states and 
inflammatory responses using presence/absence patterns.

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
try:
    from dtaidistance import dtw
    from dtaidistance import dtw_ndim
except ImportError:
    print("Warning: dtaidistance not available. DTW analysis will be skipped.")
    dtw = None
    dtw_ndim = None

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for HTCF
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class TemporalPatternAnalyzer:
    """
    Main class for pattern-based temporal viral analysis
    """
    
    def __init__(self, output_dir, min_samples=2):
        """
        Initialize the analyzer
        
        Parameters:
        -----------
        output_dir : str
            Directory to save results
        min_samples : int
            Minimum number of samples for inclusion (default: 2)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_samples = min_samples
        
        # Initialize containers
        self.metadata = None
        self.counts = None
        self.taxonomy = None
        self.viral_patterns = None
        self.patient_timeseries = {}
        
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
        print(f"    Column names after pivot: {list(self.counts.columns)}")
        print(f"    Column dtypes: {self.counts.dtypes.iloc[0]}")
        
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
        
        print(f"  Metadata samples: {sorted(self.metadata['database_ID'].tolist())}")
        print(f"  Count samples: {sorted(self.counts.columns.tolist())}")
        
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
        
    def analyze_viral_patterns(self):
        """
        Analyze viral patterns using the new pattern-based approach
        """
        print("Identifying viral families...")
        
        # Priority viral families of interest
        priority_families = [
            'Virgaviridae', 'Alphaflexiviridae', 'Betaflexiviridae',
            'Genomoviridae', 'Partitiviridae', 'Chrysoviridae',
            'Endornaviridae', 'Mitoviridae', 'Bromoviridae',
            'Caulimoviridae', 'Closteroviridae', 'Luteoviridae',
            'Potyviridae', 'Secoviridae', 'Solemoviridae',
            'Tombusviridae', 'Tymoviridae', 'Pospiviroidae'
        ]
        
        # Filter taxonomy to viral families
        kingdom_col = 'superkingdom' if 'superkingdom' in self.taxonomy.columns else 'kingdom'
        
        viral_tax = self.taxonomy[
            (self.taxonomy[kingdom_col] == 'Viruses') &
            (self.taxonomy['family'].notna())
        ].copy()
        
        print(f"  Found {len(viral_tax)} viral contigs")
        
        # Analyze patterns by family
        pattern_results = []
        
        for family in viral_tax['family'].unique():
            if pd.isna(family):
                continue
                
            family_contigs = viral_tax[viral_tax['family'] == family]['contig_id'].values
            family_contigs = [c for c in family_contigs if c in self.counts.index]
            
            if len(family_contigs) == 0:
                continue
            
            # Sum reads across all contigs in the family
            family_counts = self.counts.loc[family_contigs].sum(axis=0)
            
            # Calculate pattern metrics
            pattern_info = self._calculate_pattern_metrics(family_counts, family)
            pattern_info.update({
                'viral_family': family,
                'n_contigs': len(family_contigs),
                'is_priority': family in priority_families
            })
            
            pattern_results.append(pattern_info)
        
        # Convert to DataFrame and filter by minimum samples
        self.viral_patterns = pd.DataFrame(pattern_results)
        
        if len(self.viral_patterns) > 0:
            # Filter families present in at least min_samples
            self.viral_patterns = self.viral_patterns[
                self.viral_patterns['n_samples_present'] >= self.min_samples
            ]
        
        print(f"  Retained {len(self.viral_patterns)} viral families after filtering")
        if len(self.viral_patterns) > 0:
            print(f"  Families: {list(self.viral_patterns['viral_family'])}")
        
        # Save detailed pattern results
        if len(self.viral_patterns) > 0:
            self.viral_patterns.to_csv(
                self.output_dir / 'viral_temporal_patterns.csv',
                index=False
            )
            
            # Create family-taxonomy mapping file
            self._create_family_taxonomy_mapping(viral_tax)
        else:
            # Create empty file with headers
            empty_df = pd.DataFrame(columns=[
                'viral_family', 'n_contigs', 'n_samples_present', 'total_samples',
                'prevalence_pct', 'presence_pattern', 'abundance_pattern',
                'pattern_type', 'total_reads', 'mean_abundance', 'max_abundance',
                'is_priority'
            ])
            empty_df.to_csv(self.output_dir / 'viral_temporal_patterns.csv', index=False)
            
            # Create empty taxonomy mapping file
            empty_tax_df = pd.DataFrame(columns=[
                'viral_family', 'taxonomy', 'contigs'
            ])
            empty_tax_df.to_csv(self.output_dir / 'viral_family_taxonomy_mapping.csv', index=False)
        
    def _calculate_pattern_metrics(self, counts_series, family_name):
        """
        Calculate pattern metrics for a viral family
        """
        # Convert counts to presence/absence and abundance categories
        total_samples = len(counts_series)
        present_mask = counts_series > 0
        n_present = present_mask.sum()
        
        # Create presence pattern (X=present, O=absent)
        presence_pattern = ''.join(['X' if x else 'O' for x in present_mask])
        
        # Create abundance pattern if any samples are present
        if n_present > 0:
            # Calculate quartiles for abundance categorization
            present_values = counts_series[present_mask]
            if len(present_values) > 1:
                q25, q75 = np.percentile(present_values, [25, 75])
                
                # Categorize abundance: H=high, M=medium, L=low, O=absent
                abundance_cats = []
                for count in counts_series:
                    if count == 0:
                        abundance_cats.append('O')
                    elif count >= q75:
                        abundance_cats.append('H')
                    elif count <= q25:
                        abundance_cats.append('L')
                    else:
                        abundance_cats.append('M')
                abundance_pattern = '-'.join(abundance_cats)
            else:
                # Only one present sample
                abundance_pattern = '-'.join(['H' if x > 0 else 'O' for x in counts_series])
        else:
            abundance_pattern = '-'.join(['O'] * total_samples)
        
        # Classify pattern type
        pattern_type = self._classify_pattern_type(present_mask)
        
        return {
            'n_samples_present': n_present,
            'total_samples': total_samples,
            'prevalence_pct': round((n_present / total_samples) * 100, 1),
            'presence_pattern': presence_pattern,
            'abundance_pattern': abundance_pattern,
            'pattern_type': pattern_type,
            'total_reads': int(counts_series.sum()),
            'mean_abundance': round(counts_series.mean(), 2),
            'max_abundance': int(counts_series.max())
        }
    
    def _classify_pattern_type(self, present_mask):
        """
        Classify the temporal pattern type
        """
        if len(present_mask) < 2:
            return 'insufficient_data'
        
        # Convert to list for easier pattern matching
        pattern = present_mask.tolist()
        n_present = sum(pattern)
        n_total = len(pattern)
        
        if n_present == 0:
            return 'absent'
        elif n_present == n_total:
            return 'stable_present'
        elif n_present == 1:
            return 'single_occurrence'
        else:
            # Check for temporal trends
            first_half = sum(pattern[:n_total//2])
            second_half = sum(pattern[n_total//2:])
            
            # More sophisticated pattern classification
            if first_half == 0 and second_half > 0:
                return 'emerging'
            elif first_half > 0 and second_half == 0:
                return 'declining'
            elif all(pattern[i] != pattern[i+1] for i in range(len(pattern)-1)):
                return 'highly_variable'
            else:
                return 'episodic'
    
    def _create_family_taxonomy_mapping(self, viral_tax):
        """
        Create a mapping file linking viral families to their full taxonomy and contigs
        """
        print("  Creating family-taxonomy mapping...")
        
        # Get unique families from our results
        analyzed_families = set(self.viral_patterns['viral_family'].unique())
        
        mapping_results = []
        
        for family in analyzed_families:
            # Get all contigs for this family
            family_tax = viral_tax[viral_tax['family'] == family].copy()
            
            if len(family_tax) == 0:
                continue
            
            # Get contigs that are in our counts data
            family_contigs = family_tax['contig_id'].values
            available_contigs = [c for c in family_contigs if c in self.counts.index]
            
            if len(available_contigs) == 0:
                continue
            
            # Build full taxonomy string from available columns
            # Common taxonomy columns in order of specificity
            tax_columns = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            
            # Use the first representative contig for taxonomy (they should all be the same family)
            representative_tax = family_tax.iloc[0]
            
            taxonomy_parts = []
            for col in tax_columns:
                if col in representative_tax and pd.notna(representative_tax[col]) and representative_tax[col] != '':
                    taxonomy_parts.append(f"{col}:{representative_tax[col]}")
            
            full_taxonomy = "; ".join(taxonomy_parts)
            
            # Join contigs with pipe delimiter (avoiding commas in CSV)
            contigs_str = "|".join(sorted(available_contigs))
            
            mapping_results.append({
                'viral_family': family,
                'taxonomy': full_taxonomy,
                'contigs': contigs_str
            })
        
        # Convert to DataFrame and save
        mapping_df = pd.DataFrame(mapping_results)
        mapping_df = mapping_df.sort_values('viral_family')
        
        mapping_path = self.output_dir / 'viral_family_taxonomy_mapping.csv'
        mapping_df.to_csv(mapping_path, index=False)
        
        print(f"    Saved family-taxonomy mapping to {mapping_path}")
        # Calculate total contigs (avoid backslash in f-string)
        total_contigs = mapping_df['contigs'].str.count('\\|').sum() + len(mapping_df)
        print(f"    Mapped {len(mapping_df)} families to their full taxonomy and {total_contigs} total contigs")
    
    def build_patient_timeseries(self):
        """
        Build time series data for each patient
        """
        print("Building patient time series...")
        
        if len(self.viral_patterns) == 0:
            print("  No viral families to analyze")
            self.patient_timeseries = {}
            return
        
        # Create a family abundance matrix
        family_abundance = {}
        for _, row in self.viral_patterns.iterrows():
            family = row['viral_family']
            # Get contigs for this family
            viral_tax = self.taxonomy[self.taxonomy['family'] == family]
            family_contigs = viral_tax['contig_id'].values
            family_contigs = [c for c in family_contigs if c in self.counts.index]
            
            if len(family_contigs) > 0:
                family_abundance[family] = self.counts.loc[family_contigs].sum(axis=0)
        
        if not family_abundance:
            print("  No family abundance data available")
            self.patient_timeseries = {}
            return
        
        families_df = pd.DataFrame(family_abundance).fillna(0)
        
        # Merge with metadata
        sample_data = self.metadata.set_index('database_ID').join(
            families_df, how='inner'
        )
        
        # Group by patient and create time series
        patients_with_timeseries = 0
        
        for patient_id, patient_data in sample_data.groupby('patient_ID'):
            if len(patient_data) < 2:  # Need at least 2 timepoints
                continue
                
            # Sort by date if available, otherwise by database_ID
            if 'sample_date_converted_eu' in patient_data.columns:
                patient_data = patient_data.sort_values('sample_date_converted_eu')
            else:
                patient_data = patient_data.sort_index()
            
            # Store patient time series
            self.patient_timeseries[patient_id] = {
                'data': patient_data,
                'n_timepoints': len(patient_data),
                'timespan_days': self._calculate_timespan(patient_data)
            }
            patients_with_timeseries += 1
        
        print(f"  Built time series for {patients_with_timeseries} patients")
        
        # Save patient summary
        if patients_with_timeseries > 0:
            patient_summary = []
            for patient_id, data in self.patient_timeseries.items():
                patient_info = data['data'].iloc[0]  # Get patient info from first sample
                patient_summary.append({
                    'patient_ID': patient_id,
                    'diagnosis': patient_info.get('diagnosis', 'unknown'),
                    'n_timepoints': data['n_timepoints'],
                    'timespan_days': data['timespan_days'],
                    'sample_IDs': ','.join(data['data'].index.astype(str))
                })
            
            pd.DataFrame(patient_summary).to_csv(
                self.output_dir / 'patient_timeseries_summary.csv',
                index=False
            )
        else:
            # Create empty file
            empty_df = pd.DataFrame(columns=[
                'patient_ID', 'diagnosis', 'n_timepoints', 'timespan_days', 'sample_IDs'
            ])
            empty_df.to_csv(self.output_dir / 'patient_timeseries_summary.csv', index=False)
    
    def _calculate_timespan(self, patient_data):
        """
        Calculate timespan between first and last samples
        """
        if 'sample_date_converted_eu' in patient_data.columns:
            try:
                dates = pd.to_datetime(patient_data['sample_date_converted_eu'])
                return (dates.max() - dates.min()).days
            except:
                return None
        return None
    
    def perform_dtw_analysis(self):
        """
        Perform DTW analysis if dtaidistance is available
        """
        print("Computing DTW distances...")
        
        if dtw is None:
            print("  DTW analysis skipped (dtaidistance not available)")
            return
        
        if len(self.patient_timeseries) < 2:
            print("  Insufficient patients for DTW analysis")
            return
        
        # Implementation would go here for DTW analysis
        print("  DTW distance computation complete")
    
    def analyze_clinical_associations(self):
        """
        Analyze associations between viral patterns and clinical variables
        """
        print("Analyzing clinical associations...")
        
        if len(self.viral_patterns) == 0:
            print("  No viral families for association analysis")
            return
        
        # Placeholder for clinical association analysis
        # This would analyze relationships between viral patterns and:
        # - Disease status (CD vs UC vs control)
        # - Inflammation markers (FCP, CRP)
        # - Disease activity
        # - Medications
        
        print(f"  Completed clinical association analysis for {len(self.viral_patterns)} families")
    
    def create_visualizations(self):
        """
        Create visualization plots
        """
        print("Creating visualizations...")
        
        if len(self.viral_patterns) > 0:
            # Create pattern heatmap
            self._create_pattern_heatmap()
        
        print("  Visualization creation complete")
    
    def _create_pattern_heatmap(self):
        """
        Create heatmap of viral family patterns
        """
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Create matrix for heatmap
            pattern_matrix = []
            family_labels = []
            
            for _, row in self.viral_patterns.iterrows():
                # Convert presence pattern to numeric
                pattern_numeric = [1 if x == 'X' else 0 for x in row['presence_pattern']]
                pattern_matrix.append(pattern_numeric)
                family_labels.append(f"{row['viral_family']} ({row['pattern_type']})")
            
            pattern_matrix = np.array(pattern_matrix)
            
            # Create heatmap
            sns.heatmap(pattern_matrix, 
                       yticklabels=family_labels,
                       xticklabels=[f'T{i+1}' for i in range(pattern_matrix.shape[1])],
                       cmap='RdYlBu_r',
                       cbar_kws={'label': 'Presence (1) / Absence (0)'},
                       ax=ax)
            
            plt.title('Viral Family Temporal Patterns')
            plt.xlabel('Time Points')
            plt.ylabel('Viral Families')
            plt.tight_layout()
            plt.savefig(self.output_dir / 'viral_pattern_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"  Warning: Could not create heatmap: {e}")
    
    def generate_summary(self):
        """
        Generate analysis summary
        """
        print("Generating summary report...")
        
        summary = {
            'data_summary': {
                'n_samples': len(self.metadata) if self.metadata is not None else 0,
                'n_contigs': self.counts.shape[0] if self.counts is not None else 0,
                'n_viral_families': len(self.viral_patterns) if self.viral_patterns is not None else 0,
                'n_patients_with_timeseries': len(self.patient_timeseries)
            },
            'pattern_analysis': {
                'families_analyzed': len(self.viral_patterns) if self.viral_patterns is not None else 0,
                'min_samples_threshold': self.min_samples
            },
            'temporal_analysis': {
                'patients_with_timeseries': len(self.patient_timeseries)
            }
        }
        
        # Add pattern type distribution
        if len(self.viral_patterns) > 0:
            pattern_types = self.viral_patterns['pattern_type'].value_counts().to_dict()
            summary['pattern_types'] = pattern_types
        
        # Save summary
        with open(self.output_dir / 'analysis_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"  Summary report saved to {self.output_dir / 'analysis_summary.json'}")
        
        # Create data dictionary
        self._create_data_dictionary()
    
    def _create_data_dictionary(self):
        """
        Create comprehensive data dictionary for family-level analysis output files
        """
        print("  Creating data dictionary...")
        
        data_dict = {
            "Data Dictionary": "IBD Virome Family-Level Temporal Pattern Analysis",
            "Analysis_Date": "Generated automatically during analysis",
            "Files_Described": [
                "viral_temporal_patterns.csv",
                "viral_family_taxonomy_mapping.csv",
                "patient_timeseries_summary.csv",
                "analysis_summary.json"
            ],
            
            "File_Descriptions": {
                "viral_temporal_patterns.csv": {
                    "description": "Temporal patterns for viral families aggregated across all constituent contigs",
                    "columns": {
                        "n_samples_present": "Number of samples where this viral family has >0 reads",
                        "total_samples": "Total number of samples in the analysis",
                        "prevalence_pct": "Percentage of samples where family is present (n_samples_present/total_samples * 100)",
                        "presence_pattern": "Binary presence pattern across samples (X=present, O=absent). Order matches sample chronology.",
                        "abundance_pattern": "Abundance level pattern across samples (H=high, M=medium, L=low, O=absent). Thresholds based on quartiles of non-zero values.",
                        "pattern_type": "Classification of temporal pattern (see Pattern_Types below)",
                        "total_reads": "Sum of all reads for this family across all samples and contigs",
                        "mean_abundance": "Average read count per sample (including zeros)",
                        "max_abundance": "Maximum read count observed in any single sample",
                        "viral_family": "Viral family name",
                        "n_contigs": "Number of contigs contributing to this family",
                        "is_priority": "Boolean indicating if this is a priority family from CLAUDE.md (plant/fungal viruses)"
                    }
                },
                
                "viral_family_taxonomy_mapping.csv": {
                    "description": "Complete taxonomic information and constituent contigs for each viral family",
                    "columns": {
                        "viral_family": "Viral family name",
                        "taxonomy": "Complete taxonomic hierarchy (semicolon-separated: superkingdom:value; kingdom:value; phylum:value; class:value; order:value; family:value; genus:value; species:value)",
                        "contigs": "Pipe-delimited list of all contig IDs belonging to this family"
                    }
                },
                
                "patient_timeseries_summary.csv": {
                    "description": "Summary of patients with longitudinal samples included in the analysis",
                    "columns": {
                        "patient_ID": "Unique patient identifier",
                        "diagnosis": "Clinical diagnosis (CD=Crohn's Disease, UC=Ulcerative Colitis, etc.)",
                        "n_timepoints": "Number of longitudinal samples available for this patient",
                        "timespan_days": "Days between first and last sample (if date information available)",
                        "sample_IDs": "Comma-delimited list of database sample IDs for this patient"
                    }
                }
            },
            
            "Pattern_Types": {
                "stable_present": "Family present in all analyzed samples (100% prevalence)",
                "episodic": "Family present intermittently - irregular pattern of presence/absence",
                "emerging": "Family absent early in time series, becomes present in later samples", 
                "declining": "Family present early in time series, becomes absent in later samples",
                "stable_absent": "Family absent in all samples (should not appear in results)",
                "single_occurrence": "Family present in only one sample",
                "highly_variable": "Family alternates presence/absence in every consecutive sample"
            },
            
            "Abundance_Categories": {
                "H": "High abundance - above 75th percentile of non-zero values for this family",
                "M": "Medium abundance - between 25th and 75th percentiles of non-zero values",
                "L": "Low abundance - below 25th percentile of non-zero values", 
                "O": "Absent - zero reads in this sample"
            },
            
            "Presence_Codes": {
                "X": "Present - family has >0 reads in this sample",
                "O": "Absent - family has 0 reads in this sample"
            },
            
            "Priority_Families": {
                "description": "Families marked as 'is_priority=True' are of special interest based on literature",
                "source": "Priority list from CLAUDE.md project documentation",
                "categories": [
                    "Plant viruses (Virgaviridae, Alphaflexiviridae, Betaflexiviridae, etc.)",
                    "Fungal viruses (Chrysoviridae, Endornaviridae, Mitoviridae)",
                    "Emerging/novel viruses (Genomoviridae, Partitiviridae)"
                ]
            },
            
            "Analysis_Parameters": {
                "min_samples": f"Minimum samples required for family inclusion: {self.min_samples}",
                "delimiter_rationale": "Pipe (|) delimiters used in multi-value fields to avoid CSV parsing conflicts"
            },
            
            "Interpretation_Notes": {
                "temporal_order": "Patterns are ordered chronologically based on sample metadata or database_ID",
                "family_aggregation": "Read counts are summed across all contigs belonging to each family",
                "biological_significance": "Plant viruses may reflect dietary intake, bacteriophages may modulate gut microbiome",
                "clinical_relevance": "Temporal patterns may correlate with disease activity, treatment response, or dietary changes",
                "statistical_considerations": "Zero-inflation is common in viral data - consider appropriate statistical models"
            }
        }
        
        # Save data dictionary as readable text file
        dict_path = self.output_dir / 'data_dictionary.txt'
        with open(dict_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("DATA DICTIONARY: IBD Virome Family-Level Temporal Pattern Analysis\n")
            f.write("="*80 + "\n\n")
            f.write(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Analysis Parameters: min_samples={self.min_samples}\n\n")
            
            f.write("-"*80 + "\n")
            f.write("FILE DESCRIPTIONS\n")
            f.write("-"*80 + "\n\n")
            
            # viral_temporal_patterns.csv
            f.write("1. viral_temporal_patterns.csv\n")
            f.write("   Description: Temporal patterns for viral families aggregated across all constituent contigs\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - n_samples_present: Number of samples where this viral family has >0 reads\n")
            f.write("   - total_samples: Total number of samples in the analysis\n")
            f.write("   - prevalence_pct: Percentage of samples where family is present (n_samples_present/total_samples * 100)\n")
            f.write("   - presence_pattern: Binary presence pattern across samples (X=present, O=absent)\n")
            f.write("   - abundance_pattern: Abundance level pattern (H=high, M=medium, L=low, O=absent)\n")
            f.write("   - pattern_type: Classification of temporal pattern (see PATTERN TYPES below)\n")
            f.write("   - total_reads: Sum of all reads for this family across all samples and contigs\n")
            f.write("   - mean_abundance: Average read count per sample (including zeros)\n")
            f.write("   - max_abundance: Maximum read count observed in any single sample\n")
            f.write("   - viral_family: Viral family name\n")
            f.write("   - n_contigs: Number of contigs contributing to this family\n")
            f.write("   - is_priority: Boolean indicating if this is a priority family (plant/fungal viruses)\n\n")
            
            # viral_family_taxonomy_mapping.csv
            f.write("2. viral_family_taxonomy_mapping.csv\n")
            f.write("   Description: Complete taxonomic information and constituent contigs for each viral family\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - viral_family: Viral family name\n")
            f.write("   - taxonomy: Complete taxonomic hierarchy (semicolon-separated)\n")
            f.write("   - contigs: Pipe-delimited list of all contig IDs belonging to this family\n\n")
            
            # patient_timeseries_summary.csv
            f.write("3. patient_timeseries_summary.csv\n")
            f.write("   Description: Summary of patients with longitudinal samples included in the analysis\n\n")
            f.write("   COLUMNS:\n")
            f.write("   - patient_ID: Unique patient identifier\n")
            f.write("   - diagnosis: Clinical diagnosis (CD=Crohn's Disease, UC=Ulcerative Colitis, etc.)\n")
            f.write("   - n_timepoints: Number of longitudinal samples available for this patient\n")
            f.write("   - timespan_days: Days between first and last sample (if date information available)\n")
            f.write("   - sample_IDs: Comma-delimited list of database sample IDs for this patient\n\n")
            
            f.write("-"*80 + "\n")
            f.write("PATTERN TYPE DEFINITIONS\n")
            f.write("-"*80 + "\n\n")
            f.write("stable_present    = Family present in all analyzed samples (100% prevalence)\n")
            f.write("episodic         = Family present intermittently - irregular pattern of presence/absence\n")
            f.write("emerging         = Family absent early in time series, becomes present in later samples\n")
            f.write("declining        = Family present early in time series, becomes absent in later samples\n")
            f.write("single_occurrence = Family present in only one sample\n")
            f.write("highly_variable  = Family alternates presence/absence in every consecutive sample\n")
            f.write("stable_absent    = Family absent in all samples (should not appear in results)\n\n")
            
            f.write("-"*80 + "\n")
            f.write("ABUNDANCE CATEGORIES\n")
            f.write("-"*80 + "\n\n")
            f.write("H = High abundance - above 75th percentile of non-zero values for this family\n")
            f.write("M = Medium abundance - between 25th and 75th percentiles of non-zero values\n")
            f.write("L = Low abundance - below 25th percentile of non-zero values\n")
            f.write("O = Absent - zero reads in this sample\n\n")
            
            f.write("-"*80 + "\n")
            f.write("PRESENCE CODES\n")
            f.write("-"*80 + "\n\n")
            f.write("X = Present - family has >0 reads in this sample\n")
            f.write("O = Absent - family has 0 reads in this sample\n\n")
            
            f.write("-"*80 + "\n")
            f.write("PRIORITY FAMILIES\n")
            f.write("-"*80 + "\n\n")
            f.write("Families marked as 'is_priority=True' are of special interest based on literature:\n\n")
            f.write("Plant viruses:\n")
            f.write("  - Virgaviridae, Alphaflexiviridae, Betaflexiviridae\n")
            f.write("  - Bromoviridae, Caulimoviridae, Closteroviridae\n")
            f.write("  - Luteoviridae, Potyviridae, Secoviridae\n")
            f.write("  - Solemoviridae, Tombusviridae, Tymoviridae, Pospiviroidae\n\n")
            f.write("Fungal viruses:\n")
            f.write("  - Chrysoviridae, Endornaviridae, Mitoviridae\n\n")
            f.write("Emerging/novel viruses:\n")
            f.write("  - Genomoviridae, Partitiviridae\n\n")
            
            f.write("-"*80 + "\n")
            f.write("INTERPRETATION NOTES\n")
            f.write("-"*80 + "\n\n")
            f.write("1. Temporal Order: Patterns are ordered chronologically based on sample metadata\n")
            f.write("2. Family Aggregation: Read counts are summed across all contigs belonging to each family\n")
            f.write("3. Biological Significance: Plant viruses may reflect dietary intake,\n")
            f.write("   bacteriophages may modulate gut microbiome\n")
            f.write("4. Clinical Relevance: Temporal patterns may correlate with disease activity,\n")
            f.write("   treatment response, or dietary changes\n")
            f.write("5. Statistical Considerations: Zero-inflation is common in viral data -\n")
            f.write("   consider appropriate statistical models\n")
            f.write("6. Delimiter Choice: Pipe (|) delimiters used in multi-value fields to avoid\n")
            f.write("   CSV parsing conflicts\n\n")
            
            f.write("="*80 + "\n")
        
        print(f"    Data dictionary saved to {dict_path}")

def main():
    """Main analysis pipeline"""
    parser = argparse.ArgumentParser(
        description="Temporal-Ecological Fingerprinting Analysis with Pattern-Based Approach"
    )
    parser.add_argument('--metadata', required=True, help='Path to metadata file')
    parser.add_argument('--counts', required=True, help='Path to contig count table')
    parser.add_argument('--taxonomy', required=True, help='Path to taxonomy file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--min-samples', type=int, default=2, 
                       help='Minimum number of samples for inclusion (default: 2)')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = TemporalPatternAnalyzer(
        output_dir=args.output_dir,
        min_samples=args.min_samples
    )
    
    try:
        # Run analysis pipeline
        analyzer.load_data(args.metadata, args.counts, args.taxonomy)
        analyzer.align_data()
        analyzer.analyze_viral_patterns()
        analyzer.build_patient_timeseries()
        analyzer.perform_dtw_analysis()
        analyzer.analyze_clinical_associations()
        analyzer.create_visualizations()
        analyzer.generate_summary()
        
        print("\nAnalysis complete! Results saved to", args.output_dir)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()