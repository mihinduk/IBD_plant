#!/usr/bin/env python3
"""
Result Validation and Comparison Module

Ensures that refactored modular code produces identical results to the 
proven working analysis scripts.

Key validation targets:
- Temporal analysis: 11,694 contigs, 57 viral families
- Host response: 13 standardized output files
- Pipeline integrity: No data loss in format conversion
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any

class ResultValidator:
    """Validates analysis outputs against known working results"""
    
    def __init__(self, validation_mode='strict'):
        self.validation_mode = validation_mode
        self.known_results = {
            'temporal_analysis': {
                'total_contigs_analyzed': 11694,
                'viral_families_identified': 57,
                'pattern_types': {
                    'emerging': 7681,
                    'episodic': 3056,
                    'stable_present': 670,
                    'declining': 287
                }
            }
        }
    
    def validate_temporal_analysis(self, results_dir: Path, dataset_name: str = "UC_temporal") -> bool:
        """
        Validate temporal analysis results against known good outputs
        
        Parameters:
        -----------
        results_dir : Path
            Directory containing temporal analysis results
        dataset_name : str
            Name of the dataset being validated
            
        Returns:
        --------
        bool
            True if validation passes
        """
        print(f"🔍 Validating temporal analysis results in {results_dir}")
        
        # Check for expected output files
        expected_files = [
            'viral_temporal_patterns.csv',
            'contig_temporal_patterns.csv',
            'temporal_cluster_associations.csv',
            'family_temporal_patterns.csv'
        ]
        
        missing_files = []
        for file_name in expected_files:
            file_path = results_dir / file_name
            if not file_path.exists():
                missing_files.append(file_name)
        
        if missing_files:
            print(f"❌ Missing temporal analysis files: {missing_files}")
            return False
        
        # Validate contig counts (if this is the UC temporal dataset)
        if dataset_name == "UC_temporal":
            contig_file = results_dir / 'contig_temporal_patterns.csv'
            if contig_file.exists():
                df = pd.read_csv(contig_file)
                actual_contigs = len(df)
                expected_contigs = self.known_results['temporal_analysis']['total_contigs_analyzed']
                
                if actual_contigs != expected_contigs:
                    print(f"❌ Contig count mismatch: expected {expected_contigs}, got {actual_contigs}")
                    return False
                else:
                    print(f"✅ Contig count validated: {actual_contigs}")
        
        # Validate viral families (if family file exists)
        family_file = results_dir / 'family_temporal_patterns.csv'
        if family_file.exists() and dataset_name == "UC_temporal":
            df = pd.read_csv(family_file)
            actual_families = len(df['family'].unique()) if 'family' in df.columns else len(df)
            expected_families = self.known_results['temporal_analysis']['viral_families_identified']
            
            if actual_families != expected_families:
                print(f"❌ Viral family count mismatch: expected {expected_families}, got {actual_families}")
                return False
            else:
                print(f"✅ Viral family count validated: {actual_families}")
        
        print(f"✅ Temporal analysis validation passed for {dataset_name}")
        return True
    
    def validate_host_response_outputs(self, output_dir: Path) -> bool:
        """
        Validate host response analysis outputs
        
        Parameters:
        -----------
        output_dir : Path
            Directory containing host response results
            
        Returns:
        --------
        bool
            True if all expected files exist and are non-empty
        """
        print(f"🔍 Validating host response outputs in {output_dir}")
        
        expected_files = [
            'host_response_correlation_results.json',
            'host_response_analysis_summary.json',
            'viral_family_correlations.tsv',
            'viral_species_correlations.tsv',
            'host_response_networks_family.gml',
            'host_response_networks_species.gml',
            'correlation_heatmap_family.png',
            'correlation_heatmap_species.png',
            'network_plot_family.png',
            'network_plot_species.png',
            'host_response_signatures.tsv',
            'analysis.log',
            'job_summary.txt'
        ]
        
        missing_files = []
        empty_files = []
        
        for file_name in expected_files:
            file_path = output_dir / file_name
            
            if not file_path.exists():
                missing_files.append(file_name)
            elif file_path.stat().st_size == 0:
                empty_files.append(file_name)
        
        # Report results
        if missing_files:
            print(f"❌ Missing host response files: {missing_files}")
            return False
        
        if empty_files:
            print(f"⚠️  Empty host response files: {empty_files}")
            if self.validation_mode == 'strict':
                return False
        
        # Validate JSON files can be parsed
        json_files = [f for f in expected_files if f.endswith('.json')]
        for json_file in json_files:
            file_path = output_dir / json_file
            try:
                with open(file_path, 'r') as f:
                    json.load(f)
                print(f"✅ {json_file} is valid JSON")
            except json.JSONDecodeError as e:
                print(f"❌ {json_file} is invalid JSON: {e}")
                return False
        
        # Validate TSV files can be read
        tsv_files = [f for f in expected_files if f.endswith('.tsv')]
        for tsv_file in tsv_files:
            file_path = output_dir / tsv_file
            try:
                df = pd.read_csv(file_path, sep='\t')
                print(f"✅ {tsv_file} loaded successfully ({df.shape})")
            except Exception as e:
                print(f"❌ {tsv_file} cannot be loaded: {e}")
                return False
        
        print(f"✅ All 13 host response output files validated")
        return True
    
    def validate_format_conversion(self, original_df: pd.DataFrame, converted_df: pd.DataFrame) -> bool:
        """
        Validate that long→wide format conversion preserves data integrity
        
        Parameters:
        -----------
        original_df : pd.DataFrame
            Original long format data
        converted_df : pd.DataFrame
            Converted wide format data
            
        Returns:
        --------
        bool
            True if conversion is valid
        """
        print("🔍 Validating format conversion integrity")
        
        # Check total read counts are preserved
        if 'Reads' in original_df.columns:
            original_total = original_df['Reads'].sum()
            # Wide format: exclude contig_id column for sum
            sample_cols = [col for col in converted_df.columns if col != 'contig_id']
            converted_total = converted_df[sample_cols].sum().sum()
            
            if abs(original_total - converted_total) > 0.001:
                print(f"❌ Data loss in conversion: {original_total} → {converted_total}")
                return False
            else:
                print(f"✅ Read counts preserved: {original_total} total reads")
        
        # Check dimensions
        if 'Contig' in original_df.columns and 'Sample' in original_df.columns:
            expected_contigs = original_df['Contig'].nunique()
            expected_samples = original_df['Sample'].nunique()
            
            actual_contigs = len(converted_df)
            actual_samples = len([col for col in converted_df.columns if col != 'contig_id'])
            
            if expected_contigs != actual_contigs:
                print(f"❌ Contig count mismatch: {expected_contigs} → {actual_contigs}")
                return False
            
            if expected_samples != actual_samples:
                print(f"❌ Sample count mismatch: {expected_samples} → {actual_samples}")
                return False
            
            print(f"✅ Dimensions preserved: {actual_contigs} contigs × {actual_samples} samples")
        
        return True
    
    def validate_pipeline_stage(self, stage_name: str, input_data: Any, output_data: Any, 
                               stage_results: Dict) -> bool:
        """
        Validate a specific pipeline stage
        
        Parameters:
        -----------
        stage_name : str
            Name of the pipeline stage
        input_data : Any
            Input data to the stage
        output_data : Any
            Output data from the stage
        stage_results : Dict
            Stage-specific results and metadata
            
        Returns:
        --------
        bool
            True if stage validation passes
        """
        print(f"🔍 Validating pipeline stage: {stage_name}")
        
        # Stage-specific validation
        if stage_name == "stage1_data_prep":
            # Validate format conversion
            return self.validate_format_conversion(input_data, output_data)
        
        elif stage_name == "temporal_analysis":
            # Validate temporal analysis results
            if isinstance(stage_results, dict) and 'output_dir' in stage_results:
                return self.validate_temporal_analysis(Path(stage_results['output_dir']))
        
        elif stage_name == "host_response_analysis":
            # Validate host response results
            if isinstance(stage_results, dict) and 'output_dir' in stage_results:
                return self.validate_host_response_outputs(Path(stage_results['output_dir']))
        
        print(f"⚠️  No specific validation implemented for {stage_name}")
        return True
    
    def compare_analysis_results(self, old_results_dir: Path, new_results_dir: Path) -> bool:
        """
        Compare results from old vs new pipeline implementations
        
        Parameters:
        -----------
        old_results_dir : Path
            Directory with results from working legacy code
        new_results_dir : Path
            Directory with results from new modular code
            
        Returns:
        --------
        bool
            True if results are equivalent
        """
        print(f"🔍 Comparing old vs new results")
        print(f"   Old: {old_results_dir}")
        print(f"   New: {new_results_dir}")
        
        # Compare file lists
        old_files = set(f.name for f in old_results_dir.glob('*') if f.is_file())
        new_files = set(f.name for f in new_results_dir.glob('*') if f.is_file())
        
        missing_in_new = old_files - new_files
        extra_in_new = new_files - old_files
        
        if missing_in_new:
            print(f"❌ Files missing in new results: {missing_in_new}")
            return False
        
        if extra_in_new:
            print(f"ℹ️  Additional files in new results: {extra_in_new}")
        
        # Compare common CSV files
        common_csv_files = [f for f in old_files & new_files if f.endswith('.csv')]
        
        for csv_file in common_csv_files:
            old_df = pd.read_csv(old_results_dir / csv_file)
            new_df = pd.read_csv(new_results_dir / csv_file)
            
            if not old_df.equals(new_df):
                print(f"❌ Difference in {csv_file}")
                print(f"   Old shape: {old_df.shape}, New shape: {new_df.shape}")
                return False
            else:
                print(f"✅ {csv_file} identical")
        
        print("✅ All comparable results are identical")
        return True

def validate_all_outputs(results_dir: Path, dataset_name: str = "test") -> bool:
    """
    Convenience function to validate all analysis outputs
    
    Parameters:
    -----------
    results_dir : Path
        Directory containing all analysis results
    dataset_name : str
        Name of the dataset for context-specific validation
        
    Returns:
    --------
    bool
        True if all validations pass
    """
    validator = ResultValidator()
    
    all_passed = True
    
    # Check for temporal analysis results
    temporal_dirs = list(results_dir.glob("*temporal*"))
    for temporal_dir in temporal_dirs:
        if not validator.validate_temporal_analysis(temporal_dir, dataset_name):
            all_passed = False
    
    # Check for host response results
    host_response_dirs = list(results_dir.glob("*host_response*"))
    for hr_dir in host_response_dirs:
        if not validator.validate_host_response_outputs(hr_dir):
            all_passed = False
    
    return all_passed

if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        results_path = Path(sys.argv[1])
        dataset_name = sys.argv[2] if len(sys.argv) > 2 else "test"
        
        success = validate_all_outputs(results_path, dataset_name)
        sys.exit(0 if success else 1)
    else:
        print("Usage: python result_comparator.py <results_dir> [dataset_name]")
        sys.exit(1)