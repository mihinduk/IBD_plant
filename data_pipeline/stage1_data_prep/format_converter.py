#!/usr/bin/env python3
"""
Abundance Data Format Converter

Handles conversion between long and wide format abundance data for the IBD plant virome pipeline.
This module specifically addresses the Stage 1 format conversion issue.

Long format (input):
    Sample  Contig         Reads
    7562    contig_100018    2
    7562    contig_100067    5
    7583    contig_100018    0
    
Wide format (output):
    contig_id        7562  7583  7584  ...
    contig_100018      2     0     1
    contig_100067      5     3     0
"""

import pandas as pd
import numpy as np
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import sys
import os

# Add config directory to path to import constants
sys.path.append(str(Path(__file__).parent.parent.parent / 'config'))
try:
    from pipeline_constants import AbundanceColumns, MetadataColumns
except ImportError:
    print("Warning: Could not import pipeline_constants. Using fallback column names.")
    
    class AbundanceColumns:
        SAMPLE_ALTERNATIVES = ["Sample", "sample_id", "database_ID"]
        CONTIG_ALTERNATIVES = ["Contig", "contig_id", "contig_name"]
        READS_ALTERNATIVES = ["Reads", "count", "abundance"]

class AbundanceFormatConverter:
    """Converts abundance data between long and wide formats with validation"""
    
    def __init__(self, validation_mode: bool = True):
        """
        Initialize the format converter
        
        Parameters:
        -----------
        validation_mode : bool
            If True, perform strict validation of conversion results
        """
        self.validation_mode = validation_mode
        
    def _find_column(self, df: pd.DataFrame, alternatives: List[str]) -> str:
        """
        Find the actual column name from list of alternatives
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame to search in
        alternatives : List[str]
            List of possible column names
            
        Returns:
        --------
        str
            The actual column name found
            
        Raises:
        -------
        ValueError
            If none of the alternatives are found
        """
        for alt in alternatives:
            if alt in df.columns:
                return alt
        raise ValueError(f"None of {alternatives} found in columns: {list(df.columns)}")
    
    def _detect_format(self, df: pd.DataFrame) -> str:
        """
        Detect whether data is in long or wide format
        
        Parameters:
        -----------
        df : pd.DataFrame
            Input dataframe
            
        Returns:
        --------
        str
            'long' or 'wide'
        """
        # Check for long format indicators
        has_sample_col = any(col in df.columns for col in AbundanceColumns.SAMPLE_ALTERNATIVES)
        has_contig_col = any(col in df.columns for col in AbundanceColumns.CONTIG_ALTERNATIVES)
        has_reads_col = any(col in df.columns for col in AbundanceColumns.READS_ALTERNATIVES)
        
        # Check for wide format indicators (contig_id + many sample columns)
        has_contig_id = 'contig_id' in df.columns
        sample_like_cols = [col for col in df.columns if col.startswith(('S', 'sample', 'P'))]
        
        if has_sample_col and has_contig_col and has_reads_col:
            print(f"Detected long format: sample={has_sample_col}, contig={has_contig_col}, reads={has_reads_col}")
            return "long"
        elif has_contig_id and len(sample_like_cols) > 5:
            print(f"Detected wide format: contig_id={has_contig_id}, {len(sample_like_cols)} sample columns")
            return "wide"
        else:
            raise ValueError(
                f"Unable to determine format. Columns: {list(df.columns)[:10]}\n"
                f"Long format indicators: sample={has_sample_col}, contig={has_contig_col}, reads={has_reads_col}\n"
                f"Wide format indicators: contig_id={has_contig_id}, sample_cols={len(sample_like_cols)}"
            )
    
    def _pivot_data(self, abundance_df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert long format to wide format using pivot operation
        
        Parameters:
        -----------
        abundance_df : pd.DataFrame
            Long format abundance data
            
        Returns:
        --------
        pd.DataFrame
            Wide format abundance data
        """
        print(f"Converting long format data: {abundance_df.shape}")
        print(f"Columns: {list(abundance_df.columns)}")
        
        # Identify columns using the alternative names
        sample_col = self._find_column(abundance_df, AbundanceColumns.SAMPLE_ALTERNATIVES)
        contig_col = self._find_column(abundance_df, AbundanceColumns.CONTIG_ALTERNATIVES)
        reads_col = self._find_column(abundance_df, AbundanceColumns.READS_ALTERNATIVES)
        
        print(f"Detected columns: Sample='{sample_col}', Contig='{contig_col}', Reads='{reads_col}'")
        
        # Show sample data
        print(f"Sample values: {sorted(abundance_df[sample_col].unique())[:5]}...")
        print(f"Contig count: {abundance_df[contig_col].nunique()}")
        print(f"Total reads: {abundance_df[reads_col].sum():,}")
        
        # Perform the pivot operation
        print("Performing pivot operation...")
        wide_df = abundance_df.pivot(
            index=contig_col,      # Rows = contigs
            columns=sample_col,    # Columns = samples  
            values=reads_col       # Values = read counts
        )
        
        # Handle missing values (fill with 0 - samples where contig wasn't detected)
        wide_df = wide_df.fillna(0)
        
        # Convert to integers (read counts should be whole numbers)
        wide_df = wide_df.astype(int)
        
        # Reset index to make contig column a regular column
        wide_df = wide_df.reset_index()
        
        # Remove the column name from the columns index (clean up pandas formatting)
        wide_df.columns.name = None
        
        # Rename the contig column to standard name
        wide_df = wide_df.rename(columns={contig_col: 'contig_id'})
        
        print(f"Converted to wide format: {wide_df.shape}")
        sample_columns = [col for col in wide_df.columns if col != 'contig_id']
        print(f"Sample columns: {sample_columns[:5]}... (showing first 5 of {len(sample_columns)})")
        
        return wide_df
    
    def _validate_conversion(self, original_df: pd.DataFrame, converted_df: pd.DataFrame) -> None:
        """
        Validate that conversion preserved the data correctly
        
        Parameters:
        -----------
        original_df : pd.DataFrame
            Original long format data
        converted_df : pd.DataFrame
            Converted wide format data
            
        Raises:
        -------
        AssertionError
            If validation fails
        """
        if not self.validation_mode:
            return
        
        print("🔍 Validating conversion...")
        
        # Find the reads column in original data
        reads_col = self._find_column(original_df, AbundanceColumns.READS_ALTERNATIVES)
        
        # Check that no data was lost
        original_total = original_df[reads_col].sum()
        sample_cols = [col for col in converted_df.columns if col != 'contig_id']
        converted_total = converted_df[sample_cols].sum().sum()
        
        assert abs(original_total - converted_total) < 0.001, \
            f"Data loss in conversion: {original_total} → {converted_total}"
        
        # Check dimensions make sense
        sample_col = self._find_column(original_df, AbundanceColumns.SAMPLE_ALTERNATIVES)
        contig_col = self._find_column(original_df, AbundanceColumns.CONTIG_ALTERNATIVES)
        
        n_contigs = original_df[contig_col].nunique()
        n_samples = original_df[sample_col].nunique()
        
        assert converted_df.shape[0] == n_contigs, \
            f"Wrong number of contigs: expected {n_contigs}, got {converted_df.shape[0]}"
        
        assert converted_df.shape[1] == n_samples + 1, \
            f"Wrong number of columns: expected {n_samples + 1}, got {converted_df.shape[1]}"
        
        # Check that all original samples are represented
        original_samples = set(original_df[sample_col].unique())
        converted_samples = set(col for col in converted_df.columns if col != 'contig_id')
        
        missing_samples = original_samples - converted_samples
        extra_samples = converted_samples - original_samples
        
        assert not missing_samples, f"Missing samples in conversion: {missing_samples}"
        
        if extra_samples:
            print(f"⚠️  Extra samples in conversion: {extra_samples}")
        
        print(f"✅ Validation passed: {original_total:,} total reads preserved")
        print(f"✅ Dimensions correct: {n_contigs} contigs × {n_samples} samples")
    
    def convert_long_to_wide(self, abundance_df: pd.DataFrame) -> pd.DataFrame:
        """
        Main conversion method with format detection and validation
        
        Parameters:
        -----------
        abundance_df : pd.DataFrame
            Input abundance data (auto-detects format)
            
        Returns:
        --------
        pd.DataFrame
            Wide format abundance data
        """
        print("🔄 Starting format conversion...")
        
        # Detect current format
        current_format = self._detect_format(abundance_df)
        
        if current_format == "wide":
            print("✅ Data is already in wide format")
            return abundance_df
        
        # Convert from long to wide
        print("Converting from long to wide format...")
        converted_df = self._pivot_data(abundance_df)
        
        # Validate conversion
        if self.validation_mode:
            self._validate_conversion(abundance_df, converted_df)
        
        print("✅ Format conversion completed successfully")
        return converted_df
    
    def load_and_convert(self, file_path: Path) -> pd.DataFrame:
        """
        Load abundance data from file and convert to wide format
        
        Parameters:
        -----------
        file_path : Path
            Path to abundance data file (.tsv or .tsv.gz)
            
        Returns:
        --------
        pd.DataFrame
            Wide format abundance data
        """
        print(f"📁 Loading abundance data from: {file_path.name}")
        
        # Handle compressed files
        if file_path.suffix == '.gz':
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            df = pd.read_csv(file_path, sep='\t')
        
        print(f"Loaded data: {df.shape}")
        
        # Convert to wide format
        wide_df = self.convert_long_to_wide(df)
        
        return wide_df
    
    def save_wide_format(self, wide_df: pd.DataFrame, output_path: Path) -> None:
        """
        Save wide format data to file
        
        Parameters:
        -----------
        wide_df : pd.DataFrame
            Wide format abundance data
        output_path : Path
            Output file path
        """
        print(f"💾 Saving wide format data to: {output_path}")
        
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save as TSV
        wide_df.to_csv(output_path, sep='\t', index=False)
        
        print(f"✅ Saved {wide_df.shape} wide format data")

def test_with_uc_temporal_data(data_path: Optional[Path] = None) -> pd.DataFrame:
    """
    Test the converter with UC temporal data
    
    Parameters:
    -----------
    data_path : Optional[Path]
        Path to test data. If None, uses default test data location.
        
    Returns:
    --------
    pd.DataFrame
        Converted wide format data
    """
    if data_path is None:
        # Default test data location
        data_path = Path("/path/to/test_data/contig_count_table_UC_230_4.tsv.gz")
    
    if not data_path.exists():
        print(f"❌ Test data not found at {data_path}")
        print("Please provide path to UC temporal test data")
        return pd.DataFrame()
    
    # Initialize converter
    converter = AbundanceFormatConverter(validation_mode=True)
    
    # Load and convert
    try:
        wide_df = converter.load_and_convert(data_path)
        print(f"🎉 Test successful! Converted data shape: {wide_df.shape}")
        return wide_df
    except Exception as e:
        print(f"❌ Test failed: {e}")
        raise

if __name__ == "__main__":
    # Command line interface
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert abundance data formats")
    parser.add_argument("input_file", help="Input abundance data file")
    parser.add_argument("output_file", help="Output wide format file")
    parser.add_argument("--no-validation", action="store_true", 
                       help="Skip validation (faster but less safe)")
    
    args = parser.parse_args()
    
    # Convert the data
    converter = AbundanceFormatConverter(validation_mode=not args.no_validation)
    
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    
    wide_data = converter.load_and_convert(input_path)
    converter.save_wide_format(wide_data, output_path)
    
    print("🎉 Conversion completed successfully!")