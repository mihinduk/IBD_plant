#!/usr/bin/env python3
"""
Pipeline Constants and Naming Conventions
=========================================

Centralized filename patterns and constants for the Staged Hierarchical 
Temporal Fingerprinting Pipeline. This ensures consistent naming across
all pipeline stages and enables scalability.

Author: Claude Code for Kathie Mihindukulasuriya
Date: 2025-07-02
"""

from pathlib import Path
from typing import List, Optional

class PipelineFiles:
    """Centralized file naming conventions for pipeline stages"""
    
    # Input file patterns (for discovery)
    INPUT_ABUNDANCE_PATTERNS = [
        "contig_count_table*.tsv.gz",
        "contig_count_table*.tsv", 
        "*contig*abundance*.tsv",
        "UC_temporal_contig_abundance_with_taxonomy.tsv"  # Legacy
    ]
    
    INPUT_METADATA_PATTERNS = [
        "KM_resources/21032025_Metadata_freeze4_imputed_fin.txt",
        "*metadata*.txt",
        "*metadata*.tsv",
        "filtered_metadata_uc_temporal.tsv"  # Legacy
    ]
    
    INPUT_TAXONOMY_PATTERNS = [
        "KM_resources/05_00_assembly_1KB_tophit_linage_with_header.tsv",
        "*taxonomy*.tsv",
        "*lineage*.tsv"
    ]
    
    # Stage output files - enhanced with clinical metadata
    STAGE1_FILTERED_ABUNDANCE = "filtered_abundance_{diagnosis}_{inflammation}.tsv"
    STAGE1_FILTERED_METADATA = "filtered_metadata_{diagnosis}_{inflammation}.tsv" 
    STAGE1_PATIENT_GROUPS = "patient_groups_{diagnosis}_{inflammation}.json"
    STAGE1_DIAGNOSTICS = "stage1_diagnostics_{job_id}.json"
    
    STAGE2_INDIVIDUAL_PATTERNS = "individual_patterns_{diagnosis}_{inflammation}_{job_id}.json"
    STAGE2_PATIENT_RESULTS = "patient_results_{diagnosis}_{inflammation}_{patient_id}.json"
    STAGE2_DIAGNOSTICS = "stage2_diagnostics_{job_id}.json"
    
    STAGE3_GROUP_RESULTS = "group_aggregation_{diagnosis}_{inflammation}_{job_id}.json"
    STAGE3_TEMPORAL_SIGNATURES = "temporal_signatures_{diagnosis}_{inflammation}_{job_id}.json" 
    STAGE3_DIAGNOSTICS = "stage3_diagnostics_{job_id}.json"
    
    STAGE4_CROSS_GROUP_ANALYSIS = "cross_group_analysis_{job_id}.json"
    STAGE4_FINAL_FINGERPRINTS = "final_fingerprints_{diagnosis}_{inflammation}_{job_id}.json"
    STAGE4_DIAGNOSTICS = "stage4_diagnostics_{job_id}.json"
    
    @staticmethod
    def discover_files(directory: Path, patterns: List[str]) -> Optional[Path]:
        """
        Discover files using patterns in order of preference
        
        Parameters:
        -----------
        directory : Path
            Directory to search in
        patterns : List[str]
            List of glob patterns, in order of preference
            
        Returns:
        --------
        Optional[Path]
            First matching file found, or None
        """
        for pattern in patterns:
            # Handle subdirectory patterns (e.g., "KM_resources/file.txt")
            if "/" in pattern:
                candidate = directory / pattern
                if candidate.exists():
                    return candidate
            else:
                # Handle glob patterns
                matches = list(directory.glob(pattern))
                if matches:
                    # Return the most recently modified file
                    return max(matches, key=lambda p: p.stat().st_mtime)
        
        return None
    
    @staticmethod
    def format_filename(template: str, **kwargs) -> str:
        """
        Format filename template with variables
        
        Parameters:
        -----------
        template : str
            Filename template with {variable} placeholders
        **kwargs
            Variables to substitute
            
        Returns:
        --------
        str
            Formatted filename
        """
        return template.format(**kwargs)

class PipelineDirectories:
    """Standard directory structure for pipeline stages"""
    
    STAGE1_DIR = "stage1_filtering"
    STAGE2_DIR = "stage2_individual" 
    STAGE3_DIR = "stage3_group_aggregation"
    STAGE4_DIR = "stage4_cross_group_analysis"
    
    # Results subdirectories
    DIAGNOSTICS_DIR = "diagnostics"
    FIGURES_DIR = "figures"
    INTERMEDIATE_DIR = "intermediate"

class PipelineDefaults:
    """Default parameters and thresholds"""
    
    # Filtering thresholds
    MIN_PREVALENCE = 0.05  # 5% of samples
    MIN_ABUNDANCE = 1.0    # log10 scale
    MIN_COMBINED_SCORE = 0.1
    
    # Temporal analysis
    MIN_SAMPLES_PER_PATIENT = 2
    MIN_TOTAL_READS = 10
    
    # Patient grouping
    SAMPLES_PER_GROUP = 4  # Default timepoints
    MAX_PARALLEL_JOBS = 4
    
    # Memory and performance
    CHUNK_SIZE = 5000  # For processing large datasets
    LOG_INTERVAL = 1000  # Progress logging interval

class MetadataColumns:
    """Standardized metadata column names"""
    
    # Primary identifiers
    SAMPLE_ID = "database_ID"
    PATIENT_ID = "patient_ID"
    
    # Clinical classification
    DIAGNOSIS = "diagnosis_assoc_fin"
    INFLAMMATION = "inflammation"
    
    # Other important columns
    TIMEPOINT = "timepoint"
    DISEASE_ACTIVITY = "disease_activity"
    
    # Required columns for pipeline
    REQUIRED_COLUMNS = [SAMPLE_ID, PATIENT_ID, DIAGNOSIS]
    RECOMMENDED_COLUMNS = [SAMPLE_ID, PATIENT_ID, DIAGNOSIS, INFLAMMATION]
    
    @staticmethod
    def validate_metadata(df, logger=None):
        """Validate that metadata has required columns"""
        missing_required = [col for col in MetadataColumns.REQUIRED_COLUMNS if col not in df.columns]
        missing_recommended = [col for col in MetadataColumns.RECOMMENDED_COLUMNS if col not in df.columns]
        
        if missing_required:
            raise ValueError(f"Missing required metadata columns: {missing_required}")
        
        if missing_recommended and logger:
            logger.warning(f"Missing recommended metadata columns: {missing_recommended}")
        
        return True

class AbundanceColumns:
    """Standardized abundance data column names"""
    
    # Long format (Sample, Contig, Reads)
    LONG_FORMAT_SAMPLE = "Sample"
    LONG_FORMAT_CONTIG = "Contig" 
    LONG_FORMAT_READS = "Reads"
    
    # Wide format (contig_id as first column, samples as other columns)
    WIDE_FORMAT_CONTIG = "contig_id"
    
    # Alternative names found in different data sources
    CONTIG_ALTERNATIVES = ["Contig", "contig_id", "contig_name", "id"]
    SAMPLE_ALTERNATIVES = ["Sample", "sample_id", "database_ID"]
    READS_ALTERNATIVES = ["Reads", "count", "abundance", "reads"]
    
    @staticmethod
    def detect_format(df, logger=None):
        """Detect whether data is in long or wide format"""
        # Check for long format indicators
        has_sample_col = any(col in df.columns for col in AbundanceColumns.SAMPLE_ALTERNATIVES)
        has_contig_col = any(col in df.columns for col in AbundanceColumns.CONTIG_ALTERNATIVES) 
        has_reads_col = any(col in df.columns for col in AbundanceColumns.READS_ALTERNATIVES)
        
        # Check for wide format indicators (many columns starting with 'S')
        sample_cols = [col for col in df.columns if col.startswith('S')]
        
        if has_sample_col and has_contig_col and has_reads_col:
            if logger:
                logger.info("Detected long format abundance data")
            return "long"
        elif len(sample_cols) > 5:  # Arbitrary threshold
            if logger:
                logger.info("Detected wide format abundance data") 
            return "wide"
        else:
            raise ValueError("Unable to determine abundance data format")
    
    @staticmethod
    def get_column_mapping(df):
        """Get actual column names for standard fields"""
        mapping = {}
        
        # Find contig column
        for alt in AbundanceColumns.CONTIG_ALTERNATIVES:
            if alt in df.columns:
                mapping['contig'] = alt
                break
        
        # Find sample column
        for alt in AbundanceColumns.SAMPLE_ALTERNATIVES:
            if alt in df.columns:
                mapping['sample'] = alt
                break
                
        # Find reads column
        for alt in AbundanceColumns.READS_ALTERNATIVES:
            if alt in df.columns:
                mapping['reads'] = alt
                break
        
        return mapping