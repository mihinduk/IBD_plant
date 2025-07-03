#!/bin/bash
#SBATCH --job-name=stage1_filtering
#SBATCH --output=/scratch/sahlab/kathie/contig_linkage/results/stage1_%j.out
#SBATCH --error=/scratch/sahlab/kathie/contig_linkage/results/stage1_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --account=sahlab

# Save ALL arguments before conda activation (CRITICAL!)
SAVED_INPUT_DIR="$1"
SAVED_OUTPUT_DIR="$2"

# Clear positional parameters to prevent conda from seeing them (CRITICAL!)
set --
source /ref/sahlab/software/anaconda3/bin/activate

# Restore our variables AFTER conda activation
INPUT_DIR="$SAVED_INPUT_DIR"
OUTPUT_DIR="$SAVED_OUTPUT_DIR"

# Export as environment variables for Python script
export INPUT_DIR
export OUTPUT_DIR

# Set PYTHONUSERBASE to avoid conflicts
export PYTHONUSERBASE="/scratch/sahlab/kathie/python_user_${SLURM_JOB_ID}"

# Install required packages
python -m pip install --user pandas numpy scipy psutil openpyxl --quiet

# Change to scripts directory to ensure pipeline_constants can be imported
cd /scratch/sahlab/kathie/contig_linkage/scripts

# Run the Python script
python3 << 'EOF'
"""
Stage 1: Filtering & Preprocessing with Diagnostics
Part of the Staged Hierarchical Temporal Fingerprinting Pipeline

This stage handles:
- Data loading and validation
- Tier 1: Smart filtering (diagnosis-specific prevalence)
- Tier 2: Temporal correlation pre-screening
- Tier 3: Smart patient grouping
- Tier 4: Representative sampling
- Comprehensive diagnostics and timing

Input: Raw contig abundance data
Output: Filtered data ready for individual patient analysis
"""

import pandas as pd
import numpy as np
import json
import os
import sys
import time
import psutil
import gzip
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from pipeline_constants import PipelineFiles, PipelineDirectories, PipelineDefaults, MetadataColumns, AbundanceColumns

class DiagnosticTimer:
    """Context manager for timing operations with memory tracking"""
    
    def __init__(self, operation_name: str, logger: logging.Logger):
        self.operation_name = operation_name
        self.logger = logger
        self.start_time = None
        self.start_memory = None
        
    def __enter__(self):
        self.start_time = time.time()
        self.start_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
        self.logger.info(f"🚀 Starting: {self.operation_name}")
        self.logger.info(f"   Memory at start: {self.start_memory:.1f} MB")
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        end_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
        duration = end_time - self.start_time
        memory_delta = end_memory - self.start_memory
        
        self.logger.info(f"✅ Completed: {self.operation_name}")
        self.logger.info(f"   Duration: {duration:.2f} seconds ({duration/60:.2f} minutes)")
        self.logger.info(f"   Memory delta: {memory_delta:+.1f} MB (now {end_memory:.1f} MB)")
        
        if exc_type is not None:
            self.logger.error(f"❌ Failed: {self.operation_name} - {exc_val}")

class Stage1FilteringPreprocessing:
    def __init__(self, input_dir: str, output_dir: str, job_id: str):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.job_id = job_id
        
        # Create output structure
        self.stage1_dir = self.output_dir / "stage1_filtering"
        self.stage1_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        
        # Diagnostic data
        self.diagnostics = {
            "job_id": job_id,
            "start_time": datetime.now().isoformat(),
            "stage": "stage1_filtering_preprocessing",
            "operations": {}
        }
        
    def setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.stage1_dir / f"stage1_diagnostics_{self.job_id}.log"
        
        # Create logger
        self.logger = logging.getLogger('stage1')
        self.logger.setLevel(logging.INFO)
        
        # Create formatters
        detailed_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(detailed_formatter)
        self.logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(detailed_formatter)
        self.logger.addHandler(console_handler)
        
    def discover_input_files(self) -> Tuple[Path, Optional[Path]]:
        """Discover input files using established naming conventions"""
        # Discover abundance file
        abundance_file = PipelineFiles.discover_files(
            self.input_dir, 
            PipelineFiles.INPUT_ABUNDANCE_PATTERNS
        )
        
        if not abundance_file:
            raise FileNotFoundError(
                f"No abundance file found in {self.input_dir}. "
                f"Expected patterns: {PipelineFiles.INPUT_ABUNDANCE_PATTERNS}"
            )
        
        # Discover metadata file
        metadata_file = PipelineFiles.discover_files(
            self.input_dir,
            PipelineFiles.INPUT_METADATA_PATTERNS
        )
        
        if not metadata_file:
            self.logger.warning(
                f"No metadata file found in {self.input_dir}. "
                f"Will create minimal metadata from sample columns."
            )
        
        return abundance_file, metadata_file
    
    def load_and_validate_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Load contig abundance and metadata with validation"""
        with DiagnosticTimer("Data Loading and Validation", self.logger):
            # Discover input files
            abundance_file, metadata_file = self.discover_input_files()
            
            # Load abundance data (handle compression)
            self.logger.info(f"Loading abundance data from: {abundance_file.name}")
            if abundance_file.suffix == '.gz':
                import gzip
                with gzip.open(abundance_file, 'rt') as f:
                    abundance_df = pd.read_csv(f, sep='\t')
            else:
                abundance_df = pd.read_csv(abundance_file, sep='\t')
            self.logger.info(f"Loaded abundance data: {abundance_df.shape}")
            
            # Detect and convert data format
            data_format = AbundanceColumns.detect_format(abundance_df, self.logger)
            
            if data_format == "long":
                self.logger.info("Converting long format to wide format")
                # Get column mapping
                col_mapping = AbundanceColumns.get_column_mapping(abundance_df)
                
                # Pivot the data
                abundance_df = abundance_df.pivot(
                    index=col_mapping['contig'], 
                    columns=col_mapping['sample'], 
                    values=col_mapping['reads']
                ).fillna(0)
                
                # Reset index to make contig_id a column
                abundance_df = abundance_df.reset_index()
                abundance_df.columns.name = None  # Remove the columns name
                
                # Rename the contig column to contig_id
                abundance_df = abundance_df.rename(columns={col_mapping['contig']: 'contig_id'})
                
                self.logger.info(f"Converted to wide format: {abundance_df.shape}")
            
            # Load or create metadata
            if metadata_file:
                self.logger.info(f"Loading metadata from: {metadata_file.name}")
                # Handle different separators
                try:
                    metadata_df = pd.read_csv(metadata_file, sep='\t')
                except:
                    # Try comma separator for .txt files
                    metadata_df = pd.read_csv(metadata_file, sep=',')
                self.logger.info(f"Loaded metadata: {metadata_df.shape}")
            else:
                # Create minimal metadata from sample columns
                self.logger.info("Creating minimal metadata from sample columns")
                sample_cols = [col for col in abundance_df.columns if col != 'contig_id']
                metadata_df = pd.DataFrame({
                    MetadataColumns.SAMPLE_ID: sample_cols,
                    'database_ID': sample_cols,  # Use sample name as database ID
                    MetadataColumns.DIAGNOSIS: 'UC',  # Default diagnosis
                    'patient_ID': [f"patient_{i//4}" for i in range(len(sample_cols))],  # Group by 4
                    'timepoint': [i % 4 for i in range(len(sample_cols))]  # 0-3 timepoints
                })
                self.logger.info(f"Created minimal metadata for {len(sample_cols)} samples")
            
            # Validation
            sample_cols = [col for col in abundance_df.columns if col != 'contig_id']
            metadata_samples = set(metadata_df[MetadataColumns.SAMPLE_ID])
            abundance_samples = set(sample_cols)
            
            overlap = abundance_samples.intersection(metadata_samples)
            self.logger.info(f"Sample overlap: {len(overlap)} / {len(abundance_samples)} abundance samples")
            
            if len(overlap) < len(abundance_samples) * 0.8:
                self.logger.warning("Low sample overlap between abundance and metadata!")
                
            # Store diagnostics
            self.diagnostics["operations"]["data_loading"] = {
                "abundance_shape": abundance_df.shape,
                "metadata_shape": metadata_df.shape,
                "sample_overlap": len(overlap),
                "abundance_samples": len(abundance_samples),
                "metadata_samples": len(metadata_samples)
            }
            
            return abundance_df, metadata_df
    
    def tier1_smart_filtering(self, abundance_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
        """Tier 1: Smart filtering with diagnosis-specific prevalence"""
        with DiagnosticTimer("Tier 1: Smart Filtering", self.logger):
            # Get sample columns
            sample_cols = [col for col in abundance_df.columns if col != 'contig_id']
            
            # Create diagnosis mapping
            diagnosis_map = dict(zip(metadata_df[MetadataColumns.SAMPLE_ID], metadata_df[MetadataColumns.DIAGNOSIS]))
            
            # Separate samples by diagnosis
            uc_samples = [col for col in sample_cols if diagnosis_map.get(col) == 'UC']
            
            self.logger.info(f"UC samples for filtering: {len(uc_samples)}")
            
            # Calculate prevalence and abundance metrics
            results = []
            total_contigs = len(abundance_df)
            
            for idx, row in abundance_df.iterrows():
                if idx % 5000 == 0:
                    self.logger.info(f"Processing contig {idx}/{total_contigs}")
                
                contig_id = row['contig_id']
                uc_values = [row[col] for col in uc_samples if col in row.index]
                
                if not uc_values:
                    continue
                    
                # Prevalence (proportion of samples with reads)
                uc_prevalence = sum(1 for val in uc_values if val > 0) / len(uc_values)
                
                # Mean abundance (log scale)
                uc_abundance = np.mean([np.log10(val + 1) for val in uc_values])
                
                # Combined score (prevalence * abundance)
                combined_score = uc_prevalence * uc_abundance
                
                # Smart filtering criteria
                keep_contig = (
                    uc_prevalence >= 0.05 or  # Present in ≥5% of UC samples
                    uc_abundance >= 1.0 or    # High abundance
                    combined_score >= 0.1     # Good combined score
                )
                
                if keep_contig:
                    results.append(idx)
            
            # Filter dataframe
            filtered_df = abundance_df.iloc[results].copy()
            
            self.logger.info(f"Tier 1 filtering: {len(abundance_df)} → {len(filtered_df)} contigs")
            
            # Store diagnostics
            reduction_factor = len(abundance_df) / len(filtered_df) if len(filtered_df) > 0 else float('inf')
            self.diagnostics["operations"]["tier1_filtering"] = {
                "original_contigs": len(abundance_df),
                "filtered_contigs": len(filtered_df),
                "reduction_factor": reduction_factor,
                "uc_samples": len(uc_samples)
            }
            
            return filtered_df
    
    def tier2_temporal_prescreening(self, abundance_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
        """Tier 2: Temporal correlation pre-screening with distance weighting"""
        with DiagnosticTimer("Tier 2: Temporal Pre-screening", self.logger):
            # Get sample columns and create patient-timepoint mapping
            sample_cols = [col for col in abundance_df.columns if col.startswith('S')]
            
            # Create patient-timepoint structure
            patient_timepoints = {}
            for _, row in metadata_df.iterrows():
                sample = row[MetadataColumns.SAMPLE_ID]
                if sample in sample_cols:
                    patient_id = row['patient_ID']
                    timepoint = row['timepoint_numeric']
                    
                    if patient_id not in patient_timepoints:
                        patient_timepoints[patient_id] = []
                    patient_timepoints[patient_id].append((sample, timepoint))
            
            # Filter to patients with multiple timepoints
            multi_timepoint_patients = {
                pid: tps for pid, tps in patient_timepoints.items() 
                if len(tps) >= 2
            }
            
            self.logger.info(f"Patients with multiple timepoints: {len(multi_timepoint_patients)}")
            
            # Quick temporal screening
            temporal_scores = []
            total_contigs = len(abundance_df)
            
            for idx, row in abundance_df.iterrows():
                if idx % 2000 == 0:
                    self.logger.info(f"Temporal screening: {idx}/{total_contigs}")
                
                # Calculate temporal variation for each patient
                patient_variations = []
                
                for patient_id, timepoints in multi_timepoint_patients.items():
                    if len(timepoints) < 2:
                        continue
                        
                    # Get values for this patient's timepoints
                    values = []
                    times = []
                    for sample, timepoint in timepoints:
                        if sample in row.index:
                            values.append(row[sample])
                            times.append(timepoint)
                    
                    if len(values) >= 2:
                        # Calculate coefficient of variation
                        if np.mean(values) > 0:
                            cv = np.std(values) / np.mean(values)
                            patient_variations.append(cv)
                
                # Overall temporal score
                if patient_variations:
                    temporal_score = np.mean(patient_variations)
                else:
                    temporal_score = 0
                
                temporal_scores.append(temporal_score)
            
            # Filter based on temporal variation
            threshold = np.percentile(temporal_scores, 25)  # Keep top 75% most variable
            keep_indices = [
                idx for idx, score in enumerate(temporal_scores) 
                if score >= threshold
            ]
            
            filtered_df = abundance_df.iloc[keep_indices].copy()
            
            self.logger.info(f"Tier 2 temporal filtering: {len(abundance_df)} → {len(filtered_df)} contigs")
            
            # Store diagnostics
            self.diagnostics["operations"]["tier2_temporal"] = {
                "input_contigs": len(abundance_df),
                "output_contigs": len(filtered_df),
                "multi_timepoint_patients": len(multi_timepoint_patients),
                "temporal_threshold": threshold,
                "mean_temporal_score": np.mean(temporal_scores)
            }
            
            return filtered_df
    
    def tier3_smart_grouping(self, metadata_df: pd.DataFrame) -> Dict:
        """Tier 3: Smart patient grouping"""
        with DiagnosticTimer("Tier 3: Smart Patient Grouping", self.logger):
            # Create patient-level summary
            patient_summary = metadata_df.groupby('patient_ID').agg({
                MetadataColumns.DIAGNOSIS: 'first',
                'disease_activity': lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else x.iloc[0],
                'timepoint_numeric': ['min', 'max', 'count'],
                MetadataColumns.SAMPLE_ID: list
            }).reset_index()
            
            # Flatten column names
            patient_summary.columns = [
                'patient_ID', MetadataColumns.DIAGNOSIS, 'disease_activity', 
                'timepoint_min', 'timepoint_max', 'timepoint_count', 'samples'
            ]
            
            # Create smart groups
            groups = {}
            
            # UC patients by activity
            uc_patients = patient_summary[patient_summary[MetadataColumns.DIAGNOSIS] == 'UC']
            
            groups['UC_active'] = uc_patients[
                uc_patients['disease_activity'] == 'active'
            ]['patient_ID'].tolist()
            
            groups['UC_inactive'] = uc_patients[
                uc_patients['disease_activity'] == 'inactive'
            ]['patient_ID'].tolist()
            
            groups['UC_mixed'] = uc_patients[
                ~uc_patients['disease_activity'].isin(['active', 'inactive'])
            ]['patient_ID'].tolist()
            
            # High temporal resolution patients
            groups['high_temporal'] = patient_summary[
                patient_summary['timepoint_count'] >= 4
            ]['patient_ID'].tolist()
            
            # Transition patients (those with changing activity)
            transition_patients = []
            for patient_id in patient_summary['patient_ID']:
                patient_data = metadata_df[metadata_df['patient_ID'] == patient_id]
                activities = patient_data['disease_activity'].unique()
                if len(activities) > 1:
                    transition_patients.append(patient_id)
            
            groups['transition'] = transition_patients
            
            self.logger.info(f"Created {len(groups)} patient groups:")
            for group_name, patients in groups.items():
                self.logger.info(f"  {group_name}: {len(patients)} patients")
            
            # Store diagnostics
            self.diagnostics["operations"]["tier3_grouping"] = {
                "total_patients": len(patient_summary),
                "groups": {name: len(patients) for name, patients in groups.items()},
                "uc_patients": len(uc_patients)
            }
            
            return groups
    
    def tier4_representative_sampling(self, abundance_df: pd.DataFrame, metadata_df: pd.DataFrame, 
                                    patient_groups: Dict) -> Tuple[pd.DataFrame, Dict]:
        """Tier 4: Representative sampling for computational efficiency"""
        with DiagnosticTimer("Tier 4: Representative Sampling", self.logger):
            from collections import defaultdict
            
            # Score patients for representativeness
            patient_scores = defaultdict(dict)
            sample_cols = [col for col in abundance_df.columns if col.startswith('S')]
            
            # Calculate scores for each patient
            for patient_id in metadata_df['patient_ID'].unique():
                patient_samples = metadata_df[
                    metadata_df['patient_ID'] == patient_id
                ][MetadataColumns.SAMPLE_ID].tolist()
                
                # Filter to available samples
                available_samples = [s for s in patient_samples if s in sample_cols]
                
                if not available_samples:
                    continue
                
                # Temporal diversity score
                timepoints = metadata_df[
                    metadata_df[MetadataColumns.SAMPLE_ID].isin(available_samples)
                ]['timepoint_numeric'].values
                temporal_span = np.max(timepoints) - np.min(timepoints) if len(timepoints) > 1 else 0
                temporal_density = len(timepoints) / (temporal_span + 1)
                
                # Viral diversity score
                patient_abundance = abundance_df[available_samples].sum(axis=1)
                viral_richness = (patient_abundance > 0).sum()
                viral_load = np.log10(patient_abundance.sum() + 1)
                
                # Transition priority (patients with changing disease activity)
                is_transition = patient_id in patient_groups.get('transition', [])
                
                # Combined score
                combined_score = (
                    temporal_density * 2 +
                    viral_richness * 0.01 +
                    viral_load * 0.5 +
                    (10 if is_transition else 0)
                )
                
                patient_scores[patient_id] = {
                    'temporal_span': temporal_span,
                    'temporal_density': temporal_density,
                    'viral_richness': viral_richness,
                    'viral_load': viral_load,
                    'is_transition': is_transition,
                    'combined_score': combined_score,
                    'samples': available_samples
                }
            
            # Select representative patients from each group
            selected_patients = set()
            selection_summary = {}
            
            for group_name, group_patients in patient_groups.items():
                # Get scores for this group
                group_scores = {
                    pid: patient_scores[pid] for pid in group_patients 
                    if pid in patient_scores
                }
                
                if not group_scores:
                    selection_summary[group_name] = {'selected': 0, 'available': 0}
                    continue
                
                # Select top patients (aim for ~50% reduction)
                n_select = max(1, len(group_scores) // 2)
                
                # Prioritize transition patients and high scorers
                sorted_patients = sorted(
                    group_scores.items(), 
                    key=lambda x: x[1]['combined_score'], 
                    reverse=True
                )
                
                selected_for_group = [pid for pid, _ in sorted_patients[:n_select]]
                selected_patients.update(selected_for_group)
                
                selection_summary[group_name] = {
                    'selected': len(selected_for_group),
                    'available': len(group_scores)
                }
                
                self.logger.info(f"Group {group_name}: selected {len(selected_for_group)}/{len(group_scores)} patients")
            
            # Create filtered metadata for selected patients
            selected_metadata = metadata_df[
                metadata_df['patient_ID'].isin(selected_patients)
            ].copy()
            
            # Filter abundance data to selected samples
            selected_samples = selected_metadata[MetadataColumns.SAMPLE_ID].tolist()
            abundance_cols = [col for col in abundance_df.columns if not col.startswith('S')]
            sample_cols = [col for col in abundance_df.columns if col.startswith('S') and col in selected_samples]
            
            filtered_abundance = abundance_df[abundance_cols + sample_cols].copy()
            
            self.logger.info(f"Tier 4 sampling: {len(metadata_df['patient_ID'].unique())} → {len(selected_patients)} patients")
            self.logger.info(f"Samples: {len([c for c in abundance_df.columns if c.startswith('S')])} → {len(sample_cols)}")
            
            # Store diagnostics
            self.diagnostics["operations"]["tier4_sampling"] = {
                "original_patients": len(metadata_df['patient_ID'].unique()),
                "selected_patients": len(selected_patients),
                "original_samples": len([c for c in abundance_df.columns if c.startswith('S')]),
                "selected_samples": len(sample_cols),
                "reduction_factor": len(metadata_df['patient_ID'].unique()) / len(selected_patients),
                "group_selection": selection_summary
            }
            
            return filtered_abundance, selected_metadata
    
    def save_outputs(self, filtered_abundance: pd.DataFrame, filtered_metadata: pd.DataFrame, 
                    patient_groups: Dict):
        """Save all outputs for next stages"""
        with DiagnosticTimer("Saving Outputs", self.logger):
            # Save filtered abundance data
            abundance_file = self.stage1_dir / PipelineFiles.STAGE1_FILTERED_ABUNDANCE
            filtered_abundance.to_csv(abundance_file, sep='\t', index=False)
            self.logger.info(f"Saved filtered abundance: {abundance_file}")
            
            # Save filtered metadata
            metadata_file = self.stage1_dir / PipelineFiles.STAGE1_FILTERED_METADATA
            filtered_metadata.to_csv(metadata_file, sep='\t', index=False)
            self.logger.info(f"Saved filtered metadata: {metadata_file}")
            
            # Save patient groups
            groups_file = self.stage1_dir / PipelineFiles.STAGE1_PATIENT_GROUPS
            with open(groups_file, 'w') as f:
                json.dump(patient_groups, f, indent=2)
            self.logger.info(f"Saved patient groups: {groups_file}")
            
            # Save diagnostics
            self.diagnostics["end_time"] = datetime.now().isoformat()
            self.diagnostics["total_duration_minutes"] = (
                datetime.fromisoformat(self.diagnostics["end_time"]) - 
                datetime.fromisoformat(self.diagnostics["start_time"])
            ).total_seconds() / 60
            
            diagnostics_file = self.stage1_dir / PipelineFiles.format_filename(
                PipelineFiles.STAGE1_DIAGNOSTICS, job_id=self.job_id
            )
            with open(diagnostics_file, 'w') as f:
                json.dump(self.diagnostics, f, indent=2)
            self.logger.info(f"Saved diagnostics: {diagnostics_file}")
    
    def run(self):
        """Execute complete Stage 1 pipeline"""
        self.logger.info("🚀 Starting Stage 1: Filtering & Preprocessing")
        self.logger.info(f"Job ID: {self.job_id}")
        self.logger.info(f"Input directory: {self.input_dir}")
        self.logger.info(f"Output directory: {self.output_dir}")
        
        try:
            # Load data
            abundance_df, metadata_df = self.load_and_validate_data()
            
            # Tier 1: Smart filtering
            abundance_df = self.tier1_smart_filtering(abundance_df, metadata_df)
            
            # Tier 2: Temporal pre-screening
            abundance_df = self.tier2_temporal_prescreening(abundance_df, metadata_df)
            
            # Tier 3: Smart patient grouping
            patient_groups = self.tier3_smart_grouping(metadata_df)
            
            # Tier 4: Representative sampling
            abundance_df, metadata_df = self.tier4_representative_sampling(
                abundance_df, metadata_df, patient_groups
            )
            
            # Save outputs
            self.save_outputs(abundance_df, metadata_df, patient_groups)
            
            self.logger.info("✅ Stage 1 completed successfully!")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ Stage 1 failed: {str(e)}")
            self.diagnostics["error"] = str(e)
            self.diagnostics["end_time"] = datetime.now().isoformat()
            
            # Save diagnostics even on failure
            diagnostics_file = self.stage1_dir / PipelineFiles.format_filename(
                PipelineFiles.STAGE1_DIAGNOSTICS.replace('.json', '_FAILED.json'),
                job_id=self.job_id
            )
            with open(diagnostics_file, 'w') as f:
                json.dump(self.diagnostics, f, indent=2)
            
            raise

def main():
    import os
    
    # Get arguments from environment variables (set by bash wrapper)
    input_dir = os.environ.get('INPUT_DIR')
    output_dir = os.environ.get('OUTPUT_DIR')
    job_id = os.environ.get('SLURM_JOB_ID', 'manual')
    
    if not input_dir or not output_dir:
        print("Error: INPUT_DIR and OUTPUT_DIR environment variables must be set")
        sys.exit(1)
    
    # Create and run Stage 1
    stage1 = Stage1FilteringPreprocessing(input_dir, output_dir, job_id)
    stage1.run()

if __name__ == "__main__":
    main()
EOF