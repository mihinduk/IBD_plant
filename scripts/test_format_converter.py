#!/usr/bin/env python3
"""
Test script for the format converter using real UC temporal data

This script tests the format converter on actual IBD dataset samples
to ensure it correctly handles the Stage 1 conversion issue.
"""

import sys
from pathlib import Path

# Add data_pipeline to path
sys.path.append(str(Path(__file__).parent.parent / 'data_pipeline' / 'stage1_data_prep'))
sys.path.append(str(Path(__file__).parent.parent / 'validation'))

from format_converter import AbundanceFormatConverter
from result_comparator import ResultValidator

def test_format_converter_basic():
    """Basic test of format converter functionality"""
    print("🧪 Testing Format Converter - Basic Functionality")
    print("=" * 50)
    
    converter = AbundanceFormatConverter(validation_mode=True)
    
    # Test with sample data structure
    import pandas as pd
    
    # Create sample long format data (mimics real UC temporal structure)
    sample_data = pd.DataFrame({
        'Sample': ['7562', '7562', '7583', '7583', '7584', '7584'],
        'Contig': ['contig_100018', 'contig_100067', 'contig_100018', 'contig_100067', 'contig_100018', 'contig_100067'],
        'Reads': [2, 5, 0, 3, 1, 0]
    })
    
    print("Input data (long format):")
    print(sample_data)
    print()
    
    # Convert to wide format
    try:
        wide_data = converter.convert_long_to_wide(sample_data)
        print("Output data (wide format):")
        print(wide_data)
        print()
        
        # Verify the conversion
        validator = ResultValidator()
        if validator.validate_format_conversion(sample_data, wide_data):
            print("✅ Basic format conversion test PASSED")
            return True
        else:
            print("❌ Basic format conversion test FAILED")
            return False
            
    except Exception as e:
        print(f"❌ Format conversion failed: {e}")
        return False

def test_with_real_data_path(data_path):
    """Test with actual UC temporal data file"""
    print(f"🧪 Testing Format Converter - Real Data")
    print(f"Data path: {data_path}")
    print("=" * 50)
    
    data_file = Path(data_path)
    if not data_file.exists():
        print(f"❌ Test data file not found: {data_file}")
        print("Please provide path to UC temporal data file")
        return False
    
    converter = AbundanceFormatConverter(validation_mode=True)
    
    try:
        # Load and convert real data
        wide_data = converter.load_and_convert(data_file)
        
        print(f"✅ Real data conversion PASSED")
        print(f"   Converted shape: {wide_data.shape}")
        print(f"   Sample columns: {len([col for col in wide_data.columns if col != 'contig_id'])}")
        print(f"   Contigs: {len(wide_data)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Real data conversion FAILED: {e}")
        return False

def main():
    """Main test function"""
    print("🚀 Format Converter Test Suite")
    print("Testing the Stage 1 format conversion fix")
    print("=" * 60)
    
    # Test 1: Basic functionality
    basic_test_passed = test_format_converter_basic()
    
    print()
    
    # Test 2: Real data (if path provided)
    if len(sys.argv) > 1:
        real_data_path = sys.argv[1]
        real_test_passed = test_with_real_data_path(real_data_path)
    else:
        print("ℹ️  To test with real data, provide file path as argument:")
        print("   python test_format_converter.py /path/to/contig_count_table.tsv.gz")
        real_test_passed = True  # Skip real data test
    
    # Summary
    print()
    print("=" * 60)
    print("🏁 Test Results Summary")
    print(f"   Basic functionality: {'✅ PASSED' if basic_test_passed else '❌ FAILED'}")
    
    if len(sys.argv) > 1:
        print(f"   Real data test: {'✅ PASSED' if real_test_passed else '❌ FAILED'}")
        all_passed = basic_test_passed and real_test_passed
    else:
        all_passed = basic_test_passed
    
    if all_passed:
        print("🎉 All tests PASSED - Format converter is ready!")
        return 0
    else:
        print("❌ Some tests FAILED - Check implementation")
        return 1

if __name__ == "__main__":
    sys.exit(main())