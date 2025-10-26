#!/usr/bin/env python3
"""
Example script demonstrating gtf2utr usage
"""

import os
import sys
from pathlib import Path

# Add the parent directory to the path to import gtf2utr
sys.path.insert(0, str(Path(__file__).parent.parent))

from gtf2utr import GTFProcessor, UTRExtractor

def main():
    """Run example pipeline"""
    
    # Get the directory containing this script
    example_dir = Path(__file__).parent
    data_dir = example_dir / "data"
    
    # Input files
    gtf_file = data_dir / "sample.gtf"
    fasta_file = data_dir / "sample.fa"
    
    # Output files
    processed_gtf = example_dir / "sample_processed.gtf"
    utr_sequences = example_dir / "sample_utr_sequences.fa"
    
    print("Running gtf2utr example...")
    print(f"Input GTF: {gtf_file}")
    print(f"Input FASTA: {fasta_file}")
    print(f"Output processed GTF: {processed_gtf}")
    print(f"Output UTR sequences: {utr_sequences}")
    
    # Step 1: Process GTF file to classify UTRs
    print("\nStep 1: Processing GTF file...")
    processor = GTFProcessor(str(gtf_file), str(processed_gtf))
    processor.load_gtf()
    processor.classify_utrs()
    processor.write_output_gtf()
    processor.print_statistics()
    
    # Step 2: Extract UTR sequences
    print("\nStep 2: Extracting UTR sequences...")
    extractor = UTRExtractor()
    extractor.load_fasta(str(fasta_file))
    extractor.load_gtf(str(processed_gtf))
    extractor.extract_all_utrs(str(utr_sequences))
    
    print(f"\nExample completed successfully!")
    print(f"Check the output files:")
    print(f"  - {processed_gtf}")
    print(f"  - {utr_sequences}")

if __name__ == "__main__":
    main()