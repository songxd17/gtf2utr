#!/usr/bin/env python3
"""
Command Line Interface Module
Provides command line tools for gtf2utr package
"""

import argparse
import sys
import os
from .gtf_processor import GTFProcessor
from .utr_extractor import UTRExtractor


def process_gtf():
    """Command line interface for processing GTF files"""
    parser = argparse.ArgumentParser(
        description='GTF UTR Processor - Extract and re-classify UTRs',
        prog='gtf2utr-process'
    )
    parser.add_argument('input_gtf', help='Input GTF file path (supports .gz compressed files)')
    parser.add_argument('output_gtf', help='Output GTF file path')
    
    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.input_gtf):
        print(f"Error: Input file does not exist: {args.input_gtf}")
        sys.exit(1)
    
    # Create processor and run
    processor = GTFProcessor(args.input_gtf, args.output_gtf)
    processor.process()


def extract_utrs():
    """Command line interface for extracting UTR sequences"""
    parser = argparse.ArgumentParser(
        description='Extract UTR sequences from reference genome',
        prog='gtf2utr-extract'
    )
    parser.add_argument('gtf_file', help='Input GTF file (with classified UTRs)')
    parser.add_argument('fasta_file', help='Reference genome FASTA file')
    parser.add_argument('output_file', help='Output FASTA file')
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.gtf_file):
        print(f"Error: GTF file does not exist: {args.gtf_file}")
        sys.exit(1)
    
    if not os.path.exists(args.fasta_file):
        print(f"Error: FASTA file does not exist: {args.fasta_file}")
        sys.exit(1)
    
    # Create extractor and run
    extractor = UTRExtractor()
    extractor.process(args.gtf_file, args.fasta_file, args.output_file)


def main():
    """Main command line interface"""
    parser = argparse.ArgumentParser(
        description='gtf2utr - GTF to UTR sequence extraction tool',
        prog='gtf2utr'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Process GTF subcommand
    process_parser = subparsers.add_parser(
        'process', 
        help='Process GTF file and classify UTR regions'
    )
    process_parser.add_argument('input_gtf', help='Input GTF file path (supports .gz compressed files)')
    process_parser.add_argument('output_gtf', help='Output GTF file path')
    
    # Extract UTR subcommand
    extract_parser = subparsers.add_parser(
        'extract', 
        help='Extract UTR sequences from reference genome'
    )
    extract_parser.add_argument('gtf_file', help='Input GTF file (with classified UTRs)')
    extract_parser.add_argument('fasta_file', help='Reference genome FASTA file')
    extract_parser.add_argument('output_file', help='Output FASTA file')
    
    # Complete pipeline subcommand
    pipeline_parser = subparsers.add_parser(
        'pipeline', 
        help='Run complete GTF processing and UTR extraction pipeline'
    )
    pipeline_parser.add_argument('input_gtf', help='Input GTF file path (supports .gz compressed files)')
    pipeline_parser.add_argument('fasta_file', help='Reference genome FASTA file')
    pipeline_parser.add_argument('output_fasta', help='Output UTR sequences FASTA file')
    pipeline_parser.add_argument('--temp-gtf', help='Temporary GTF file path (optional)', 
                                default='temp_processed.gtf')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.command == 'process':
            # Check input file
            if not os.path.exists(args.input_gtf):
                print(f"Error: Input file does not exist: {args.input_gtf}")
                sys.exit(1)
            
            # Create processor and run
            processor = GTFProcessor(args.input_gtf, args.output_gtf)
            processor.process()
            
        elif args.command == 'extract':
            # Check input files
            if not os.path.exists(args.gtf_file):
                print(f"Error: GTF file does not exist: {args.gtf_file}")
                sys.exit(1)
            
            if not os.path.exists(args.fasta_file):
                print(f"Error: FASTA file does not exist: {args.fasta_file}")
                sys.exit(1)
            
            # Create extractor and run
            extractor = UTRExtractor()
            extractor.process(args.gtf_file, args.fasta_file, args.output_file)
            
        elif args.command == 'pipeline':
            # Check input files
            if not os.path.exists(args.input_gtf):
                print(f"Error: Input GTF file does not exist: {args.input_gtf}")
                sys.exit(1)
            
            if not os.path.exists(args.fasta_file):
                print(f"Error: FASTA file does not exist: {args.fasta_file}")
                sys.exit(1)
            
            print("=== Running complete GTF2UTR pipeline ===")
            
            # Step 1: Process GTF file
            print("\nStep 1: Processing GTF file...")
            processor = GTFProcessor(args.input_gtf, args.temp_gtf)
            processor.process()
            
            # Step 2: Extract UTR sequences
            print("\nStep 2: Extracting UTR sequences...")
            extractor = UTRExtractor()
            extractor.process(args.temp_gtf, args.fasta_file, args.output_fasta)
            
            # Clean up temporary files (optional)
            if args.temp_gtf == 'temp_processed.gtf' and os.path.exists(args.temp_gtf):
                print(f"\nCleaning up temporary file: {args.temp_gtf}")
                os.remove(args.temp_gtf)
            
            print("\n=== GTF2UTR pipeline completed ===")
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()