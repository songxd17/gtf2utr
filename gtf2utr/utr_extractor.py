#!/usr/bin/env python3
"""
UTR Sequence Extractor Module

Extract UTR sequences from reference genome FASTA files based on UTR annotations in processed GTF files.

Features:
1. Load reference genome FASTA files (supports gzip compression)
2. Parse UTR annotations from GTF files
3. Extract 5'UTR and 3'UTR sequences
4. Handle strand direction and sequence concatenation
5. Generate FASTA output with detailed metadata
"""

import gzip
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


class UTRExtractor:
    """UTR sequence extractor for extracting UTR sequences from reference genome"""
    
    def __init__(self):
        """Initialize UTR extractor"""
        self.chromosomes = {}
        self.utr_data = defaultdict(lambda: {'5utr': [], '3utr': []})
        self.stats = {
            'transcripts_processed': 0,
            'five_utr_regions': 0,
            'three_utr_regions': 0,
            'sequences_extracted': 0
        }
    
    def parse_gtf_attributes(self, attributes_str: str) -> Dict[str, str]:
        """Parse GTF attributes string into dictionary"""
        attributes = {}
        for attr in attributes_str.strip().split(';'):
            attr = attr.strip()
            if attr and ' ' in attr:
                key, value = attr.split(' ', 1)
                # Remove quotes from values
                value = value.strip('"')
                attributes[key] = value
        return attributes
    
    def load_fasta(self, fasta_file: str) -> None:
        """Load reference genome FASTA file"""
        print(f"Loading FASTA file: {fasta_file}")
        
        open_func = gzip.open if fasta_file.endswith('.gz') else open
        
        current_chr = None
        sequence_lines = []
        
        with open_func(fasta_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous chromosome
                    if current_chr and sequence_lines:
                        self.chromosomes[current_chr] = ''.join(sequence_lines)
                        print(f"Loaded chromosome: {current_chr} ({len(self.chromosomes[current_chr])} bp)")
                    
                    # Start new chromosome
                    current_chr = line[1:].split()[0]  # Take first part of header
                    sequence_lines = []
                else:
                    sequence_lines.append(line)
            
            # Save last chromosome
            if current_chr and sequence_lines:
                self.chromosomes[current_chr] = ''.join(sequence_lines)
                print(f"Loaded chromosome: {current_chr} ({len(self.chromosomes[current_chr])} bp)")
        
        print(f"Total chromosomes loaded: {len(self.chromosomes)}")
    
    def load_gtf(self, gtf_file: str) -> None:
        """Load GTF file and extract UTR information"""
        print(f"Loading GTF file: {gtf_file}")
        
        open_func = gzip.open if gtf_file.endswith('.gz') else open
        
        with open_func(gtf_file, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 100000 == 0:
                    print(f"Processed {line_num} lines...")
                
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) != 9:
                    continue
                
                seqname, source, feature, start, end, score, strand, frame, attributes = fields
                
                # Only process UTR features
                if feature not in ['five_prime_utr', 'three_prime_utr']:
                    continue
                
                # Parse attributes
                attr_dict = self.parse_gtf_attributes(attributes)
                transcript_id = attr_dict.get('transcript_id')
                gene_id = attr_dict.get('gene_id', '')
                gene_name = attr_dict.get('gene_name', '')
                
                if not transcript_id:
                    continue
                
                # Store UTR information and additional metadata
                utr_info = {
                    'chr': seqname,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'gene_id': gene_id,
                    'gene_name': gene_name
                }
                
                if feature == 'five_prime_utr':
                    self.utr_data[transcript_id]['5utr'].append(utr_info)
                    self.stats['five_utr_regions'] += 1
                elif feature == 'three_prime_utr':
                    self.utr_data[transcript_id]['3utr'].append(utr_info)
                    self.stats['three_utr_regions'] += 1
        
        print(f"Loaded UTR data for {len(self.utr_data)} transcripts")
        print(f"5'UTR regions: {self.stats['five_utr_regions']}")
        print(f"3'UTR regions: {self.stats['three_utr_regions']}")
    
    def reverse_complement(self, sequence: str) -> str:
        """Return reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base.upper(), 'N') for base in reversed(sequence))
    
    def extract_sequence(self, chr_name: str, start: int, end: int, strand: str) -> str:
        """Extract sequence from chromosome"""
        if chr_name not in self.chromosomes:
            return ""
        
        # Convert to 0-based indexing
        start_idx = start - 1
        end_idx = end
        
        # Check boundaries
        if start_idx < 0 or end_idx > len(self.chromosomes[chr_name]):
            return ""
        
        sequence = self.chromosomes[chr_name][start_idx:end_idx]
        
        # Negative strand requires reverse complement
        if strand == '-':
            sequence = self.reverse_complement(sequence)
        
        return sequence
    
    def concatenate_utr_sequences(self, utr_regions: List[Dict], transcript_id: str) -> Tuple[str, Dict]:
        """Concatenate UTR sequences for transcript and return sequence and metadata"""
        if not utr_regions:
            return "", {}
        
        # Sort by genomic position
        strand = utr_regions[0]['strand']
        
        if strand == '+':
            # Positive strand, sort by start position
            sorted_regions = sorted(utr_regions, key=lambda x: x['start'])
        else:
            # Negative strand, sort by start position in reverse
            sorted_regions = sorted(utr_regions, key=lambda x: x['start'], reverse=True)
        
        sequences = []
        ranges = []
        
        for region in sorted_regions:
            seq = self.extract_sequence(
                region['chr'], 
                region['start'], 
                region['end'], 
                region['strand']
            )
            if seq:
                sequences.append(seq)
                ranges.append(f"{region['chr']}:{region['start']}-{region['end']}")
        
        # Prepare metadata
        metadata = {
            'strand': strand,
            'gene_id': utr_regions[0].get('gene_id', ''),
            'gene_name': utr_regions[0].get('gene_name', ''),
            'ranges': ';'.join(ranges),
            'length': sum(len(seq) for seq in sequences)
        }
        
        return ''.join(sequences), metadata
    
    def extract_all_utrs(self, output_file: str) -> None:
        """Extract all UTR sequences and write to FASTA file"""
        print(f"Extracting UTR sequences to: {output_file}")
        
        with open(output_file, 'w') as f:
            for transcript_id, utr_info in self.utr_data.items():
                self.stats['transcripts_processed'] += 1
                
                # Extract 5'UTR
                if utr_info['5utr']:
                    five_utr_seq, metadata = self.concatenate_utr_sequences(utr_info['5utr'], transcript_id)
                    if five_utr_seq:
                        # Create detailed FASTA header
                        gene_info = f"{metadata['gene_id']}"
                        if metadata['gene_name']:
                            gene_info += f"|{metadata['gene_name']}"
                        
                        header = (f">{transcript_id}_5UTR "
                                f"length={metadata['length']} "
                                f"strand={metadata['strand']} "
                                f"gene={gene_info} "
                                f"range={metadata['ranges']}")
                        
                        f.write(f"{header}\n")
                        f.write(f"{five_utr_seq}\n")
                        self.stats['sequences_extracted'] += 1
                
                # Extract 3'UTR
                if utr_info['3utr']:
                    three_utr_seq, metadata = self.concatenate_utr_sequences(utr_info['3utr'], transcript_id)
                    if three_utr_seq:
                        # Create detailed FASTA header
                        gene_info = f"{metadata['gene_id']}"
                        if metadata['gene_name']:
                            gene_info += f"|{metadata['gene_name']}"
                        
                        header = (f">{transcript_id}_3UTR "
                                f"length={metadata['length']} "
                                f"strand={metadata['strand']} "
                                f"gene={gene_info} "
                                f"range={metadata['ranges']}")
                        
                        f.write(f"{header}\n")
                        f.write(f"{three_utr_seq}\n")
                        self.stats['sequences_extracted'] += 1
                
                if self.stats['transcripts_processed'] % 1000 == 0:
                    print(f"Processed {self.stats['transcripts_processed']} transcripts...")
        
        print("\nExtraction completed!")
        print(f"Transcripts processed: {self.stats['transcripts_processed']}")
        print(f"Sequences extracted: {self.stats['sequences_extracted']}")
    
    def get_statistics(self) -> Dict:
        """Get statistics information"""
        return self.stats.copy()
    
    def process(self, gtf_file: str, fasta_file: str, output_file: str) -> None:
        """Complete processing workflow"""
        try:
            # Load reference genome
            self.load_fasta(fasta_file)
            
            # Load GTF file
            self.load_gtf(gtf_file)
            
            # Extract UTR sequences
            self.extract_all_utrs(output_file)
            
        except Exception as e:
            print(f"Error: {e}")
            raise