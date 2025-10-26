#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GTF Processor Module
Extract UTR and CDS information from GTF files and correctly annotate 5'UTR and 3'UTR based on CDS position and strand direction

Features:
1. Parse GTF files (supports gzip compression)
2. Extract exon, UTR, CDS information from protein_coding transcripts
3. Re-annotate UTRs as 5'UTR or 3'UTR based on CDS position and strand direction
4. Output processed GTF file
"""

import gzip
import os
from collections import defaultdict
import re


class GTFProcessor:
    """GTF file processor for extracting and classifying UTR regions"""
    
    def __init__(self, input_gtf, output_gtf):
        """
        Initialize GTF processor
        
        Args:
            input_gtf (str): Input GTF file path
            output_gtf (str): Output GTF file path
        """
        self.input_gtf = input_gtf
        self.output_gtf = output_gtf
        
        # Store transcript information
        self.transcripts = defaultdict(lambda: {
            'gene_id': '',
            'gene_name': '',
            'chrom': '',
            'strand': '',
            'exons': [],
            'utrs': [],
            'cds': [],
            'gene_type': '',
            'transcript_type': ''
        })
        
        # Statistics information
        self.stats = {
            'total_transcripts': 0,
            'protein_coding_transcripts': 0,
            'transcripts_with_cds': 0,
            'transcripts_with_utrs': 0,
            'output_transcripts': 0
        }

    def parse_attributes(self, attr_string):
        """Parse GTF attributes string"""
        attributes = {}
        # Match key "value" format
        pattern = r'(\w+)\s+"([^"]*)"'
        matches = re.findall(pattern, attr_string)
        for key, value in matches:
            attributes[key] = value
        return attributes

    def parse_gtf_line(self, line):
        """Parse GTF line"""
        if line.startswith('#') or not line.strip():
            return None
        
        fields = line.strip().split('\t')
        if len(fields) != 9:
            return None
        
        chrom, source, feature, start, end, score, strand, frame, attributes = fields
        
        return {
            'chrom': chrom,
            'source': source,
            'feature': feature,
            'start': int(start),
            'end': int(end),
            'score': score,
            'strand': strand,
            'frame': frame,
            'attributes': self.parse_attributes(attributes)
        }

    def load_gtf(self):
        """Load GTF file"""
        print(f"Loading GTF file: {self.input_gtf}")
        
        # Check if it's a gzip file
        open_func = gzip.open if self.input_gtf.endswith('.gz') else open
        mode = 'rt' if self.input_gtf.endswith('.gz') else 'r'
        
        with open_func(self.input_gtf, mode, encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 100000 == 0:
                    print(f"Processed {line_num} lines...")
                
                record = self.parse_gtf_line(line)
                if not record:
                    continue
                
                # Only process protein_coding genes
                if record['attributes'].get('gene_type') != 'protein_coding':
                    continue
                
                transcript_id = record['attributes'].get('transcript_id')
                if not transcript_id:
                    continue
                
                # Only process protein_coding transcripts
                if record['attributes'].get('transcript_type') != 'protein_coding':
                    continue
                
                # Store basic transcript information
                transcript = self.transcripts[transcript_id]
                transcript['gene_id'] = record['attributes'].get('gene_id', '')
                transcript['gene_name'] = record['attributes'].get('gene_name', '')
                transcript['chrom'] = record['chrom']
                transcript['strand'] = record['strand']
                transcript['gene_type'] = record['attributes'].get('gene_type', '')
                transcript['transcript_type'] = record['attributes'].get('transcript_type', '')
                
                # Store information based on feature type
                if record['feature'] == 'exon':
                    transcript['exons'].append({
                        'start': record['start'],
                        'end': record['end'],
                        'exon_number': record['attributes'].get('exon_number', '')
                    })
                elif record['feature'] in ['UTR', 'five_prime_utr', 'three_prime_utr']:
                    transcript['utrs'].append({
                        'start': record['start'],
                        'end': record['end'],
                        'original_type': record['feature']
                    })
                elif record['feature'] == 'CDS':
                    transcript['cds'].append({
                        'start': record['start'],
                        'end': record['end'],
                        'frame': record['frame']
                    })
        
        print(f"GTF file loading completed, processed {len(self.transcripts)} protein_coding transcripts")

    def classify_utrs(self):
        """Re-classify UTRs based on CDS position and strand direction"""
        print("Re-classifying UTRs...")
        
        for transcript_id, transcript in self.transcripts.items():
            if not transcript['cds'] or not transcript['utrs']:
                continue
            
            # Get CDS start and end positions
            cds_starts = [cds['start'] for cds in transcript['cds']]
            cds_ends = [cds['end'] for cds in transcript['cds']]
            cds_min = min(cds_starts)
            cds_max = max(cds_ends)
            
            strand = transcript['strand']
            
            # Re-classify UTRs
            classified_utrs = []
            for utr in transcript['utrs']:
                utr_start = utr['start']
                utr_end = utr['end']
                
                if strand == '+':
                    # Positive strand: 5'UTR upstream of CDS, 3'UTR downstream of CDS
                    if utr_end < cds_min:
                        utr_type = 'five_prime_utr'
                    elif utr_start > cds_max:
                        utr_type = 'three_prime_utr'
                    elif utr_start < cds_min and utr_end >= cds_min:
                        # UTR spans CDS start position, split as 5'UTR
                        utr_type = 'five_prime_utr'
                        utr['end'] = cds_min - 1  # Adjust end position
                    elif utr_start <= cds_max and utr_end > cds_max:
                        # UTR spans CDS end position, split as 3'UTR
                        utr_type = 'three_prime_utr'
                        utr['start'] = cds_max + 1  # Adjust start position
                    else:
                        # UTR completely within CDS, skip
                        continue
                else:  # strand == '-'
                    # Negative strand: 5'UTR downstream of CDS, 3'UTR upstream of CDS
                    if utr_start > cds_max:
                        utr_type = 'five_prime_utr'
                    elif utr_end < cds_min:
                        utr_type = 'three_prime_utr'
                    elif utr_start <= cds_max and utr_end > cds_max:
                        # UTR spans CDS end position, split as 5'UTR
                        utr_type = 'five_prime_utr'
                        utr['start'] = cds_max + 1  # Adjust start position
                    elif utr_start < cds_min and utr_end >= cds_min:
                        # UTR spans CDS start position, split as 3'UTR
                        utr_type = 'three_prime_utr'
                        utr['end'] = cds_min - 1  # Adjust end position
                    else:
                        # UTR completely within CDS, skip
                        continue
                
                # Ensure adjusted coordinates are valid
                if utr['start'] <= utr['end']:
                    utr['classified_type'] = utr_type
                    classified_utrs.append(utr)
            
            transcript['classified_utrs'] = classified_utrs

    def write_output_gtf(self):
        """Output processed GTF file"""
        print(f"Writing output file: {self.output_gtf}")
        
        with open(self.output_gtf, 'w', encoding='utf-8') as f:
            # Write header comments
            f.write("##gff-version 2\n")
            f.write("##source: GTF UTR Processor\n")
            f.write("##description: Processed GTF with classified UTRs\n")
            
            for transcript_id, transcript in self.transcripts.items():
                if not transcript['cds']:
                    continue  # Skip transcripts without CDS
                
                chrom = transcript['chrom']
                strand = transcript['strand']
                gene_id = transcript['gene_id']
                gene_name = transcript['gene_name']
                gene_type = transcript['gene_type']
                transcript_type = transcript['transcript_type']
                
                # Build basic attributes string
                base_attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_type "{gene_type}"; gene_name "{gene_name}"; transcript_type "{transcript_type}";'
                
                # Write CDS
                for i, cds in enumerate(sorted(transcript['cds'], key=lambda x: x['start'])):
                    f.write(f"{chrom}\tprocessed\tCDS\t{cds['start']}\t{cds['end']}\t.\t{strand}\t{cds['frame']}\t{base_attrs}\n")
                
                # Write classified UTRs
                if 'classified_utrs' in transcript:
                    for utr in sorted(transcript['classified_utrs'], key=lambda x: x['start']):
                        utr_type = utr['classified_type']
                        f.write(f"{chrom}\tprocessed\t{utr_type}\t{utr['start']}\t{utr['end']}\t.\t{strand}\t.\t{base_attrs}\n")
                
                self.stats['output_transcripts'] += 1

    def print_statistics(self):
        """Print statistics information"""
        print("\n=== Processing Statistics ===")
        print(f"Total input transcripts: {len(self.transcripts)}")
        
        transcripts_with_cds = sum(1 for t in self.transcripts.values() if t['cds'])
        transcripts_with_utrs = sum(1 for t in self.transcripts.values() if t['utrs'])
        transcripts_with_classified_utrs = sum(1 for t in self.transcripts.values() if 'classified_utrs' in t and t['classified_utrs'])
        
        print(f"Transcripts with CDS: {transcripts_with_cds}")
        print(f"Transcripts with UTRs: {transcripts_with_utrs}")
        print(f"Transcripts with classified UTRs: {transcripts_with_classified_utrs}")
        print(f"Output transcripts: {self.stats['output_transcripts']}")

    def process(self):
        """Run processing workflow"""
        print("Starting GTF UTR processing...")
        self.load_gtf()
        self.classify_utrs()
        self.write_output_gtf()
        self.print_statistics()
        print("GTF UTR processing completed!")