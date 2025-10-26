#!/usr/bin/env python3
"""
Unit tests for UTRExtractor class
"""

import unittest
import tempfile
import os
from pathlib import Path
import sys

# Add the parent directory to the path to import gtf2utr
sys.path.insert(0, str(Path(__file__).parent.parent))

from gtf2utr.utr_extractor import UTRExtractor

class TestUTRExtractor(unittest.TestCase):
    """Test cases for UTRExtractor"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.extractor = UTRExtractor()
        
        # Sample FASTA content
        self.sample_fasta_content = """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
        
        # Sample GTF content with UTR annotations
        self.sample_gtf_content = """##description: test GTF file
chr1	HAVANA	five_prime_utr	1000	1099	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1";
chr1	HAVANA	three_prime_utr	2801	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1";
"""
    
    def test_parse_gtf_attributes(self):
        """Test GTF attribute parsing"""
        attr_string = 'gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1";'
        attributes = self.extractor.parse_gtf_attributes(attr_string)
        
        self.assertEqual(attributes['gene_id'], 'ENSG00000001')
        self.assertEqual(attributes['transcript_id'], 'ENST00000001')
        self.assertEqual(attributes['gene_name'], 'GENE1')
    
    def test_load_fasta(self):
        """Test FASTA file loading"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(self.sample_fasta_content)
            temp_file = f.name
        
        try:
            self.extractor.load_fasta(temp_file)
            
            # Check that chromosome was loaded
            self.assertIn('chr1', self.extractor.chromosomes)
            self.assertGreater(len(self.extractor.chromosomes['chr1']), 0)
            
        finally:
            os.unlink(temp_file)
    
    def test_load_gtf(self):
        """Test GTF file loading"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(self.sample_gtf_content)
            temp_file = f.name
        
        try:
            self.extractor.load_gtf(temp_file)
            
            # Check that UTR data was loaded
            self.assertIn('ENST00000001', self.extractor.utr_data)
            self.assertGreater(len(self.extractor.utr_data['ENST00000001']['5utr']), 0)
            self.assertGreater(len(self.extractor.utr_data['ENST00000001']['3utr']), 0)
            
        finally:
            os.unlink(temp_file)
    
    def test_reverse_complement(self):
        """Test reverse complement function"""
        sequence = "ATCG"
        expected = "CGAT"
        result = self.extractor.reverse_complement(sequence)
        self.assertEqual(result, expected)
        
        # Test with longer sequence
        sequence = "ATCGATCG"
        expected = "CGATCGAT"
        result = self.extractor.reverse_complement(sequence)
        self.assertEqual(result, expected)
    
    def test_extract_sequence(self):
        """Test sequence extraction"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(self.sample_fasta_content)
            temp_file = f.name
        
        try:
            self.extractor.load_fasta(temp_file)
            
            # Test forward strand extraction
            sequence = self.extractor.extract_sequence('chr1', 1, 10, '+')
            self.assertEqual(len(sequence), 10)
            
            # Test reverse strand extraction
            sequence = self.extractor.extract_sequence('chr1', 1, 10, '-')
            self.assertEqual(len(sequence), 10)
            
        finally:
            os.unlink(temp_file)
    
    def test_get_statistics(self):
        """Test statistics retrieval"""
        stats = self.extractor.get_statistics()
        
        # Check that all expected keys are present
        expected_keys = ['transcripts_processed', 'five_utr_regions', 'three_utr_regions', 'sequences_extracted']
        for key in expected_keys:
            self.assertIn(key, stats)
            self.assertIsInstance(stats[key], int)

if __name__ == '__main__':
    unittest.main()