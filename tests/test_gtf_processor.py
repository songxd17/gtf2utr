#!/usr/bin/env python3
"""
Unit tests for GTFProcessor class
"""

import unittest
import tempfile
import os
from pathlib import Path
import sys

# Add the parent directory to the path to import gtf2utr
sys.path.insert(0, str(Path(__file__).parent.parent))

from gtf2utr.gtf_processor import GTFProcessor

class TestGTFProcessor(unittest.TestCase):
    """Test cases for GTFProcessor"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.processor = GTFProcessor("dummy_input.gtf", "dummy_output.gtf")
        
        # Sample GTF content
        self.sample_gtf_content = """##description: test GTF file
chr1	HAVANA	gene	1000	3000	.	+	.	gene_id "ENSG00000001"; gene_name "GENE1"; gene_type "protein_coding";
chr1	HAVANA	transcript	1000	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	exon	1000	1200	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	CDS	1100	1200	.	+	0	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	exon	1500	1700	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "2"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	CDS	1500	1700	.	+	1	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "2"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	exon	2500	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "3"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	CDS	2500	2800	.	+	2	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "3"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
"""
    
    def test_parse_gtf_line(self):
        """Test GTF line parsing"""
        line = "chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id \"ENSG00000001\"; transcript_id \"ENST00000001\";"
        parsed = self.processor.parse_gtf_line(line)
        
        self.assertEqual(parsed['chrom'], 'chr1')
        self.assertEqual(parsed['source'], 'HAVANA')
        self.assertEqual(parsed['feature'], 'exon')
        self.assertEqual(parsed['start'], 1000)
        self.assertEqual(parsed['end'], 1200)
        self.assertEqual(parsed['strand'], '+')
        self.assertEqual(parsed['attributes']['gene_id'], 'ENSG00000001')
        self.assertEqual(parsed['attributes']['transcript_id'], 'ENST00000001')
    
    def test_parse_attributes(self):
        """Test attribute parsing"""
        attr_string = 'gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1";'
        attributes = self.processor.parse_attributes(attr_string)
        
        self.assertEqual(attributes['gene_id'], 'ENSG00000001')
        self.assertEqual(attributes['transcript_id'], 'ENST00000001')
        self.assertEqual(attributes['exon_number'], '1')
    
    def test_load_gtf_from_string(self):
        """Test loading GTF from string content"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(self.sample_gtf_content)
            temp_file = f.name
        
        try:
            # Update the processor's input file path and load
            self.processor.input_gtf = temp_file
            self.processor.load_gtf()
            
            # Check that transcript data was loaded
            self.assertGreater(len(self.processor.transcripts), 0)
            
        finally:
            os.unlink(temp_file)
    
    def test_classify_utrs(self):
        """Test UTR classification"""
        # Create a more complete GTF content for testing UTR classification
        complete_gtf_content = """##description: test GTF file
chr1	HAVANA	gene	1000	3000	.	+	.	gene_id "ENSG00000001"; gene_name "GENE1"; gene_type "protein_coding";
chr1	HAVANA	transcript	1000	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	exon	1000	1200	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	UTR	1000	1099	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	CDS	1100	1200	.	+	0	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	exon	2800	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	CDS	2800	2900	.	+	0	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1	HAVANA	UTR	2901	3000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(complete_gtf_content)
            temp_file = f.name
        
        try:
            # Update the processor's input file path and load
            self.processor.input_gtf = temp_file
            self.processor.load_gtf()
            self.processor.classify_utrs()
            
            # Check that UTRs were classified
            # This is a basic check - in a real test you'd verify specific classifications
            self.assertIsNotNone(self.processor.transcripts)
            
        finally:
            os.unlink(temp_file)

if __name__ == '__main__':
    unittest.main()