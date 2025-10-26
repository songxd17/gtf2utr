"""
gtf2utr: A Python package for extracting UTR sequences from GTF and FASTA files.

This package provides tools to:
1. Process GTF files to classify UTR regions
2. Extract UTR sequences from reference genome FASTA files
3. Generate FASTA files with UTR sequences and detailed metadata

Author: gtf2utr team
License: MIT
"""

__version__ = "1.0.0"
__author__ = "gtf2utr team"
__email__ = "contact@gtf2utr.org"

from .gtf_processor import GTFProcessor
from .utr_extractor import UTRExtractor

__all__ = ["GTFProcessor", "UTRExtractor"]