# gtf2utr

[![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/gtf2utr.svg)](https://badge.fury.io/py/gtf2utr)

A Python package for extracting UTR sequences from GTF and FASTA files.

## Features

- 🧬 **GTF File Processing**: Parse GTF files and correctly classify 5'UTR and 3'UTR based on CDS positions
- 🔍 **Sequence Extraction**: Extract UTR sequences from reference genome FASTA files
- 🧵 **Strand Handling**: Correctly handle sequence orientation for both forward and reverse strands
- 📊 **Detailed Metadata**: Output FASTA headers containing length, strand, gene information, and coordinate ranges
- 🗜️ **Compressed File Support**: Support for gzip-compressed GTF and FASTA files
- ⚡ **Efficient Processing**: Optimized algorithms for handling large genomic datasets

## Installation

### Install from PyPI (Recommended)

```bash
pip install gtf2utr
```

### Install from Source

```bash
git clone https://github.com/songxd17/gtf2utr.git
cd gtf2utr
pip install -e .
```

## Quick Start

### Command Line Usage

#### 1. Complete Pipeline (Recommended)

```bash
# Run the complete GTF processing and UTR extraction pipeline
gtf2utr pipeline input.gtf.gz reference.fa.gz output_utrs.fa
```

#### 2. Step-by-Step Execution

```bash
# Step 1: Process GTF file and classify UTR regions
gtf2utr process input.gtf.gz processed_utrs.gtf

# Step 2: Extract UTR sequences from reference genome
gtf2utr extract processed_utrs.gtf reference.fa.gz output_utrs.fa
```

#### 3. Individual Commands

```bash
# Process GTF file only
gtf2utr-process input.gtf.gz processed_utrs.gtf

# Extract UTR sequences only
gtf2utr-extract processed_utrs.gtf reference.fa.gz output_utrs.fa
```

### Python API Usage

```python
from gtf2utr import GTFProcessor, UTRExtractor

# Process GTF file
processor = GTFProcessor('input.gtf.gz', 'processed_utrs.gtf')
processor.process()

# Extract UTR sequences
extractor = UTRExtractor()
extractor.process('processed_utrs.gtf', 'reference.fa.gz', 'output_utrs.fa')

# Get statistics
stats = extractor.get_statistics()
print(f"Transcripts processed: {stats['transcripts_processed']}")
print(f"Sequences extracted: {stats['sequences_extracted']}")
```

## Input File Formats

### GTF Files

Supports standard GTF format files containing the following features:
- `exon`: Exon regions
- `CDS`: Coding sequence regions
- `UTR`, `five_prime_utr`, `three_prime_utr`: UTR regions

Example:
```
chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chr1    HAVANA  UTR     11869   12227   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328";
```

### FASTA Files

Standard genomic FASTA files supporting multiple chromosomes:
```
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>chr2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

## Output Format

The output FASTA file contains detailed sequence information:

```
>ENST00000456328_5UTR length=358 strand=+ gene=ENSG00000223972|DDX11L1 range=chr1:11869-12227
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>ENST00000456328_3UTR length=1021 strand=+ gene=ENSG00000223972|DDX11L1 range=chr1:12613-13221;chr1:13453-13670
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

Header information includes:
- `transcript_id_UTRtype`: Transcript ID and UTR type (5UTR or 3UTR)
- `length`: Total UTR sequence length
- `strand`: Strand direction (+ or -)
- `gene`: Gene ID and gene name
- `range`: Genomic coordinate ranges (multiple regions separated by semicolons)


## Example Data

Example data and usage scripts are provided in the `examples/` directory:

```bash
# Run example
cd examples
python run_example.py
```

## Testing

Run the test suite:

```bash
# Install test dependencies
pip install -e .[dev]

# Run tests
pytest tests/

# Run tests with coverage report
pytest --cov=gtf2utr tests/
```


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
- Email: songxiaodong@lglab.ac.cn

