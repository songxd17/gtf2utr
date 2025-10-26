#!/usr/bin/env python3
"""
Setup script for gtf2utr package
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="gtf2utr",
    version="1.0.0",
    author="gtf2utr team",
    author_email="contact@gtf2utr.org",
    description="A Python package for extracting UTR sequences from GTF and FASTA files",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/gtf2utr",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "gtf2utr-process=gtf2utr.cli:process_gtf",
            "gtf2utr-extract=gtf2utr.cli:extract_utrs",
            "gtf2utr=gtf2utr.cli:main",
        ],
    },
    keywords="bioinformatics genomics UTR GTF FASTA sequence-analysis",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/gtf2utr/issues",
        "Source": "https://github.com/yourusername/gtf2utr",
        "Documentation": "https://github.com/yourusername/gtf2utr/blob/main/README.md",
    },
    include_package_data=True,
    zip_safe=False,
)