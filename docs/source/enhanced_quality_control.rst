==========================
Enhanced Quality Control
==========================

.. versionadded:: 2.0
   Enhanced quality control pipeline with aDNA artifact removal and advanced filtering.

The enhanced quality control pipeline provides advanced tools for processing ancient DNA sequences,
focusing on artifact removal, quality filtering, and genetic diversity analysis to optimize 
haplogroup classification results.

Overview
========

The enhanced pipeline consists of three main components:

1. **aDNA Sequence Cleaner** - Removes ancient DNA artifacts and ambiguous nucleotides
2. **Improved FASTA to HSD Converter** - Advanced conversion with quality filtering
3. **HSD Diversity Analyzer** - Comprehensive genetic diversity assessment

Features
========

aDNA Sequence Cleaning
---------------------

The :class:`~sanger_pipeline.utils.adna_sequence_cleaner.aDNASequenceCleaner` provides:

* **Artifact Removal**: Eliminates common ancient DNA artifacts
* **Ambiguous Nucleotide Resolution**: Converts ambiguous nucleotides to preferred bases
* **Quality Assessment**: Evaluates sequence quality based on valid nucleotide proportion
* **Length Filtering**: Removes sequences below minimum length thresholds

Quality Filtering
----------------

The :class:`~sanger_pipeline.utils.improved_fasta_to_hsd_converter.ImprovedFastaToHSDConverter` offers:

* **Quality Thresholds**: Configurable quality cutoffs (default: 70%)
* **Artifact Detection**: Advanced detection of alignment and sequencing artifacts
* **HVS Region Processing**: Independent processing of HVS1, HVS2, and HVS3 regions
* **Variant Optimization**: Produces optimal variant counts for haplogroup analysis

Diversity Analysis
-----------------

The :class:`~sanger_pipeline.utils.hsd_diversity_analyzer.HSDDiversityAnalyzer` provides:

* **Sample Comparison**: Pairwise similarity analysis between samples
* **Variant Statistics**: Comprehensive statistics on variant distributions
* **Position Diversity**: Analysis of diversity at specific genomic positions
* **Quality Flags**: Automatic detection of potential quality issues

Usage
=====

Quick Start
----------

Run the enhanced quality control pipeline on your processed samples:

.. code-block:: bash

   # Run enhanced quality control
   python enhanced_hsd_converter.py

This will automatically:

1. Combine consensus sequences from the pipeline output
2. Clean sequences using aDNA-specific algorithms
3. Convert to high-quality HSD format with filtering
4. Perform comprehensive diversity analysis

Manual Usage
-----------

You can also run each component individually:

.. code-block:: python

   from sanger_pipeline.utils.adna_sequence_cleaner import aDNASequenceCleaner
   from sanger_pipeline.utils.improved_fasta_to_hsd_converter import ImprovedFastaToHSDConverter
   from sanger_pipeline.utils.hsd_diversity_analyzer import HSDDiversityAnalyzer

   # Clean sequences
   cleaner = aDNASequenceCleaner(min_length=50, min_quality=0.6)
   cleaned_sequences = cleaner.clean_fasta_file("input.fasta", "cleaned.fasta")

   # Convert to HSD with quality filtering
   converter = ImprovedFastaToHSDConverter(min_quality_threshold=0.7)
   converter.convert_fasta_to_hsd("cleaned.fasta", "output.hsd")

   # Analyze diversity
   analyzer = HSDDiversityAnalyzer()
   samples = analyzer.parse_hsd_file("output.hsd")
   diversity_report = analyzer.analyze_diversity(samples)

Configuration
============

aDNA Sequence Cleaner
--------------------

.. code-block:: python

   cleaner = aDNASequenceCleaner(
       min_length=50,        # Minimum sequence length
       min_quality=0.6       # Minimum proportion of valid nucleotides
   )

Improved HSD Converter
---------------------

.. code-block:: python

   converter = ImprovedFastaToHSDConverter(
       reference_path="ref/rCRS.fasta",    # Reference sequence path
       min_quality_threshold=0.7           # Quality threshold (0-1)
   )

Output Files
===========

The enhanced quality control pipeline produces several output files:

* **{output}_final_cleaned.fasta**: Cleaned consensus sequences after artifact removal
* **{output}_final_high_quality.hsd**: High-quality HSD file with filtered variants
* **Diversity Analysis Report**: Printed statistics on genetic diversity and quality

Quality Metrics
==============

The pipeline provides comprehensive quality metrics:

Variant Statistics
-----------------

* **Total Samples**: Number of samples processed
* **Variant Range**: Minimum and maximum variants per sample
* **Mean Variants**: Average number of variants across samples
* **Unique Positions**: Total number of variant positions detected

Sample Similarity
-----------------

* **Pairwise Similarities**: Jaccard similarity between all sample pairs
* **Mean Similarity**: Average genetic similarity across samples
* **Similarity Range**: Minimum and maximum similarity values

Quality Flags
------------

* **Low Variant Count**: Samples with unusually low variant numbers
* **High Similarity**: Samples with potential contamination or duplication
* **Quality Issues**: Other potential quality concerns

Best Practices
=============

Quality Thresholds
-----------------

* Use **70% quality threshold** for optimal results (default)
* Lower thresholds (60-65%) for degraded samples if needed
* Higher thresholds (75-80%) for very high-quality requirements

Sample Selection
---------------

* Review diversity analysis report before final analysis
* Consider removing samples flagged for quality issues
* Prioritize samples with variant counts within expected ranges

Troubleshooting
==============

Common Issues
------------

**Low Sample Retention**

If too few samples pass quality filtering:

1. Lower the quality threshold (e.g., from 0.7 to 0.65)
2. Check input sequence quality
3. Review cleaning parameters

**High Sample Similarity**

If samples show unexpectedly high similarity:

1. Check for contamination or sample mix-ups
2. Review primer specificity
3. Consider additional quality filtering

**Unusual Variant Counts**

If variant counts are outside expected ranges:

1. Verify reference sequence alignment
2. Check HVS region detection
3. Review input sequence quality

API Reference
============

.. automodule:: sanger_pipeline.utils.adna_sequence_cleaner
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: sanger_pipeline.utils.improved_fasta_to_hsd_converter
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: sanger_pipeline.utils.hsd_diversity_analyzer
   :members:
   :undoc-members:
   :show-inheritance:
