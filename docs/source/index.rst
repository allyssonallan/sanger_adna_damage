=====================================
Sanger DNA Damage Analysis Pipeline
=====================================

.. image:: https://img.shields.io/badge/Python-3.8+-blue.svg
   :target: https://python.org
   :alt: Python Version

.. image:: https://img.shields.io/badge/License-MIT-green.svg
   :target: https://opensource.org/licenses/MIT
   :alt: MIT License

.. image:: https://img.shields.io/badge/Documentation-Sphinx-orange.svg
   :target: #
   :alt: Documentation

A comprehensive, modular pipeline for processing Sanger sequencing AB1 files, including quality control, 
alignment, consensus building, and ancient DNA damage analysis.

.. important::
   **IMPORTANT DISCLAIMER - Tool Purpose & Limitations**
   
   This pipeline is **NOT** a tool for authenticating ancient DNA samples. It is designed for:
   
   * **Prioritizing haplogroups** for follow-up analysis
   * **Evaluating sample quality** based on insert size and damage patterns
   * **Providing surrogate bootstrapped damage indicators** 
   * **Assisting in haplogroup origin assessment**
   * **Guiding selection** of promising samples for NGS sequencing
   
   **⚠️ All ancient DNA authentication must be performed using NGS-based methods with appropriate controls, contamination assessment, and phylogenetic analysis.**
   
   This tool provides preliminary screening to help researchers prioritize samples and resources before proceeding to more comprehensive NGS-based ancient DNA authentication workflows.

🚀 Features
===========

* **Modular Architecture**: Well-organized codebase with clear separation of concerns
* **Quality Control**: Convert AB1 files with Phred quality filtering and visualization
* **Sequence Processing**: Align forward/reverse reads and build consensus sequences
* **HVS Region Processing**: Independent processing of HVS1, HVS2, and HVS3 regions with intelligent merging
* **Ancient DNA Analysis**: Comprehensive aDNA damage pattern detection and assessment
* **Statistical Validation**: Bootstrap analysis for damage assessment
* **Beautiful QC Reports**: Interactive HTML reports with charts, tables, and analysis summaries
* **Command Line Interface**: Easy-to-use CLI for all pipeline operations
* **Extensible Design**: Easy to add new analysis modules and features

📚 Documentation Contents
=========================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   
   installation
   quickstart
   configuration

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   tutorials/index
   howto/index
   cli_reference

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics
   
   understanding_damage_analysis
   customization
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api/index
   api/core
   api/utils
   api/cli

.. toctree::
   :maxdepth: 1
   :caption: Development
   
   contributing
   changelog
   license

🎯 Quick Start
==============

Installation
------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/allyssonallan/sanger_adna_damage.git
   cd sanger_adna_damage

   # Create and activate virtual environment
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\\Scripts\\activate

   # Install dependencies
   pip install -r requirements.txt

   # Install the package in development mode
   pip install -e .

Basic Usage
-----------

.. code-block:: bash

   # Run the complete pipeline
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --input-dir ./input \\
       --output-dir ./output \\
       --config ./config/default_config.yaml

   # Generate comprehensive QC report
   python -m src.sanger_pipeline.cli.main generate-report \\
       --output-dir ./output \\
       --open-browser

📊 Pipeline Overview
===================

The Sanger DNA Damage Analysis Pipeline processes AB1 sequencing files through several stages:

1. **AB1 Conversion**: Convert binary AB1 files to FASTA format with quality filtering
2. **Sequence Alignment**: Align forward and reverse reads using MAFFT
3. **Consensus Building**: Generate consensus sequences for each HVS region
4. **HVS Region Merging**: Intelligently combine available HVS regions
5. **Damage Analysis**: Detect and assess ancient DNA damage patterns
6. **Report Generation**: Create comprehensive QC reports with visualizations

.. mermaid::

   graph TD
       A[AB1 Files] --> B[AB1 Conversion]
       B --> C[Quality Filtering]
       C --> D[Sequence Alignment]
       D --> E[Consensus Building]
       E --> F[HVS Region Merging]
       F --> G[Damage Analysis]
       G --> H[QC Report Generation]
       H --> I[Final Results]

🔬 Ancient DNA Analysis
======================

The pipeline includes sophisticated ancient DNA damage analysis:

* **Damage Pattern Detection**: Identifies characteristic C→T and G→A transitions
* **Statistical Validation**: Bootstrap analysis with 10,000 iterations
* **Damage Assessment**: Quantitative scoring of damage patterns
* **Visual Reports**: Damage plots and interactive visualizations

📈 Output Structure
==================

The pipeline creates organized output directories:

.. code-block:: text

   output/
   ├── fasta/              # Raw FASTA files converted from AB1
   ├── filtered/           # Quality-filtered sequences
   ├── consensus/          # Consensus sequences for each HVS region
   ├── aligned/            # Intermediate alignment files
   ├── final/              # Merged HVS region sequences
   ├── damage_analysis/    # aDNA damage analysis results
   ├── plots/              # Quality score plots
   └── reports/            # Interactive HTML QC reports

📞 Support and Contributing
===========================

* **Issues**: Report bugs and request features on GitHub
* **Discussions**: Join community discussions for help and ideas
* **Contributing**: See our contributing guide for development setup
* **Documentation**: This documentation is built with Sphinx

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
