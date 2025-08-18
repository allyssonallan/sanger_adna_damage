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
alignment, consensus building, ancient DNA damage analysis, and enhanced quality control for optimal 
haplogroup classification.

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

Core Pipeline
-------------

* **Modular Architecture**: Well-organized codebase with clear separation of concerns
* **Quality Control**: Convert AB1 files with Phred quality filtering and visualization
* **Sequence Processing**: Align forward/reverse reads and build consensus sequences
* **HVS Region Processing**: Independent processing of HVS1, HVS2, and HVS3 regions with intelligent merging
* **Ancient DNA Analysis**: Comprehensive aDNA damage pattern detection and assessment
* **Statistical Validation**: Bootstrap analysis for damage assessment
* **Beautiful QC Reports**: Interactive HTML reports with charts, tables, and analysis summaries
* **Command Line Interface**: Easy-to-use CLI for all pipeline operations
* **Extensible Design**: Easy to add new analysis modules and features

Enhanced Quality Control (NEW!)
-------------------------------

* **aDNA Sequence Cleaning**: Advanced removal of ancient DNA artifacts and ambiguous nucleotides
* **Quality Filtering**: Configurable quality thresholds with 70% default for optimal results
* **Diversity Analysis**: Comprehensive genetic diversity assessment and sample comparison
* **Sample Prioritization**: Automated identification of highest-quality samples for downstream analysis
* **Quality Metrics**: Detailed reports on variant counts, sample similarity, and potential quality issues
* **Artifact Detection**: Advanced detection and removal of alignment and sequencing artifacts

HSD Conversion Methods
---------------------

* **Regional Hybrid Method**: Optimal approach with 52.4 average variants per sample (recommended)
* **Direct Method**: Alternative approach with 66.0 average variants per sample  
* **Enhanced Converter**: Improved quality control with artifact detection and filtering

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
   enhanced_quality_control
   pipeline_workflow
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

The pipeline processes Sanger sequencing data through multiple quality-controlled stages with comprehensive branching for different quality control approaches and output formats:

.. mermaid::

   graph TB
       subgraph "Input Stage"
           A[📁 AB1 Files<br/>- Forward reads<br/>- Reverse reads<br/>- HVS1/2/3 regions]
       end
       
       subgraph "Core Processing"
           B[🔄 AB1 Conversion<br/>Quality filtering<br/>Phred scores ≥ Q20/Q30]
           C[🧹 Quality Control<br/>Length filtering<br/>Base quality assessment]
           D[🔗 Consensus Building<br/>Forward/reverse alignment<br/>Per HVS region]
           E[🧩 Region Merging<br/>Combine HVS regions<br/>Sample consolidation]
       end
       
       subgraph "Analysis & QC"
           F[🧬 Damage Analysis<br/>C→T, G→A transitions<br/>Bootstrap statistics<br/>P-value calculation]
           G[📊 Interactive Reports<br/>HTML dashboard<br/>Quality visualizations<br/>Statistical summaries]
       end
       
       subgraph "Enhanced Quality Control (v2.0+)"
           H[✨ Enhanced Pipeline Entry]
           I[🧪 aDNA Sequence Cleaner<br/>- Remove artifacts<br/>- Resolve ambiguous bases<br/>- Filter poly-N regions<br/>- Quality scoring]
           J[📝 Improved HSD Converter<br/>- Reference alignment<br/>- Quality metrics<br/>- Variant filtering<br/>- Statistical validation]
           K[📈 Diversity Analyzer<br/>- Haplogroup diversity<br/>- Sample comparison<br/>- Quality assessment<br/>- Priority ranking]
       end
       
       subgraph "Output Options"
           L1[📋 Standard HSD<br/>Basic variant calling<br/>Regional/Direct methods]
           L2[🎯 Enhanced HSD<br/>Quality-filtered variants<br/>Statistical confidence<br/>Diversity metrics]
           L3[📊 Quality Reports<br/>Interactive dashboards<br/>Damage plots<br/>Statistical summaries]
           L4[📁 Processed Sequences<br/>FASTA files<br/>Consensus sequences<br/>Quality scores]
       end
       
       subgraph "Configuration Variables"
           V1[⚙️ Quality Thresholds<br/>--min-quality: 15-30<br/>--min-length: 30-100bp<br/>--quality-filter: 0.6-0.8]
           V2[🔧 Pipeline Options<br/>--alignment-tool: mafft/muscle<br/>--damage-threshold: 0.02<br/>--bootstrap-iterations: 1000]
           V3[📂 I/O Directories<br/>--input-dir: AB1 files<br/>--output-dir: Results<br/>--config: YAML settings]
       end
       
       %% Main workflow
       A --> B
       B --> C
       C --> D
       D --> E
       E --> F
       F --> G
       
       %% Enhanced workflow branch
       E -.-> H
       H --> I
       I --> J
       J --> K
       
       %% Output generation
       E --> L1
       F --> L3
       G --> L3
       J --> L2
       K --> L2
       D --> L4
       E --> L4
       
       %% Configuration influences
       V1 -.-> B
       V1 -.-> C
       V1 -.-> I
       V2 -.-> D
       V2 -.-> F
       V2 -.-> J
       V3 -.-> A
       V3 -.-> L1
       V3 -.-> L2
       V3 -.-> L3
       V3 -.-> L4
       
       %% Alternative paths
       G -.-> L1
       L3 -.-> H
       
       %% Styling
       style A fill:#e3f2fd
       style H fill:#fff3e0
       style L2 fill:#e8f5e8
       style F fill:#fce4ec
       style K fill:#f3e5f5
       
       classDef inputNode fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
       classDef coreNode fill:#e8f5e8,stroke:#388e3c,stroke-width:2px
       classDef enhancedNode fill:#fff3e0,stroke:#f57c00,stroke-width:2px
       classDef outputNode fill:#fce4ec,stroke:#c2185b,stroke-width:2px
       classDef configNode fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
       
       class A inputNode
       class B,C,D,E coreNode  
       class H,I,J,K enhancedNode
       class L1,L2,L3,L4 outputNode
       class V1,V2,V3 configNode

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
