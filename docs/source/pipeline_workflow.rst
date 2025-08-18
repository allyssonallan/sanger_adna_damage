=======================================
Complete Pipeline Workflow Reference
=======================================

This document provides a comprehensive overview of the entire Sanger aDNA damage analysis pipeline, including all possible pathways, quality control options, configuration variables, and output types.

Pipeline Architecture Overview
==============================

The pipeline is designed with a modular architecture that supports multiple processing pathways:

1. **Standard Pipeline**: Core processing workflow for routine analysis
2. **Enhanced Quality Control**: Advanced aDNA-specific processing (v2.0+)
3. **Manual Tools**: Individual components for custom workflows
4. **Reporting System**: Comprehensive quality assessment and visualization

Complete Workflow Diagram
==========================

.. mermaid::

   graph TB
       subgraph "📁 Input Data"
           A1[AB1 Files<br/>Forward Reads]
           A2[AB1 Files<br/>Reverse Reads]  
           A3[Reference Sequences<br/>rCRS, HVS regions]
           A4[Configuration Files<br/>YAML settings]
       end
       
       subgraph "🔧 Configuration Variables"
           V1[Quality Settings<br/>--min-quality: 15-30<br/>--min-length: 30-100bp<br/>--quality-threshold: 0.6-0.8]
           V2[Pipeline Parameters<br/>--alignment-tool: mafft/muscle<br/>--alignment-params: --auto<br/>--damage-threshold: 0.02]
           V3[I/O Configuration<br/>--input-dir: source path<br/>--output-dir: results path<br/>--config: settings file]
           V4[Enhanced QC Options<br/>--aggressive-cleaning: bool<br/>--reference-aware: bool<br/>--bootstrap-iterations: 1000]
       end
       
       subgraph "🔄 Stage 1: File Conversion & Initial QC"
           B1[AB1 Converter<br/>Extract sequences & quality scores]
           B2[Quality Filtering<br/>Phred score filtering<br/>Length requirements]
           B3[Format Conversion<br/>AB1 → FASTA/FASTQ<br/>Quality score preservation]
           B4[Initial QC Check<br/>File integrity<br/>Sequence validity]
       end
       
       subgraph "🧬 Stage 2: Sequence Processing"
           C1[Forward Sequence Processing<br/>Quality trimming<br/>Artifact removal]
           C2[Reverse Sequence Processing<br/>Quality trimming<br/>Artifact removal]
           C3[Sequence Alignment<br/>MAFFT/MUSCLE alignment<br/>Parameter optimization]
           C4[Consensus Generation<br/>Forward/reverse merging<br/>Conflict resolution]
       end
       
       subgraph "🧩 Stage 3: Regional Analysis"
           D1[HVS1 Processing<br/>16024-16365 bp<br/>Regional alignment]
           D2[HVS2 Processing<br/>57-372 bp<br/>Regional alignment]
           D3[HVS3 Processing<br/>438-574 bp<br/>Regional alignment]
           D4[Regional Merging<br/>Combine available regions<br/>Sample consolidation]
       end
       
       subgraph "🔬 Stage 4: Damage Analysis"
           E1[Damage Pattern Detection<br/>C→T transitions (5')<br/>G→A transitions (3')]
           E2[Statistical Analysis<br/>Bootstrap validation<br/>P-value calculation<br/>Confidence intervals]
           E3[Background Comparison<br/>Modern DNA controls<br/>Significance testing]
           E4[Damage Scoring<br/>Composite damage scores<br/>Assessment categories]
       end
       
       subgraph "✨ Enhanced Quality Control Branch (v2.0+)"
           F1[Pipeline Entry Point<br/>Enhanced mode trigger]
           F2[aDNA Sequence Cleaner<br/>- Artifact removal<br/>- Ambiguous base resolution<br/>- Poly-N filtering<br/>- Quality rescoring]
           F3[Improved HSD Converter<br/>- Reference alignment<br/>- Quality-based filtering<br/>- Statistical validation<br/>- Enhanced variant calling]
           F4[Diversity Analyzer<br/>- Haplogroup diversity<br/>- Sample comparison<br/>- Quality ranking<br/>- Priority assessment]
       end
       
       subgraph "📊 Stage 5: Quality Control & Reporting"
           G1[Quality Metrics Calculation<br/>Sequence quality scores<br/>Coverage statistics<br/>Processing success rates]
           G2[Statistical Summaries<br/>Sample-level statistics<br/>Batch-level summaries<br/>Comparative analysis]
           G3[Visualization Generation<br/>Quality plots<br/>Damage profiles<br/>Interactive charts]
           G4[Report Compilation<br/>HTML dashboard<br/>PDF summaries<br/>CSV exports]
       end
       
       subgraph "📝 Stage 6: Output Generation"
           H1[Standard HSD Output<br/>Basic variant calling<br/>Regional method<br/>Direct method]
           H2[Enhanced HSD Output<br/>Quality-filtered variants<br/>Statistical confidence<br/>Reference-aligned calls]
           H3[FASTA Sequences<br/>Raw conversions<br/>Filtered sequences<br/>Consensus sequences<br/>Final merged sequences]
           H4[Quality Reports<br/>Interactive HTML<br/>Statistical summaries<br/>Processing logs<br/>Error reports]
           H5[Diagnostic Files<br/>Alignment files<br/>Intermediate outputs<br/>Debug information]
       end
       
       subgraph "🎯 Alternative Processing Paths"
           I1[Manual Tool Access<br/>Individual component usage<br/>Custom parameter sets]
           I2[Batch Processing<br/>Multiple sample handling<br/>Parallel execution]
           I3[Reprocessing Options<br/>Parameter adjustment<br/>Selective re-running]
           I4[Integration Endpoints<br/>External tool compatibility<br/>Pipeline chaining]
       end
       
       %% Main workflow connections
       A1 --> B1
       A2 --> B1
       A3 --> C3
       A4 --> V3
       
       B1 --> B2
       B2 --> B3
       B3 --> B4
       B4 --> C1
       B4 --> C2
       
       C1 --> C3
       C2 --> C3
       C3 --> C4
       C4 --> D1
       C4 --> D2
       C4 --> D3
       
       D1 --> D4
       D2 --> D4
       D3 --> D4
       D4 --> E1
       
       E1 --> E2
       E2 --> E3
       E3 --> E4
       E4 --> G1
       
       %% Enhanced QC branch
       D4 -.-> F1
       F1 --> F2
       F2 --> F3
       F3 --> F4
       F4 --> H2
       
       %% Reporting and outputs
       G1 --> G2
       G2 --> G3
       G3 --> G4
       G4 --> H4
       
       D4 --> H1
       E4 --> H1
       C4 --> H3
       D4 --> H3
       F3 --> H3
       
       %% Alternative paths
       H4 -.-> I1
       B1 -.-> I2
       G4 -.-> I3
       H1 -.-> I4
       H2 -.-> I4
       
       %% Configuration influences
       V1 -.-> B2
       V1 -.-> F2
       V2 -.-> C3
       V2 -.-> E2
       V3 -.-> B1
       V3 -.-> H1
       V3 -.-> H2
       V4 -.-> F2
       V4 -.-> F3
       V4 -.-> F4
       
       %% Processing logs and diagnostics
       B1 --> H5
       C3 --> H5
       E2 --> H5
       F3 --> H5
       
       %% Styling
       style F1 fill:#fff3e0,stroke:#f57c00,stroke-width:3px
       style H2 fill:#e8f5e8,stroke:#388e3c,stroke-width:3px
       style E4 fill:#fce4ec,stroke:#c2185b,stroke-width:2px
       style G4 fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
       
       classDef inputNode fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
       classDef configNode fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
       classDef coreNode fill:#e8f5e8,stroke:#388e3c,stroke-width:2px
       classDef enhancedNode fill:#fff3e0,stroke:#f57c00,stroke-width:2px
       classDef outputNode fill:#fce4ec,stroke:#c2185b,stroke-width:2px
       classDef altNode fill:#f5f5f5,stroke:#616161,stroke-width:1px
       
       class A1,A2,A3,A4 inputNode
       class V1,V2,V3,V4 configNode
       class B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,G1,G2,G3,G4 coreNode
       class F1,F2,F3,F4 enhancedNode
       class H1,H2,H3,H4,H5 outputNode
       class I1,I2,I3,I4 altNode

Pipeline Stages Detailed
=========================

Stage 1: File Conversion & Initial QC
--------------------------------------

**Purpose**: Convert proprietary AB1 files to standard formats with initial quality assessment.

**Key Components**:

- **AB1 Converter**: Extracts DNA sequences and quality scores from ABI format
- **Quality Filtering**: Applies Phred score thresholds and length requirements
- **Format Conversion**: Produces FASTA/FASTQ outputs with preserved quality information
- **Initial QC**: Validates file integrity and sequence completeness

**Configuration Variables**:

- ``--min-quality``: Phred score threshold (15-30, default: 20)
- ``--min-length``: Minimum sequence length (30-100bp, default: 30)
- ``--quality-window``: Quality assessment window size

**Outputs**:

- Raw FASTA files (``fasta/`` directory)
- Quality score plots (``plots/`` directory)
- Processing logs

Stage 2: Sequence Processing
----------------------------

**Purpose**: Process forward and reverse sequences, generate alignments, and build consensus sequences.

**Key Components**:

- **Forward/Reverse Processing**: Independent quality trimming and artifact removal
- **Sequence Alignment**: MAFFT or MUSCLE alignment with parameter optimization
- **Consensus Generation**: Intelligent merging with conflict resolution

**Configuration Variables**:

- ``--alignment-tool``: Alignment software (mafft/muscle)
- ``--alignment-params``: Tool-specific parameters
- ``--consensus-threshold``: Minimum agreement for consensus calls

**Outputs**:

- Filtered sequences (``filtered/`` directory)
- Alignment files (intermediate)
- Consensus sequences per region (``consensus/`` directory)

Stage 3: Regional Analysis
---------------------------

**Purpose**: Process specific HVS regions and merge available regions per sample.

**Key Components**:

- **HVS1 Processing**: Mitochondrial positions 16024-16365
- **HVS2 Processing**: Mitochondrial positions 57-372  
- **HVS3 Processing**: Mitochondrial positions 438-574
- **Regional Merging**: Combines available regions into final sequences

**Configuration Variables**:

- ``--hvs-regions``: Specify which regions to process
- ``--region-overlap``: Handling of overlapping regions
- ``--merge-strategy``: Approach for combining regions

**Outputs**:

- Regional consensus files
- Merged sequences (``final/`` directory)
- Region coverage statistics

Stage 4: Damage Analysis
-------------------------

**Purpose**: Analyze ancient DNA damage patterns with statistical validation.

**Key Components**:

- **Damage Pattern Detection**: Identifies C→T and G→A transitions
- **Statistical Analysis**: Bootstrap validation with confidence intervals
- **Background Comparison**: Compares against modern DNA controls
- **Damage Scoring**: Generates composite scores and assessments

**Configuration Variables**:

- ``--damage-threshold``: Minimum damage level for significance
- ``--bootstrap-iterations``: Number of bootstrap samples (default: 1000)
- ``--modern-controls``: Reference modern DNA datasets

**Outputs**:

- Damage analysis results (``damage_analysis/`` directory)
- Statistical summaries (JSON format)
- Damage profile plots

Enhanced Quality Control (v2.0+)
---------------------------------

**Purpose**: Advanced aDNA-specific processing with enhanced quality control.

**Key Components**:

- **aDNA Sequence Cleaner**: Removes artifacts, resolves ambiguous bases
- **Improved HSD Converter**: Reference-aware variant calling with quality metrics
- **Diversity Analyzer**: Comprehensive haplogroup diversity assessment

**Configuration Variables**:

- ``--aggressive-cleaning``: Enable intensive artifact removal
- ``--reference-aware``: Use reference-guided processing
- ``--quality-filter``: Enhanced quality threshold (0.6-0.8)

**Outputs**:

- Cleaned sequences (``*_cleaned.fasta``)
- High-quality HSD files (``*_high_quality.hsd``)
- Diversity analysis reports

Stage 5: Quality Control & Reporting
-------------------------------------

**Purpose**: Generate comprehensive quality assessments and interactive reports.

**Key Components**:

- **Quality Metrics**: Sequence quality, coverage, success rates
- **Statistical Summaries**: Sample and batch-level statistics
- **Visualization**: Quality plots, damage profiles, interactive charts
- **Report Compilation**: HTML dashboards, PDF summaries, CSV exports

**Outputs**:

- Interactive HTML reports (``reports/`` directory)
- Quality visualization plots (``plots/`` directory)
- Statistical summary files (CSV/JSON)

Stage 6: Output Generation
---------------------------

**Purpose**: Produce final analysis outputs in multiple formats.

**Available Outputs**:

1. **Standard HSD Files**: Basic variant calling using regional or direct methods
2. **Enhanced HSD Files**: Quality-filtered variants with statistical confidence
3. **FASTA Sequences**: Raw, filtered, consensus, and final merged sequences
4. **Quality Reports**: Interactive dashboards and summary statistics
5. **Diagnostic Files**: Alignment files, logs, and debug information

Alternative Processing Paths
=============================

Manual Tool Access
-------------------

Access individual pipeline components for custom workflows:

.. code-block:: bash

   # Individual AB1 conversion
   python -m src.sanger_pipeline.cli.main convert-ab1 input.ab1 output.fasta
   
   # Manual damage analysis
   python -m src.sanger_pipeline.cli.main analyze-damage sequences/ results/
   
   # Standalone HSD conversion
   python -m src.sanger_pipeline.cli.main convert-to-hsd consensus/ output.hsd

Batch Processing
----------------

Process multiple samples efficiently:

.. code-block:: bash

   # Batch pipeline execution
   python scripts/batch_processor.py \
       --input-root ./samples/ \
       --output-root ./results/ \
       --parallel 4

Reprocessing Options
--------------------

Adjust parameters and reprocess selectively:

.. code-block:: bash

   # Reprocess with different quality threshold
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output_q25 \
       --min-quality 25
   
   # Regenerate reports only
   python generate_report.py ./existing_output/

Integration Endpoints
---------------------

Pipeline outputs compatible with external tools:

- **HaploGrep**: Direct HSD file upload
- **BEAST**: Sequence alignment formats
- **Custom Analysis**: CSV/JSON data exports
- **Database Systems**: Structured output formats

Configuration Reference
========================

Complete configuration file example:

.. code-block:: yaml

   # Complete pipeline configuration
   quality:
     min_phred_score: 20
     min_sequence_length: 30
     quality_window: 15
     quality_threshold: 0.7
   
   alignment:
     tool: "mafft"
     parameters: "--auto"
     consensus_threshold: 0.6
   
   hvs_regions:
     HVS1: {start: 16024, end: 16365}
     HVS2: {start: 57, end: 372}
     HVS3: {start: 438, end: 574}
   
   damage:
     damage_threshold: 0.02
     bootstrap_iterations: 1000
     significance_level: 0.05
   
   enhanced_qc:
     enabled: true
     aggressive_cleaning: false
     reference_aware: true
     quality_filter: 0.7
   
   output:
     generate_plots: true
     interactive_reports: true
     export_formats: ["hsd", "fasta", "csv"]

Performance Considerations
==========================

**Resource Requirements**:

- **Memory**: 2-8GB depending on dataset size
- **CPU**: Multi-core recommended for alignment steps
- **Storage**: 2-5x input size for intermediate files
- **Network**: Optional for reference downloads

**Optimization Strategies**:

- Use parallel processing for large datasets
- Adjust quality thresholds based on sample quality
- Enable caching for repeated analyses
- Configure temporary directory for large datasets

**Troubleshooting**:

- Monitor memory usage during alignment steps
- Check disk space for intermediate files
- Validate input file integrity
- Review log files for processing errors

Error Handling & Recovery
==========================

**Common Issues**:

1. **Input File Problems**: Corrupted AB1 files, missing files
2. **Quality Issues**: Low-quality sequences, insufficient coverage
3. **Alignment Failures**: Reference mismatches, parameter issues
4. **Resource Limitations**: Memory exhaustion, disk space

**Recovery Strategies**:

- Automatic retry with relaxed parameters
- Graceful degradation to available data
- Detailed error logging and reporting
- Checkpoint-based resumption

**Support Resources**:

- Comprehensive log analysis
- Interactive troubleshooting guide
- Community support forums
- Developer contact information

This comprehensive workflow reference provides complete coverage of all pipeline capabilities, configuration options, and processing pathways available in the Sanger aDNA damage analysis system.
