=============
CLI Reference
=============

Complete reference for the Sanger DNA Damage Analysis Pipeline command-line interface.

🖥️ Overview
===========

The pipeline provides a comprehensive command-line interface (CLI) for all operations. All commands are accessed through the main CLI module:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main [COMMAND] [OPTIONS]

📋 Available Commands
====================

.. contents::
   :local:
   :depth: 2

run-pipeline
============

Run the complete analysis pipeline from AB1 files to final results.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main run-pipeline [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--input-dir``
     - PATH
     - **Required**. Directory containing AB1 files
   * - ``--output-dir`` 
     - PATH
     - **Required**. Directory for pipeline outputs
   * - ``--config``
     - PATH
     - Configuration file path (default: config/default_config.yaml)
   * - ``--quality-threshold``
     - INTEGER
     - Override quality threshold from config
   * - ``--min-length``
     - INTEGER
     - Override minimum sequence length from config
   * - ``--force``
     - FLAG
     - Overwrite existing output directory
   * - ``--dry-run``
     - FLAG
     - Show what would be processed without running
   * - ``--verbose``
     - FLAG
     - Enable verbose logging
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Basic usage
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./ab1_files \
       --output-dir ./results

   # With custom configuration
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./ab1_files \
       --output-dir ./results \
       --config ./custom_config.yaml

   # Override quality threshold
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./ab1_files \
       --output-dir ./results \
       --quality-threshold 25

   # Dry run to check what will be processed
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./ab1_files \
       --output-dir ./results \
       --dry-run

**Return Codes**:

* ``0``: Success
* ``1``: General error
* ``2``: Invalid arguments
* ``3``: Input files not found
* ``4``: Configuration error

generate-report
===============

Generate interactive HTML QC reports from pipeline outputs.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main generate-report [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--output-dir``
     - PATH
     - **Required**. Directory containing pipeline outputs
   * - ``--report-dir``
     - PATH
     - Directory for report files (default: output-dir/reports)
   * - ``--open-browser``
     - FLAG
     - Open report in browser after generation
   * - ``--title``
     - TEXT
     - Custom report title
   * - ``--template``
     - PATH
     - Custom HTML template file
   * - ``--format``
     - CHOICE
     - Report format: html, json (default: html)
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Generate report and open in browser
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./results \
       --open-browser

   # Custom report location and title
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./results \
       --report-dir ./custom_reports \
       --title "Ancient DNA Analysis Report"

   # JSON format for programmatic access
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./results \
       --format json

analyze-damage
==============

Perform ancient DNA damage analysis on sequences.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main analyze-damage [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--input-dir``
     - PATH
     - **Required**. Directory containing FASTA sequences
   * - ``--output-dir``
     - PATH
     - **Required**. Directory for damage analysis results
   * - ``--config``
     - PATH
     - Configuration file path
   * - ``--threshold``
     - FLOAT
     - P-value threshold for significance (0.0-1.0)
   * - ``--iterations``
     - INTEGER
     - Bootstrap iterations (1000-100000)
   * - ``--reference``
     - PATH
     - Reference sequence file
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Basic damage analysis
   python -m src.sanger_pipeline.cli.main analyze-damage \
       --input-dir ./results/final \
       --output-dir ./damage_results

   # With custom parameters
   python -m src.sanger_pipeline.cli.main analyze-damage \
       --input-dir ./results/final \
       --output-dir ./damage_results \
       --threshold 0.01 \
       --iterations 50000

status
======

Check pipeline status and output summary.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main status [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--output-dir``
     - PATH
     - Pipeline output directory to check
   * - ``--input-dir``
     - PATH
     - Original input directory
   * - ``--config``
     - PATH
     - Configuration file used
   * - ``--detailed``
     - FLAG
     - Show detailed per-file status
   * - ``--json``
     - FLAG
     - Output status in JSON format
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Basic status check
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./results

   # Detailed status with original inputs
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./results \
       --input-dir ./ab1_files \
       --detailed

   # JSON output for scripts
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./results \
       --json

validate
========

Validate configuration files and check system requirements.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main validate [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--config``
     - PATH
     - Configuration file to validate
   * - ``--check-deps``
     - FLAG
     - Check external dependencies (MAFFT, etc.)
   * - ``--check-input``
     - PATH
     - Validate input directory
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Validate configuration
   python -m src.sanger_pipeline.cli.main validate \
       --config ./my_config.yaml

   # Check all dependencies
   python -m src.sanger_pipeline.cli.main validate \
       --check-deps

   # Validate input directory
   python -m src.sanger_pipeline.cli.main validate \
       --check-input ./ab1_files

convert
=======

Convert AB1 files to FASTA format only.

**Syntax**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main convert [OPTIONS]

**Options**:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--input-dir``
     - PATH
     - **Required**. Directory containing AB1 files
   * - ``--output-dir``
     - PATH
     - **Required**. Directory for FASTA outputs
   * - ``--quality-filter``
     - FLAG
     - Apply quality filtering during conversion
   * - ``--quality-threshold``
     - INTEGER
     - Quality threshold for filtering (default: 20)
   * - ``--help``
     - FLAG
     - Show help message

**Examples**:

.. code-block:: bash

   # Simple conversion
   python -m src.sanger_pipeline.cli.main convert \
       --input-dir ./ab1_files \
       --output-dir ./fasta_files

   # With quality filtering
   python -m src.sanger_pipeline.cli.main convert \
       --input-dir ./ab1_files \
       --output-dir ./fasta_files \
       --quality-filter \
       --quality-threshold 25

🔧 Global Options
================

These options work with all commands:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``--version``
     - FLAG
     - Show pipeline version
   * - ``--help``
     - FLAG
     - Show help for command
   * - ``--verbose``
     - FLAG
     - Enable verbose output
   * - ``--quiet``
     - FLAG
     - Suppress non-error output
   * - ``--log-file``
     - PATH
     - Write logs to file
   * - ``--config-help``
     - FLAG
     - Show configuration help

**Examples**:

.. code-block:: bash

   # Check version
   python -m src.sanger_pipeline.cli.main --version

   # Get help for any command
   python -m src.sanger_pipeline.cli.main run-pipeline --help

   # Verbose logging to file
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --verbose \
       --log-file ./pipeline.log \
       --input-dir ./input \
       --output-dir ./output

📝 Configuration via CLI
========================

Many configuration parameters can be overridden via command line:

**Quality Control Overrides**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main run-pipeline \
       --quality-threshold 25 \
       --min-length 75 \
       --input-dir ./input \
       --output-dir ./output

**Damage Analysis Overrides**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main analyze-damage \
       --threshold 0.01 \
       --iterations 50000 \
       --input-dir ./sequences \
       --output-dir ./damage

🔄 Chaining Commands
===================

Commands can be chained for custom workflows:

.. code-block:: bash

   # Step-by-step processing
   
   # 1. Convert AB1 to FASTA
   python -m src.sanger_pipeline.cli.main convert \
       --input-dir ./ab1_files \
       --output-dir ./fasta_files
   
   # 2. Run full pipeline on converted files
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./ab1_files \
       --output-dir ./results
   
   # 3. Generate report
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./results \
       --open-browser
   
   # 4. Check status
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./results

📊 Exit Codes
=============

All commands return standard exit codes:

.. list-table::
   :widths: 10 90
   :header-rows: 1

   * - Code
     - Meaning
   * - ``0``
     - Success - command completed without errors
   * - ``1``
     - General error - something went wrong during execution
   * - ``2``
     - Invalid arguments - check command syntax and options
   * - ``3``
     - Input error - files not found or invalid input data
   * - ``4``
     - Configuration error - invalid or missing configuration
   * - ``5``
     - Dependency error - external tools not found or not working
   * - ``6``
     - Output error - cannot write to output directory

**Using exit codes in scripts**:

.. code-block:: bash

   #!/bin/bash
   
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output
   
   if [ $? -eq 0 ]; then
       echo "Pipeline completed successfully"
       python -m src.sanger_pipeline.cli.main generate-report \
           --output-dir ./output \
           --open-browser
   else
       echo "Pipeline failed with exit code $?"
       exit 1
   fi

🌍 Environment Variables
=======================

The CLI respects several environment variables:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Variable
     - Description
   * - ``SANGER_CONFIG``
     - Default configuration file path
   * - ``SANGER_OUTPUT_DIR``
     - Default output directory
   * - ``SANGER_QUALITY_THRESHOLD``
     - Default quality threshold
   * - ``TMPDIR``
     - Temporary directory for processing
   * - ``MAFFT_BINDIR``
     - MAFFT installation directory

**Using environment variables**:

.. code-block:: bash

   # Set default configuration
   export SANGER_CONFIG=/path/to/my/config.yaml
   
   # Set default output location
   export SANGER_OUTPUT_DIR=/data/sanger_results
   
   # Run with environment defaults
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input

🔍 Debugging and Troubleshooting
================================

**Enable verbose output**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main run-pipeline \
       --verbose \
       --input-dir ./input \
       --output-dir ./output

**Save logs to file**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main run-pipeline \
       --log-file ./debug.log \
       --input-dir ./input \
       --output-dir ./output

**Dry run to check inputs**:

.. code-block:: bash

   python -m src.sanger_pipeline.cli.main run-pipeline \
       --dry-run \
       --input-dir ./input \
       --output-dir ./output

**Validate before running**:

.. code-block:: bash

   # Check configuration
   python -m src.sanger_pipeline.cli.main validate \
       --config ./my_config.yaml \
       --check-deps \
       --check-input ./input

📝 Scripting Examples
====================

**Bash script for automated processing**:

.. code-block:: bash

   #!/bin/bash
   # automated_analysis.sh
   
   INPUT_DIR="$1"
   OUTPUT_DIR="$2"
   CONFIG_FILE="${3:-config/default_config.yaml}"
   
   # Validate inputs
   python -m src.sanger_pipeline.cli.main validate \
       --config "$CONFIG_FILE" \
       --check-input "$INPUT_DIR"
   
   if [ $? -ne 0 ]; then
       echo "Validation failed"
       exit 1
   fi
   
   # Run pipeline
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir "$INPUT_DIR" \
       --output-dir "$OUTPUT_DIR" \
       --config "$CONFIG_FILE" \
       --verbose
   
   # Generate report if pipeline succeeded
   if [ $? -eq 0 ]; then
       python -m src.sanger_pipeline.cli.main generate-report \
           --output-dir "$OUTPUT_DIR" \
           --open-browser
       
       # Show final status
       python -m src.sanger_pipeline.cli.main status \
           --output-dir "$OUTPUT_DIR" \
           --detailed
   fi

**Python script for batch processing**:

.. code-block:: python

   #!/usr/bin/env python3
   # batch_process.py
   
   import subprocess
   import sys
   from pathlib import Path
   
   def run_pipeline(input_dir, output_dir, config_file):
       """Run pipeline with error handling"""
       cmd = [
           sys.executable, "-m", "src.sanger_pipeline.cli.main", 
           "run-pipeline",
           "--input-dir", str(input_dir),
           "--output-dir", str(output_dir),
           "--config", str(config_file)
       ]
       
       result = subprocess.run(cmd, capture_output=True, text=True)
       
       if result.returncode == 0:
           print(f"✓ Successfully processed {input_dir}")
           return True
       else:
           print(f"✗ Failed to process {input_dir}: {result.stderr}")
           return False
   
   # Process multiple directories
   base_dir = Path("./samples")
   output_base = Path("./results")
   config = Path("./config/default_config.yaml")
   
   for sample_dir in base_dir.iterdir():
       if sample_dir.is_dir():
           output_dir = output_base / sample_dir.name
           run_pipeline(sample_dir, output_dir, config)

This comprehensive CLI reference covers all available commands and options for the Sanger DNA Damage Analysis Pipeline. Use it as a quick reference while working with the pipeline, or for developing automated workflows and scripts.
