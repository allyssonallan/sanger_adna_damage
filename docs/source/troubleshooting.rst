===============
Troubleshooting
===============

Common issues, solutions, and debugging strategies for the Sanger DNA Damage Analysis Pipeline.

🚨 Quick Diagnostic Commands
============================

Before diving into specific issues, run these diagnostic commands:

.. code-block:: bash

   # Check pipeline installation
   python -c "from src.sanger_pipeline.core.pipeline import SangerPipeline; print('✓ Pipeline imported successfully')"
   
   # Check external dependencies
   mafft --version
   
   # Validate your configuration
   python -m src.sanger_pipeline.cli.main validate --config your_config.yaml --check-deps
   
   # Check input files
   python -m src.sanger_pipeline.cli.main validate --check-input ./your_input_dir

🔧 Installation and Setup Issues
=================================

Pipeline Won't Import
---------------------

**Error**: ``ModuleNotFoundError: No module named 'src.sanger_pipeline'``

**Solutions**:

1. **Check Installation**:

   .. code-block:: bash

      # Verify you're in the right directory
      ls -la src/sanger_pipeline/
      
      # Install in development mode
      pip install -e .

2. **Python Path Issues**:

   .. code-block:: bash

      # Add to Python path temporarily
      export PYTHONPATH="${PYTHONPATH}:$(pwd)"
      
      # Or use absolute imports
      python -m src.sanger_pipeline.cli.main --help

3. **Virtual Environment**:

   .. code-block:: bash

      # Make sure virtual environment is activated
      which python
      pip list | grep -i bio  # Check if BioPython is installed

MAFFT Not Found
---------------

**Error**: ``MAFFT executable not found`` or ``Command 'mafft' not found``

**Solutions**:

**macOS**:

.. code-block:: bash

   # Install via Homebrew
   brew install mafft
   
   # Verify installation
   mafft --version

**Ubuntu/Debian**:

.. code-block:: bash

   # Update package list and install
   sudo apt update
   sudo apt install mafft
   
   # Verify installation
   which mafft

**Manual Installation**:

.. code-block:: bash

   # Download and install manually
   wget https://mafft.cbrc.jp/alignment/software/mafft-7.490-without-extensions-src.tgz
   tar -xzf mafft-7.490-without-extensions-src.tgz
   cd mafft-7.490-without-extensions/core/
   make clean
   make
   sudo make install

**Path Issues**:

.. code-block:: bash

   # Check if MAFFT is in PATH
   echo $PATH
   
   # Add MAFFT directory to PATH
   export PATH="/usr/local/bin:$PATH"
   
   # Make permanent in ~/.bashrc or ~/.zshrc
   echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.bashrc

BioPython Issues
---------------

**Error**: ``ImportError: No module named 'Bio'``

**Solutions**:

.. code-block:: bash

   # Install BioPython
   pip install biopython
   
   # If installation fails, try with conda
   conda install -c conda-forge biopython
   
   # Verify installation
   python -c "import Bio; print(f'BioPython version: {Bio.__version__}')"

**Common BioPython Installation Issues**:

.. code-block:: bash

   # If you get compiler errors on macOS
   xcode-select --install
   
   # If you get permission errors
   pip install --user biopython

📁 File and Data Issues
=======================

No AB1 Files Found
------------------

**Error**: ``No AB1 files found in input directory``

**Diagnostic Steps**:

.. code-block:: bash

   # Check file extensions
   ls -la input/ | grep -i ab1
   
   # Check for hidden characters or wrong extensions
   file input/*

**Solutions**:

1. **Check File Extensions**:

   .. code-block:: bash

      # Rename files with wrong extensions
      for file in input/*.AB1; do
          mv "$file" "${file%.AB1}.ab1"
      done

2. **Check File Permissions**:

   .. code-block:: bash

      # Fix permissions
      chmod 644 input/*.ab1

3. **Verify File Format**:

   .. code-block:: bash

      # Check if files are actually AB1 format
      python -c "
      from Bio import SeqIO
      try:
          record = SeqIO.read('input/sample.ab1', 'abi')
          print('✓ File is valid AB1 format')
      except:
          print('✗ File is not valid AB1 format')
      "

Corrupted AB1 Files
------------------

**Error**: ``Error reading AB1 file`` or ``Invalid trace file format``

**Diagnostic**:

.. code-block:: bash

   # Check file size and integrity
   ls -lh input/*.ab1
   
   # Test with BioPython
   python -c "
   from Bio import SeqIO
   import sys
   try:
       record = SeqIO.read(sys.argv[1], 'abi')
       print(f'Sequence length: {len(record.seq)}')
       print(f'Quality scores available: {hasattr(record, \"letter_annotations\")}')
   except Exception as e:
       print(f'Error: {e}')
   " input/sample.ab1

**Solutions**:

1. **Re-download Files**: Get fresh copies from sequencing provider
2. **Check Transfer Method**: Ensure files weren't corrupted during transfer
3. **Alternative Conversion**: Try manual conversion with other tools

File Pairing Issues
-------------------

**Error**: ``Cannot pair forward and reverse reads``

**Diagnostic**:

.. code-block:: bash

   # Check naming patterns
   ls -la input/ | grep -E "_[FR]\.ab1|_[12]\.ab1|_(forward|reverse)\.ab1"

**Solutions**:

1. **Fix Naming Convention**:

   .. code-block:: bash

      # Rename to standard pattern
      mv sample_forward.ab1 sample_F.ab1
      mv sample_reverse.ab1 sample_R.ab1

2. **Manual Pairing**: Edit configuration to specify custom pairing rules

⚙️ Processing and Analysis Issues
=================================

Quality Filtering Removes All Sequences
---------------------------------------

**Error**: ``No sequences passed quality filtering``

**Diagnostic**:

.. code-block:: bash

   # Check original sequence quality
   python -c "
   from Bio import SeqIO
   record = SeqIO.read('input/sample_F.ab1', 'abi')
   qualities = record.letter_annotations['phred_quality']
   print(f'Average quality: {sum(qualities)/len(qualities):.1f}')
   print(f'Min quality: {min(qualities)}')
   print(f'Max quality: {max(qualities)}')
   "

**Solutions**:

1. **Lower Quality Threshold**:

   .. code-block:: bash

      # Try with lower threshold
      python -m src.sanger_pipeline.cli.main run-pipeline \
          --input-dir ./input \
          --output-dir ./output_lowq \
          --quality-threshold 10

2. **Check Input Quality**:

   .. code-block:: bash

      # Generate quality plots first
      python -m src.sanger_pipeline.cli.main convert \
          --input-dir ./input \
          --output-dir ./converted \
          --quality-filter

3. **Adjust Minimum Length**:

   .. code-block:: yaml

      # In configuration file
      quality_threshold: 15
      min_sequence_length: 25

Alignment Failures
-----------------

**Error**: ``MAFFT alignment failed`` or ``Empty alignment result``

**Diagnostic**:

.. code-block:: bash

   # Test MAFFT manually
   mafft --version
   mafft --auto test_sequences.fasta

**Solutions**:

1. **Check Sequence Compatibility**:

   .. code-block:: bash

      # Ensure sequences are from same organism/region
      head -20 output/filtered/sample*_filtered.fasta

2. **Manual Alignment Test**:

   .. code-block:: bash

      # Test alignment manually
      cat output/filtered/sample_F_filtered.fasta output/filtered/sample_R_filtered.fasta > test_align.fasta
      mafft test_align.fasta > aligned_test.fasta

3. **Alternative Alignment Parameters**:

   Modify pipeline to use different MAFFT parameters for difficult sequences.

No HVS Regions Detected
-----------------------

**Error**: ``No HVS regions found in consensus sequence``

**Diagnostic**:

.. code-block:: bash

   # Check sequence content and length
   python -c "
   from Bio import SeqIO
   for record in SeqIO.parse('output/consensus/sample_consensus.fasta', 'fasta'):
       print(f'Sequence ID: {record.id}')
       print(f'Length: {len(record.seq)}')
       print(f'First 100 bp: {record.seq[:100]}')
   "

**Solutions**:

1. **Check Sequence Type**:
   - Ensure sequences are mitochondrial DNA
   - Verify they contain hypervariable regions

2. **Adjust HVS Coordinates**:

   .. code-block:: yaml

      # Modify in configuration
      hvs_regions:
        HVS1:
          start: 16000  # More permissive coordinates
          end: 16400

3. **Manual Region Identification**:
   Use BLAST or other tools to identify actual sequence content.

🧬 Ancient DNA Analysis Issues
=============================

Bootstrap Analysis Fails
------------------------

**Error**: ``Bootstrap analysis failed`` or ``Insufficient data for bootstrap``

**Diagnostic**:

.. code-block:: bash

   # Check sequence count and length
   python -c "
   from Bio import SeqIO
   sequences = list(SeqIO.parse('output/final/sample_final.fasta', 'fasta'))
   print(f'Number of sequences: {len(sequences)}')
   for i, seq in enumerate(sequences):
       print(f'Sequence {i+1}: {len(seq.seq)} bp')
   "

**Solutions**:

1. **Reduce Bootstrap Iterations**:

   .. code-block:: yaml

      # In configuration
      bootstrap_iterations: 1000  # Reduced from 10000

2. **Check Minimum Requirements**:
   - Sequences should be >50bp
   - Need sufficient sequence data for statistical analysis

3. **Alternative Analysis**:

   .. code-block:: bash

      # Run damage analysis with different parameters
      python -m src.sanger_pipeline.cli.main analyze-damage \
          --input-dir ./output/final \
          --output-dir ./damage_simple \
          --iterations 1000

Unrealistic Damage Scores
-------------------------

**Issue**: Damage scores don't match expected values

**Diagnostic**:

.. code-block:: bash

   # Examine damage analysis details
   python -c "
   import json
   with open('output/damage_analysis/sample_damage_analysis.json') as f:
       data = json.load(f)
   print('Damage score:', data['damage_score'])
   print('C→T rate:', data['c_to_t_rate'])
   print('G→A rate:', data['g_to_a_rate'])
   print('Background rate:', data['background_rate'])
   "

**Solutions**:

1. **Check Reference Expectations**:
   - Modern DNA: damage score <0.2
   - Ancient DNA: damage score >0.3
   - Consider sample age and preservation

2. **Validate with Controls**:
   - Include known modern samples
   - Include known ancient samples
   - Compare patterns across samples

3. **Review Sample History**:
   - Check extraction methods
   - Review storage conditions
   - Consider contamination sources

📊 Report Generation Issues
===========================

Report Generation Fails
-----------------------

**Error**: ``Failed to generate QC report``

**Diagnostic**:

.. code-block:: bash

   # Check output directory structure
   tree output/
   
   # Test report generation with verbose output
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./output \
       --verbose

**Solutions**:

1. **Check Dependencies**:

   .. code-block:: bash

      # Ensure all required files exist
      ls -la output/damage_analysis/
      ls -la output/final/

2. **Permissions**:

   .. code-block:: bash

      # Check write permissions
      touch output/reports/test_file.html
      rm output/reports/test_file.html

3. **Manual Report Generation**:

   .. code-block:: python

      # Test report generation in Python
      from src.sanger_pipeline.utils.report_generator import ReportGenerator
      
      generator = ReportGenerator()
      report_path = generator.generate_report("./output")

Browser Won't Open Report
-------------------------

**Issue**: Report generates but browser doesn't open

**Solutions**:

.. code-block:: bash

   # Open report manually
   open output/reports/qc_report_*.html  # macOS
   xdg-open output/reports/qc_report_*.html  # Linux
   
   # Or specify browser
   firefox output/reports/qc_report_*.html

🐛 Performance and Memory Issues
===============================

Slow Processing
--------------

**Issue**: Pipeline takes very long to run

**Diagnostic**:

.. code-block:: bash

   # Monitor resource usage
   top -p $(pgrep python)  # Linux
   # or
   Activity Monitor  # macOS

**Solutions**:

1. **Reduce Bootstrap Iterations**:

   .. code-block:: yaml

      bootstrap_iterations: 1000  # Instead of 10000

2. **Process Smaller Batches**:

   .. code-block:: bash

      # Split large datasets
      mkdir batch1 batch2
      mv input/sample1*.ab1 batch1/
      mv input/sample2*.ab1 batch2/

3. **Increase Quality Threshold**:

   .. code-block:: yaml

      quality_threshold: 25  # Higher threshold = less data to process

Memory Issues
------------

**Error**: ``MemoryError`` or system becomes unresponsive

**Solutions**:

1. **Process Individual Samples**:

   .. code-block:: bash

      # Process one sample at a time
      for sample in input/*_F.ab1; do
          base=$(basename "$sample" _F.ab1)
          mkdir "temp_$base"
          cp "input/${base}_F.ab1" "input/${base}_R.ab1" "temp_$base/"
          python -m src.sanger_pipeline.cli.main run-pipeline \
              --input-dir "temp_$base" \
              --output-dir "output_$base"
      done

2. **Adjust System Limits**:

   .. code-block:: bash

      # Increase memory limits (Linux)
      ulimit -v 4000000  # 4GB virtual memory limit

🔍 Advanced Debugging
====================

Enable Debug Logging
--------------------

.. code-block:: bash

   # Run with maximum verbosity
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output \
       --verbose \
       --log-file debug.log

   # Examine log file
   tail -f debug.log

Python Debugging
----------------

.. code-block:: python

   # Debug pipeline programmatically
   import logging
   logging.basicConfig(level=logging.DEBUG)
   
   from src.sanger_pipeline.core.pipeline import SangerPipeline
   
   pipeline = SangerPipeline(debug=True)
   try:
       results = pipeline.run("./input", "./output")
   except Exception as e:
       import traceback
       traceback.print_exc()

Isolate Problem Stage
--------------------

.. code-block:: bash

   # Test each stage individually
   
   # 1. Test AB1 conversion only
   python -m src.sanger_pipeline.cli.main convert \
       --input-dir ./input \
       --output-dir ./test_convert
   
   # 2. Test damage analysis only
   python -m src.sanger_pipeline.cli.main analyze-damage \
       --input-dir ./existing_sequences \
       --output-dir ./test_damage

📞 Getting Help
===============

Before Asking for Help
----------------------

1. **Check this troubleshooting guide** for your specific issue
2. **Run diagnostic commands** to gather information
3. **Try simple solutions** first (restart, reinstall, etc.)
4. **Prepare detailed information** about your problem

Information to Include
---------------------

When reporting issues, include:

1. **Error messages** (complete text, not screenshots)
2. **Command that failed** (exact command line)
3. **Configuration file** (if using custom config)
4. **System information**:

   .. code-block:: bash

      # Gather system info
      python --version
      pip list | grep -i bio
      mafft --version
      uname -a  # Linux/macOS
      
5. **Sample data characteristics**:
   - Number of AB1 files
   - Approximate file sizes
   - Expected sample types (modern/ancient)

Where to Get Help
----------------

1. **GitHub Issues**: https://github.com/allyssonallan/sanger_adna_damage/issues
2. **Documentation**: Check all relevant documentation sections
3. **Community Discussions**: GitHub Discussions for general questions

Creating Good Bug Reports
-------------------------

.. code-block:: markdown

   ## Bug Report Template
   
   **Problem Description**: Brief description of what's wrong
   
   **Expected Behavior**: What should happen
   
   **Actual Behavior**: What actually happens
   
   **Error Message**: 
   ```
   Paste exact error message here
   ```
   
   **Command Used**:
   ```bash
   python -m src.sanger_pipeline.cli.main run-pipeline --input-dir ./input --output-dir ./output
   ```
   
   **Environment**:
   - OS: macOS 12.0 / Ubuntu 20.04 / Windows 11
   - Python version: 3.9.7
   - Pipeline version: 1.0.0
   - MAFFT version: 7.490
   
   **Configuration** (if using custom config):
   ```yaml
   quality_threshold: 20
   # ... other relevant settings
   ```
   
   **Sample Data**: 
   - Number of AB1 files: 4
   - Sample types: Ancient DNA from bone
   
   **Additional Context**: Any other relevant information

This comprehensive troubleshooting guide should help you diagnose and solve most common issues with the Sanger DNA Damage Analysis Pipeline.
