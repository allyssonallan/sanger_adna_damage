===============
Your First Analysis
===============

This tutorial walks you through running your first complete analysis with the Sanger DNA Damage Analysis Pipeline. By the end, you'll have processed AB1 files and generated a comprehensive QC report.

.. attention::
   **🎯 Understanding Tool Purpose**
   
   This analysis provides **preliminary screening results** to help you:
   
   * Identify samples with promising damage patterns for NGS follow-up
   * Assess sequence quality and insert sizes
   * Prioritize haplogroups for further investigation
   * Make informed decisions about resource allocation
   
   **This is NOT a substitute for proper NGS-based aDNA authentication.** Use these results to guide your NGS sample selection and experimental design.

🎯 Tutorial Goals
================

By completing this tutorial, you will:

* ✅ Run a complete pipeline analysis from start to finish
* ✅ Understand the output directory structure
* ✅ Generate and interpret an interactive QC report
* ✅ Assess ancient DNA damage patterns
* ✅ Know how to troubleshoot common issues

⏱️ **Estimated Time**: 15-20 minutes

📋 Prerequisites
===============

Before starting this tutorial:

* ✅ Pipeline installed and tested (see :doc:`../installation`)
* ✅ At least 2-4 AB1 files (forward and reverse reads for 1-2 samples)
* ✅ Basic familiarity with command line
* ✅ MAFFT installed and accessible

💾 Sample Data Setup
===================

If you don't have AB1 files, you can create a test scenario:

.. code-block:: bash

   # Create tutorial workspace
   mkdir sanger_tutorial
   cd sanger_tutorial
   
   # Create input directory
   mkdir input
   
   # Copy your AB1 files to input/ directory
   # Files should be named like: sample1_F.ab1, sample1_R.ab1

Expected File Structure
----------------------

Your input directory should look like:

.. code-block:: text

   input/
   ├── sample1_F.ab1      # Forward read for sample 1
   ├── sample1_R.ab1      # Reverse read for sample 1
   ├── sample2_F.ab1      # Forward read for sample 2 (optional)
   └── sample2_R.ab1      # Reverse read for sample 2 (optional)

.. note::
   The pipeline automatically detects forward (F) and reverse (R) reads based on filename patterns.

🚀 Step 1: Configuration Setup
==============================

Copy the default configuration to your working directory:

.. code-block:: bash

   # Copy default configuration
   cp /path/to/sanger_adna_damage/config/default_config.yaml ./tutorial_config.yaml

View and understand the configuration:

.. code-block:: bash

   # View configuration
   cat tutorial_config.yaml

You should see something like:

.. code-block:: yaml

   quality_threshold: 20
   min_sequence_length: 50
   damage_threshold: 0.05
   bootstrap_iterations: 10000
   
   hvs_regions:
     HVS1:
       start: 16024
       end: 16365
     HVS2:
       start: 57
       end: 372
     HVS3:
       start: 438
       end: 574

**Configuration Explanation**:

* ``quality_threshold: 20`` - Keep bases with Q20+ quality (99% accuracy)
* ``min_sequence_length: 50`` - Sequences must be at least 50bp after filtering
* ``damage_threshold: 0.05`` - P-value threshold for damage significance
* ``bootstrap_iterations: 10000`` - Number of statistical iterations

🔧 Step 2: Run the Complete Pipeline
====================================

Now let's run the complete analysis:

.. code-block:: bash

   # Run the complete pipeline
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output \
       --config ./tutorial_config.yaml

**What happens during processing**:

1. **AB1 Conversion**: Binary AB1 files converted to FASTA format
2. **Quality Filtering**: Low-quality bases and sequences removed
3. **Sequence Alignment**: Forward and reverse reads aligned using MAFFT
4. **Consensus Building**: Consensus sequences created for each HVS region
5. **HVS Merging**: Available HVS regions combined into final sequences
6. **Damage Analysis**: Ancient DNA damage patterns analyzed

**Expected Output**:

.. code-block:: text

   Starting Sanger pipeline...
   Processing AB1 files...
   ✓ Converted sample1_F.ab1 to FASTA
   ✓ Converted sample1_R.ab1 to FASTA
   ✓ Quality filtering completed
   ✓ Sequence alignment completed
   ✓ Consensus sequences generated
   ✓ HVS regions merged
   ✓ Damage analysis completed
   Pipeline completed successfully!

🗂️ Step 3: Explore the Output Structure
=======================================

After successful completion, examine your output directory:

.. code-block:: bash

   # Explore the output structure
   tree output/

You should see:

.. code-block:: text

   output/
   ├── fasta/                    # Raw FASTA conversions
   │   ├── sample1_F.fasta
   │   ├── sample1_R.fasta
   │   └── ...
   ├── filtered/                 # Quality-filtered sequences
   │   ├── sample1_F_filtered.fasta
   │   ├── sample1_R_filtered.fasta
   │   └── ...
   ├── consensus/                # Consensus by HVS region
   │   ├── sample1_HVS1_consensus.fasta
   │   ├── sample1_HVS2_consensus.fasta
   │   ├── sample1_HVS3_consensus.fasta
   │   └── ...
   ├── aligned/                  # Intermediate alignments
   │   ├── sample1_aligned.fasta
   │   └── ...
   ├── final/                    # Final merged sequences
   │   ├── sample1_final.fasta
   │   └── ...
   ├── damage_analysis/          # Ancient DNA damage analysis
   │   ├── sample1_damage_analysis.json
   │   └── ...
   └── plots/                    # Quality visualizations
       ├── sample1_F_quality.png
       └── ...

**Directory Explanations**:

* **fasta/**: Raw conversions from AB1 format
* **filtered/**: Quality-filtered sequences ready for analysis
* **consensus/**: Consensus sequences for each HVS region independently
* **final/**: Your final processed sequences (main results)
* **damage_analysis/**: Ancient DNA damage assessment results
* **plots/**: Quality score visualizations

📊 Step 4: Generate Interactive QC Report
=========================================

Create a comprehensive QC report with visualizations:

.. code-block:: bash

   # Generate interactive QC report
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./output \
       --open-browser

This command will:

1. Analyze all pipeline outputs
2. Generate statistical summaries
3. Create interactive visualizations
4. Open the report in your default browser

**Expected Browser Output**:

The report opens with several tabs:

* **Overview**: Processing summary and key metrics
* **Damage Analysis**: Ancient DNA damage assessment
* **Quality Control**: Sequence quality distributions
* **Sample Details**: Per-sample detailed results

🔍 Step 5: Interpret Your Results
=================================

Overview Tab Analysis
--------------------

Look for these key metrics:

.. code-block:: text

   ✓ Samples Processed: 2/2 (100%)
   ✓ Average Quality Score: 28.5
   ✓ HVS Regions Detected: HVS1, HVS2, HVS3
   ✓ Total Sequences: 4 (2 samples × 2 reads)

**Good indicators**:
* High success rate (>90%)
* Quality scores >20
* Multiple HVS regions detected

Damage Analysis Tab
------------------

Key damage metrics to examine:

.. code-block:: text

   Damage Score: 0.23 (Low-Moderate)
   P-value: 0.045 (Significant)
   C→T Transitions: 12%
   G→A Transitions: 8%

**Interpretation**:
* **Damage Score 0-0.3**: Low damage (modern DNA or well-preserved)
* **Damage Score 0.3-0.7**: Moderate damage (possible ancient DNA)
* **Damage Score >0.7**: High damage (likely ancient DNA)
* **P-value <0.05**: Statistically significant damage pattern

Quality Control Tab
------------------

Examine quality distributions:

* **Mean Quality**: Should be >20 for reliable analysis
* **Length Distribution**: Check if sequences meet minimum length
* **Processing Efficiency**: High success rates indicate good data

Sample Details Tab
------------------

Per-sample breakdown shows:

* Individual quality metrics
* HVS region coverage
* Damage assessment for each sample
* File processing status

📈 Step 6: Examine Individual Results
====================================

Look at specific output files:

Final Sequences
--------------

.. code-block:: bash

   # View your final processed sequences
   cat output/final/sample1_final.fasta

Example output:

.. code-block:: text

   >sample1_HVS1_HVS2_final
   GATTTCACGGAGGATGGTGGTCAAGGGACCCCCCCTCCCCCATGCTTACAAGCAAGTACA...

Damage Analysis Results
----------------------

.. code-block:: bash

   # View damage analysis (formatted JSON)
   python -c "
   import json
   with open('output/damage_analysis/sample1_damage_analysis.json') as f:
       data = json.load(f)
       print(json.dumps(data, indent=2))
   "

Example output:

.. code-block:: json

   {
     "sample_id": "sample1",
     "damage_score": 0.23,
     "p_value": 0.045,
     "c_to_t_rate": 0.12,
     "g_to_a_rate": 0.08,
     "assessment": "Low-Moderate damage detected",
     "significance": "statistically_significant"
   }

Quality Plots
------------

.. code-block:: bash

   # View quality plots (if you have image viewer)
   open output/plots/sample1_F_quality.png  # macOS
   # or
   xdg-open output/plots/sample1_F_quality.png  # Linux

🎯 Step 7: Understanding Your Results
====================================

Success Indicators
-----------------

Your analysis was successful if:

* ✅ All AB1 files were converted to FASTA
* ✅ Quality filtering produced sequences >50bp
* ✅ At least one HVS region was successfully processed
* ✅ Final sequences were generated
* ✅ QC report generated without errors

Ancient DNA Assessment
---------------------

Based on damage analysis results:

**Modern DNA Indicators**:
* Damage score <0.2
* P-value >0.05 (not significant)
* Low C→T and G→A rates (<5%)

**Ancient DNA Indicators**:
* Damage score >0.3
* P-value <0.05 (significant)
* Elevated C→T and G→A rates (>10%)

**Borderline Cases**:
* Damage score 0.2-0.3
* May need additional validation
* Consider sample preservation conditions

🔧 Step 8: Check Pipeline Status
===============================

Get a summary of pipeline results:

.. code-block:: bash

   # Check overall pipeline status
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./output

Expected output:

.. code-block:: text

   Pipeline Status Report
   =====================
   
   Input Files: 4 AB1 files detected
   ✓ FASTA Conversion: 4/4 successful
   ✓ Quality Filtering: 4/4 passed
   ✓ Consensus Building: 6/6 HVS regions processed
   ✓ Final Sequences: 2/2 samples completed
   ✓ Damage Analysis: 2/2 samples analyzed
   
   HVS Region Coverage:
   - HVS1: 100% (2/2 samples)
   - HVS2: 100% (2/2 samples)  
   - HVS3: 50% (1/2 samples)
   
   Quality Summary:
   - Average Quality Score: 28.5
   - Average Sequence Length: 245bp
   - Overall Success Rate: 100%

🚨 Troubleshooting Common Issues
===============================

Issue 1: No AB1 Files Found
---------------------------

**Error**: ``No AB1 files found in input directory``

**Solution**:

.. code-block:: bash

   # Check file extensions and naming
   ls -la input/
   
   # Ensure files have .ab1 extension
   # Rename if necessary:
   mv sample1_forward.AB1 sample1_F.ab1

Issue 2: Quality Filtering Removes All Sequences
-----------------------------------------------

**Error**: ``No sequences passed quality filtering``

**Solution**:

.. code-block:: bash

   # Lower quality threshold temporarily
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output_lowq \
       --quality-threshold 15

Issue 3: MAFFT Not Found
------------------------

**Error**: ``MAFFT executable not found``

**Solution**:

.. code-block:: bash

   # Check MAFFT installation
   mafft --version
   
   # Install if missing
   # macOS: brew install mafft
   # Ubuntu: sudo apt install mafft

Issue 4: Empty Final Sequences
------------------------------

**Problem**: Final sequences are very short or missing

**Diagnosis**:

.. code-block:: bash

   # Check intermediate results
   ls -la output/filtered/
   cat output/filtered/sample1_F_filtered.fasta

**Solution**: Lower quality threshold or check input data quality

🎉 Congratulations!
==================

You've successfully completed your first Sanger DNA analysis! You now have:

* ✅ Processed AB1 files into high-quality sequences
* ✅ Generated consensus sequences for HVS regions
* ✅ Assessed ancient DNA damage patterns
* ✅ Created a comprehensive QC report
* ✅ Learned to interpret key results

🎯 Next Steps
=============

Now that you've completed your first analysis:

1. **Explore Advanced Features**: Try :doc:`ancient_dna_workflow` for specialized ancient DNA analysis
2. **Customize Configuration**: Learn about :doc:`../configuration` options
3. **Batch Processing**: Process multiple samples with :doc:`batch_processing`
4. **Understand Damage Analysis**: Deep dive into :doc:`damage_assessment`

🔄 Practice Exercises
====================

To reinforce your learning:

1. **Try Different Quality Thresholds**: Re-run with quality thresholds of 15, 25, and 30
2. **Analyze Different Sample Types**: Process both modern and potentially ancient samples
3. **Compare Results**: Run the same samples with different configurations
4. **Explore Command Options**: Try different CLI commands and options

💡 Key Takeaways
================

* The pipeline processes AB1 files through multiple quality-controlled stages
* Interactive QC reports provide comprehensive analysis summaries
* Damage analysis helps identify ancient DNA patterns
* Configuration files allow customization for different sample types
* Each processing stage has specific outputs that can be examined individually

You're now ready to tackle more complex analyses and explore the advanced features of the Sanger DNA Damage Analysis Pipeline!
