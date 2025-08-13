===================
Process Single Sample
===================

**Goal**: Process a single AB1 sample pair (forward and reverse reads) through the complete pipeline to generate final consensus sequences and damage assessment.

**Level**: 🟢 Beginner

**Time**: 10-15 minutes

**Prerequisites**:
* Pipeline installed and configured
* One sample with forward and reverse AB1 files
* Basic command line familiarity

📂 Sample Preparation
====================

Expected File Structure
-----------------------

For this guide, we'll process a single sample called "sample001":

.. code-block:: text

   project/
   ├── input/
   │   ├── sample001_F.ab1    # Forward read
   │   └── sample001_R.ab1    # Reverse read
   └── output/                # Will be created

File Naming Requirements
-----------------------

The pipeline automatically detects paired reads based on filename patterns:

* **Forward reads**: Must contain ``_F``, ``_forward``, or ``_1`` 
* **Reverse reads**: Must contain ``_R``, ``_reverse``, or ``_2``

**Valid naming examples**:

.. code-block:: text

   sample001_F.ab1 / sample001_R.ab1
   sample001_forward.ab1 / sample001_reverse.ab1
   sample001_1.ab1 / sample001_2.ab1
   MySample_F.ab1 / MySample_R.ab1

🔧 Setup and Configuration
=========================

1. **Create Working Directory**

   .. code-block:: bash

      # Create project directory
      mkdir single_sample_analysis
      cd single_sample_analysis
      
      # Create input directory
      mkdir input
      
      # Copy your AB1 files
      cp /path/to/your/sample001_F.ab1 input/
      cp /path/to/your/sample001_R.ab1 input/

2. **Copy Configuration File**

   .. code-block:: bash

      # Copy default configuration
      cp /path/to/sanger_adna_damage/config/default_config.yaml ./sample_config.yaml

3. **Verify Input Files**

   .. code-block:: bash

      # Check your input files
      ls -la input/
      
      # Should show:
      # sample001_F.ab1
      # sample001_R.ab1

🚀 Step-by-Step Processing
==========================

Step 1: Validate Setup
----------------------

Before running the analysis, validate your setup:

.. code-block:: bash

   # Check if pipeline can find your files
   python -m src.sanger_pipeline.cli.main validate \
       --check-input ./input \
       --config ./sample_config.yaml

**Expected output**:

.. code-block:: text

   ✓ Configuration file is valid
   ✓ Input directory exists
   ✓ Found 2 AB1 files (1 sample pair)
   ✓ External dependencies available
   Validation passed!

Step 2: Run the Complete Pipeline
---------------------------------

Process your sample through all pipeline stages:

.. code-block:: bash

   # Run complete analysis
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output \
       --config ./sample_config.yaml \
       --verbose

**Processing stages** (watch for these in the output):

.. code-block:: text

   Starting Sanger pipeline...
   [1/6] Converting AB1 files to FASTA...
   ✓ Converted sample001_F.ab1 (285 bases, avg quality: 32.4)
   ✓ Converted sample001_R.ab1 (298 bases, avg quality: 29.8)
   
   [2/6] Applying quality filtering...
   ✓ sample001_F: 267 bases retained (93.7%)
   ✓ sample001_R: 275 bases retained (92.3%)
   
   [3/6] Aligning forward and reverse reads...
   ✓ sample001: Successfully aligned using MAFFT
   
   [4/6] Building consensus sequences...
   ✓ sample001_HVS1: 142 bases
   ✓ sample001_HVS2: 198 bases
   ✓ sample001_HVS3: 89 bases
   
   [5/6] Merging HVS regions...
   ✓ sample001: Combined HVS1+HVS2+HVS3 (429 bases total)
   
   [6/6] Analyzing damage patterns...
   ✓ sample001: Damage score = 0.15 (p-value = 0.23)
   
   Pipeline completed successfully!

Step 3: Examine Output Structure
--------------------------------

Explore what the pipeline created:

.. code-block:: bash

   # View the complete output structure
   tree output/

**Output explanation**:

.. code-block:: text

   output/
   ├── fasta/                           # Raw FASTA conversions
   │   ├── sample001_F.fasta           # Forward read as FASTA
   │   └── sample001_R.fasta           # Reverse read as FASTA
   ├── filtered/                        # Quality-filtered sequences
   │   ├── sample001_F_filtered.fasta
   │   └── sample001_R_filtered.fasta
   ├── aligned/                         # Aligned forward+reverse reads
   │   └── sample001_aligned.fasta
   ├── consensus/                       # Consensus by HVS region
   │   ├── sample001_HVS1_consensus.fasta
   │   ├── sample001_HVS2_consensus.fasta
   │   └── sample001_HVS3_consensus.fasta
   ├── final/                          # ⭐ Your main result
   │   └── sample001_final.fasta
   ├── damage_analysis/                # Ancient DNA assessment
   │   └── sample001_damage_analysis.json
   └── plots/                          # Quality visualizations
       ├── sample001_F_quality.png
       └── sample001_R_quality.png

📊 Examine Your Results
=======================

Step 4: View Final Sequence
---------------------------

Your main result is the final consensus sequence:

.. code-block:: bash

   # View your final processed sequence
   cat output/final/sample001_final.fasta

**Example output**:

.. code-block:: text

   >sample001_HVS1_HVS2_HVS3_final
   GATTTCACGGAGGATGGTGGTCAAGGGACCCCCCCTCCCCCATGCTTACAAGCAAGTACA
   TGTTTGTTTGAGATGCTTTGCTCACCCCCTCTCTTTGTTTGCTTTGGAGCACTTGGAACC
   GATGGTGCTGGTTCCGGAGCCCTGTTTATCCACCTTGTTTCCCCTGTATTCCATCTCTAC
   CTTCCAACCCATTCCCACCCCACTCGTTGGTGAATCTTATTTTTCGGTTAGAGTCCCACC
   CTGTGTGACCCTGCTTGTGATGCCGTTAGAGATGGTAACAGAGGTTATCATGCTTCCCTA
   GGCTACTACTGTGCAAGGCCCCCATTTGTTCAATGGAAAGATTTCGTTGATCCGTGTGAC
   CTGGAAACAGGCAAAGATGGGGATGATGGCGCCTCTAGGATAATAGGGCGTGTTTCACGG
   AGGATGGTGGTCAAGGGACCCCCCCTCCCCCATGCTTACAAGCAAGTACATG

Step 5: Check Damage Analysis
-----------------------------

Examine the ancient DNA damage assessment:

.. code-block:: bash

   # View damage analysis results (formatted)
   python -c "
   import json
   with open('output/damage_analysis/sample001_damage_analysis.json') as f:
       data = json.load(f)
   print('Sample:', data['sample_id'])
   print('Damage Score:', data['damage_score'])
   print('P-value:', data['p_value'])
   print('Assessment:', data['assessment'])
   print('C→T rate:', f\"{data['c_to_t_rate']:.2%}\")
   print('G→A rate:', f\"{data['g_to_a_rate']:.2%}\")
   "

**Example output**:

.. code-block:: text

   Sample: sample001
   Damage Score: 0.15
   P-value: 0.23
   Assessment: Low damage detected
   C→T rate: 8.50%
   G→A rate: 6.20%

**Interpretation**:
* **Damage Score 0.15**: Low damage (consistent with modern DNA)
* **P-value 0.23**: Not statistically significant (p > 0.05)
* **Assessment**: Low damage suggests modern DNA or well-preserved sample

Step 6: Generate Interactive Report
----------------------------------

Create a comprehensive QC report:

.. code-block:: bash

   # Generate interactive HTML report
   python -m src.sanger_pipeline.cli.main generate-report \
       --output-dir ./output \
       --title "Sample001 Analysis Report" \
       --open-browser

This opens a detailed report in your browser with:

* **Overview**: Processing summary and quality metrics
* **Damage Analysis**: Detailed damage assessment with plots
* **Quality Control**: Sequence quality distributions
* **Sample Details**: Per-file processing results

🔍 Understanding Your Results
============================

Quality Metrics
---------------

Check these key indicators in your report:

.. code-block:: bash

   # Get quick status summary
   python -m src.sanger_pipeline.cli.main status \
       --output-dir ./output \
       --detailed

**Good quality indicators**:
* ✅ Both forward and reverse reads processed successfully
* ✅ Average quality scores >20
* ✅ Sequence lengths >50bp after filtering
* ✅ Multiple HVS regions detected

HVS Region Coverage
------------------

Your sample should ideally cover multiple HVS regions:

* **HVS1 only**: Partial coverage, adequate for basic analysis
* **HVS1 + HVS2**: Good coverage for most applications  
* **HVS1 + HVS2 + HVS3**: Excellent coverage for comprehensive analysis

Damage Assessment Interpretation
--------------------------------

Based on your damage analysis results:

**Modern DNA Pattern** (damage score <0.3, p>0.05):
* Low C→T and G→A transition rates
* Even damage distribution
* High confidence in sequence authenticity

**Potential Ancient DNA** (damage score >0.3, p<0.05):
* Elevated C→T transitions at 5' ends
* Elevated G→A transitions at 3' ends
* Characteristic ancient DNA damage pattern

**Borderline Cases** (damage score 0.2-0.4):
* May need additional validation
* Consider sample age and preservation conditions
* Additional quality controls recommended

📈 View Quality Plots
====================

Examine quality score distributions:

.. code-block:: bash

   # View quality plots (if you have image viewer)
   open output/plots/sample001_F_quality.png    # macOS
   # or
   xdg-open output/plots/sample001_F_quality.png  # Linux

These plots show:
* Quality score distribution along sequence length
* Areas of high/low quality
* Regions that were filtered out

🔧 Troubleshooting Common Issues
===============================

Issue 1: Low Quality Sequences
------------------------------

**Problem**: Sequences too short after quality filtering

**Solution**: Lower quality threshold temporarily

.. code-block:: bash

   # Re-run with lower quality threshold
   python -m src.sanger_pipeline.cli.main run-pipeline \
       --input-dir ./input \
       --output-dir ./output_lowq \
       --quality-threshold 15

Issue 2: No HVS Regions Detected
--------------------------------

**Problem**: No consensus sequences in HVS regions

**Check sequence content**:

.. code-block:: bash

   # Examine filtered sequences
   cat output/filtered/sample001_F_filtered.fasta
   cat output/filtered/sample001_R_filtered.fasta

**Solutions**:
* Check if sequences are mitochondrial DNA
* Verify HVS region coordinates in configuration
* Consider if sample covers different regions

Issue 3: Alignment Failures
---------------------------

**Problem**: MAFFT alignment fails

**Check**:

.. code-block:: bash

   # Verify MAFFT installation
   mafft --version
   
   # Check sequence compatibility
   head -n 20 output/filtered/sample001_*_filtered.fasta

**Solutions**:
* Ensure MAFFT is properly installed
* Check that sequences are from same organism
* Verify sequences have sufficient length

🎯 Next Steps
=============

After successfully processing your single sample:

1. **Analyze More Samples**: Use :doc:`batch_processing` for multiple samples
2. **Customize Analysis**: Try :doc:`create_custom_config` for different settings
3. **Ancient DNA Focus**: Explore :doc:`assess_damage_patterns` for ancient samples
4. **Publication Reports**: Create publication-ready outputs with :doc:`generate_publication_reports`

🏆 Success Checklist
====================

Your single sample analysis is successful if:

* ✅ Both AB1 files were converted to FASTA
* ✅ Quality filtering retained reasonable sequence lengths
* ✅ Forward and reverse reads were successfully aligned
* ✅ At least one HVS region consensus was generated
* ✅ Final merged sequence was created
* ✅ Damage analysis completed without errors
* ✅ Interactive QC report generated successfully

📊 Typical Results Summary
=========================

For a successful single sample analysis, expect:

* **Input**: 2 AB1 files (forward + reverse)
* **Output**: 1 final consensus sequence
* **Processing time**: 2-5 minutes
* **Quality retention**: 80-95% of original bases
* **HVS coverage**: 1-3 regions depending on sample
* **Final sequence length**: 100-500bp typically

**File sizes** (approximate):
* AB1 files: 50-200KB each
* Final FASTA: 1-2KB
* Damage analysis JSON: 2-5KB
* QC report HTML: 500KB-2MB

This guide provides a complete workflow for processing individual samples, giving you the foundation to understand the pipeline before scaling up to larger analyses.
