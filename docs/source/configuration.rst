=============
Configuration
=============

The Sanger DNA Damage Analysis Pipeline uses YAML configuration files to control its behavior. This guide explains all configuration options and how to customize them for your specific needs.

.. tip::
   **Configuration for Screening Workflows**
   
   Remember that this pipeline is designed for **sample screening and prioritization**. Configure parameters to optimize for identifying promising samples for NGS follow-up, not for definitive authentication.

📁 Configuration Files
======================

Default Configuration
---------------------

The pipeline comes with a default configuration file at ``config/default_config.yaml``. This contains sensible defaults for most use cases:

.. code-block:: yaml

   # config/default_config.yaml
   quality_threshold: 20
   min_sequence_length: 30
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

Custom Configuration
-------------------

Create your own configuration file for specific projects:

.. code-block:: bash

   # Copy default configuration
   cp config/default_config.yaml my_project_config.yaml
   
   # Edit with your preferred editor
   nano my_project_config.yaml
   
   # Use in pipeline
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --config my_project_config.yaml \\
       --input-dir ./input \\
       --output-dir ./output

🔧 Configuration Parameters
===========================

Quality Control Settings
------------------------

quality_threshold
~~~~~~~~~~~~~~~~~

**Type**: Integer  
**Default**: 20  
**Range**: 0-40+  
**Description**: Minimum Phred quality score for sequence filtering

.. code-block:: yaml

   quality_threshold: 20  # Keep bases with Q20+ (99% accuracy)

**Usage Guidelines**:

* **Modern DNA (Q30+)**: Use 25-30 for high-quality modern samples
* **Ancient DNA (Q15-20)**: Use 15-20 for degraded ancient samples  
* **Exploratory (Q10-15)**: Use 10-15 for initial assessment of poor samples

min_sequence_length
~~~~~~~~~~~~~~~~~~~

**Type**: Integer  
**Default**: 50  
**Range**: 10-1000+  
**Description**: Minimum sequence length after quality filtering

.. code-block:: yaml

   min_sequence_length: 50  # Sequences must be at least 50bp

**Usage Guidelines**:

* **Standard Analysis**: 50-100bp minimum for reliable analysis
* **Ancient DNA**: 30-50bp for highly degraded samples
* **High Quality**: 100-200bp for modern, high-quality samples

Ancient DNA Analysis Settings
-----------------------------

damage_threshold
~~~~~~~~~~~~~~~~

**Type**: Float  
**Default**: 0.05  
**Range**: 0.001-0.1  
**Description**: P-value threshold for damage assessment significance

.. code-block:: yaml

   damage_threshold: 0.05  # 5% significance level

**Usage Guidelines**:

* **Conservative**: 0.01 (1%) for strict damage assessment
* **Standard**: 0.05 (5%) for typical analysis
* **Liberal**: 0.1 (10%) for exploratory analysis

bootstrap_iterations
~~~~~~~~~~~~~~~~~~~~

**Type**: Integer  
**Default**: 10000  
**Range**: 1000-100000  
**Description**: Number of bootstrap iterations for damage analysis

.. code-block:: yaml

   bootstrap_iterations: 10000  # 10,000 iterations

**Usage Guidelines**:

* **Quick Testing**: 1000-5000 iterations
* **Standard Analysis**: 10000 iterations
* **High Precision**: 50000-100000 iterations (slower)

HVS Region Definitions
---------------------

hvs_regions
~~~~~~~~~~~

**Type**: Dictionary  
**Description**: Coordinates for hypervariable regions of mitochondrial DNA

.. code-block:: yaml

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

**Customization**:

You can modify these coordinates or add new regions:

.. code-block:: yaml

   hvs_regions:
     HVS1:
       start: 16000  # Extended HVS1 region
       end: 16400
     HVS2:
       start: 50     # Extended HVS2 region
       end: 400
     CUSTOM_REGION:  # Add custom region
       start: 1000
       end: 1500

🎯 Configuration Templates
==========================

Ancient DNA Configuration
-------------------------

Optimized for degraded ancient DNA samples:

.. code-block:: yaml

   # ancient_dna_config.yaml
   
   # Relaxed quality filtering for degraded samples
   quality_threshold: 15
   min_sequence_length: 30
   
   # Sensitive damage detection
   damage_threshold: 0.1
   bootstrap_iterations: 50000
   
   # Standard HVS regions
   hvs_regions:
     HVS1:
       start: 16024
       end: 16365
     HVS2:
       start: 57
       end: 372

Modern DNA Configuration
-----------------------

Optimized for high-quality modern samples:

.. code-block:: yaml

   # modern_dna_config.yaml
   
   # Strict quality filtering
   quality_threshold: 30
   min_sequence_length: 100
   
   # Conservative damage detection (expecting no damage)
   damage_threshold: 0.01
   bootstrap_iterations: 10000
   
   # Standard HVS regions
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

Exploratory Configuration
------------------------

For initial assessment of unknown samples:

.. code-block:: yaml

   # exploratory_config.yaml
   
   # Permissive quality filtering
   quality_threshold: 10
   min_sequence_length: 25
   
   # Liberal damage detection
   damage_threshold: 0.1
   bootstrap_iterations: 5000
   
   # All HVS regions
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

High-Throughput Configuration
----------------------------

For processing large numbers of samples quickly:

.. code-block:: yaml

   # high_throughput_config.yaml
   
   # Balanced quality filtering
   quality_threshold: 20
   min_sequence_length: 50
   
   # Fast damage analysis
   damage_threshold: 0.05
   bootstrap_iterations: 5000  # Reduced for speed
   
   # Focus on most informative regions
   hvs_regions:
     HVS1:
       start: 16024
       end: 16365
     HVS2:
       start: 57
       end: 372

🧪 Validation and Testing
=========================

Configuration Validation
------------------------

Test your configuration before running large analyses:

.. code-block:: bash

   # Validate configuration syntax
   python -c "import yaml; yaml.safe_load(open('my_config.yaml'))"
   
   # Test with pipeline status command
   python -m src.sanger_pipeline.cli.main status --config my_config.yaml
   
   # Run on small test dataset
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --input-dir ./test_data \\
       --output-dir ./test_output \\
       --config my_config.yaml

Parameter Testing
----------------

Test different parameter values systematically:

.. code-block:: bash

   # Test different quality thresholds
   for threshold in 15 20 25 30; do
       echo "Testing quality threshold: $threshold"
       python -m src.sanger_pipeline.cli.main run-pipeline \\
           --input-dir ./test_data \\
           --output-dir ./output_q${threshold} \\
           --quality-threshold $threshold
   done

🔍 Advanced Configuration
========================

Environment Variables
---------------------

Some settings can be controlled via environment variables:

.. code-block:: bash

   # Override configuration file location
   export SANGER_CONFIG=/path/to/my/config.yaml
   
   # Set temporary directory
   export TMPDIR=/path/to/large/temp/space
   
   # Control memory usage
   export MAX_MEMORY_GB=8

Command Line Overrides
----------------------

You can override configuration values from the command line:

.. code-block:: bash

   # Override quality threshold
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --config my_config.yaml \\
       --quality-threshold 25 \\
       --input-dir ./input \\
       --output-dir ./output

Configuration Validation Schema
------------------------------

The pipeline validates configuration files against a schema. Required fields:

.. code-block:: yaml

   # Minimum required configuration
   quality_threshold: 20
   damage_threshold: 0.05
   
   hvs_regions:
     HVS1:
       start: 16024
       end: 16365

🔄 Configuration Management
===========================

Version Control
---------------

Track your configuration files in version control:

.. code-block:: bash

   # Add configuration to git
   git add my_project_config.yaml
   git commit -m "Add project-specific configuration"

Multiple Configurations
-----------------------

Organize configurations by project or sample type:

.. code-block:: text

   configs/
   ├── default_config.yaml
   ├── ancient_dna/
   │   ├── permafrost_samples.yaml
   │   └── cave_samples.yaml
   ├── modern_dna/
   │   ├── reference_samples.yaml
   │   └── population_study.yaml
   └── exploratory/
       └── unknown_samples.yaml

Configuration Documentation
---------------------------

Document your custom configurations:

.. code-block:: yaml

   # ancient_permafrost_config.yaml
   # Configuration for ancient DNA from permafrost samples
   # Created: 2024-01-15
   # Author: Research Team
   # Purpose: Optimized for highly degraded permafrost samples
   
   quality_threshold: 12  # Very permissive due to degradation
   min_sequence_length: 25  # Short fragments expected
   damage_threshold: 0.1   # Liberal due to expected damage
   bootstrap_iterations: 50000  # High precision for publication

⚠️ Common Configuration Issues
==============================

YAML Syntax Errors
------------------

.. code-block:: yaml

   # ❌ Incorrect indentation
   hvs_regions:
   HVS1:
     start: 16024
   
   # ✅ Correct indentation
   hvs_regions:
     HVS1:
       start: 16024

Invalid Parameter Values
-----------------------

.. code-block:: yaml

   # ❌ Invalid quality threshold
   quality_threshold: 45  # Too high
   
   # ✅ Valid quality threshold
   quality_threshold: 25  # Reasonable for high-quality samples

Missing Required Fields
----------------------

.. code-block:: yaml

   # ❌ Missing required fields
   quality_threshold: 20
   # Missing hvs_regions!
   
   # ✅ All required fields
   quality_threshold: 20
   damage_threshold: 0.05
   hvs_regions:
     HVS1:
       start: 16024
       end: 16365

🎯 Best Practices
=================

1. **Start with Defaults**: Use the default configuration as a starting point
2. **Document Changes**: Comment your modifications and reasoning
3. **Test Thoroughly**: Validate configurations on small datasets first
4. **Version Control**: Track configuration changes alongside code
5. **Project-Specific**: Create separate configurations for different projects
6. **Parameter Testing**: Systematically test different parameter values
7. **Backup Configs**: Keep copies of working configurations

📊 Performance Tuning
=====================

Quality vs. Speed Trade-offs
---------------------------

.. code-block:: yaml

   # Fast processing (lower quality)
   quality_threshold: 15
   bootstrap_iterations: 5000
   
   # High quality (slower processing)
   quality_threshold: 25
   bootstrap_iterations: 50000

Memory Optimization
------------------

.. code-block:: yaml

   # For large datasets, reduce memory usage
   min_sequence_length: 100  # Filter short sequences early
   quality_threshold: 25     # Strict filtering reduces data volume

🔗 Integration with Other Tools
==============================

Export Configuration
--------------------

.. code-block:: bash

   # Convert to other formats for external tools
   python -c "
   import yaml, json
   with open('my_config.yaml') as f:
       config = yaml.safe_load(f)
   with open('my_config.json', 'w') as f:
       json.dump(config, f, indent=2)
   "

Pipeline Integration
-------------------

.. code-block:: bash

   # Use in automated pipelines
   CONFIG_FILE="configs/production_config.yaml"
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --config "$CONFIG_FILE" \\
       --input-dir "$INPUT_DIR" \\
       --output-dir "$OUTPUT_DIR"

This comprehensive configuration system allows you to fine-tune the pipeline for your specific research needs, from quick exploratory analyses to publication-ready ancient DNA assessments.
