===========
Core Module
===========

The core module contains the main pipeline orchestration and coordination logic.

.. automodule:: sanger_pipeline.core
   :members:
   :undoc-members:
   :show-inheritance:

pipeline module
===============

.. automodule:: sanger_pipeline.core.pipeline
   :members:
   :undoc-members:
   :show-inheritance:

SangerPipeline Class
-------------------

.. autoclass:: sanger_pipeline.core.pipeline.SangerPipeline
   :members:
   :undoc-members:
   :show-inheritance:

   The main pipeline orchestrator that coordinates all processing steps from AB1 conversion to final damage analysis.

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.pipeline import SangerPipeline

      # Initialize with default configuration
      pipeline = SangerPipeline()

      # Run complete analysis
      results = pipeline.run(
          input_dir="./ab1_files",
          output_dir="./results"
      )

      # Check results
      print(f"Processed {len(results.final_files)} samples")
      print(f"Success rate: {results.processing_stats.success_rate:.1%}")

   **Configuration**:

   The pipeline accepts configuration through YAML files or direct parameter setting:

   .. code-block:: python

      # Load from configuration file
      pipeline = SangerPipeline(config_path="my_config.yaml")

      # Or set parameters directly
      pipeline = SangerPipeline()
      pipeline.config.quality_threshold = 25
      pipeline.config.damage_threshold = 0.01

processor module
================

.. automodule:: sanger_pipeline.core.processor
   :members:
   :undoc-members:
   :show-inheritance:

AB1Processor Class
------------------

.. autoclass:: sanger_pipeline.core.processor.AB1Processor
   :members:
   :undoc-members:
   :show-inheritance:

   Handles AB1 file processing including conversion, quality filtering, and alignment.

   **Processing Steps**:

   1. **AB1 Conversion**: Convert binary AB1 files to FASTA format
   2. **Quality Filtering**: Remove low-quality bases and sequences
   3. **Read Pairing**: Match forward and reverse reads
   4. **Alignment**: Align paired reads using MAFFT
   5. **Consensus Building**: Generate consensus sequences

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.processor import AB1Processor

      processor = AB1Processor(config)
      
      # Process all AB1 files in directory
      results = processor.process_directory(
          input_dir="./ab1_files",
          output_dir="./processed"
      )

      # Process individual file pair
      consensus = processor.process_pair(
          forward_file="sample_F.ab1",
          reverse_file="sample_R.ab1"
      )

HVSProcessor Class
------------------

.. autoclass:: sanger_pipeline.core.processor.HVSProcessor
   :members:
   :undoc-members:
   :show-inheritance:

   Specialized processor for hypervariable sequence (HVS) region analysis.

   **HVS Region Processing**:

   - Identifies HVS1, HVS2, and HVS3 regions in sequences
   - Builds consensus sequences for each region independently
   - Merges available regions into final sequences
   - Handles missing or partial region coverage

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.processor import HVSProcessor

      hvs_processor = HVSProcessor(config)
      
      # Process HVS regions
      hvs_results = hvs_processor.process_hvs_regions(
          consensus_sequence,
          sample_id="sample001"
      )

      # Merge available regions
      final_sequence = hvs_processor.merge_regions(hvs_results)

analyzer module
===============

.. automodule:: sanger_pipeline.core.analyzer
   :members:
   :undoc-members:
   :show-inheritance:

DamageAnalyzer Class
--------------------

.. autoclass:: sanger_pipeline.core.analyzer.DamageAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   Core engine for ancient DNA damage analysis and statistical validation.

   **Analysis Methods**:

   - **Position-based damage detection**: Analyzes C→T and G→A transitions by position
   - **Statistical validation**: Bootstrap analysis with configurable iterations
   - **Significance testing**: P-value calculation and assessment
   - **Pattern recognition**: Identifies authentic ancient DNA patterns

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.analyzer import DamageAnalyzer

      analyzer = DamageAnalyzer(config)
      
      # Analyze damage in sequences
      damage_results = analyzer.analyze_directory("./final_sequences")

      # Get detailed results for specific sample
      sample_result = analyzer.analyze_sample("sample001.fasta")
      
      print(f"Damage score: {sample_result.damage_score:.3f}")
      print(f"P-value: {sample_result.p_value:.4f}")
      print(f"Assessment: {sample_result.assessment}")

QualityAnalyzer Class
---------------------

.. autoclass:: sanger_pipeline.core.analyzer.QualityAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   Comprehensive quality assessment for sequences and processing results.

   **Quality Metrics**:

   - Sequence quality score distributions
   - Length statistics and filtering efficiency
   - Processing success rates by stage
   - Coverage assessment for HVS regions

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.analyzer import QualityAnalyzer

      quality_analyzer = QualityAnalyzer()
      
      # Analyze processing results
      quality_report = quality_analyzer.analyze_pipeline_output("./results")
      
      # Get quality metrics for specific files
      file_metrics = quality_analyzer.analyze_file("sequence.fasta")

Core Data Structures
====================

PipelineResults
---------------

.. autoclass:: sanger_pipeline.core.pipeline.PipelineResults
   :members:
   :undoc-members:
   :show-inheritance:

   Container for complete pipeline execution results.

   **Attributes**:

   - ``input_files``: List of original AB1 input files
   - ``converted_files``: Successfully converted FASTA files
   - ``filtered_files``: Quality-filtered sequence files
   - ``consensus_files``: Consensus sequences by HVS region
   - ``final_files``: Final merged sequence files
   - ``damage_results``: Damage analysis results by sample
   - ``processing_stats``: Overall processing statistics
   - ``errors``: List of processing errors encountered

ProcessingStats
---------------

.. autoclass:: sanger_pipeline.core.pipeline.ProcessingStats
   :members:
   :undoc-members:
   :show-inheritance:

   Statistical summary of pipeline processing efficiency.

   **Metrics**:

   - Total files processed and success rates
   - Processing time by stage
   - Quality filtering statistics
   - HVS region coverage rates

DamageResult
------------

.. autoclass:: sanger_pipeline.core.analyzer.DamageResult
   :members:
   :undoc-members:
   :show-inheritance:

   Comprehensive ancient DNA damage analysis results for a single sample.

   **Core Metrics**:

   - ``damage_score``: Overall damage assessment (0-1 scale)
   - ``p_value``: Statistical significance of damage pattern
   - ``c_to_t_rate``: C→T transition rate at 5' ends
   - ``g_to_a_rate``: G→A transition rate at 3' ends
   - ``assessment``: Human-readable interpretation
   - ``significance``: Statistical significance category

Configuration Management
=========================

ConfigManager
-------------

.. autoclass:: sanger_pipeline.core.config.ConfigManager
   :members:
   :undoc-members:
   :show-inheritance:

   Handles configuration loading, validation, and management.

   **Configuration Sources**:

   - YAML configuration files
   - Environment variables
   - Command-line parameter overrides
   - Default fallback values

   **Example Usage**:

   .. code-block:: python

      from sanger_pipeline.core.config import ConfigManager

      # Load configuration
      config_manager = ConfigManager()
      config = config_manager.load("my_config.yaml")

      # Validate configuration
      is_valid, errors = config_manager.validate(config)

      # Override parameters
      config.quality_threshold = 25
      config.bootstrap_iterations = 50000

Error Handling
==============

Core Exceptions
---------------

.. autoexception:: sanger_pipeline.core.exceptions.PipelineError
   :members:
   :show-inheritance:

   Base exception for all pipeline-related errors.

.. autoexception:: sanger_pipeline.core.exceptions.ProcessingError
   :members:
   :show-inheritance:

   Raised when file processing fails.

.. autoexception:: sanger_pipeline.core.exceptions.AnalysisError
   :members:
   :show-inheritance:

   Raised when damage analysis fails.

.. autoexception:: sanger_pipeline.core.exceptions.ConfigurationError
   :members:
   :show-inheritance:

   Raised when configuration is invalid or missing.

**Error Handling Example**:

.. code-block:: python

   from sanger_pipeline.core.exceptions import ProcessingError, AnalysisError

   try:
       results = pipeline.run(input_dir, output_dir)
   except ProcessingError as e:
       print(f"Processing failed: {e}")
       # Handle processing error
   except AnalysisError as e:
       print(f"Analysis failed: {e}")
       # Handle analysis error

Performance Considerations
=========================

**Memory Management**:

The core modules are designed for efficient memory usage:

- Streaming processing for large files
- Lazy loading of sequence data
- Automatic cleanup of temporary files

**Parallel Processing**:

Some operations support parallel execution:

.. code-block:: python

   # Enable parallel processing
   pipeline = SangerPipeline(n_jobs=4)
   results = pipeline.run(input_dir, output_dir)

**Performance Monitoring**:

.. code-block:: python

   # Enable performance monitoring
   pipeline = SangerPipeline(monitor_performance=True)
   results = pipeline.run(input_dir, output_dir)
   
   # Access performance metrics
   perf_report = pipeline.get_performance_report()

This core module documentation provides comprehensive coverage of the main pipeline orchestration and analysis components.
