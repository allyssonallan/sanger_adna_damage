=============
API Reference
=============

Complete API documentation for the Sanger DNA Damage Analysis Pipeline modules, classes, and functions.

📚 Overview
===========

The pipeline is organized into several main modules:

* **Core Pipeline** (:doc:`core`) - Main pipeline orchestration and coordination
* **Utilities** (:doc:`utils`) - Helper functions for file processing, analysis, and reporting
* **CLI** (:doc:`cli`) - Command-line interface implementation

Each module provides both high-level interfaces for common use cases and low-level functions for custom workflows.

🏗️ Architecture Overview
========================

The pipeline follows a modular architecture:

.. code-block:: text

   src/sanger_pipeline/
   ├── core/                 # Core pipeline logic
   │   ├── pipeline.py      # Main pipeline orchestrator
   │   ├── processor.py     # AB1 file processing
   │   └── analyzer.py      # Damage analysis engine
   ├── utils/               # Utility modules
   │   ├── ab1_converter.py # AB1 to FASTA conversion
   │   ├── quality_filter.py # Quality filtering
   │   ├── consensus.py     # Consensus sequence building
   │   ├── damage_analyzer.py # Ancient DNA damage analysis
   │   └── report_generator.py # QC report generation
   └── cli/                 # Command-line interface
       └── main.py         # CLI entry point

🎯 Quick Start
==============

**Basic Pipeline Usage**:

.. code-block:: python

   from src.sanger_pipeline.core.pipeline import SangerPipeline
   
   # Initialize pipeline
   pipeline = SangerPipeline(config_path="config/default_config.yaml")
   
   # Run complete analysis
   results = pipeline.run(
       input_dir="./ab1_files",
       output_dir="./results"
   )

**Direct Module Usage**:

.. code-block:: python

   from src.sanger_pipeline.utils.ab1_converter import AB1Converter
   from src.sanger_pipeline.utils.damage_analyzer import DamageAnalyzer
   
   # Convert AB1 files
   converter = AB1Converter()
   fasta_files = converter.convert_directory("./ab1_files", "./fasta_output")
   
   # Analyze damage patterns
   analyzer = DamageAnalyzer()
   damage_results = analyzer.analyze_sequences("./final_sequences")

📖 Module Documentation
=======================

.. toctree::
   :maxdepth: 2
   
   core
   utils
   cli

🔧 Configuration API
====================

**Configuration Loading**:

.. code-block:: python

   from src.sanger_pipeline.utils.config import load_config, validate_config
   
   # Load configuration
   config = load_config("my_config.yaml")
   
   # Validate configuration
   is_valid, errors = validate_config(config)

**Dynamic Configuration**:

.. code-block:: python

   # Override configuration at runtime
   pipeline = SangerPipeline()
   pipeline.config.quality_threshold = 25
   pipeline.config.bootstrap_iterations = 50000

📊 Data Structures
==================

**Pipeline Results**:

.. code-block:: python

   class PipelineResults:
       """Results from complete pipeline execution"""
       
       def __init__(self):
           self.input_files: List[str] = []
           self.converted_files: List[str] = []
           self.filtered_files: List[str] = []
           self.consensus_files: List[str] = []
           self.final_files: List[str] = []
           self.damage_results: Dict[str, DamageResult] = {}
           self.processing_stats: ProcessingStats = ProcessingStats()
           self.errors: List[ProcessingError] = []

**Damage Analysis Results**:

.. code-block:: python

   class DamageResult:
       """Results from ancient DNA damage analysis"""
       
       def __init__(self):
           self.sample_id: str
           self.damage_score: float
           self.p_value: float
           self.c_to_t_rate: float
           self.g_to_a_rate: float
           self.background_rate: float
           self.assessment: str
           self.significance: str
           self.position_data: Dict[str, List[float]]
           self.bootstrap_stats: BootstrapStats

🛠️ Utility Functions
====================

**File Processing**:

.. code-block:: python

   from src.sanger_pipeline.utils.file_utils import (
       find_ab1_files, pair_reads, validate_input_dir
   )
   
   # Find and pair AB1 files
   ab1_files = find_ab1_files("./input")
   paired_reads = pair_reads(ab1_files)
   
   # Validate input directory
   is_valid, message = validate_input_dir("./input")

**Quality Assessment**:

.. code-block:: python

   from src.sanger_pipeline.utils.quality_utils import (
       calculate_quality_stats, assess_sequence_quality
   )
   
   # Calculate quality statistics
   stats = calculate_quality_stats("sequence.fasta")
   
   # Assess overall quality
   quality_report = assess_sequence_quality(stats, threshold=20)

📈 Analysis Functions
====================

**Damage Analysis**:

.. code-block:: python

   from src.sanger_pipeline.utils.damage_analyzer import (
       calculate_damage_score, bootstrap_analysis, assess_significance
   )
   
   # Calculate damage score
   score = calculate_damage_score(sequences)
   
   # Perform bootstrap analysis
   p_value = bootstrap_analysis(sequences, iterations=10000)
   
   # Assess statistical significance
   significance = assess_significance(p_value, threshold=0.05)

**Consensus Building**:

.. code-block:: python

   from src.sanger_pipeline.utils.consensus import (
       build_consensus, merge_hvs_regions, validate_consensus
   )
   
   # Build consensus sequence
   consensus = build_consensus(forward_seq, reverse_seq)
   
   # Merge HVS regions
   final_seq = merge_hvs_regions(hvs_sequences)
   
   # Validate consensus quality
   is_valid = validate_consensus(consensus, min_length=50)

🎨 Visualization API
===================

**Report Generation**:

.. code-block:: python

   from src.sanger_pipeline.utils.report_generator import ReportGenerator
   
   # Generate interactive HTML report
   generator = ReportGenerator()
   report_path = generator.generate_report(
       output_dir="./results",
       title="My Analysis Report"
   )

**Custom Plots**:

.. code-block:: python

   from src.sanger_pipeline.utils.plotting import (
       plot_quality_scores, plot_damage_profile, plot_length_distribution
   )
   
   # Create quality score plot
   plot_quality_scores("sequence.fasta", "quality_plot.png")
   
   # Create damage profile plot
   plot_damage_profile(damage_results, "damage_plot.png")

🔌 Extension Points
==================

**Custom Processors**:

.. code-block:: python

   from src.sanger_pipeline.core.processor import BaseProcessor
   
   class CustomProcessor(BaseProcessor):
       """Custom processing step"""
       
       def process(self, input_data):
           # Custom processing logic
           return processed_data

**Custom Analyzers**:

.. code-block:: python

   from src.sanger_pipeline.utils.damage_analyzer import BaseDamageAnalyzer
   
   class CustomDamageAnalyzer(BaseDamageAnalyzer):
       """Custom damage analysis method"""
       
       def analyze(self, sequences):
           # Custom analysis logic
           return damage_results

🚨 Error Handling
================

**Exception Types**:

.. code-block:: python

   from src.sanger_pipeline.exceptions import (
       PipelineError, ConfigurationError, FileProcessingError,
       AnalysisError, ValidationError
   )
   
   try:
       pipeline.run(input_dir, output_dir)
   except FileProcessingError as e:
       print(f"File processing failed: {e}")
   except AnalysisError as e:
       print(f"Analysis failed: {e}")

**Error Recovery**:

.. code-block:: python

   # Robust processing with error handling
   results = pipeline.run_with_recovery(
       input_dir="./input",
       output_dir="./output",
       continue_on_error=True
   )
   
   # Check for errors
   if results.errors:
       for error in results.errors:
           print(f"Error in {error.file}: {error.message}")

🧪 Testing Utilities
====================

**Test Data Generation**:

.. code-block:: python

   from src.sanger_pipeline.testing.generators import (
       generate_test_ab1, generate_test_sequences, generate_damage_pattern
   )
   
   # Generate test AB1 file
   test_ab1 = generate_test_ab1(length=300, quality_mean=25)
   
   # Generate test sequences with damage
   damaged_seqs = generate_damage_pattern(sequences, damage_rate=0.3)

**Mock Objects**:

.. code-block:: python

   from src.sanger_pipeline.testing.mocks import MockPipeline, MockAnalyzer
   
   # Use mock pipeline for testing
   mock_pipeline = MockPipeline()
   results = mock_pipeline.run("./test_input", "./test_output")

🔍 Debugging and Logging
========================

**Logging Configuration**:

.. code-block:: python

   import logging
   from src.sanger_pipeline.utils.logging import setup_logging
   
   # Set up logging
   setup_logging(level=logging.DEBUG, log_file="pipeline.log")
   
   # Use logger in your code
   logger = logging.getLogger(__name__)
   logger.info("Starting analysis")

**Debug Mode**:

.. code-block:: python

   # Enable debug mode
   pipeline = SangerPipeline(debug=True)
   pipeline.run(input_dir, output_dir)
   
   # Access debug information
   debug_info = pipeline.get_debug_info()

📊 Performance Monitoring
=========================

**Performance Metrics**:

.. code-block:: python

   from src.sanger_pipeline.utils.performance import PerformanceMonitor
   
   # Monitor performance
   monitor = PerformanceMonitor()
   
   with monitor.time_block("conversion"):
       converter.convert_files(ab1_files)
   
   # Get performance report
   report = monitor.get_report()

**Memory Management**:

.. code-block:: python

   from src.sanger_pipeline.utils.memory import MemoryManager
   
   # Monitor memory usage
   memory_manager = MemoryManager()
   memory_manager.start_monitoring()
   
   # Process data
   results = pipeline.run(input_dir, output_dir)
   
   # Get memory report
   memory_report = memory_manager.get_report()

This API reference provides comprehensive documentation for all public interfaces, making it easy to integrate the pipeline into custom workflows or extend its functionality.
