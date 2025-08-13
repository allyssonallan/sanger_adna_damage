# Sanger DNA Damage Analysis Pipeline

A comprehensive, modular pipeline for processing Sanger sequencing AB1 files, including quality control, alignment, consensus building, and damage analysis.

## 🚀 Features

- **Modular Architecture**: Well-organized codebase with clear separation of concerns
- **Quality Control**: Convert AB1 files with Phred quality filtering and visualization
- **Sequence Processing**: Align forward/reverse reads and build consensus sequences
- **HVS Region Processing**: Independent processing of HVS1, HVS2, and HVS3 regions with intelligent merging
- **Ancient DNA Analysis**: Comprehensive aDNA damage pattern detection and assessment
- **Statistical Validation**: Bootstrap analysis for damage assessment
- **Beautiful QC Reports**: Interactive HTML reports with charts, tables, and analysis summaries
- **Command Line Interface**: Easy-to-use CLI for all pipeline operations
- **Extensible Design**: Easy to add new analysis modules and features

## 📈 Output Structure

The pipeline creates organized output directories:

```bash
results/
├── fasta/              # Raw FASTA files converted from AB1
├── filtered/           # Quality-filtered sequences with N-substitution  
├── consensus/          # Consensus sequences from paired F/R reads for each HVS region
├── aligned/            # Intermediate alignment files
├── final/              # Merged HVS region sequences (HVS1_HVS2_HVS3, HVS1_HVS2, etc.)
├── damage_analysis/    # aDNA damage analysis results
│   ├── *_damage_results.json  # Detailed damage statistics
│   └── *_damage_plots.png     # Damage pattern visualizations
├── plots/              # Quality score plots for each AB1 file
└── reports/            # Comprehensive QC HTML reports with interactive visualizations
```

## 📊 Pipeline Workflow

1. **AB1 Conversion** (`src/sanger_pipeline/core/ab1_converter.py`)

```bash
├── src/                          # Source code
│   └── sanger_pipeline/          # Main package
│       ├── core/                 # Core functionality
│       │   ├── ab1_converter.py  # AB1 conversion and quality filtering
│       │   ├── consensus_builder.py # Sequence alignment and consensus
│       │   ├── quality_filter.py # Quality filtering utilities
│       │   ├── adna_damage_analyzer.py # Ancient DNA damage analysis
│       │   └── pipeline.py       # Main pipeline orchestrator
│       ├── io/                   # Input/Output handling
│       ├── plotting/             # Visualization utilities
│       ├── utils/                # Helper functions and constants
│       └── cli/                  # Command-line interface
├── scripts/                      # Standalone scripts
│   └── run_pipeline.py           # Main pipeline runner
├── config/                       # Configuration files
│   └── default_config.yaml       # Default settings
├── tests/                        # Test suite
├── input/                        # Input data directory
├── output/                       # Pipeline output directory
├── ref/                          # Reference sequences (e.g., rCRS.fasta)
├── venv/                         # Python virtual environment
├── requirements.txt              # Python dependencies
├── setup.py                     # Package installation
├── LICENSE                      # MIT License
└── README.md                    # This file
```

## 🛠️ Installation

### Prerequisites

- **Python 3.8+**
- **MAFFT** (for sequence alignment)
- **Quarto** (for report generation)

### System Dependencies

```bash
# On macOS with Homebrew
brew install mafft quarto

# On Ubuntu/Debian
sudo apt-get install mafft
# Install Quarto from https://quarto.org/docs/get-started/

# On CentOS/RHEL
sudo yum install mafft
```

### Python Environment

```bash
# Clone the repository
git clone https://github.com/allyssonallan/sanger_adna_damage.git
cd sanger_adna_damage

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install the package in development mode
pip install -e .
```

## 🎯 Quick Start

### Basic Usage

```bash
# Run the complete pipeline
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output \
    --config ./config/default_config.yaml

# Convert a single AB1 file
python -m src.sanger_pipeline.cli.main convert-ab1 \
    --input-file sample.ab1 \
    --output-file sample.fasta \
    --generate-plot

# Analyze aDNA damage patterns for a single sequence
python -m src.sanger_pipeline.cli.main analyze-damage \
    --input-file consensus.fasta \
    --reference ref/rCRS.fasta \
    --output-dir damage_results/

# Generate damage analysis summary
python -m src.sanger_pipeline.cli.main damage-summary \
    --results-dir ./output/damage_analysis/

# Generate comprehensive QC report with beautiful visualizations
python -m src.sanger_pipeline.cli.main generate-report \
    --output-dir ./output/ \
    --open-browser

# Check processing status
python -m src.sanger_pipeline.cli.main status \
    --input-dir ./output/
```

### 📊 QC Report Generation

Generate beautiful, interactive HTML reports with comprehensive analysis summaries:

```bash
# Generate QC report with all analysis results
python -m src.sanger_pipeline.cli.main generate-report \
    --output-dir ./output/ \
    --open-browser

# Or use the standalone script
python generate_report.py ./output/
```

The QC report includes:

- **📈 Overview Dashboard**: Sample counts, processing statistics, and key metrics
- **📁 Directory Analysis**: File counts, sizes, and types for each output directory
- **🧬 Sample Processing**: Individual sample status with HVS region breakdown
- **☢️ Damage Analysis**: aDNA damage patterns, damage assessment scores, and quality metrics
- **🔗 HVS Combinations**: Distribution of merged HVS region combinations
- **📊 Interactive Charts**: Beautiful visualizations powered by Chart.js
- **⏰ Timestamps**: Report generation time and processing timeline

**Report Features:**
- 🎨 Modern, responsive design with Bootstrap 5
- 📱 Mobile-friendly tabbed interface
- 🎯 Color-coded status indicators
- 📈 Interactive charts and progress bars
- 🔍 Searchable and sortable data tables
- ⬇️ Self-contained HTML (no internet required)

### Advanced Configuration

Create a custom configuration file:

```yaml
# my_config.yaml
quality:
  min_phred_score: 25
  terminal_length: 15

alignment:
  tool: "mafft"
  parameters: "--auto --maxiterate 1000"

bootstrap:
  iterations: 50000

plotting:
  figure_size: [15, 6]
  dpi: 300
```

Run with custom configuration:

```bash
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir data/ \
    --output-dir results/ \
    --config my_config.yaml
```

## � Processing Steps

1. **AB1 Conversion** (`src/sanger_pipeline/core/ab1_converter.py`)

   - Convert AB1 files to FASTA format
   - Filter low-quality bases (replace with 'N')
   - Generate quality score plots

2. **Sequence Alignment** (`src/sanger_pipeline/core/consensus_builder.py`)

   - Reverse-complement reverse reads
   - Align forward and reverse sequences using MAFFT
   - Build consensus sequences from alignments for each HVS region

3. **HVS Region Processing**

   - Parse HVS region information from filenames (HVS1, HVS2, HVS3)
   - Process each HVS region independently
   - Generate consensus sequences for available regions
   - Merge regions in proper order (HVS1, HVS2, HVS3) based on availability

4. **Ancient DNA Damage Analysis** (`src/sanger_pipeline/core/adna_damage_analyzer.py`)

   - Detect C→T transitions at 5' terminus
   - Analyze G→A transitions at 3' terminus
   - Perform bootstrap statistical validation (10,000 iterations)
   - Generate damage assessment scores
   - Create damage pattern visualizations ("smile plots")

5. **Quality Control Report**

   - Generate interactive HTML reports with Quarto
   - Include quality statistics, damage analysis, plots, and recommendations

## 🧬 Ancient DNA Analysis

The pipeline includes comprehensive ancient DNA damage analysis capabilities:

### HVS Region Processing

The pipeline expects AB1 files with HVS region information in the filename:

- Format: `SAMPLE_NAME_HVS#-Direction.ab1`
- Examples: 
  - `1_SJT_P1_PB_A_B1_HVS1-F.ab1` (HVS1 forward)
  - `1_SJT_P1_PB_A_B1_HVS1-R.ab1` (HVS1 reverse)
  - `11_SJC_FO_PB_C2_B1_HVS2-F.ab1` (HVS2 forward)

The pipeline will:

1. Group files by sample name and HVS region
2. Process each HVS region independently (align F/R, build consensus)
3. Merge available HVS regions per sample with appropriate naming:
   - `sample_HVS1_HVS2_HVS3_merged.fasta` (all three regions)
   - `sample_HVS1_HVS2_merged.fasta` (two regions)
   - `sample_HVS2_merged.fasta` (single region)

### Damage Analysis Features

- **Damage Pattern Detection**: Identifies C→T and G→A transitions characteristic of aDNA
- **Terminal Analysis**: Focuses on 5' and 3' termini where damage is most prevalent
- **Bootstrap Validation**: Statistical damage assessment using 10,000 bootstrap iterations
- **Damage Scoring**: Confidence assessment for ancient vs. modern DNA
- **Visualization**: "Smile plots" showing damage patterns across sequence positions

### Usage Examples

```bash
# Run full pipeline with aDNA analysis
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir data/ab1_files/ \
    --output-dir results/

# Analyze single sequence for damage patterns
python -m src.sanger_pipeline.cli.main analyze-damage \
    --input-file consensus.fasta \
    --reference ref/rCRS.fasta \
    --output-dir damage_results/

# Generate summary of multiple damage analyses
python -m src.sanger_pipeline.cli.main damage-summary \
    --results-dir results/damage_analysis/
```

### Output

The aDNA analysis generates:

- **JSON Results**: Detailed damage statistics and damage assessment scores
- **Visual Plots**: Characteristic damage plots and pattern visualizations
- **Summary Reports**: Batch analysis summaries with damage assessments

### Configuration

```yaml
# aDNA-specific configuration
quality:
  terminal_length: 10        # Terminal bases to analyze (5-25bp)
  
bootstrap:
  iterations: 10000          # Bootstrap iterations for validation
  
damage:
  ct_threshold: 0.05         # C→T transition threshold
  ga_threshold: 0.05         # G→A transition threshold
  damage_threshold: 0.5      # Threshold for damage pattern classification
```

## 🔧 Module Documentation

### Core Modules

- **`AB1Converter`**: Handles AB1 file conversion and quality filtering
- **`ConsensusBuilder`**: Manages sequence alignment and consensus building
- **`QualityFilter`**: Provides quality filtering utilities
- **`ADNADamageAnalyzer`**: Ancient DNA damage pattern detection and assessment
- **`SangerPipeline`**: Main pipeline orchestrator with integrated aDNA analysis

### Configuration Options

The pipeline uses YAML configuration files for flexible parameter management:

```yaml
quality:
  min_phred_score: 20        # Minimum quality threshold for bases
  terminal_length: 10        # Length for terminal damage analysis (bp)

alignment:
  tool: "mafft"             # Alignment tool (currently supports MAFFT)
  parameters: "--auto"       # Tool-specific parameters

bootstrap:
  iterations: 10000          # Bootstrap iterations for statistical validation

damage:
  ct_threshold: 0.05         # C→T transition threshold for damage assessment
  ga_threshold: 0.05         # G→A transition threshold for damage assessment
  damage_threshold: 0.5      # Threshold for damage pattern classification

output:
  directories:              # Output directory structure
    fasta: "fasta"
    filtered: "filtered"
    consensus: "consensus"
    aligned: "aligned"
    final: "final"
    damage_analysis: "damage_analysis"
    plots: "plots"
    reports: "reports"
```

### Adding New Modules

To add new functionality:

1. Create new module in appropriate `src/sanger_pipeline/` subdirectory
1. Add imports to `src/sanger_pipeline/__init__.py`
1. Update pipeline orchestrator if needed
1. Add tests in `tests/` directory
1. Update documentation

## 🧪 Testing

```bash
# Run all tests
python -m pytest tests/

# Run with coverage
python -m pytest tests/ --cov=sanger_pipeline

# Run specific test file
python -m pytest tests/test_core/test_ab1_converter.py
```

## � Output Structure

The pipeline creates organized output directories:

```bash
results/
├── fasta/              # Raw FASTA files converted from AB1
├── filtered/           # Quality-filtered sequences with N-substitution  
├── consensus/          # Consensus sequences from paired F/R reads
├── aligned/            # Intermediate alignment files
├── final/              # Merged HVS region sequences
├── damage_analysis/    # aDNA damage analysis results
│   ├── *_damage_results.json  # Detailed damage statistics
│   └── *_damage_plots.png     # Damage pattern visualizations
├── plots/              # Quality score plots for each AB1 file
└── reports/            # Comprehensive QC HTML reports with interactive visualizations
```

## 🔍 Quality Control

The pipeline includes comprehensive quality control:

- **Quality Score Analysis**: Visual plots of Phred scores
- **N Content Tracking**: Monitor low-quality base replacement
- **Length Statistics**: Track sequence lengths through processing
- **Processing Status**: Monitor success/failure of each step
- **Interactive Reports**: HTML reports with summary statistics

## 🤝 Contributing

1. Fork the repository
1. Create a feature branch (`git checkout -b feature/amazing-feature`)
1. Commit changes (`git commit -m 'Add amazing feature'`)
1. Push to branch (`git push origin feature/amazing-feature`)
1. Open a Pull Request

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Support

For questions, issues, or contributions:

- **Issues**: Open an issue on GitHub
- **Documentation**: See inline code documentation
- **Examples**: Check the `examples/` directory (coming soon)

## 🏗️ Development

### Setting up Development Environment

```bash
# Install development dependencies
pip install -r requirements.txt

# Install pre-commit hooks (optional)
pip install pre-commit
pre-commit install

# Run code formatting
black src/ tests/
flake8 src/ tests/

# Type checking
mypy src/
```
