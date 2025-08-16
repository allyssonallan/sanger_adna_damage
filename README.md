# 🧬 Sanger aDNA Damage Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-sphinx-blue.svg)](docs/)

A comprehensive pipeline for processing Sanger sequencing data from ancient DNA (aDNA) samples, with automatic damage pattern analysis and quality assessment.

> [!IMPORTANT]
> **🚨 IMPORTANT DISCLAIMER - Tool Purpose & Limitations**
>
> This pipeline is **NOT** a tool for authenticating ancient DNA samples. It is designed for:
>
> - **Prioritizing haplogroups** for follow-up analysis
> - **Evaluating sample quality** based on insert size and damage patterns
> - **Providing surrogate bootstrapped damage indicators**
> - **Assisting in haplogroup origin assessment**
> - **Guiding selection** of promising samples for NGS sequencing
>
> **⚠️ All ancient DNA authentication must be performed using NGS-based methods** with appropriate controls, contamination assessment, and phylogenetic analysis.
>
> This tool provides preliminary screening to help researchers prioritize samples and resources before proceeding to comprehensive NGS-based ancient DNA authentication workflows.

## 🚀 Quick Start

### Installation

```bash
# Clone and install
git clone https://github.com/yourusername/sanger_adna_damage.git
cd sanger_adna_damage
pip install -r requirements.txt
```

### Basic Usage

```bash
# Run the complete pipeline
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output

# Convert single AB1 file
python -m src.sanger_pipeline.cli.main convert-ab1 \
    sample.ab1 output.fasta \
    --min-quality 20 \
    --min-sequence-length 30

# Generate QC report
python generate_report.py
```

## 📚 Documentation

📖 **[Complete Documentation](https://allysson.dev.br/sanger_adna_damage/)** - Comprehensive guides, tutorials, and API reference

**Quick Links:**

- [Installation Guide](https://allysson.dev.br/sanger_adna_damage/installation.html) - Detailed setup instructions
- [Usage Tutorial](https://allysson.dev.br/sanger_adna_damage/quickstart.html) - Step-by-step workflow guide  
- [Configuration](https://allysson.dev.br/sanger_adna_damage/configuration.html) - Customization options
- [API Reference](https://allysson.dev.br/sanger_adna_damage/api/) - Function and parameter documentation
- [Examples](https://allysson.dev.br/sanger_adna_damage/tutorials/) - Real-world use cases
- [Troubleshooting](https://allysson.dev.br/sanger_adna_damage/troubleshooting.html) - Common issues and solutions

## ✨ Key Features

- **📊 AB1 Processing** - Convert Sanger files to FASTQ with quality scores
- **🔗 Consensus Generation** - Merge forward/reverse reads intelligently  
- **📈 Quality Control** - Comprehensive QC reports with visualizations
- **🧪 Damage Analysis** - Bootstrap analysis of aDNA damage patterns
- **🔧 Automated Pipeline** - Single-command execution

## 🏗️ Development

```bash
# Development setup
git clone https://github.com/yourusername/sanger_adna_damage.git
cd sanger_adna_damage

# Set up virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests
python -m pytest tests/

# Generate documentation
cd docs && python -m sphinx.cmd.build source _build
```

See the [Contributing Guide](https://allysson.dev.br/sanger_adna_damage/contributing.html) for development workflows.

## 📝 License

MIT License - see [LICENSE](LICENSE) for details.

## 🆘 Support

- 📖 [Documentation](https://allysson.dev.br/sanger_adna_damage/) - Complete guides and references
- 🐛 [Issues](https://github.com/yourusername/sanger_adna_damage/issues) - Bug reports and feature requests
- 💬 [Discussions](https://github.com/yourusername/sanger_adna_damage/discussions) - Questions and community
