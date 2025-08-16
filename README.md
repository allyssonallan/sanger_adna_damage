# 🧬 Sanger aDNA Damage Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-sphinx-blue.svg)](docs/)

A comprehensive pipeline for processing Sanger sequencing data from ancient DNA (aDNA) samples, with automatic damage pattern analysis and quality assessment.

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
# Process your AB1 files in 4 simple steps
./00_process_ab1.sh                    # Process AB1 files
python 01_convert_ab1_quality.py       # Convert to FASTQ
python 02_make_consensus.py            # Generate consensus
quarto render 03_qc_report.qmd         # Create QC report
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
pip install -r requirements.txt
Rscript renv_manager.R
```

See the [Contributing Guide](https://allysson.dev.br/sanger_adna_damage/contributing.html) for development workflows.

## 📝 License

MIT License - see [LICENSE](LICENSE) for details.

## 🆘 Support

- 📖 [Documentation](https://allysson.dev.br/sanger_adna_damage/) - Complete guides and references
- 🐛 [Issues](https://github.com/yourusername/sanger_adna_damage/issues) - Bug reports and feature requests
- 💬 [Discussions](https://github.com/yourusername/sanger_adna_damage/discussions) - Questions and community
