<!-- README.md for AB1 Sequencing QC Pipeline -->
# AB1 Sanger Sequencing Quality Control Pipeline

This repository contains a complete pipeline for processing Sanger sequencing AB1 files, including base-calling with Phred quality filtering, alignment and consensus building, HVSI/HVSII merging, and generation of an interactive HTML QC report.

## Features

- Convert raw `.ab1` chromatograms to FASTA sequences
- Filter low-quality bases (Phred < 20) and plot quality profiles
- Align forward and reverse reads and build consensus sequences
- Merge HVSI and HVSII regions into a single sequence per sample
- Generate an HTML QC report with summary tables and plots via Quarto

## Requirements

### System

- Python 3.6 or higher
- R 4.0 or higher
- MAFFT (for sequence alignment)
- Quarto (for report rendering)
- Git

### Python Dependencies

Install and activate a virtual environment, then install with pip:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### R Dependencies

Manage R packages using `renv`. On first run, `renv_manager.R` will initialize the environment and create a `renv.lock` file to capture package versions; on subsequent runs, it will restore the environment from `renv.lock`.

```bash
Rscript renv_manager.R
```

## Usage

Clone the repository and run the main pipeline script:

```bash
git clone <repository_url>
cd <repository_dir>

# Python setup (create and activate a virtual environment; venv/ is ignored by .gitignore)
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# R setup (initialize or restore the R environment; only renv.lock and renv/settings are tracked)
Rscript renv_manager.R

# Run the main AB1 processing pipeline
bash 00_process_ab1.sh

# View QC report
open qc_report.html  # or browse in your web browser
```

## Damage Analysis (optional)
After generating consensus sequences, perform bootstrap-based damage estimation:
```bash
Rscript 04_bootstrap_damage.R
```

## Pipeline Steps

1. **Conversion & Filtering** (`01_convert_ab1_quality.py`)
   - Extract sequences and Phred quality from `.ab1` files
   - Produce raw and filtered FASTA, quality plots saved in `plots/`
2. **Alignment & Consensus** (`02_make_consensus.py`)
   - Reverse-complement reverse reads, align with MAFFT
   - Generate consensus FASTA per sample
3. **Merge HVSI/HVSII**
   - Concatenate HVSI and HVSII consensus into `final/`
4. **QC Report** (`03_qc_report.qmd`)
   - Render HTML report summarizing sample metrics and plots
5. **Damage Analysis** (`04_bootstrap_damage.R`)
   - Perform bootstrap-based damage estimation on consensus sequences
   - Outputs summary statistics and plots in `damage_results/`

## Directory Structure

```{bash}
.
├── 00_process_ab1.sh           # Main Bash pipeline script
├── 01_convert_ab1_quality.py   # Convert .ab1 to FASTA and filter by quality
├── 02_make_consensus.py        # Build consensus from alignments
├── 03_qc_report.qmd            # Quarto notebook for QC report
├── 04_bootstrap_damage.R       # R script for damage bootstrapping (optional)
├── requirements.txt            # Python package requirements
├── renv_manager.R              # R environment manager via renv
├── renv.lock                   # Lockfile for R package versions
├── venv/                       # Python virtual environment (ignored by .gitignore)
├── fasta/                      # Raw FASTA outputs
├── filtered/                   # Quality-filtered FASTA
├── consensus/                  # Per-sample consensus FASTA
├── aligned/                    # MAFFT alignment files
├── final/                      # Merged HVSI/HVSII sequences
├── damage_results/             # Output directory for damage analysis (bootstrapping)
├── logs/                       # Pipeline log files
└── plots/                      # Quality and summary plots
```

## Contributing

- Please submit bug reports and pull requests via GitHub issues.
- Follow the existing code style (Bash scripts and Python 3).

## License

This project is released under the [MIT License](LICENSE).

## GitHub Publishing Best Practices

- Use a `.gitignore` file to prevent committing unnecessary or sensitive files (see the provided `.gitignore`).
- Never commit credentials or environment files (e.g., `.env`, `.Renviron`, keys, certificates).
- Write clear, descriptive commit messages and pull request titles/descriptions.
- Protect the main branch with branch protection rules; require reviews before merging.
- Automate testing and linting using GitHub Actions or other CI services.
- Keep dependencies up-to-date; consider enabling Dependabot or similar tools.
- Include a `SECURITY.md` or specify security policies for handling issues and vulnerabilities.
- Tag releases and maintain a `CHANGELOG.md` for project history.
- Consider adding a `CODEOWNERS` file to define responsible maintainers.
