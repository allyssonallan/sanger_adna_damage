🧬 Complete Sanger aDNA Pipeline Execution Guide
================================================

This guide shows you how to run the complete Sanger ancient DNA damage analysis pipeline from start to finish.

## 📋 Prerequisites

1. **Environment Setup:**
   ```bash
   cd /Users/allyssonallan/sanger_adna_damage
   # Ensure Python environment is activated
   # Install dependencies if needed: pip install -r requirements.txt
   ```

2. **Input Data:**
   - AB1 files in the `input/` directory ✅
   - Reference sequence at `ref/rCRS.fasta` ✅
   - Configuration file at `config/default_config.yaml` ✅

## 🚀 Method 1: Complete Pipeline with CLI

### Step 1: Run the Complete Pipeline

```bash
sanger-pipeline run \
    --input-dir ./input \
    --output-dir ./output_q30 \
    --min-quality 30 \
    --config config/default_config.yaml
```

**What this does:**
- Converts all AB1 files to FASTA format
- Applies quality filtering (Q30)
- Builds consensus sequences for each HVS region
- Merges regions per sample
- Generates initial QC reports

### Step 2: Generate Comprehensive HTML Report

```bash
sanger-pipeline generate-report --output-dir ./output_q30
```

**Output:** `output_q30/reports/sanger_qc_report_YYYYMMDD_HHMMSS.html`

### Step 3: Convert to HSD Format for HaploGrep

```bash
# BWA-MEM based converter (recommended - uses proper alignment mapper)
sanger-pipeline hsd bwa \
    --consensus-dir ./output_q30/consensus/ \
    --output ./output_q30/haplogroups_bwa.hsd
```

## 🚀 Method 2: Alternative CLI Entry Points

### Using Python Module Direct Execution

If you prefer the module execution method:

```bash
python -m sanger_pipeline.cli.main run \
    --input-dir ./input \
    --output-dir ./output_q30 \
    --min-quality 30 \
    --config config/default_config.yaml
```

### Generate Reports

```bash
python -m sanger_pipeline.cli.main generate-report \
    --output-dir ./output_q30
```

### Convert to HSD

```bash
# BWA method (recommended)
sanger-pipeline hsd bwa \
    --consensus-dir ./output_q30/consensus/ \
    --output ./output_q30/haplogroups_bwa.hsd
```

## 🚀 Method 3: One-Command Complete Workflow

Create a simple shell script to run everything:

```bash
#!/bin/bash
# complete_pipeline.sh

echo "🧬 Starting Complete Sanger aDNA Pipeline..."

# Step 1: Run main pipeline
echo "📊 Step 1: Processing AB1 files..."
sanger-pipeline run \
    --input-dir ./input \
    --output-dir ./output_q30 \
    --min-quality 30 \
    --config config/default_config.yaml

# Step 2: Generate HTML report
echo "📋 Step 2: Generating HTML report..."
sanger-pipeline generate-report --output-dir ./output_q30

# Step 3: Convert to HSD with BWA method
echo "🧬 Step 3: Converting to HSD (BWA method)..."
sanger-pipeline hsd bwa \
    --consensus-dir ./output_q30/consensus/ \
    --output ./output_q30/haplogroups_bwa.hsd

echo "✅ Pipeline complete! Check outputs in ./output_q30/"
echo "📊 HTML Report: ./output_q30/reports/"
echo "🧬 HSD Files: ./output_q30/haplogroups_bwa.hsd"
```

## 📁 Expected Output Structure

After running the complete pipeline, you'll have:

```text
output_q30/
├── aligned/                    # Aligned sequences
├── consensus/                  # Consensus FASTA files
│   ├── *_HVS1_consensus.fasta
│   ├── *_HVS2_consensus.fasta
│   └── *_HVS3_consensus.fasta
├── damage_analysis/            # Damage pattern analysis
├── fasta/                      # Converted FASTA files
├── filtered/                   # Quality-filtered sequences
├── final/                      # Final merged sequences
├── plots/                      # Damage analysis plots
├── reports/                    # HTML QC reports
│   └── sanger_qc_report_*.html
└── haplogroups_bwa.hsd         # HSD file (BWA method)
```

## 🎯 Recommended HSD Converter Choice

**Use BWA-MEM method for all samples:**
- Uses proper BWA-MEM alignment mapping for accurate variant calling
- Handles complex indels and degraded sequences common in aDNA
- Conservative, high-confidence variant identification
- Optimized for both modern DNA and ancient DNA samples
- Single recommended approach after removal of less reliable methods

## 🔧 Additional Commands

### Check Pipeline Status

```bash
sanger-pipeline status --input-dir ./input
```

### Convert Single AB1 File

```bash
sanger-pipeline convert-ab1 \
    input/sample.ab1 \
    output/sample.fasta \
    --min-quality 30
```

### Validate Primer Pairs

```bash
sanger-pipeline validate-primers \
    input/sample.ab1 \
    --verbose
```

### Analyze Damage Patterns

```bash
sanger-pipeline analyze-damage \
    ./output_q30/consensus/ \
    ./output_q30/damage_analysis/
```

## 🚨 Troubleshooting

### Common Issues:

1. **Missing dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **BWA not found:**
   ```bash
   # BWA is available ✅
   which bwa
   ```

3. **Permission errors:**
   ```bash
   chmod +x scripts/run_pipeline.py
   ```

4. **Python path issues:**
   ```bash
   export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"
   ```

## 🎯 Next Steps

After pipeline completion:

1. **Review HTML QC Report** - Check quality metrics and damage patterns ✅
2. **Upload HSD to HaploGrep** - Use the .hsd file for haplogroup classification
3. **Analyze Results** - Review consensus sequences and variant calls
4. **Select Samples** - Use quality metrics to prioritize samples for NGS

## 📊 Recent Improvements

- ✅ **Dashboard/Report Generation Fixed** - The QC report system is now working correctly
- ✅ **BWA-Only Architecture** - Removed less reliable regional alignment methods
- ✅ **Enhanced CLI Interface** - Simplified commands with `sanger-pipeline` entry point
- ✅ **Proper Error Handling** - Better error messages and troubleshooting

The pipeline is now ready to process your ancient DNA samples with comprehensive quality control and reliable BWA-MEM based analysis!
