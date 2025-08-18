# 🚀 Quick Start: Complete Pipeline with HTML + HSD Outputs

This guide shows you how to run the complete Sanger aDNA pipeline to get both comprehensive HTML reports and HSD files for haplogroup analysis.

## 📋 Prerequisites

1. **Input data**: AB1 files in a directory (e.g., `./input/`)
2. **File naming**: Files should follow pattern `*_HVS*-F.ab1` and `*_HVS*-R.ab1`
3. **Dependencies**: Python 3.8+, BioPython, MAFFT installed

## 🔧 Complete Workflow (3 Steps)

### Step 1: Process AB1 Files → Consensus Sequences

```bash
# Run the main pipeline
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output_min30_q30 \
    --min-quality 30 \
    --min-sequence-length 30
```

**What this does:**
- Converts AB1 → FASTA with quality scores
- Filters low-quality sequences 
- Aligns forward/reverse reads (F+R only, no reference)
- Builds consensus sequences for each HVS region
- Performs damage analysis with bootstrap statistics

**Output directories created:**
```
output_min30_q30/
├── fasta/           # Raw FASTA conversions
├── filtered/        # Quality-filtered sequences  
├── aligned/         # F+R read alignments
├── consensus/       # HVS consensus sequences (49 files)
├── damage_analysis/ # Individual damage analysis JSONs
└── plots/          # Quality plots for each sample
```

### Step 2: Generate HTML QC Report

```bash
# Create comprehensive HTML report
python generate_report.py output_min30_q30
```

**What this creates:**
- **File**: `output_min30_q30/reports/sanger_qc_report_YYYYMMDD_HHMMSS.html`
- **Size**: ~8.6 MB (includes embedded plots and data)
- **Contents**: 
  - Sample processing statistics
  - Quality metrics dashboard
  - Interactive damage analysis plots
  - Bootstrap confidence intervals
  - Sample prioritization recommendations

**View the report:**
```bash
# Open in browser (macOS)
open output_min30_q30/reports/sanger_qc_report_*.html

# Open in browser (Linux)
xdg-open output_min30_q30/reports/sanger_qc_report_*.html
```

### Step 3: Create HSD File for Haplogroup Analysis

```bash
# Convert consensus sequences to HSD format
python convert_hvs_consensus_to_hsd.py output_min30_q30/consensus/ my_samples.hsd
```

**What this creates:**
- **File**: `my_samples.hsd` (7KB)
- **Format**: HaploGrep-compatible HSD format
- **Contents**: 
  - Sample IDs and variant positions
  - Regional variant mapping (HVS1, HVS2, HVS3)
  - Conservative variant calling (10-35 variants per region)

**Upload to HaploGrep:**
1. Go to https://haplogrep.i-med.ac.at/
2. Upload your `my_samples.hsd` file
3. Get haplogroup classifications

## 📊 Expected Results

After running all three steps, you'll have:

### ✅ HTML Report (`output_min30_q30/reports/`)
- **Sample statistics**: Processing success rates, quality metrics
- **Damage analysis**: C→T and G→A transition patterns  
- **Quality assessment**: Sequence lengths, Phred scores
- **Visualization**: Interactive plots with damage patterns
- **Recommendations**: Which samples to prioritize for NGS

### ✅ HSD File (`my_samples.hsd`)
- **Format**: Ready for HaploGrep analysis
- **Content**: ~24 samples with regional variants
- **Quality**: Scientifically accurate variant calls
- **Usage**: Direct upload to haplogroup classification tools

## 🎯 One-Line Complete Workflow

```bash
# Complete pipeline: AB1 → HTML + HSD
python -m src.sanger_pipeline.cli.main run-pipeline --input-dir ./input --output-dir ./output_min30_q30 --min-quality 30 --min-sequence-length 30 && python generate_report.py output_min30_q30 && python convert_hvs_consensus_to_hsd.py output_min30_q30/consensus/ my_samples.hsd && echo "✅ Complete! Check output_min30_q30/reports/ for HTML and my_samples.hsd for HaploGrep"
```

## 🔧 Customization Options

### Adjust Quality Thresholds

```bash
# More stringent quality (fewer but higher quality sequences)
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output_min40_q35 \
    --min-quality 35 \
    --min-sequence-length 40

# More permissive quality (more sequences, potentially lower quality)  
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output_min20_q25 \
    --min-quality 25 \
    --min-sequence-length 20
```

### Use Custom Configuration

```bash
# Copy and modify default config
cp config/default_config.yaml my_config.yaml
# Edit my_config.yaml as needed

# Run with custom config
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output \
    --config my_config.yaml
```

## 📁 File Organization

```
your_project/
├── input/                           # Your AB1 files
│   ├── sample1_HVS1-F.ab1
│   ├── sample1_HVS1-R.ab1
│   └── ...
├── output_min30_q30/                # Pipeline outputs
│   ├── consensus/                   # 49 consensus FASTA files
│   ├── reports/                     # HTML QC report
│   └── ...
├── my_samples.hsd                   # HSD file for HaploGrep
└── analysis_notes.txt               # Your analysis notes
```

## 🆘 Troubleshooting

### No consensus files generated
- Check AB1 file naming (must include HVS and -F/-R)
- Verify quality thresholds aren't too strict
- Check `output_*/logs/` for error messages

### HTML report is empty
- Ensure pipeline completed successfully first
- Check that consensus files exist in `output_*/consensus/`

### HSD file has no samples
- Verify consensus files have sufficient length (>50bp)
- Check for proper HVS region identification in filenames

## 🎉 Success Indicators

✅ **Pipeline Success**: Console shows "✓" marks for each step  
✅ **HTML Report**: Large file (~8MB) with interactive content  
✅ **HSD File**: Contains sample IDs and variant lists  
✅ **Ready for Analysis**: Both outputs ready for downstream analysis

Your pipeline is now set up with the optimal regionalized approach for ancient DNA analysis!
