# Minimum Sequence Length Implementation

## Overview

This document describes the implementation of minimum sequence length filtering (≥30bp) in the Sanger DNA damage analysis pipeline.

## Key Implementation Details

### 1. Sequence Length Filtering Location

**WHERE**: Length filtering occurs **AFTER** quality filtering but **BEFORE** consensus building
**WHY**: This ensures we exclude sequences that become too short after quality degradation

### 2. Strand Alignment Process

✅ **CORRECT IMPLEMENTATION**: The reverse complement is created **BEFORE** alignment

```python
# In consensus_builder.py - process_paired_reads()
def process_paired_reads(self, forward_file, reverse_file, alignment_output, consensus_output, sample_name):
    # 1. Create reverse complement of reverse read FIRST
    reverse_rc_file = reverse_file.parent / f"{reverse_file.stem}_rc.fasta"
    self.reverse_complement_sequence(reverse_file, reverse_rc_file)
    
    # 2. THEN align forward + reverse_complement
    self.align_sequences([forward_file, reverse_rc_file], alignment_output)
    
    # 3. Build consensus from aligned sequences
    consensus_record = self.build_consensus(alignment_output, consensus_output, ...)
```

### 3. Processing Workflow with Length Filtering

```
AB1 Files
    ↓
1. AB1 → Raw FASTA (ab1_converter.py)
    ↓
2. Quality Filtering (replace Q<threshold with 'N')
    ↓
3. Length Validation (count A,T,C,G ≥ min_sequence_length)
    ↓ (EXCLUDE if too short)
4. Paired Read Processing:
   - Forward FASTA + Reverse FASTA
   - Reverse complement creation
   - MAFFT alignment
   - Consensus building
    ↓
5. HVS Region Merging
    ↓
6. Ancient DNA Damage Analysis
```

## Files Modified

### 1. `src/sanger_pipeline/utils/constants.py`
- Added `DEFAULT_MIN_SEQUENCE_LENGTH = 30`

### 2. `config/default_config.yaml`
- Added `min_sequence_length: 30` under quality section

### 3. `src/sanger_pipeline/core/quality_filter.py`
- Added `min_sequence_length` parameter to `__init__()`
- Added `validate_sequence_length()` method
- Added `filter_sequence_with_length_check()` method

### 4. `src/sanger_pipeline/core/ab1_converter.py`
- Added `min_sequence_length` parameter to `__init__()`
- Modified `filter_by_quality()` to return `Optional[SeqRecord]`
- Added length validation in filtering process
- Updated `process_ab1_file()` return type

### 5. `src/sanger_pipeline/core/pipeline.py`
- Added `min_sequence_length` parameter to `__init__()`
- Updated initialization to use config values
- Modified `_step_1_convert_ab1_files()` to handle excluded sequences
- Added statistics reporting for filtered vs excluded sequences

### 6. `src/sanger_pipeline/cli/main.py`
- Added `--min-sequence-length` CLI option
- Updated both `run_pipeline` and `convert_ab1` commands

## Usage Examples

### Command Line Usage

```bash
# Run pipeline with default 30bp minimum
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output

# Run with custom minimum length
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ./input \
    --output-dir ./output \
    --min-sequence-length 50

# Convert single file with length checking
python -m src.sanger_pipeline.cli.main convert-ab1 \
    sample.ab1 output.fasta \
    --min-sequence-length 30
```

### Configuration File Usage

```yaml
# config/default_config.yaml
quality:
  min_phred_score: 15
  min_sequence_length: 30  # New setting
  terminal_length: 10
```

### Programmatic Usage

```python
from sanger_pipeline.core.pipeline import SangerPipeline

pipeline = SangerPipeline(
    input_dir="./input",
    output_dir="./output",
    min_quality=20,
    min_sequence_length=30  # New parameter
)

pipeline.run()
```

## Length Validation Logic

The length validation counts **only valid nucleotides**:

```python
def validate_sequence_length(self, sequence: str) -> bool:
    # Count only A, T, C, G (not N's or gaps)
    valid_bases = sum(1 for base in sequence if base in 'ATCG')
    return valid_bases >= self.min_sequence_length
```

## Impact on Pipeline

### Sequences That Pass Filtering
- Continue through normal pipeline workflow
- Generate consensus sequences
- Included in final results

### Sequences That Are Excluded
- Raw FASTA files are still created (for diagnostics)
- Filtered FASTA files are NOT created (or deleted if empty)
- Quality plots are still generated (for quality assessment)
- Excluded from consensus building
- Pipeline reports statistics: "X passed filtering, Y excluded (too short)"

## Benefits

1. **Quality Control**: Prevents unreliable consensus from very short sequences
2. **Ancient DNA Compatibility**: 30bp minimum suitable for degraded aDNA
3. **Configurable**: Can be adjusted via CLI or config file
4. **Early Filtering**: Removes short sequences before expensive alignment operations
5. **Preserves Raw Data**: Original data still available for analysis
6. **Clear Reporting**: Users know exactly how many sequences were excluded

## Validation

The implementation has been tested with:
- ✅ Sequences above minimum length (40bp) → Pass
- ✅ Sequences below minimum length (20bp) → Excluded
- ✅ Sequences that become too short after quality filtering → Excluded
- ✅ Configuration loading from YAML
- ✅ CLI parameter handling
- ✅ Statistics reporting

## Reference Alignment

**IMPORTANT**: The reference sequence is NOT reversed. The alignment works as follows:

1. **Forward read**: Used as-is
2. **Reverse read**: Converted to reverse complement BEFORE alignment
3. **Alignment**: MAFFT aligns Forward + Reverse-complement
4. **Reference**: Remains in standard orientation for damage analysis

This ensures proper alignment regardless of sequencing direction.
