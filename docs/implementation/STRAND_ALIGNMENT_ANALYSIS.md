# 🧬 Strand Alignment Analysis: Complete Workflow

## 📋 **CRITICAL FINDINGS**

### ✅ **1. Strand Alignment Timing: BEFORE Consensus Merge**

**ALIGNMENT HAPPENS FIRST, THEN CONSENSUS, THEN MERGING**

```
Step 1: AB1 → Raw FASTA (quality filtering)
Step 2: Strand Alignment + Consensus Building  ← ALIGNMENT HAPPENS HERE
Step 3: HVS Region Merging                     ← MERGING HAPPENS HERE
Step 4: Damage Analysis
```

### ✅ **2. Reverse Complement Timing: BEFORE Alignment**

**Location**: `consensus_builder.py` lines 235-250

```python
def process_paired_reads(self, forward_file, reverse_file, alignment_output, consensus_output, sample_name):
    # STEP 1: Create reverse complement BEFORE alignment
    reverse_rc_file = reverse_file.parent / f"{reverse_file.stem}_rc.fasta"
    self.reverse_complement_sequence(reverse_file, reverse_rc_file)
    
    # STEP 2: Align Forward + Reverse-complement
    self.align_sequences([forward_file, reverse_rc_file], alignment_output)
    
    # STEP 3: Build consensus from alignment
    consensus_record = self.build_consensus(alignment_output, consensus_output, ...)
```

### ✅ **3. Reference Sequence: NOT REVERSED**

**Location**: `adna_damage_analyzer.py` lines 46-55

```python
def analyze_sequence_damage(self, sequence_file: Path, reference_file: Path):
    # Read sequences AS-IS (reference NOT reversed)
    query_seq = SeqIO.read(sequence_file, "fasta")
    ref_seq = SeqIO.read(reference_file, "fasta")
    
    # Align query TO reference (reference stays in standard orientation)
    alignment = self._align_sequences(str(ref_seq.seq), str(query_seq.seq))
```

## 🔄 **Complete Workflow Analysis**

### **Step-by-Step Process:**

#### **Step 1: AB1 Conversion**
```
sample_HVS1-F.ab1 → sample_HVS1-F.fasta
sample_HVS1-R.ab1 → sample_HVS1-R.fasta
```

#### **Step 2: Alignment & Consensus (PER HVS REGION)**
```
FOR EACH HVS REGION (HVS1, HVS2, HVS3):
  1. Forward read: sample_HVS1-F.fasta (as-is)
  2. Reverse read: sample_HVS1-R.fasta → REVERSE COMPLEMENT
  3. Alignment: MAFFT(Forward + RC_Reverse) → aligned.fasta
  4. Consensus: aligned.fasta → sample_HVS1_consensus.fasta
```

#### **Step 3: Region Merging**
```
sample_HVS1_consensus.fasta
sample_HVS2_consensus.fasta  → CONCATENATE → sample_HVS1_HVS2_HVS3_merged.fasta
sample_HVS3_consensus.fasta
```

#### **Step 4: Damage Analysis**
```
sample_consensus.fasta + reference.fasta → damage_patterns
(Both sequences in standard orientation)
```

## 🎯 **Key Technical Details**

### **Reverse Complement Implementation**
**Location**: `consensus_builder.py` lines 35-55

```python
def reverse_complement_sequence(self, input_file: Path, output_file: Path):
    record = SeqIO.read(input_file, "fasta")
    record.seq = record.seq.reverse_complement()  # BioPython reverse complement
    record.id += "_RC"
    record.description = "Reverse-complemented"
    SeqIO.write(record, output_file, "fasta")
```

### **Alignment Tool Configuration**
**Default**: MAFFT with `--auto` parameters
**Location**: `consensus_builder.py` lines 57-95

```python
def align_sequences(self, sequence_files: List[Path], output_file: Path):
    # Combine Forward + Reverse-complement into temp file
    # Run: mafft --auto temp_file.fasta > alignment_output.fasta
    cmd = [self.alignment_tool] + self.alignment_params.split() + [str(temp_file)]
```

### **Reference Alignment for Damage Analysis**
**Location**: `adna_damage_analyzer.py` lines 68-89

```python
def _align_sequences(self, reference: str, query: str):
    # Global alignment: reference (standard) vs query (consensus)
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(reference, query)  # ref FIRST, query SECOND
```

## 📊 **Alignment Strategy Summary**

| Stage | Forward Read | Reverse Read | Reference | Notes |
|-------|-------------|--------------|-----------|-------|
| **Input** | As-is | As-is | Standard orientation | Raw AB1 data |
| **Pre-alignment** | As-is | **REVERSE COMPLEMENT** | Not involved | Prepare for alignment |
| **Consensus alignment** | Forward | RC_Reverse | Not involved | MAFFT alignment |
| **Damage analysis** | Not used | Not used | **Standard orientation** | Query vs Reference |

## ✅ **Validation of Approach**

### **Why Reverse Complement Before Alignment?**
1. **Sanger sequencing**: Forward and reverse reads are from opposite strands
2. **Proper orientation**: RC makes both reads face same direction
3. **Accurate consensus**: Alignment can properly overlap sequences
4. **Standard practice**: Matches established bioinformatics workflows

### **Why Reference NOT Reversed?**
1. **Standard orientation**: References (like rCRS) are in canonical orientation
2. **Damage patterns**: C→T transitions are directional and position-specific
3. **Consistency**: All analyses use same reference orientation
4. **Comparison**: Allows comparison with published damage patterns

## 🎯 **ANSWERS TO YOUR QUESTIONS**

### **Q1: Is strand alignment before or after consensus merge?**
**A1**: ✅ **BEFORE** - Alignment happens in Step 2, merging happens in Step 3

### **Q2: When is reverse complement created - before or after alignment?**
**A2**: ✅ **BEFORE** - RC is created first, then aligned with forward read

### **Q3: Is the reference reversed for alignment?**
**A3**: ✅ **NOT REVERSED** - Reference stays in standard orientation

## 🔍 **Code Evidence Summary**

- **Pipeline order**: `run()` calls `_step_2_align_and_consensus()` then `_step_3_merge_regions()`
- **RC timing**: `process_paired_reads()` calls `reverse_complement_sequence()` before `align_sequences()`
- **Reference orientation**: `analyze_sequence_damage()` reads reference as-is, no reversal

**The implementation is scientifically sound and follows best practices for paired-end sequence processing and ancient DNA damage analysis.**
