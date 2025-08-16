# Implementation Documentation

This directory contains detailed technical documentation about the implementation of various features in the Sanger aDNA Damage Pipeline.

## Files

### Core Implementation Analysis

- **[`STRAND_ALIGNMENT_ANALYSIS.md`](STRAND_ALIGNMENT_ANALYSIS.md)** - Complete analysis of strand alignment workflow, reverse complement timing, and reference sequence handling
- **[`SEQUENCE_LENGTH_IMPLEMENTATION.md`](SEQUENCE_LENGTH_IMPLEMENTATION.md)** - Detailed documentation of minimum sequence length filtering implementation (≥30bp)
- **[`IMPLEMENTATION_SUMMARY.md`](IMPLEMENTATION_SUMMARY.md)** - High-level summary of all implementation changes and validation results

## Purpose

These documents provide:

1. **Technical Details** - In-depth analysis of how the pipeline processes sequences
2. **Implementation Evidence** - Code locations and validation of processing logic
3. **Workflow Documentation** - Step-by-step breakdown of the complete pipeline
4. **Developer Reference** - Technical specifications for contributors and maintainers

## Target Audience

- Pipeline developers and contributors
- Researchers wanting to understand the technical details
- Users implementing custom modifications
- Code reviewers and auditors

## Key Findings Summary

### Strand Processing (✅ Verified Correct)
- Reverse complement is created **BEFORE** alignment
- Alignment occurs **BEFORE** consensus merging  
- Reference sequence is **NOT** reversed for damage analysis

### Sequence Length Filtering (✅ Implemented)
- Minimum 30bp filter implemented after quality filtering
- Configurable via CLI (`--min-sequence-length`) and config file
- Early exclusion saves computation on unusable sequences

### Workflow Order (✅ Confirmed)
```
AB1 → Quality Filter → Length Check → Alignment → Consensus → Region Merge → Damage Analysis
```

## Related Documentation

For user-focused documentation, see the main [docs/source/](../source/) directory.
