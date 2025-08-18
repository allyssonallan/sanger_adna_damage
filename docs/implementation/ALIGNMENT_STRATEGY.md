"""
Reference Alignment Strategy Analysis for Sanger aDNA Pipeline

This document analyzes the optimal integration of reference alignment in the workflow.
"""

# STRATEGIC ANALYSIS: When to align with reference

## Option 1: Early Reference Alignment (Pre-Consensus)
### When: After quality filtering, before consensus generation
### Process: Raw reads → Reference alignment → Consensus from aligned reads → HSD

PROS:
✅ More accurate consensus sequences
✅ Better handling of indels and structural variants
✅ Improved positioning of sequences within HVS regions
✅ Can detect and correct orientation issues early
✅ Reference-guided consensus reduces errors

CONS:
❌ More complex pipeline
❌ Requires additional alignment step
❌ May introduce reference bias in ancient samples
❌ Could mask authentic ancient variants

BEST FOR:
- High-quality modern samples
- Cases where precise positioning is critical
- When reference bias is acceptable


## Option 2: Late Reference Alignment (Post-Consensus)
### When: After consensus generation, before HSD conversion
### Process: Raw reads → Consensus → Reference alignment → HSD

PROS:
✅ Preserves authentic sequence variants in consensus
✅ Simpler pipeline (current approach)
✅ Less reference bias
✅ Maintains independence of consensus step

CONS:
❌ May have positioning errors in consensus
❌ Harder to handle structural variants
❌ Less precise HVS region boundaries

BEST FOR:
- Ancient DNA samples (current use case)
- When preserving authentic variants is priority
- Simpler workflow requirements


## Option 3: Dual Strategy (Conditional Alignment)
### When: Different strategies based on sample quality/type
### Process: Quality assessment → Choose alignment strategy → Proceed accordingly

PROS:
✅ Optimal for each sample type
✅ Flexible approach
✅ Can handle diverse sample qualities

CONS:
❌ Complex decision logic
❌ Harder to maintain and test
❌ More parameters to configure

BEST FOR:
- Mixed sample types
- Research environments with varying data quality


## Option 4: Hybrid Approach (Multiple Reference Points)
### When: Use reference at multiple stages
### Process: Pre-alignment guidance + Post-consensus validation

PROS:
✅ Best of both worlds
✅ Quality control at multiple stages
✅ Can detect and correct errors

CONS:
❌ Most complex approach
❌ Computational overhead
❌ May be overkill for current needs


# RECOMMENDATION FOR CURRENT PIPELINE

## Primary Recommendation: Enhanced Option 2 (Smart Post-Consensus Alignment)

### Implementation Strategy:

1. **Keep current consensus generation** (preserves authentic ancient DNA variants)

2. **Add reference-aware HSD conversion** with:
   - Sequence alignment before variant calling
   - Smart positioning within HVS regions
   - Quality-based alignment parameters

3. **Add alignment quality assessment** to detect problematic sequences

### Technical Implementation:

```python
class ReferenceAwareHSDConverter:
    def __init__(self, reference_file, alignment_tool="mafft"):
        self.reference_seq = self.load_reference(reference_file)
        self.alignment_tool = alignment_tool
        self.hvs_regions = {...}
    
    def find_variants_with_alignment(self, sample_seq, hvs_region):
        # 1. Extract reference HVS region
        ref_region = self.get_hvs_reference(hvs_region)
        
        # 2. Align sample sequence to reference region
        alignment = self.align_to_reference(sample_seq, ref_region)
        
        # 3. Extract variants from alignment
        variants = self.extract_variants_from_alignment(alignment, hvs_region)
        
        return variants
    
    def align_to_reference(self, sample_seq, ref_seq):
        # Create temporary files for alignment
        # Run MAFFT/similar to align sample to reference
        # Return alignment object
        pass
    
    def extract_variants_from_alignment(self, alignment, hvs_region):
        # Parse alignment to find true variants
        # Handle indels properly
        # Map to correct mitochondrial positions
        pass
```

### Benefits of This Approach:

1. **Preserves authenticity**: Consensus sequences maintain original ancient DNA signatures
2. **Improves accuracy**: Reference alignment catches positioning errors in HSD conversion
3. **Handles indels**: Proper alignment can detect insertions/deletions
4. **Maintains simplicity**: Adds complexity only where needed (HSD stage)
5. **Backward compatible**: Existing consensus files still work

### When to Use Each Approach:

**Use Current Method (No Reference Alignment) When:**
- Ancient DNA samples with expected damage
- Quick screening/prioritization
- Low computational resources
- Simple variant detection sufficient

**Use Enhanced Method (Reference-Aware HSD) When:**
- Publishing results
- Detailed haplogroup analysis
- Samples with complex variants
- Maximum accuracy required

**Future Enhancement (Pre-Consensus Alignment) When:**
- Modern DNA samples
- Known high-quality sequences
- Structural variant detection needed
- Reference genomes are highly trusted

### Implementation Priority:

1. **Phase 1**: Enhance current HSD converter with optional reference alignment
2. **Phase 2**: Add alignment quality metrics and reporting
3. **Phase 3**: Consider pre-consensus alignment for specific use cases

This strategy maintains the pipeline's focus on ancient DNA while providing enhanced accuracy when needed.
