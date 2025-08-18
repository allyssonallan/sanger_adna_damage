# Pipeline Integration Strategy for Regional HSD Conversion

## Summary of Analysis

After testing three different approaches for HVS consensus to HSD conversion:

1. **❌ Aligned Method**: Too many artifacts (90+ variants/sample)
2. **⚠️ Direct Method**: Good but concatenates regions (66 variants/sample) 
3. **✅ Hybrid Regional**: Best balance (52 variants/sample)**

The **Hybrid Regional approach** provides the optimal solution by:
- Processing each HVS region independently  
- Using direct comparison (no alignment artifacts)
- Maintaining regional specificity
- Producing scientifically reasonable variant counts

## Pipeline Integration Plan

### Phase 1: Add Regional HSD Converter to Pipeline

#### 1. Move hybrid converter to pipeline utils
```bash
cp hybrid_regional_hsd_converter.py src/sanger_pipeline/utils/regional_hsd_converter.py
```

#### 2. Update pipeline to include HSD conversion step
- Add step 7: Regional HSD conversion  
- Process consensus files by region
- Generate HSD output for HaploGrep

#### 3. Add CLI command
```bash
python -m src.sanger_pipeline.cli.main convert-to-hsd \
    --consensus-dir output/consensus/ \
    --output samples.hsd
```

### Phase 2: Integration with Report Generator

#### 1. Add HSD conversion results to HTML report
- Summary of HSD conversion
- Number of variants per region  
- Link to download HSD file
- Instructions for HaploGrep usage

#### 2. Add HSD conversion button to dashboard
- "Convert to HSD" button
- Real-time conversion progress
- Download link for results

### Phase 3: Automatic Pipeline Integration

#### 1. Make HSD conversion optional pipeline step
```yaml
# config.yaml
hsd_conversion:
  enabled: true
  method: "regional"  # regional, direct, or aligned
  max_variants_per_region: 35
```

#### 2. Add to main pipeline flow
```python
def _step_7_hsd_conversion(self) -> None:
    """Step 7: Convert consensus sequences to HSD format for haplogroup analysis."""
    if self.config.get('hsd_conversion', {}).get('enabled', False):
        # Run regional HSD conversion
        # Add results to pipeline output
```

## Implementation Priority

### **High Priority (Immediate)**
1. ✅ **Regional HSD converter implemented and tested**
2. 🔄 **Move to pipeline utils directory** 
3. 🔄 **Add CLI command for standalone usage**
4. 🔄 **Update documentation with usage examples**

### **Medium Priority (Next Sprint)**  
1. 🔍 **Integration with HTML report generator**
2. 🔍 **Add to main pipeline as optional step**
3. 🔍 **Configuration file integration**

### **Low Priority (Future Enhancement)**
1. 🔮 **Dashboard integration with real-time conversion**
2. 🔮 **Batch processing capabilities**
3. 🔮 **Integration with HaploGrep API (if available)**

## Technical Benefits

### **Regional Processing Advantages**
- **Reduced alignment artifacts**: Each region aligned independently
- **Better variant quality**: Direct comparison eliminates spurious indels
- **Biological relevance**: Maintains HVS region boundaries
- **Performance**: Faster processing than full-genome alignment

### **Pipeline Integration Benefits** 
- **Seamless workflow**: One command from AB1 to HSD
- **Consistent results**: Standardized processing across all samples
- **Quality control**: Built-in validation and error handling
- **Documentation**: Automatic reporting of conversion results

## Usage Examples

### **Standalone Usage**
```bash
# Convert existing consensus files
python hybrid_regional_hsd_converter.py output/consensus/ samples.hsd

# With custom reference
python hybrid_regional_hsd_converter.py data/consensus/ results.hsd --reference custom_ref.fasta
```

### **Pipeline Integration** (Future)
```bash
# Full pipeline with HSD conversion
python -m src.sanger_pipeline.cli.main run-pipeline \
    --input-dir ab1_files/ \
    --output-dir output/ \
    --enable-hsd-conversion

# HSD conversion only
python -m src.sanger_pipeline.cli.main convert-to-hsd \
    --consensus-dir output/consensus/ \
    --output samples.hsd
```

## Next Steps

1. **✅ Complete**: Hybrid regional converter implementation and testing
2. **🔄 Current**: Move converter to pipeline structure  
3. **📋 Next**: Add CLI command and documentation
4. **🎯 Goal**: Seamless AB1 → HSD workflow for haplogroup analysis

The regional approach provides the optimal balance of accuracy, performance, and biological relevance for ancient DNA haplogroup classification.
