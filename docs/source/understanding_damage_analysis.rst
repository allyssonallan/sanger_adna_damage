==============================
Understanding Damage Analysis
==============================

This guide provides a comprehensive understanding of ancient DNA damage analysis in the Sanger DNA Damage Analysis Pipeline, covering the scientific background, methods, and interpretation of results.

🧬 Scientific Background
========================

Ancient DNA and Degradation
---------------------------

Ancient DNA (aDNA) undergoes characteristic chemical degradation over time, primarily through:

**Depurination**
  Loss of purine bases (A and G) leaving apurinic sites that can cause strand breaks during PCR amplification.

**Cytosine Deamination**
  Spontaneous hydrolytic deamination of cytosine to uracil, which is read as thymine during PCR, resulting in C→T transitions.

**5-methylcytosine Deamination**
  Deamination of 5-methylcytosine to thymine, also causing C→T transitions, particularly common in CpG dinucleotides.

Characteristic Damage Patterns
------------------------------

Ancient DNA shows distinctive damage patterns:

1. **5' C→T Transitions**: High frequency of C→T misincorporations at the 5' ends of sequences
2. **3' G→A Transitions**: High frequency of G→A misincorporations at the 3' ends of sequences  
3. **Position-Dependent Damage**: Damage rates decrease with distance from sequence ends
4. **Strand Asymmetry**: Different damage patterns on plus and minus strands

.. note::
   These patterns are caused by depurination during DNA extraction and library preparation, where ancient templates with apurinic sites are converted to C→T and G→A misincorporations.

📊 Damage Detection Methods
===========================

Position-Based Analysis
-----------------------

The pipeline analyzes damage by examining base transitions at different positions:

**5' End Analysis** (positions 1-20):
  * Count C→T transitions in each position
  * Calculate transition frequencies
  * Compare against background rates

**3' End Analysis** (positions -20 to -1):
  * Count G→A transitions in each position  
  * Calculate transition frequencies
  * Compare against background rates

**Middle Positions** (control):
  * Calculate baseline transition rates
  * Used for comparison with end positions

Damage Score Calculation
------------------------

The damage score is computed as:

.. math::

   \text{Damage Score} = \frac{(\text{5' C→T rate} + \text{3' G→A rate})}{2} - \text{Background rate}

Where:
- **5' C→T rate**: Average C→T frequency in positions 1-5
- **3' G→A rate**: Average G→A frequency in positions -5 to -1  
- **Background rate**: Average transition rate in middle positions

**Score Interpretation**:
- **0.0-0.2**: Minimal damage (modern DNA)
- **0.2-0.4**: Low-moderate damage (recent/well-preserved)
- **0.4-0.7**: High damage (ancient DNA likely)
- **>0.7**: Very high damage (definite ancient DNA)

🎲 Statistical Validation
=========================

Bootstrap Analysis
-----------------

The pipeline uses bootstrap resampling to assess statistical significance:

1. **Resampling**: Create 10,000 bootstrap samples from the original data
2. **Recalculation**: Calculate damage score for each bootstrap sample
3. **Distribution**: Build distribution of bootstrap damage scores
4. **P-value**: Calculate probability of observing score by chance

**Bootstrap Process**:

.. code-block:: python

   for i in range(10000):
       bootstrap_sample = resample(original_sequences)
       bootstrap_score = calculate_damage_score(bootstrap_sample)
       bootstrap_scores.append(bootstrap_score)
   
   p_value = sum(score >= observed_score for score in bootstrap_scores) / 10000

Significance Testing
-------------------

**Null Hypothesis**: Observed damage patterns are due to random sequencing errors

**Alternative Hypothesis**: Observed damage patterns indicate authentic ancient DNA

**P-value Interpretation**:
- **p < 0.01**: Highly significant ancient DNA damage
- **p < 0.05**: Significant ancient DNA damage  
- **p < 0.10**: Marginally significant
- **p ≥ 0.10**: Not significant (modern DNA pattern)

🔍 Damage Pattern Recognition
============================

Authentic Ancient DNA Patterns
------------------------------

**Characteristic Features**:

1. **High 5' C→T rates** (>15% in first few positions)
2. **High 3' G→A rates** (>10% in last few positions)
3. **Exponential decay** from sequence ends toward middle
4. **Strand asymmetry** (different patterns on forward/reverse strands)
5. **Statistical significance** (p < 0.05)

**Visual Indicators** in damage plots:
- Sharp peaks at sequence ends
- Gradual decline toward sequence middle
- Clear asymmetry between 5' and 3' ends

Modern DNA Patterns
-------------------

**Characteristic Features**:

1. **Low transition rates** (<5% across all positions)
2. **Uniform distribution** (no position-dependent effects)
3. **Random error patterns** (not systematically at ends)
4. **No strand asymmetry**
5. **No statistical significance** (p > 0.05)

Contamination Patterns
---------------------

**Mixed Ancient/Modern**:
- Intermediate damage scores (0.2-0.4)
- Irregular position-dependent patterns
- Variable significance levels

**Modern Contamination**:
- Lower damage scores than expected
- Reduced statistical significance
- Inconsistent patterns across samples

📈 Interpreting Damage Analysis Results
======================================

JSON Output Structure
---------------------

The damage analysis produces detailed JSON output:

.. code-block:: json

   {
     "sample_id": "sample001",
     "total_sequences": 245,
     "total_bases": 62847,
     "damage_score": 0.34,
     "p_value": 0.0234,
     "assessment": "Moderate damage detected",
     "significance": "statistically_significant",
     "c_to_t_rate": 0.187,
     "g_to_a_rate": 0.156,
     "background_rate": 0.034,
     "position_data": {
       "5_prime": [0.21, 0.18, 0.15, 0.12, 0.09, ...],
       "3_prime": [0.17, 0.14, 0.11, 0.08, 0.06, ...]
     },
     "bootstrap_stats": {
       "mean": 0.032,
       "std": 0.018,
       "confidence_95": [0.028, 0.036]
     }
   }

**Key Fields Explained**:

- ``damage_score``: Overall damage assessment (0-1 scale)
- ``p_value``: Statistical significance of damage pattern
- ``assessment``: Human-readable damage interpretation
- ``c_to_t_rate``: C→T transition rate at 5' ends
- ``g_to_a_rate``: G→A transition rate at 3' ends
- ``position_data``: Damage rates by sequence position
- ``bootstrap_stats``: Statistical validation results

Result Classification
---------------------

**High Confidence Ancient DNA**:
.. code-block:: json

   {
     "damage_score": 0.67,
     "p_value": 0.0001,
     "assessment": "High damage detected - likely ancient DNA"
   }

**Moderate Confidence**:
.. code-block:: json

   {
     "damage_score": 0.28,
     "p_value": 0.0423,
     "assessment": "Moderate damage detected"
   }

**Low Confidence/Modern**:
.. code-block:: json

   {
     "damage_score": 0.12,
     "p_value": 0.2567,
     "assessment": "Low damage detected - consistent with modern DNA"
   }

🎨 Visualization and Plots
=========================

Damage Profile Plots
--------------------

The QC report includes interactive damage profile plots showing:

**5' End Damage** (C→T transitions):
- X-axis: Position from 5' end (1-20)
- Y-axis: C→T transition frequency (%)
- Expected: High values at position 1, exponential decay

**3' End Damage** (G→A transitions):
- X-axis: Position from 3' end (-20 to -1)
- Y-axis: G→A transition frequency (%)
- Expected: High values at position -1, exponential decay

**Combined Damage Plot**:
- Shows both 5' and 3' patterns together
- Highlights asymmetric damage patterns
- Includes confidence intervals from bootstrap analysis

Statistical Validation Plots
----------------------------

**Bootstrap Distribution**:
- Histogram of bootstrap damage scores
- Observed score marked as vertical line
- P-value calculation visualization

**Confidence Intervals**:
- 95% confidence intervals for damage estimates
- Error bars on position-specific damage rates
- Statistical significance indicators

⚙️ Configuration Parameters
===========================

Damage Analysis Settings
------------------------

Key configuration parameters that affect damage analysis:

.. code-block:: yaml

   # Damage analysis configuration
   damage_threshold: 0.05        # P-value significance threshold
   bootstrap_iterations: 10000   # Number of bootstrap samples
   min_sequence_length: 50       # Minimum length for damage analysis
   position_range: 20            # Positions to analyze from each end
   
   # Advanced settings
   background_region: [20, -20]  # Positions for background calculation
   transition_types: ["C>T", "G>A"]  # Transition types to analyze
   confidence_level: 0.95        # Confidence interval level

**Parameter Effects**:

- **Higher bootstrap_iterations**: More precise p-values (slower)
- **Lower damage_threshold**: More conservative significance testing
- **Larger position_range**: Analyzes more positions from ends
- **Higher min_sequence_length**: Excludes short, potentially unreliable sequences

🧪 Quality Control and Validation
=================================

Internal Quality Checks
-----------------------

The pipeline performs several quality control checks:

1. **Sequence Length Validation**: Ensures sequences are long enough for reliable analysis
2. **Coverage Assessment**: Checks that sufficient positions have adequate coverage
3. **Bootstrap Convergence**: Verifies bootstrap analysis has converged
4. **Outlier Detection**: Identifies and flags unusual damage patterns

**Quality Flags**:
- ``insufficient_coverage``: Too few sequences or positions
- ``short_sequences``: Average sequence length below threshold
- ``bootstrap_warning``: Bootstrap analysis may be unreliable
- ``outlier_pattern``: Unusual damage distribution

External Validation
-------------------

**Cross-Sample Consistency**:
- Compare damage patterns across samples from same context
- Look for consistent archaeological/temporal patterns
- Check for batch effects or processing artifacts

**Positive/Negative Controls**:
- Include known ancient samples (positive controls)
- Include modern DNA samples (negative controls)
- Compare results with established methods

🔬 Advanced Interpretation
=========================

Age Estimation
--------------

While damage patterns can indicate ancient DNA, they cannot precisely determine age:

**General Guidelines**:
- **Very high damage** (>0.6): Likely >1000 years old
- **High damage** (0.4-0.6): Potentially 500-1000 years
- **Moderate damage** (0.2-0.4): Recent or well-preserved
- **Low damage** (<0.2): Modern or extremely well-preserved

.. warning::
   Age estimation from damage is approximate and depends on preservation conditions, temperature, pH, and other environmental factors.

Preservation Assessment
----------------------

Damage patterns can indicate preservation quality:

**Excellent Preservation**:
- Low damage despite age
- Even coverage across regions
- High sequence quality

**Poor Preservation**:
- High damage relative to age
- Uneven damage patterns
- Low sequence quality

**Variable Preservation**:
- Inconsistent damage across samples
- Position-dependent quality variations
- May indicate heterogeneous conditions

🎯 Best Practices
================

Sample Processing
-----------------

1. **Include Controls**: Always process positive and negative controls
2. **Replicate Extractions**: Process multiple extractions when possible
3. **Document Context**: Record archaeological/environmental context
4. **Blind Analysis**: Analyze samples without knowing expected results

Data Interpretation
------------------

1. **Consider Context**: Interpret results in light of archaeological context
2. **Multiple Lines of Evidence**: Use damage analysis alongside other authenticity criteria
3. **Conservative Approach**: Be cautious with borderline results
4. **Expert Review**: Have results reviewed by experienced researchers

Reporting Standards
------------------

1. **Full Disclosure**: Report all samples, including failures
2. **Method Details**: Describe analysis parameters and settings
3. **Statistical Results**: Include p-values and confidence intervals
4. **Visual Evidence**: Provide damage profile plots
5. **Raw Data**: Make underlying data available for review

🚨 Common Pitfalls
=================

Interpretation Errors
---------------------

**Over-interpretation**:
- Calling modern DNA "ancient" based on borderline damage
- Ignoring statistical significance
- Not considering preservation context

**Under-interpretation**:
- Dismissing significant damage patterns
- Requiring unrealistically high damage levels
- Ignoring consistent patterns across samples

Technical Issues
---------------

**Sample Preparation**:
- Contamination with modern DNA
- PCR artifacts mimicking damage
- Library preparation effects

**Analysis Parameters**:
- Inappropriate quality thresholds
- Insufficient bootstrap iterations
- Wrong reference sequences

📚 Further Reading
=================

**Key Publications**:

1. Briggs et al. (2007) - Patterns of damage in genomic DNA sequences from a Neandertal
2. Skoglund et al. (2014) - Separating ancient DNA from modern contamination
3. Jónsson et al. (2013) - mapDamage2.0: fast approximate Bayesian estimates

**Related Methods**:
- mapDamage: Alternative damage analysis software
- PMDtools: Authentication based on damage patterns  
- EAGER: Ancient DNA analysis pipeline

This comprehensive guide provides the theoretical background and practical knowledge needed to understand and interpret ancient DNA damage analysis results from the Sanger DNA Damage Analysis Pipeline.
