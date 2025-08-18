"""
FASTA to HSD Format Converter

Converts concatenated mitochondrial DNA FASTA sequences to HSD format
for haplogroup analysis with tools like HaploGrep.
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq


class FastaToHSDConverter:
    """Converts FASTA sequences to HSD format for haplogroup analysis."""
    
    def __init__(self, reference_sequence: Optional[str] = None):
        """
        Initialize converter with reference sequence.
        
        Args:
            reference_sequence: Path to reference sequence file or sequence string.
                               If None, looks for reference in ref/ directory.
        """
        self.reference_seq = self._load_reference_sequence(reference_sequence)
        print(f"Reference sequence loaded: {len(self.reference_seq)} bases")
        
        # Define HVS regions (1-based coordinates)
        self.hvs_regions = {
            'HVS1': {'start': 16024, 'end': 16365},
            'HVS2': {'start': 57, 'end': 372}, 
            'HVS3': {'start': 438, 'end': 574}
        }
        
    def _load_reference_sequence(self, reference: Optional[str]) -> str:
        """Load reference sequence from file or use default location."""
        # Try provided reference first
        if reference and Path(reference).exists():
            print(f"Loading reference from: {reference}")
            record = SeqIO.read(reference, "fasta")
            return str(record.seq).upper()
        
        # Look for reference in standard locations
        possible_refs = [
            Path("ref/rCRS.fasta"),
            Path("ref/reference.fasta"),
            Path("ref/mtDNA_reference.fasta"),
            Path("reference/rCRS.fasta"),
            Path("reference.fasta")
        ]
        
        for ref_path in possible_refs:
            if ref_path.exists():
                print(f"Loading reference from: {ref_path}")
                record = SeqIO.read(ref_path, "fasta")
                return str(record.seq).upper()
        
        # If no reference found, raise error
        raise FileNotFoundError(
            "Reference sequence not found. Please provide a reference file or "
            "place it in ref/rCRS.fasta"
        )
    
    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing gaps and converting to uppercase."""
        return re.sub(r'[^ATCG]', '', sequence.upper())
    
    def _detect_hvs_regions(self, sample_name: str) -> List[str]:
        """Detect which HVS regions are present based on sample name."""
        hvs_present = []
        sample_upper = sample_name.upper()
        
        if 'HVS1' in sample_upper:
            hvs_present.append('HVS1')
        if 'HVS2' in sample_upper:
            hvs_present.append('HVS2')
        if 'HVS3' in sample_upper:
            hvs_present.append('HVS3')
            
        # If no HVS regions detected, try to infer from sequence length
        if not hvs_present:
            # This is a fallback - you may need to adjust based on your data
            hvs_present = ['HVS1', 'HVS2']  # Default assumption
            
        return hvs_present
    
    def _extract_hvs_references(self, hvs_regions: List[str]) -> str:
        """Extract and concatenate HVS reference sequences."""
        concatenated_ref = ""
        
        for region in hvs_regions:
            if region in self.hvs_regions:
                start = self.hvs_regions[region]['start'] - 1  # Convert to 0-based
                end = self.hvs_regions[region]['end']
                ref_segment = self.reference_seq[start:end]
                concatenated_ref += ref_segment
                print(f"  Added {region}: positions {start+1}-{end} ({len(ref_segment)} bp)")
        
        return concatenated_ref
    
    def find_variants(self, sample_seq: str, sample_name: str) -> List[str]:
        """
        Find variants between sample and reference sequence.
        
        Args:
            sample_seq: Sample DNA sequence
            sample_name: Sample identifier
            
        Returns:
            List of HSD format variants
        """
        variants = []
        
        # Clean sample sequence
        sample_clean = self._clean_sequence(sample_seq)
        
        print(f"Processing {sample_name}: Sample={len(sample_clean)}bp")
        
        # Detect HVS regions from sample name
        hvs_regions = self._detect_hvs_regions(sample_name)
        print(f"  Detected HVS regions: {', '.join(hvs_regions)}")
        
        if not hvs_regions:
            print(f"Warning: No HVS regions detected for {sample_name}")
            return variants
        
        # Check if we need to process individual HVS files or a concatenated file
        if len(hvs_regions) == 1:
            # Single HVS region - direct mapping
            region = hvs_regions[0]
            ref_seq = self._extract_hvs_references([region])
            variants = self._find_variants_single_region(sample_clean, region, ref_seq, sample_name)
        else:
            # Multiple HVS regions - try to detect boundaries in the concatenated sequence
            variants = self._find_variants_concatenated(sample_clean, hvs_regions, sample_name)
        
        print(f"  Found {len(variants)} variants for {sample_name}")
        return variants
    
    def _find_variants_single_region(self, sample_seq: str, region: str, ref_seq: str, sample_name: str) -> List[str]:
        """Find variants for a single HVS region."""
        variants = []
        
        if region not in self.hvs_regions:
            return variants
            
        region_start = self.hvs_regions[region]['start']
        min_length = min(len(sample_seq), len(ref_seq))
        
        if min_length < 50:  # Skip very short sequences
            print(f"    Warning: Sequence too short ({min_length}bp) for {sample_name}")
            return variants
        
        # Find variants with quality filtering
        variant_count = 0
        for pos in range(min_length):
            if pos >= len(sample_seq) or pos >= len(ref_seq):
                break
                
            ref_base = ref_seq[pos]
            sample_base = sample_seq[pos]
            
            if ref_base != sample_base:
                mt_position = region_start + pos
                variant = f"{mt_position}{sample_base}"
                variants.append(variant)
                variant_count += 1
                
                # Debug: show first few variants
                if variant_count <= 5:
                    print(f"    Variant {variant_count}: pos {mt_position} {ref_base}→{sample_base}")
        
        return variants
    
    def _find_variants_concatenated(self, sample_seq: str, hvs_regions: List[str], sample_name: str) -> List[str]:
        """Find variants for concatenated HVS regions with boundary detection."""
        variants = []
        
        # Extract reference sequences for each region
        ref_sequences = []
        region_starts = []
        
        for region in hvs_regions:
            if region in self.hvs_regions:
                region_ref = self.reference_seq[
                    self.hvs_regions[region]['start']-1:self.hvs_regions[region]['end']
                ]
                ref_clean = self._clean_sequence(region_ref)
                ref_sequences.append(ref_clean)
                region_starts.append(self.hvs_regions[region]['start'])
        
        if not ref_sequences:
            return variants
        
        # Try to identify region boundaries in the sample sequence
        # This is a heuristic approach - in practice, you might need more sophisticated alignment
        total_ref_length = sum(len(ref) for ref in ref_sequences)
        
        if len(sample_seq) < total_ref_length * 0.5:  # Sample is too short
            print(f"    Warning: Sample sequence ({len(sample_seq)}bp) much shorter than expected ({total_ref_length}bp)")
            return variants
        
        # Simple approach: assume regions are concatenated in order
        sample_pos = 0
        
        for i, (region, ref_seq, mt_start) in enumerate(zip(hvs_regions, ref_sequences, region_starts)):
            region_length = len(ref_seq)
            
            # Extract the corresponding part of the sample sequence
            sample_region = sample_seq[sample_pos:sample_pos + region_length]
            
            if len(sample_region) < region_length * 0.7:  # Region too short
                print(f"    Warning: {region} region too short in sample")
                sample_pos += region_length
                continue
            
            # Compare this region
            comparison_length = min(len(sample_region), region_length)
            region_variants = 0
            
            for pos in range(comparison_length):
                if pos >= len(sample_region) or pos >= len(ref_seq):
                    break
                    
                ref_base = ref_seq[pos]
                sample_base = sample_region[pos]
                
                if ref_base != sample_base:
                    mt_position = mt_start + pos
                    variant = f"{mt_position}{sample_base}"
                    variants.append(variant)
                    region_variants += 1
                    
                    # Debug: show first few variants per region
                    if region_variants <= 3:
                        print(f"    {region} Variant {region_variants}: pos {mt_position} {ref_base}→{sample_base}")
            
            print(f"    {region}: {region_variants} variants")
            sample_pos += region_length
        
        return variants
    
    def convert_fasta_to_hsd(self, fasta_file: Path, output_file: Path, 
                            range_notation: str = "16024-16365;57-372;438-574", 
                            haplogroup_column: str = "?") -> None:
        """
        Convert FASTA file to HSD format.
        
        Args:
            fasta_file: Input FASTA file path
            output_file: Output HSD file path
            range_notation: Range notation for sequences (e.g., "16024-16365;57-372;438-574")
            haplogroup_column: Default haplogroup assignment
        """
        hsd_lines = []
        
        print(f"Converting {fasta_file} to HSD format...")
        
        # Parse FASTA file
        sequence_count = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_count += 1
            sample_name = record.id
            sequence = str(record.seq)
            
            print(f"\nProcessing sequence {sequence_count}: {sample_name}")
            print(f"  Raw sequence length: {len(sequence)} bases")
            
            # Find variants
            variants = self.find_variants(sequence, sample_name)
            
            # Create HSD line
            if variants:
                variants_str = "\t".join(variants)
                hsd_line = f"{sample_name}\t{range_notation}\t{haplogroup_column}\t{variants_str}"
            else:
                # No variants found - matches reference perfectly
                hsd_line = f"{sample_name}\t{range_notation}\t{haplogroup_column}"
            
            hsd_lines.append(hsd_line)
            print(f"  HSD entry: {len(variants)} variants")
        
        # Write HSD file
        with open(output_file, 'w') as f:
            for line in hsd_lines:
                f.write(line + "\n")
        
        print(f"\n✅ Converted {len(hsd_lines)} sequences to HSD format: {output_file}")
    
    def convert_pipeline_output(self, pipeline_output_dir: Path, 
                               output_file: Path) -> None:
        """
        Convert Sanger pipeline merged FASTA outputs to HSD format.
        
        Args:
            pipeline_output_dir: Directory containing pipeline outputs
            output_file: Output HSD file path
        """
        # Find all merged FASTA files from pipeline
        final_dir = pipeline_output_dir / "final"
        
        if not final_dir.exists():
            print(f"Final directory not found: {final_dir}")
            print("Looking for alternative directories...")
            
            # Try other possible locations
            alternatives = [
                pipeline_output_dir / "consensus",
                pipeline_output_dir / "merged", 
                pipeline_output_dir
            ]
            
            for alt_dir in alternatives:
                if alt_dir.exists():
                    fasta_files = list(alt_dir.glob("*.fasta"))
                    if fasta_files:
                        final_dir = alt_dir
                        print(f"Found FASTA files in: {final_dir}")
                        break
        
        # Look for FASTA files
        fasta_patterns = ["*_merged.fasta", "*_final.fasta", "*_consensus.fasta", "*.fasta"]
        merged_files = []
        
        for pattern in fasta_patterns:
            files = list(final_dir.glob(pattern))
            merged_files.extend(files)
            if files:
                print(f"Found {len(files)} files matching {pattern}")
        
        # Remove duplicates
        merged_files = list(set(merged_files))
        
        if not merged_files:
            print(f"No FASTA files found in {final_dir}")
            print("Available files:")
            for f in final_dir.iterdir():
                print(f"  {f.name}")
            return
        
        print(f"\nProcessing {len(merged_files)} FASTA files...")
        hsd_lines = []
        
        for fasta_file in merged_files:
            print(f"\n📁 Processing: {fasta_file.name}")
            
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sample_name = record.id
                    sequence = str(record.seq)
                    
                    # Extract sample info from filename if needed
                    if sample_name == "consensus" or not sample_name:
                        # Use filename for sample name
                        sample_name = fasta_file.stem.replace("_merged", "").replace("_final", "").replace("_consensus", "")
                    
                    print(f"  Sample: {sample_name}")
                    print(f"  Length: {len(sequence)} bases")
                    
                    # Find variants
                    variants = self.find_variants(sequence, sample_name)
                    
                    # Create HSD line with regional ranges
                    if variants:
                        variants_str = " ".join(variants)  # Space-separated
                        hsd_line = f"{sample_name}\t16024-16365;57-372;438-574\t?\t{variants_str}"
                    else:
                        hsd_line = f"{sample_name}\t16024-16365;57-372;438-574\t?"
                    
                    hsd_lines.append(hsd_line)
                    
            except Exception as e:
                print(f"  ❌ Error processing {fasta_file.name}: {e}")
                continue
        
        # Write HSD file
        with open(output_file, 'w') as f:
            for line in hsd_lines:
                f.write(line + "\n")
        
        print(f"\n✅ Converted {len(hsd_lines)} sequences to HSD format: {output_file}")
        print(f"📄 HSD file ready for HaploGrep analysis!")


def main():
    """Command line interface for FASTA to HSD conversion."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Convert FASTA sequences to HSD format for haplogroup analysis"
    )
    parser.add_argument(
        "--input", "-i", 
        required=True,
        help="Input FASTA file or pipeline output directory"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output HSD file"
    )
    parser.add_argument(
        "--reference", "-r",
        help="Reference sequence file (auto-detects if not provided)"
    )
    parser.add_argument(
        "--pipeline-mode",
        action="store_true",
        help="Process entire pipeline output directory"
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize converter
        converter = FastaToHSDConverter(args.reference)
        
        input_path = Path(args.input)
        output_path = Path(args.output)
        
        if args.pipeline_mode:
            converter.convert_pipeline_output(input_path, output_path)
        else:
            converter.convert_fasta_to_hsd(input_path, output_path)
            
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())