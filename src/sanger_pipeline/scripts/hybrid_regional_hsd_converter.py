#!/usr/bin/env python3
"""
Hybrid HSD Converter: Regional direct comparison without alignment artifacts.

This approach:
1. Extracts reference segments for each HVS region
2. Performs DIRECT comparison (no MAFFT alignment)
3. Processes each region independently
4. Combines results per sample

This avoids the alignment artifacts while maintaining regional specificity.

Usage:
    python -m sanger_pipeline.scripts.hybrid_regional_hsd_converter <consensus_dir> <output_hsd>
"""

import sys
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    print("Error: BioPython not installed. Install with: pip install biopython")
    sys.exit(1)


class HybridRegionalHSDConverter:
    """Convert HVS consensus files using hybrid regional direct comparison."""

    def __init__(self, reference_file: str = "ref/rCRS.fasta"):
        """Initialize converter with reference file."""
        self.reference_file = reference_file
        self.reference_seq = None

        # HVS region definitions (1-based coordinates)
        self.hvs_regions = {
            "HVS1": {"start": 16024, "end": 16365},
            "HVS2": {"start": 57, "end": 372},
            "HVS3": {"start": 438, "end": 574},
        }

        # Load reference sequence
        self.load_reference()

        # Extract reference segments for each HVS region
        self.reference_segments = self._extract_reference_segments()

    def load_reference(self):
        """Load reference sequence from FASTA file."""
        try:
            with open(self.reference_file, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    self.reference_seq = str(record.seq)
                    break
            if self.reference_seq is None:
                raise Exception("No sequences found in reference file")
            print(f"✓ Reference sequence loaded: {len(self.reference_seq)} bases")
        except Exception as e:
            raise Exception(f"Error loading reference: {e}")

    def _extract_reference_segments(self) -> Dict[str, str]:
        """Extract reference segments for each HVS region."""
        if self.reference_seq is None:
            raise Exception("Reference sequence not loaded")

        segments = {}
        for region, coords in self.hvs_regions.items():
            # Convert to 0-based indexing
            start_idx = coords["start"] - 1
            end_idx = coords["end"]
            segment = self.reference_seq[start_idx:end_idx]
            segments[region] = self._clean_sequence(segment)
            print(
                f"✓ {region}: {len(segment)} bases ({coords['start']}-{coords['end']})"
            )
        return segments

    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing non-ATCG characters."""
        return re.sub(r"[^ATCG]", "", sequence.upper())

    def parse_sample_name(self, filename: str) -> Tuple[str, Optional[str]]:
        """Parse sample name and HVS region from filename."""
        basename = filename.replace(".fasta", "").replace("_consensus", "")

        # Extract HVS region
        hvs_region = None
        if "HVS1" in basename.upper():
            hvs_region = "HVS1"
        elif "HVS2" in basename.upper():
            hvs_region = "HVS2"
        elif "HVS3" in basename.upper():
            hvs_region = "HVS3"

        # Extract sample ID
        if hvs_region:
            pattern = f"_{hvs_region.lower()}"
            sample_id = basename.lower().replace(pattern, "")
        else:
            sample_id = basename.lower()

        return sample_id, hvs_region

    def find_variants_direct_regional(
        self, consensus_seq: str, hvs_region: str
    ) -> List[str]:
        """
        Find variants using direct comparison within the specific HVS region.

        Args:
            consensus_seq: Cleaned consensus sequence
            hvs_region: HVS region (HVS1, HVS2, or HVS3)

        Returns:
            List of variants in HSD format
        """
        variants = []

        if hvs_region not in self.reference_segments:
            print(f"⚠️  Unknown HVS region: {hvs_region}")
            return variants

        reference_segment = self.reference_segments[hvs_region]
        region_start = self.hvs_regions[hvs_region]["start"]

        # Direct comparison - no alignment
        min_length = min(len(consensus_seq), len(reference_segment))

        # Skip very short sequences
        if min_length < 50:
            return variants

        # Conservative approach: limit differences for ancient DNA
        max_differences = min_length // 10  # Allow max 10% differences
        differences_found = 0

        for pos in range(min_length):
            if differences_found >= max_differences:
                break

            ref_base = reference_segment[pos]
            cons_base = consensus_seq[pos]

            if ref_base != cons_base and cons_base in "ATCG":
                # Convert to mitochondrial coordinates
                mt_position = region_start + pos
                variant = f"{mt_position}{cons_base}"
                variants.append(variant)
                differences_found += 1

        return variants

    def process_consensus_file(
        self, fasta_file: Path
    ) -> Tuple[str, Optional[str], List[str]]:
        """
        Process a single consensus FASTA file using regional direct comparison.

        Returns:
            (sample_id, hvs_region, variants)
        """
        sample_id, hvs_region = self.parse_sample_name(fasta_file.name)

        if not hvs_region:
            print(f"⚠️  Could not determine HVS region for {fasta_file.name}")
            return sample_id, "unknown", []

        try:
            # Read consensus sequence
            with open(fasta_file, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    consensus_seq = self._clean_sequence(str(record.seq))

                    if len(consensus_seq) < 50:
                        print(
                            f"⚠️  {fasta_file.name}: Sequence too short ({len(consensus_seq)}bp)"
                        )
                        return sample_id, hvs_region, []

                    # Find variants using direct regional comparison
                    variants = self.find_variants_direct_regional(
                        consensus_seq, hvs_region
                    )

                    print(
                        f"📁 {fasta_file.name}: {len(variants)} variants ({hvs_region}) - hybrid regional"
                    )
                    return sample_id, hvs_region, variants

        except Exception as e:
            print(f"❌ Error processing {fasta_file.name}: {e}")
            return sample_id, hvs_region or "unknown", []

        return sample_id, hvs_region or "unknown", []

    def process_consensus_directory(self, consensus_dir: str) -> Dict[str, List[str]]:
        """
        Process all consensus files in a directory.

        Returns:
            Dictionary mapping sample_id to list of all variants
        """
        consensus_path = Path(consensus_dir)
        if not consensus_path.exists():
            raise FileNotFoundError(f"Directory not found: {consensus_dir}")

        fasta_files = list(consensus_path.glob("*_consensus.fasta"))
        if not fasta_files:
            print(f"⚠️  No consensus files found in {consensus_dir}")
            return {}

        print(f"🔍 Found {len(fasta_files)} consensus files")

        # Group results by sample
        sample_variants = defaultdict(list)

        for fasta_file in sorted(fasta_files):
            sample_id, hvs_region, variants = self.process_consensus_file(fasta_file)

            if variants:
                sample_variants[sample_id].extend(variants)

        # Sort variants for each sample
        for sample_id in sample_variants:
            sample_variants[sample_id] = sorted(set(sample_variants[sample_id]))

        return dict(sample_variants)

    def write_hsd_file(self, sample_variants: Dict[str, List[str]], output_file: str):
        """Write variants to HSD format file."""
        with open(output_file, "w") as f:
            # Write header
            f.write("SampleID\tRange\n")

            # Write each sample
            for sample_id, variants in sample_variants.items():
                if variants:
                    variants_str = " ".join(variants)
                    f.write(f"{sample_id}\t{variants_str}\n")
                else:
                    f.write(f"{sample_id}\trCRS\n")

        print(f"✅ HSD file written: {output_file}")
        print(f"📊 Processed {len(sample_variants)} samples")


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description="Convert HVS consensus files to HSD format using hybrid regional direct comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python hybrid_regional_hsd_converter.py output/consensus/ samples.hsd
  python hybrid_regional_hsd_converter.py data/consensus/ --output results.hsd
        """,
    )

    parser.add_argument(
        "consensus_dir", help="Directory containing consensus FASTA files"
    )
    parser.add_argument("output_file", nargs="?", help="Output HSD file (optional)")
    parser.add_argument(
        "--output", "-o", help="Output HSD file (alternative way to specify)"
    )
    parser.add_argument(
        "--reference",
        "-r",
        default="ref/rCRS.fasta",
        help="Reference FASTA file (default: ref/rCRS.fasta)",
    )

    args = parser.parse_args()

    # Determine output file
    if args.output_file:
        output_file = args.output_file
    elif args.output:
        output_file = args.output
    else:
        output_file = "hybrid_regional_output.hsd"

    try:
        print("🧬 Hybrid Regional HSD Converter")
        print("=" * 50)
        print(f"Input directory: {args.consensus_dir}")
        print(f"Output file: {output_file}")
        print(f"Reference: {args.reference}")
        print()

        # Create converter
        converter = HybridRegionalHSDConverter(args.reference)
        print()

        # Process consensus files
        print("🔄 Processing consensus files...")
        sample_variants = converter.process_consensus_directory(args.consensus_dir)
        print()

        # Write results
        print("💾 Writing HSD file...")
        converter.write_hsd_file(sample_variants, output_file)
        print()

        # Summary
        total_samples = len(sample_variants)
        samples_with_variants = sum(
            1 for variants in sample_variants.values() if variants
        )
        total_variants = sum(len(variants) for variants in sample_variants.values())

        print("📈 Summary:")
        print(f"   Total samples: {total_samples}")
        print(f"   Samples with variants: {samples_with_variants}")
        print(f"   Total variants found: {total_variants}")
        if total_samples > 0:
            avg_variants = total_variants / total_samples
            print(f"   Average variants per sample: {avg_variants:.1f}")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
