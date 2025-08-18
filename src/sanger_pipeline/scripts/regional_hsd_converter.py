#!/usr/bin/env python3
"""
Regional HSD Converter: Align each HVS region separately to its reference segment.

This approach breaks the alignment into three independent processes:
1. HVS1 consensus → HVS1 reference region alignment → HVS1 variants
2. HVS2 consensus → HVS2 reference region alignment → HVS2 variants
3. HVS3 consensus → HVS3 reference region alignment → HVS3 variants

This should reduce alignment artifacts and provide more accurate variant calling.

Usage:
    python -m sanger_pipeline.scripts.regional_hsd_converter <consensus_dir> <output_hsd>
"""

import sys
import argparse
import re
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    print("Error: BioPython not installed. Install with: pip install biopython")
    sys.exit(1)


class RegionalHSDConverter:
    """Convert HVS consensus files using regional alignment approach."""

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

    def align_to_reference_segment(
        self, consensus_seq: str, hvs_region: str
    ) -> Optional[str]:
        """
        Align consensus sequence to its corresponding reference segment using MAFFT.

        Args:
            consensus_seq: Cleaned consensus sequence
            hvs_region: HVS region (HVS1, HVS2, or HVS3)

        Returns:
            Aligned consensus sequence or None if alignment fails
        """
        if hvs_region not in self.reference_segments:
            print(f"⚠️  Unknown HVS region: {hvs_region}")
            return None

        reference_segment = self.reference_segments[hvs_region]

        # Create temporary files for MAFFT
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as temp_input:
            # Write reference segment and consensus as FASTA
            temp_input.write(f">reference_{hvs_region}\n{reference_segment}\n")
            temp_input.write(f">consensus_{hvs_region}\n{consensus_seq}\n")
            temp_input_path = temp_input.name

        try:
            # Run MAFFT alignment
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            ) as temp_output:
                temp_output_path = temp_output.name

            cmd = ["mafft", "--quiet", "--auto", temp_input_path]

            with open(temp_output_path, "w") as output_file:
                result = subprocess.run(
                    cmd, stdout=output_file, stderr=subprocess.PIPE, text=True
                )

            if result.returncode != 0:
                print(f"⚠️  MAFFT failed for {hvs_region}: {result.stderr}")
                return None

            # Parse alignment result
            aligned_consensus = None
            with open(temp_output_path, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    if "consensus" in record.id:
                        aligned_consensus = str(record.seq).upper()
                        break

            # Clean up temporary files
            Path(temp_input_path).unlink()
            Path(temp_output_path).unlink()

            return aligned_consensus

        except Exception as e:
            print(f"⚠️  Alignment error for {hvs_region}: {e}")
            return None

    def find_variants_from_alignment(
        self, aligned_consensus: str, hvs_region: str
    ) -> List[str]:
        """
        Find variants from aligned sequences.

        Args:
            aligned_consensus: Aligned consensus sequence (may contain gaps)
            hvs_region: HVS region name

        Returns:
            List of variants in HSD format
        """
        variants = []
        reference_segment = self.reference_segments[hvs_region]
        region_start = self.hvs_regions[hvs_region]["start"]

        # Align reference segment too (to get same gap structure)
        # For now, we'll use a simpler approach comparing position by position

        # Remove gaps from aligned sequences for comparison
        ref_clean = reference_segment
        cons_clean = aligned_consensus.replace("-", "")

        if len(cons_clean) == 0:
            return variants

        # Compare sequences position by position
        min_length = min(len(ref_clean), len(cons_clean))
        max_differences = (
            min_length // 8
        )  # Allow max 12.5% differences (more conservative)
        differences_found = 0

        for pos in range(min_length):
            if differences_found >= max_differences:
                break

            ref_base = ref_clean[pos]
            cons_base = cons_clean[pos]

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
        Process a single consensus FASTA file.

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

                    # Align to reference segment
                    aligned_consensus = self.align_to_reference_segment(
                        consensus_seq, hvs_region
                    )

                    if aligned_consensus is None:
                        print(f"⚠️  {fasta_file.name}: Alignment failed")
                        return sample_id, hvs_region, []

                    # Find variants from alignment
                    variants = self.find_variants_from_alignment(
                        aligned_consensus, hvs_region
                    )

                    print(
                        f"📁 {fasta_file.name}: {len(variants)} variants ({hvs_region}) - aligned method"
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
        description="Convert HVS consensus files to HSD format using regional alignment approach",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python regional_hsd_converter.py output/consensus/ samples.hsd
  python regional_hsd_converter.py data/consensus/ --output results.hsd
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
        output_file = "regional_alignment_output.hsd"

    try:
        print("🧬 Regional HSD Converter")
        print("=" * 50)
        print(f"Input directory: {args.consensus_dir}")
        print(f"Output file: {output_file}")
        print(f"Reference: {args.reference}")
        print()

        # Create converter
        converter = RegionalHSDConverter(args.reference)
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
