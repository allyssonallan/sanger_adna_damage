#!/usr/bin/env python3
"""
Whole Reference HSD Converter: Align consensus sequences to complete mitochondrial reference.

This approach:
1. Aligns each consensus sequence to the COMPLETE reference genome
2. Automatically determines the best-matching HVS region
3. Extracts variants based on the aligned position
4. Provides better positioning accuracy for consensus sequences

This is especially useful when:
- HVS region identification from filename is uncertain
- Consensus sequences may span multiple regions
- You want to find the best matching position automatically
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
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: BioPython not installed. Install with: pip install biopython")
    sys.exit(1)


class WholeReferenceHSDConverter:
    """Convert HVS consensus files using whole reference alignment approach."""

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

    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing non-ATCG characters."""
        return re.sub(r"[^ATCG]", "", sequence.upper())

    def parse_sample_name(self, filename: str) -> Tuple[str, Optional[str]]:
        """Parse sample name and try to identify HVS region from filename."""
        basename = filename.replace(".fasta", "").replace("_consensus", "")

        # Extract HVS region hint from filename
        hvs_region_hint = None
        if "HVS1" in basename.upper():
            hvs_region_hint = "HVS1"
        elif "HVS2" in basename.upper():
            hvs_region_hint = "HVS2"
        elif "HVS3" in basename.upper():
            hvs_region_hint = "HVS3"

        # Extract sample ID
        if hvs_region_hint:
            pattern = f"_{hvs_region_hint.lower()}"
            sample_id = basename.lower().replace(pattern, "")
        else:
            sample_id = basename.lower()

        return sample_id, hvs_region_hint

    def align_to_whole_reference(
        self, consensus_seq: str, sample_name: str
    ) -> Optional[Dict]:
        """
        Align consensus sequence to the complete reference genome.

        Args:
            consensus_seq: Cleaned consensus sequence
            sample_name: Sample identifier for logging

        Returns:
            Dictionary with alignment information or None if alignment fails
        """
        # Create temporary files for MAFFT
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as temp_input:
            # Write complete reference and consensus as FASTA
            temp_input.write(f">reference_complete\n{self.reference_seq}\n")
            temp_input.write(f">consensus_{sample_name}\n{consensus_seq}\n")
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
                print(f"⚠️  MAFFT failed for {sample_name}: {result.stderr}")
                return None

            # Parse alignment result
            alignment_info = self._parse_whole_reference_alignment(
                temp_output_path, sample_name
            )

            # Clean up temporary files
            Path(temp_input_path).unlink()
            Path(temp_output_path).unlink()

            return alignment_info

        except Exception as e:
            print(f"⚠️  Alignment error for {sample_name}: {e}")
            return None

    def _parse_whole_reference_alignment(
        self, alignment_file: str, sample_name: str
    ) -> Dict:
        """
        Parse the alignment result to find the best matching position.

        Returns:
            Dictionary with alignment position information
        """
        try:
            alignment = AlignIO.read(alignment_file, "fasta")

            # Get sequences (reference should be first, consensus second)
            ref_aligned = None
            cons_aligned = None

            for record in alignment:
                if "reference" in record.id:
                    ref_aligned = str(record.seq)
                elif "consensus" in record.id:
                    cons_aligned = str(record.seq)

            if ref_aligned is None or cons_aligned is None:
                print(f"⚠️  Could not parse alignment for {sample_name}")
                return None

            # Find the actual start position in the reference
            ref_start_pos = self._find_alignment_start_position(
                ref_aligned, cons_aligned
            )

            # Determine which HVS region(s) this sequence overlaps
            matching_regions = self._determine_hvs_regions(
                ref_start_pos, len(cons_aligned.replace("-", ""))
            )

            return {
                "ref_aligned": ref_aligned,
                "cons_aligned": cons_aligned,
                "ref_start_pos": ref_start_pos,
                "matching_regions": matching_regions,
                "alignment_length": len(cons_aligned.replace("-", "")),
            }

        except Exception as e:
            print(f"⚠️  Error parsing alignment for {sample_name}: {e}")
            return None

    def _find_alignment_start_position(
        self, ref_aligned: str, cons_aligned: str
    ) -> int:
        """
        Find the actual start position in the reference where the consensus aligns.

        Returns:
            1-based position in the reference genome
        """
        # Count bases before the first non-gap in consensus
        consensus_start_in_alignment = 0
        for i, base in enumerate(cons_aligned):
            if base != "-":
                consensus_start_in_alignment = i
                break

        # Count actual reference bases up to this position
        ref_bases_count = 0
        for i in range(consensus_start_in_alignment):
            if ref_aligned[i] != "-":
                ref_bases_count += 1

        # Convert to 1-based position
        return ref_bases_count + 1

    def _determine_hvs_regions(self, start_pos: int, sequence_length: int) -> List[str]:
        """
        Determine which HVS regions this sequence overlaps.

        Args:
            start_pos: 1-based start position in reference
            sequence_length: Length of the consensus sequence

        Returns:
            List of overlapping HVS region names
        """
        end_pos = start_pos + sequence_length - 1
        overlapping_regions = []

        for region_name, coords in self.hvs_regions.items():
            # Check for overlap
            if not (end_pos < coords["start"] or start_pos > coords["end"]):
                overlapping_regions.append(region_name)

        return overlapping_regions

    def extract_variants_from_whole_alignment(
        self, alignment_info: Dict, sample_name: str
    ) -> List[str]:
        """
        Extract variants from the whole reference alignment.

        Args:
            alignment_info: Alignment information dictionary
            sample_name: Sample identifier for logging

        Returns:
            List of variants in HSD format
        """
        variants = []
        ref_aligned = alignment_info["ref_aligned"]
        cons_aligned = alignment_info["cons_aligned"]
        ref_start_pos = alignment_info["ref_start_pos"]

        # Track current position in reference
        current_ref_pos = ref_start_pos
        max_variants = 50  # Conservative limit for ancient DNA

        for i in range(len(ref_aligned)):
            ref_base = ref_aligned[i]
            cons_base = cons_aligned[i]

            # Skip if we've found too many variants (likely alignment issue)
            if len(variants) >= max_variants:
                print(
                    f"⚠️  {sample_name}: Too many variants found ({len(variants)}), stopping"
                )
                break

            # Handle gaps and matches
            if ref_base == "-":
                # Insertion in consensus (skip)
                continue
            elif cons_base == "-":
                # Deletion in consensus (move reference position)
                current_ref_pos += 1
                continue
            elif ref_base == cons_base:
                # Match (move both positions)
                current_ref_pos += 1
                continue
            else:
                # Mismatch - this is a variant
                if (
                    cons_base in "ATCG" and current_ref_pos <= 16569
                ):  # Valid mtDNA position
                    variant = f"{current_ref_pos}{cons_base}"
                    variants.append(variant)
                current_ref_pos += 1

        print(
            f"📍 {sample_name}: Found {len(variants)} variants from whole reference alignment"
        )
        return variants

    def process_consensus_file(
        self, fasta_file: Path
    ) -> Tuple[str, List[str], List[str]]:
        """
        Process a single consensus FASTA file using whole reference alignment.

        Returns:
            (sample_id, detected_regions, variants)
        """
        sample_id, region_hint = self.parse_sample_name(fasta_file.name)

        try:
            # Read consensus sequence
            with open(fasta_file, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    consensus_seq = self._clean_sequence(str(record.seq))

                    if len(consensus_seq) < 50:
                        print(
                            f"⚠️  {fasta_file.name}: Sequence too short ({len(consensus_seq)}bp)"
                        )
                        return sample_id, [], []

                    # Align to whole reference
                    alignment_info = self.align_to_whole_reference(
                        consensus_seq, sample_id
                    )

                    if alignment_info is None:
                        print(f"⚠️  {fasta_file.name}: Alignment failed")
                        return sample_id, [], []

                    # Get detected regions
                    detected_regions = alignment_info["matching_regions"]

                    # Extract variants
                    variants = self.extract_variants_from_whole_alignment(
                        alignment_info, sample_id
                    )

                    print(
                        f"🎯 {fasta_file.name}: Detected regions {detected_regions}, {len(variants)} variants"
                    )

                    if region_hint and region_hint not in detected_regions:
                        print(
                            f"⚠️  {fasta_file.name}: Filename suggests {region_hint}, but aligned to {detected_regions}"
                        )

                    return sample_id, detected_regions, variants

        except Exception as e:
            print(f"❌ Error processing {fasta_file.name}: {e}")
            return sample_id, [], []

        return sample_id, [], []

    def convert_directory(self, consensus_dir: str, output_file: str) -> None:
        """
        Convert all consensus files in a directory to HSD format.

        Args:
            consensus_dir: Directory containing consensus FASTA files
            output_file: Output HSD file path
        """
        consensus_path = Path(consensus_dir)

        if not consensus_path.exists():
            print(f"❌ Error: Consensus directory not found: {consensus_dir}")
            return

        # Find all FASTA files
        fasta_files = list(consensus_path.glob("*.fasta")) + list(
            consensus_path.glob("*.fa")
        )

        if not fasta_files:
            print(f"❌ Error: No FASTA files found in {consensus_dir}")
            return

        print(f"🔍 Found {len(fasta_files)} consensus files")
        print(f"📊 Using whole reference alignment approach")
        print(f"🧬 Reference: {self.reference_file}")

        # Process each file and collect results by sample
        sample_data = defaultdict(lambda: {"regions": [], "variants": []})

        for fasta_file in sorted(fasta_files):
            sample_id, detected_regions, variants = self.process_consensus_file(
                fasta_file
            )

            if variants:
                sample_data[sample_id]["regions"].extend(detected_regions)
                sample_data[sample_id]["variants"].extend(variants)

        # Write HSD file
        self._write_hsd_file(sample_data, output_file)

    def _write_hsd_file(self, sample_data: Dict, output_file: str) -> None:
        """Write sample data to HSD format file."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"\n📝 Writing HSD file: {output_file}")

        with open(output_file, "w") as f:
            f.write("SampleID\\tRange\\n")

            for sample_id, data in sorted(sample_data.items()):
                if data["variants"]:
                    # Remove duplicate variants and sort numerically
                    unique_variants = list(set(data["variants"]))
                    unique_variants.sort(
                        key=lambda x: (
                            int(re.findall(r"\\d+", x)[0])
                            if re.findall(r"\\d+", x)
                            else 0
                        )
                    )

                    variants_str = "\\t".join(unique_variants)
                    regions_str = "+".join(sorted(set(data["regions"])))

                    f.write(f"{sample_id}\\t1-16569\\t{variants_str}\\n")

                    print(
                        f"✅ {sample_id}: {len(unique_variants)} variants ({regions_str})"
                    )

        print(f"\\n🎉 HSD conversion completed!")
        print(f"📄 Output file: {output_file}")
        print(f"🌐 Ready for HaploGrep analysis: https://haplogrep.i-med.ac.at/")


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description="Convert HVS consensus files to HSD using whole reference alignment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python whole_reference_hsd_converter.py consensus/ output.hsd
  python whole_reference_hsd_converter.py output_min30_q30/consensus/ haplogroups.hsd

This approach aligns each consensus sequence to the complete mitochondrial reference
to automatically determine the best matching position and extract variants.
        """,
    )

    parser.add_argument(
        "consensus_dir", help="Directory containing consensus FASTA files"
    )
    parser.add_argument("output_hsd", help="Output HSD file path")
    parser.add_argument(
        "--reference",
        "-r",
        default="ref/rCRS.fasta",
        help="Reference FASTA file (default: ref/rCRS.fasta)",
    )

    args = parser.parse_args()

    print("🧬 Whole Reference HSD Converter")
    print("=" * 50)

    try:
        converter = WholeReferenceHSDConverter(reference_file=args.reference)
        converter.convert_directory(args.consensus_dir, args.output_hsd)

    except KeyboardInterrupt:
        print("\\n⚠️  Conversion interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
