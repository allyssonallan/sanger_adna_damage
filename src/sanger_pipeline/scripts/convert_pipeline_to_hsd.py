#!/usr/bin/env python3
"""
Convert pipeline output FASTA files to HSD format for haplogroup analysis.

Usage:
    python -m sanger_pipeline.scripts.convert_pipeline_to_hsd <input_directory> <output_hsd_file> [options]
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import the converter
sys.path.append(str(Path(__file__).parent.parent))

from scripts.bwa_aligned_hsd_converter import BWAAlignedHSDConverter


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python convert_pipeline_to_hsd.py <pipeline_output_dir> <output_hsd_file> [--reference ref_file]"
        )
        print("Example: python convert_pipeline_to_hsd.py output/ samples.hsd")
        print(
            "Example: python convert_pipeline_to_hsd.py output/ samples.hsd --reference ref/rCRS.fasta"
        )
        sys.exit(1)

    pipeline_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    # Check for reference argument
    reference_file = "ref/rCRS.fasta"  # Default reference
    if len(sys.argv) > 3 and sys.argv[3] == "--reference":
        if len(sys.argv) > 4:
            reference_file = sys.argv[4]
        else:
            print("Error: --reference flag requires a file path")
            sys.exit(1)

    if not pipeline_dir.exists():
        print(f"Error: Pipeline output directory not found: {pipeline_dir}")
        sys.exit(1)

    try:
        # Initialize converter (use BWA-MEM for proper alignment)
        print("🧬 Initializing BWA HSD converter...")
        converter = BWAAlignedHSDConverter()

        # Convert pipeline output to HSD
        print("🔄 Converting pipeline output to HSD format...")
        sample_variants = converter.process_consensus_directory(str(pipeline_dir))
        converter.write_hsd_file(sample_variants, str(output_file))

        print(f"\n✅ HSD file created: {output_file}")
        print(
            "📋 You can now use this file with HaploGrep for haplogroup classification"
        )
        print("🌐 HaploGrep web interface: https://haplogrep.i-med.ac.at/")

    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        print("💡 Tips:")
        print("  - Make sure your reference sequence is available in ref/ directory")
        print("  - Check that your pipeline output contains FASTA files")
        print("  - Verify file permissions and paths")
        sys.exit(1)


if __name__ == "__main__":
    main()
