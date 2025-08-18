#!/usr/bin/env python3
"""
Convert pipeline output FASTA files to HSD format for haplogroup analysis.

Usage:
    python -m sanger_pipeline.scripts.convert_pipeline_to_hsd <input_directory> <output_hsd_file> [options]
"""

import sys
from pathlib import Path
from ...utils.fasta_to_hsd_converter import FastaToHSDConverter


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_pipeline_to_hsd.py <pipeline_output_dir> <output_hsd_file> [--reference ref_file]")
        print("Example: python convert_pipeline_to_hsd.py output/ samples.hsd")
        print("Example: python convert_pipeline_to_hsd.py output/ samples.hsd --reference ref/rCRS.fasta")
        sys.exit(1)
    
    pipeline_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    
    # Check for reference argument
    reference_file = None
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
        # Initialize converter
        print("🧬 Initializing FASTA to HSD converter...")
        converter = FastaToHSDConverter(reference_file)
        
        # Convert pipeline output to HSD
        print("🔄 Converting pipeline output to HSD format...")
        converter.convert_pipeline_output(pipeline_dir, output_file)
        
        print(f"\n✅ HSD file created: {output_file}")
        print(f"📋 You can now use this file with HaploGrep for haplogroup classification")
        print(f"🌐 HaploGrep web interface: https://haplogrep.i-med.ac.at/")
        
    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        print("💡 Tips:")
        print("  - Make sure your reference sequence is available in ref/ directory")
        print("  - Check that your pipeline output contains FASTA files")
        print("  - Verify file permissions and paths")
        sys.exit(1)


if __name__ == "__main__":
    main()
