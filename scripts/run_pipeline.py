#!/usr/bin/env python3
"""
Main pipeline runner script for Sanger DNA damage analysis with HVS region processing.

This script provides command-line interface for:
- Running the complete HVS region-aware pipeline
- Converting single AB1 files
- Checking analysis status with HVS region breakdown
"""

import sys
import argparse
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Now import our modules  # noqa: E402
from sanger_pipeline.core.pipeline import SangerPipeline  # noqa: E402
from sanger_pipeline.core.enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter  # noqa: E402
from sanger_pipeline.utils.helpers import setup_logging  # noqa: E402


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Sanger DNA damage analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Pipeline command
    pipeline_parser = subparsers.add_parser("run-pipeline", help="Run the complete pipeline")
    pipeline_parser.add_argument(
        "-i", "--input-dir", required=True, help="Directory containing AB1 files"
    )
    pipeline_parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for results"
    )
    pipeline_parser.add_argument("-c", "--config", help="Configuration file (YAML)")
    pipeline_parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=20,
        help="Minimum Phred quality score (default: 20)",
    )
    pipeline_parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    # Convert command
    convert_parser = subparsers.add_parser("convert-ab1", help="Convert single AB1 file")
    convert_parser.add_argument("ab1_file", help="Input AB1 file")
    convert_parser.add_argument("output_fasta", help="Output FASTA file")
    convert_parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=20,
        help="Minimum Phred quality score (default: 20)",
    )
    convert_parser.add_argument(
        "--generate-plot", action="store_true", help="Generate quality plot"
    )
    convert_parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    # Status command
    status_parser = subparsers.add_parser("status", help="Show analysis status")
    status_parser.add_argument("-i", "--input-dir", required=True, help="Directory to analyze")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Setup logging
    log_level = "DEBUG" if getattr(args, "verbose", False) else "INFO"
    setup_logging(log_level)

    try:
        if args.command == "run-pipeline":
            return run_pipeline(args)
        elif args.command == "convert-ab1":
            return convert_ab1(args)
        elif args.command == "status":
            return show_status(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


def run_pipeline(args):
    """Run the complete pipeline."""
    print(f"Running Sanger pipeline: {args.input_dir} -> {args.output_dir}")

    pipeline = SangerPipeline(
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        config_file=Path(args.config) if args.config else None,
        min_quality=args.min_quality,
    )

    pipeline.run()

    # Print summary
    summary = pipeline.get_summary()
    print("\nPipeline Summary:")
    for key, value in summary.items():
        print(f"  {key.replace('_', ' ').title()}: {value}")

    return 0


def convert_ab1(args):
    """Convert single AB1 file."""
    print(f"Converting {args.ab1_file} to {args.output_fasta}")

    converter = AB1Converter(min_quality=args.min_quality)
    output_path = Path(args.output_fasta)

    # Convert to FASTA
    record = converter.convert_to_fasta(Path(args.ab1_file), output_path)

    # Generate filtered/processed version using enhanced processing
    if args.min_quality > 0:
        filtered_path = output_path.with_suffix("").with_suffix("_processed.fasta")
        plot_path = output_path.with_suffix(".png") if args.generate_plot else None
        
        # Use enhanced processing
        original_record, processed_record, stats = converter.process_ab1_file_enhanced(
            Path(args.ab1_file), output_path, filtered_path, plot_path
        )
        
        if processed_record is not None:
            print(f"Generated processed FASTA: {filtered_path}")
        else:
            print("Sequence excluded due to insufficient quality")

    # Generate plot if requested and not already generated
    elif args.generate_plot:
        plot_path = output_path.with_suffix(".png")
        converter.generate_quality_plot(record, plot_path)
        print(f"Generated quality plot: {plot_path}")

    print("Conversion completed successfully")
    return 0


def show_status(args):
    """Show analysis status."""
    input_path = Path(args.input_dir)

    # Count different file types
    ab1_files = list(input_path.glob("*.ab1"))
    fasta_dir = input_path / "fasta"
    filtered_dir = input_path / "filtered"
    consensus_dir = input_path / "consensus"
    final_dir = input_path / "final"
    plots_dir = input_path / "plots"

    print(f"Analysis Status for: {args.input_dir}")
    print("=" * 50)
    print(f"AB1 files: {len(ab1_files)}")

    if fasta_dir.exists():
        fasta_files = list(fasta_dir.glob("*.fasta"))
        print(f"FASTA files: {len(fasta_files)}")

    if filtered_dir.exists():
        filtered_files = list(filtered_dir.glob("*_filtered.fasta"))
        print(f"Filtered files: {len(filtered_files)}")

    if consensus_dir.exists():
        consensus_files = list(consensus_dir.glob("*_consensus.fasta"))
        print(f"Consensus files: {len(consensus_files)}")

    if final_dir.exists():
        final_files = list(final_dir.glob("*_merged.fasta"))
        print(f"Merged files: {len(final_files)}")
        
        # Show breakdown of HVS combinations
        hvs_combinations = {}
        for f in final_files:
            # Extract HVS combination from filename
            name = f.stem
            if "_HVS" in name:
                hvs_part = name.split("_HVS", 1)[1].split("_merged")[0]
                hvs_combo = f"HVS{hvs_part.replace('_HVS', '_HVS')}"
                hvs_combinations[hvs_combo] = hvs_combinations.get(hvs_combo, 0) + 1
        
        if hvs_combinations:
            print("  HVS region combinations:")
            for combo, count in sorted(hvs_combinations.items()):
                print(f"    {combo}: {count} samples")

    if plots_dir.exists():
        plot_files = list(plots_dir.glob("*_quality.png"))
        print(f"Quality plots: {len(plot_files)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
