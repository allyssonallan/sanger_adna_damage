"""
HSD command module - HSD format conversion utilities.

This module contains CLI commands for various HSD conversion methods.
"""

import click
import sys


@click.group()
def hsd():
    """HSD converter utilities."""
    pass


@hsd.command("enhanced")
@click.option(
    "--consensus-dir",
    "-i",
    required=True,
    help="Directory containing consensus FASTA files",
)
@click.option("--output", "-o", required=True, help="Output HSD file")
@click.option(
    "--method",
    "-m",
    default="aligned",
    type=click.Choice(["aligned", "direct"]),
    help="Conversion method (default: aligned)",
)
def hsd_enhanced(consensus_dir, output, method):
    """Enhanced HSD converter using BWA-MEM alignment (recommended)."""
    try:
        from ...scripts.bwa_aligned_hsd_converter import BWAAlignedHSDConverter

        click.echo(
            f"🧬 Converting consensus files to HSD format using BWA-MEM alignment..."
        )

        converter = BWAAlignedHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output)

        click.echo(f"✅ Enhanced HSD conversion completed: {output}")

    except Exception as e:
        click.echo(f"❌ Error during enhanced HSD conversion: {e}", err=True)
        raise click.Abort()


@hsd.command("bwa")
@click.option(
    "--consensus-dir",
    "-i",
    required=True,
    help="Directory containing consensus FASTA files",
)
@click.option("--output", "-o", required=True, help="Output HSD file")
def hsd_bwa(consensus_dir, output):
    """BWA-MEM based HSD converter (recommended - uses proper alignment mapper)."""
    try:
        from ...scripts.bwa_aligned_hsd_converter import BWAAlignedHSDConverter

        click.echo(
            "🧬 Converting consensus files to HSD format using BWA-MEM alignment..."
        )

        converter = BWAAlignedHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output)

        click.echo(f"✅ BWA HSD conversion completed: {output}")

    except Exception as e:
        click.echo(f"❌ Error during BWA HSD conversion: {e}", err=True)
        raise click.Abort()


@hsd.command("pipeline")
@click.option(
    "--input-dir",
    "-i",
    required=True,
    help="Pipeline output directory containing FASTA files",
)
@click.option("--output", "-o", required=True, help="Output HSD file")
@click.option(
    "--reference", "-r", default="ref/rCRS.fasta", help="Reference sequence file"
)
def hsd_pipeline(input_dir, output, reference):
    """Convert pipeline output to HSD format."""
    try:
        from ...scripts.convert_pipeline_to_hsd import main as pipeline_main

        click.echo("🧬 Converting pipeline output to HSD format...")

        # Temporarily replace sys.argv for the script
        original_argv = sys.argv
        sys.argv = [
            "convert_pipeline_to_hsd.py",
            input_dir,
            output,
            "--reference",
            reference,
        ]

        pipeline_main()

        sys.argv = original_argv

        click.echo(f"✅ Pipeline HSD conversion completed: {output}")

    except Exception as e:
        click.echo(f"❌ Error during pipeline HSD conversion: {e}", err=True)
        raise click.Abort()
