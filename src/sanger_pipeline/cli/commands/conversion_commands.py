"""
Conversion command module - AB1 and format conversion utilities.

This module contains CLI commands for converting AB1 files and other format conversions.
"""

import click
from pathlib import Path

from ...core.ab1_converter import AB1Converter


@click.command()
@click.argument("ab1_file", type=click.Path(exists=True))
@click.argument("output_fasta", type=click.Path())
@click.option(
    "--min-quality", "-q", default=20, help="Minimum Phred quality score (default: 20)"
)
@click.option(
    "--min-sequence-length",
    "-l",
    default=30,
    help="Minimum sequence length after filtering (default: 30)",
)
@click.option("--generate-plot", is_flag=True, help="Generate quality plot")
def convert_ab1(
    ab1_file, output_fasta, min_quality, min_sequence_length, generate_plot
):
    """Convert single AB1 file to FASTA format."""
    click.echo(f"Converting {ab1_file} to {output_fasta}")

    try:
        converter = AB1Converter(
            min_quality=min_quality, min_sequence_length=min_sequence_length
        )
        output_path = Path(output_fasta)

        # Convert to FASTA
        record = converter.convert_to_fasta(Path(ab1_file), output_path)

        # Generate filtered version if quality filtering is enabled
        if min_quality > 0:
            filtered_path = output_path.with_suffix("").with_suffix("_filtered.fasta")
            filtered_record = converter.filter_by_quality(record, filtered_path)
            if filtered_record is not None:
                click.echo(f"Generated filtered FASTA: {filtered_path}")
            else:
                click.echo(
                    f"Sequence excluded due to insufficient length (minimum: {min_sequence_length} valid bases)"
                )

        # Generate plot if requested
        if generate_plot:
            plot_path = output_path.with_suffix(".png")
            converter.generate_quality_plot(record, plot_path)
            click.echo(f"Generated quality plot: {plot_path}")

        click.echo("Conversion completed successfully")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


@click.command()
@click.option(
    "--consensus-dir",
    "-i",
    required=True,
    help="Directory containing consensus FASTA files",
)
@click.option("--output", "-o", required=True, help="Output HSD file")
@click.option(
    "--reference",
    "-r",
    default="ref/rCRS.fasta",
    help="Reference sequence file (default: ref/rCRS.fasta)",
)
@click.option(
    "--method",
    "-m",
    default="regional",
    type=click.Choice(["regional", "direct", "aligned"]),
    help="Conversion method (default: regional)",
)
def convert_to_hsd(consensus_dir, output, reference, method):
    """Convert HVS consensus sequences to HSD format for haplogroup analysis."""
    from ...utils.regional_hsd_converter import HybridRegionalHSDConverter
    from pathlib import Path

    click.echo("🧬 Converting consensus sequences to HSD format")
    click.echo(f"Input directory: {consensus_dir}")
    click.echo(f"Output file: {output}")
    click.echo(f"Method: {method}")
    click.echo(f"Reference: {reference}")

    try:
        if method == "regional":
            converter = HybridRegionalHSDConverter(reference)
            sample_variants = converter.process_consensus_directory(consensus_dir)
            converter.write_hsd_file(sample_variants, output)
        else:
            # For backwards compatibility, fallback to original converter for other methods
            from ...utils.fasta_to_hsd_converter import FastaToHSDConverter

            converter = FastaToHSDConverter(reference)
            converter.convert_pipeline_output(Path(consensus_dir), Path(output))

        click.echo(f"✅ HSD conversion completed: {output}")
        click.echo("💡 Upload the HSD file to HaploGrep for haplogroup analysis:")
        click.echo("   https://haplogrep.i-med.ac.at/")

    except Exception as e:
        click.echo(f"❌ Error during HSD conversion: {e}", err=True)
        raise click.Abort()
