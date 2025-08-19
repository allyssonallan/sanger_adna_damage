"""
Core pipeline commands for the Sanger CLI.

This module contains commands related to the main pipeline execution
and core processing functionality.
"""

import click
from pathlib import Path

from ...core.pipeline import SangerPipeline
from ...core.enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter


@click.group()
def pipeline():
    """Core pipeline commands."""
    pass


@pipeline.command(name="run")
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Directory containing AB1 files",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False),
    required=True,
    help="Output directory for results",
)
@click.option(
    "--config", "-c", type=click.Path(exists=True), help="Configuration file (YAML)"
)
@click.option(
    "--min-quality", "-q", default=20, help="Minimum Phred quality score (default: 20)"
)
@click.option(
    "--min-sequence-length",
    "-l",
    default=30,
    help="Minimum sequence length after filtering (default: 30)",
)
@click.option(
    "--alignment-tool", default="mafft", help="Alignment tool (default: mafft)"
)
@click.option(
    "--alignment-params",
    default="--auto",
    help="Alignment parameters (default: --auto)",
)
def run_pipeline(
    input_dir,
    output_dir,
    config,
    min_quality,
    min_sequence_length,
    alignment_tool,
    alignment_params,
):
    """Run the complete Sanger sequencing pipeline."""
    click.echo(f"Running Sanger pipeline: {input_dir} -> {output_dir}")

    try:
        pipeline = SangerPipeline(
            input_dir=Path(input_dir),
            output_dir=Path(output_dir),
            config_file=Path(config) if config else None,
            min_quality=min_quality,
            min_sequence_length=min_sequence_length,
            alignment_tool=alignment_tool,
            alignment_params=alignment_params,
        )

        results = pipeline.run()
        click.echo(f"Pipeline completed successfully. Results: {results}")

    except Exception as e:
        click.echo(f"Pipeline failed: {e}", err=True)
        raise click.Abort()


@pipeline.command(name="convert-ab1")
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
    """Convert single AB1 file to FASTA format with quality filtering."""
    click.echo(f"Converting AB1 file: {ab1_file} -> {output_fasta}")

    try:
        converter = AB1Converter(
            min_quality=min_quality, min_sequence_length=min_sequence_length
        )

        # Convert to FASTA
        result = converter.convert_to_fasta(Path(ab1_file), Path(output_fasta))

        # Generate plot if requested
        if generate_plot:
            plot_path = Path(output_fasta).with_suffix(".png")
            converter.generate_quality_plot(result, plot_path)
            click.echo(f"Generated quality plot: {plot_path}")

        click.echo(
            f"Conversion completed successfully. Sequence length: {len(result.seq)}"
        )

    except Exception as e:
        click.echo(f"Conversion failed: {e}", err=True)
        raise click.Abort()


@click.command()
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Directory to analyze",
)
def status(input_dir):
    """Show status of analysis directory."""
    input_path = Path(input_dir)

    # Count different file types in input
    ab1_files = list(input_path.glob("*.ab1"))

    # Infer output directory structure - assume output directory is sibling to input
    output_base = input_path.parent / "output"

    fasta_dir = output_base / "fasta"
    filtered_dir = output_base / "filtered"
    consensus_dir = output_base / "consensus"
    final_dir = output_base / "final"
    plots_dir = output_base / "plots"
    damage_dir = output_base / "damage_analysis"

    click.echo(f"Analysis Status for: {input_dir}")
    click.echo("=" * 50)
    click.echo(f"AB1 files: {len(ab1_files)}")

    if fasta_dir.exists():
        fasta_files = list(fasta_dir.glob("*.fasta"))
        click.echo(f"FASTA files: {len(fasta_files)}")

    if filtered_dir.exists():
        filtered_files = list(filtered_dir.glob("*_filtered.fasta"))
        click.echo(f"Filtered files: {len(filtered_files)}")

    if consensus_dir.exists():
        consensus_files = list(consensus_dir.glob("*_consensus.fasta"))
        click.echo(f"Consensus files: {len(consensus_files)}")

    if final_dir.exists():
        final_files = list(final_dir.glob("*_merged.fasta"))
        click.echo(f"Merged files: {len(final_files)}")

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
            click.echo("  HVS region combinations:")
            for combo, count in sorted(hvs_combinations.items()):
                click.echo(f"    {combo}: {count} samples")

    if plots_dir.exists():
        plot_files = list(plots_dir.glob("*_quality.png"))
        click.echo(f"Quality plots: {len(plot_files)}")

    # Check for damage analysis results
    if damage_dir.exists():
        damage_files = list(damage_dir.glob("*_damage_results.json"))
        damage_plots = list(damage_dir.glob("*_damage_plots.png"))
        click.echo(f"Damage analysis results: {len(damage_files)}")
        click.echo(f"Damage plots: {len(damage_plots)}")
