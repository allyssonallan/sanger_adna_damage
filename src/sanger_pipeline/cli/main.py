"""
Command-line interface for the Sanger pipeline.
"""

import click
from pathlib import Path

from ..core.pipeline import SangerPipeline
from ..core.ab1_converter import AB1Converter
from ..utils.helpers import setup_logging


@click.group()
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def cli(verbose):
    """Sanger DNA damage analysis pipeline."""
    level = "DEBUG" if verbose else "INFO"
    setup_logging(level)


@cli.command()
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
@click.option("--config", "-c", type=click.Path(exists=True), help="Configuration file (YAML)")
@click.option("--min-quality", "-q", default=20, help="Minimum Phred quality score (default: 20)")
@click.option("--min-sequence-length", "-l", default=30, help="Minimum sequence length after filtering (default: 30)")
@click.option("--alignment-tool", default="mafft", help="Alignment tool (default: mafft)")
@click.option("--alignment-params", default="--auto", help="Alignment parameters (default: --auto)")
def run_pipeline(input_dir, output_dir, config, min_quality, min_sequence_length, alignment_tool, alignment_params):
    """Run the complete Sanger sequencing pipeline."""
    click.echo(f"Running Sanger pipeline: {input_dir} -> {output_dir}")

    try:
        pipeline = SangerPipeline(
            input_dir=Path(input_dir),
            output_dir=Path(output_dir),
            config_file=Path(config) if config else None,
            min_quality=min_quality,
            min_sequence_length=min_sequence_length,
            alignment={"tool": alignment_tool, "parameters": alignment_params},
        )

        pipeline.run()

        # Print summary
        summary = pipeline.get_summary()
        click.echo("\nPipeline Summary:")
        for key, value in summary.items():
            click.echo(f"  {key.replace('_', ' ').title()}: {value}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


@cli.command()
@click.argument("ab1_file", type=click.Path(exists=True))
@click.argument("output_fasta", type=click.Path())
@click.option("--min-quality", "-q", default=20, help="Minimum Phred quality score (default: 20)")
@click.option("--min-sequence-length", "-l", default=30, help="Minimum sequence length after filtering (default: 30)")
@click.option("--generate-plot", is_flag=True, help="Generate quality plot")
def convert_ab1(ab1_file, output_fasta, min_quality, min_sequence_length, generate_plot):
    """Convert single AB1 file to FASTA format."""
    click.echo(f"Converting {ab1_file} to {output_fasta}")

    try:
        converter = AB1Converter(min_quality=min_quality, min_sequence_length=min_sequence_length)
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
                click.echo(f"Sequence excluded due to insufficient length (minimum: {min_sequence_length} valid bases)")

        # Generate plot if requested
        if generate_plot:
            plot_path = output_path.with_suffix(".png")
            converter.generate_quality_plot(record, plot_path)
            click.echo(f"Generated quality plot: {plot_path}")

        click.echo("Conversion completed successfully")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


@cli.command()
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


@cli.command()
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Input FASTA file to analyze",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(exists=True),
    required=True,
    help="Reference sequence file (FASTA)",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False),
    required=True,
    help="Output directory for damage analysis results",
)
@click.option("--sample-name", "-s", help="Sample name (default: filename)")
def analyze_damage(input_file, reference, output_dir, sample_name):
    """Analyze aDNA damage patterns in a single sequence."""
    from ..core.adna_damage_analyzer import ADNADamageAnalyzer
    from pathlib import Path
    import json

    if not sample_name:
        sample_name = Path(input_file).stem

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    click.echo(f"Analyzing aDNA damage patterns for: {sample_name}")

    try:
        # Initialize damage analyzer
        analyzer = ADNADamageAnalyzer()

        # Analyze damage patterns for single file
        input_path = Path(input_file)
        ref_path = Path(reference)
        
        results = analyzer.analyze_sequence_damage(input_path, ref_path)

        # Generate plots for the single file
        analyzer.generate_damage_plots([input_path], ref_path, output_path)

        # Perform bootstrap analysis using the single file
        bootstrap_results = analyzer.bootstrap_damage_analysis([input_path], ref_path)

        # Assess damage indicators
        damage_assessment = analyzer.assess_damage_indicators(bootstrap_results)

        # Save results
        results_file = output_path / f"{sample_name}_damage_results.json"
        with open(results_file, 'w') as f:
            json.dump({
                'damage_patterns': results,
                'bootstrap_analysis': bootstrap_results,
                'damage_assessment': damage_assessment
            }, f, indent=2)

        click.echo("Damage analysis complete")
        click.echo(f"  Results saved to: {results_file}")
        click.echo(f"  Plots saved to: {output_path}")
        click.echo(f"  Status: {damage_assessment['status']}")
        if damage_assessment['status'] == "DAMAGE_INDICATED":
            click.echo("  ✓ Damage patterns are indicative of ancient DNA")
        elif damage_assessment['status'] == "PARTIAL_DAMAGE_SIGNATURE":
            click.echo("  ~ Partial damage signature detected")
        else:
            click.echo("  ✗ No significant damage signature")

    except Exception as e:
        click.echo(f"Error analyzing damage patterns: {e}", err=True)
        raise click.Abort()


@cli.command()
@click.option(
    "--results-dir",
    "-d",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Directory containing damage analysis results",
)
def damage_summary(results_dir):
    """Generate summary report of aDNA damage analysis results."""
    import json
    
    results_path = Path(results_dir)
    result_files = list(results_path.glob("*_damage_results.json"))

    if not result_files:
        click.echo("No damage analysis results found in directory")
        return

    click.echo("aDNA Damage Indicator Analysis Summary")
    click.echo("=" * 50)

    damage_indicated_samples = 0
    partial_damage_samples = 0
    total_samples = len(result_files)

    for result_file in result_files:
        try:
            with open(result_file, 'r') as f:
                data = json.load(f)
            
            sample_name = result_file.stem.replace("_damage_results", "")
            
            # Check for new or old format
            if 'damage_assessment' in data:
                assessment = data['damage_assessment']
                status = assessment.get('status', 'UNKNOWN')
            else:
                # Handle old format for backward compatibility
                assessment = data.get('authenticity_assessment', {})
                score = assessment.get('authenticity_score', 0)
                status = "DAMAGE_INDICATED" if score > 0.5 else "NO_DAMAGE_SIGNATURE"
            
            # Get sequence quality info
            damage_patterns = data.get('damage_patterns', {})
            n_percentage = damage_patterns.get('sequence_quality', {}).get('n_percentage', 0)
            valid_percentage = damage_patterns.get('sequence_quality', {}).get('valid_percentage', 100)
            
            if status == "DAMAGE_INDICATED":
                damage_indicated_samples += 1
                status_display = "✓ Damage Indicated"
            elif status == "PARTIAL_DAMAGE_SIGNATURE":
                partial_damage_samples += 1
                status_display = "~ Partial Damage"
            else:
                status_display = "✗ No Damage Signal"
            
            click.echo(f"{sample_name:20} | N: {n_percentage:4.1f}% | Valid: {valid_percentage:4.1f}% | {status_display}")
            
        except Exception as e:
            click.echo(f"Error reading {result_file}: {e}")

    click.echo("-" * 70)
    click.echo(f"Total samples: {total_samples}")
    click.echo(f"Damage indicated: {damage_indicated_samples}")
    click.echo(f"Partial damage signatures: {partial_damage_samples}")
    click.echo(f"No damage signatures: {total_samples - damage_indicated_samples - partial_damage_samples}")
    click.echo(f"Damage indication rate: {(damage_indicated_samples/total_samples)*100:.1f}%")
    
    if total_samples > 0:
        avg_valid = sum([
            json.load(open(f))['damage_patterns'].get('sequence_quality', {}).get('valid_percentage', 100)
            for f in result_files if f.exists()
        ]) / total_samples
        click.echo(f"Average sequence quality: {avg_valid:.1f}% valid bases")


@cli.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Pipeline output directory containing results",
)
@click.option(
    "--open-browser",
    is_flag=True,
    help="Automatically open the report in the default browser",
)
def generate_report(output_dir, open_browser):
    """Generate comprehensive QC report with analysis summaries."""
    click.echo(f"Generating QC report for: {output_dir}")
    
    try:
        from ..utils.report_generator import QCReportGenerator
        
        report_generator = QCReportGenerator(Path(output_dir))
        report_file = report_generator.generate_report()
        
        click.echo("✅ QC report generated successfully!")
        click.echo(f"📄 Report file: {report_file}")
        
        if open_browser:
            import webbrowser
            webbrowser.open(f"file://{report_file.absolute()}")
            click.echo("🌐 Report opened in browser")
        else:
            click.echo(f"🌐 Open in browser: file://{report_file.absolute()}")
            
    except Exception as e:
        click.echo(f"❌ Error generating report: {e}", err=True)
        raise click.ClickException(f"Report generation failed: {e}")


if __name__ == "__main__":
    cli()
