"""
Analysis command module - Damage analysis and reporting commands.

This module contains CLI commands for analyzing aDNA damage patterns and generating reports.
"""

import click
import json
from pathlib import Path


@click.command()
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
    from ...core.adna_damage_analyzer import ADNADamageAnalyzer

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


@click.command()
@click.option(
    "--results-dir",
    "-d",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Directory containing damage analysis results",
)
def damage_summary(results_dir):
    """Generate summary report of aDNA damage analysis results."""
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
