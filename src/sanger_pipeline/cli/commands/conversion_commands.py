"""
Conversion command module - AB1 and format conversion utilities.

This module contains CLI commands for converting AB1 files and other format conversions.
"""

import click
from pathlib import Path

from ...core.enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter, EnhancedAB1Converter
from ...core.primer_config import parse_cli_primers


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
            processed_path = output_path.with_suffix("").with_suffix("_processed.fasta")
            plot_path = output_path.with_suffix(".png") if generate_plot else None
            
            # Use enhanced processing for quality filtering
            original_record, processed_record, stats = converter.process_ab1_file_enhanced(
                Path(ab1_file), output_path, processed_path, plot_path
            )
            
            if processed_record is not None:
                click.echo(f"Generated processed FASTA: {processed_path}")
                click.echo(f"Processing stats: {stats.get('final_length', 'N/A')} final bases")
            else:
                click.echo(
                    f"Sequence excluded due to insufficient length (minimum: {min_sequence_length} valid bases)"
                )

        # Generate plot if requested and not already generated
        if generate_plot and min_quality == 0:
            plot_path = output_path.with_suffix(".png")
            converter.generate_quality_plot(record, plot_path)
            click.echo(f"Generated quality plot: {plot_path}")

        click.echo("Conversion completed successfully")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


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
@click.option(
    "--primer-config", 
    type=click.Path(exists=True), 
    help="YAML file containing primer configurations"
)
@click.option(
    "--primer-forward", 
    help="Custom forward primers (format: region1:sequence1,region2:sequence2)"
)
@click.option(
    "--primer-reverse", 
    help="Custom reverse primers (format: region1:sequence1,region2:sequence2)"
)
@click.option(
    "--disable-primer-removal", 
    is_flag=True, 
    help="Disable automatic primer removal"
)
@click.option(
    "--adna-mode/--modern-mode", 
    default=True, 
    help="Enable ancient DNA mode with lower primer matching thresholds"
)
@click.option(
    "--show-primer-info", 
    is_flag=True, 
    help="Display detailed primer detection information"
)
def convert_ab1_enhanced(
    ab1_file, 
    output_fasta, 
    min_quality, 
    min_sequence_length, 
    generate_plot,
    primer_config,
    primer_forward,
    primer_reverse,
    disable_primer_removal,
    adna_mode,
    show_primer_info
):
    """Enhanced AB1 to FASTA conversion with primer removal and ancient DNA support."""
    click.echo(f"Converting AB1 file with enhanced features: {ab1_file} -> {output_fasta}")
    
    try:
        # Parse custom primers from CLI
        custom_forward = parse_cli_primers(primer_forward) if primer_forward else None
        custom_reverse = parse_cli_primers(primer_reverse) if primer_reverse else None
        
        # Initialize enhanced converter
        converter = EnhancedAB1Converter(
            min_quality=min_quality,
            min_sequence_length=min_sequence_length,
            enable_primer_removal=not disable_primer_removal,
            adna_damage_mode=adna_mode,
            custom_primers_forward=custom_forward,
            custom_primers_reverse=custom_reverse,
            primer_config_file=Path(primer_config) if primer_config else None
        )
        
        # Show primer configuration if requested
        if show_primer_info:
            click.echo("\n📋 Primer Configuration:")
            for region, primer_data in converter.primers.items():
                click.echo(f"  {region}:")
                click.echo(f"    Forward:  {primer_data.get('forward', 'N/A')}")
                click.echo(f"    Reverse:  {primer_data.get('reverse_original', 'N/A')}")
                if 'description' in primer_data:
                    click.echo(f"    Description: {primer_data['description']}")
            click.echo()
        
        # Create output paths
        ab1_path = Path(ab1_file)
        output_path = Path(output_fasta)
        processed_output = output_path.with_suffix('.processed.fasta')
        plot_output = output_path.with_suffix('.png')
        
        # Convert with enhanced features
        original_record, processed_record, stats = converter.process_ab1_file_enhanced(
            ab1_path, output_path, processed_output, plot_output
        )
        
        # Display processing statistics
        click.echo("\n📊 Processing Statistics:")
        click.echo(f"  Original length: {stats.get('original_length', 'N/A')} bases")
        click.echo(f"  Final length: {stats.get('final_length', 'N/A')} bases")
        click.echo(f"  HVS region detected: {stats.get('hvs_region', 'Unknown')}")
        click.echo(f"  Primers removed: {stats.get('primers_removed', False)}")
        
        if stats.get('primers_removed'):
            primer_info = stats.get('primer_removal_info', {})
            click.echo(f"    Forward removed: {primer_info.get('forward_removed', False)}")
            click.echo(f"    Reverse removed: {primer_info.get('reverse_removed', False)}")
            click.echo(f"    Length reduction: {primer_info.get('length_reduction', 0)} bases")
        
        click.echo(f"  Quality trimmed: {stats.get('quality_trimmed', (0, 0))}")
        
        if generate_plot:
            click.echo(f"Generated quality plot: {plot_output}")
        
        click.echo("\n✅ Enhanced conversion completed successfully")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


@click.command()
@click.option(
    "--primer-config", 
    type=click.Path(exists=True), 
    help="YAML file containing primer configurations to validate"
)
@click.option(
    "--show-details", 
    is_flag=True, 
    help="Show detailed primer information"
)
def validate_primers(primer_config, show_details):
    """Validate primer configuration file and show primer information."""
    from ...core.primer_config import PrimerConfig, validate_primer_file
    
    if not primer_config:
        # Use default config
        click.echo("🔍 Validating default primer configuration...")
        config = PrimerConfig()
    else:
        click.echo(f"🔍 Validating primer configuration: {primer_config}")
        valid_file, file_issues = validate_primer_file(Path(primer_config))
        
        if not valid_file:
            click.echo("❌ Configuration file validation failed:")
            for issue in file_issues:
                click.echo(f"  - {issue}")
            return
        
        config = PrimerConfig(Path(primer_config))
    
    # Validate primers
    valid, issues = config.validate_primers()
    
    if valid:
        click.echo("✅ All primers are valid!")
    else:
        click.echo("⚠️  Primer validation issues found:")
        for issue in issues:
            click.echo(f"  - {issue}")
    
    # Show primer details
    if show_details:
        click.echo("\n📋 Primer Configuration Details:")
        click.echo("=" * 50)
        
        for region in config.get_regions():
            primer_data = config.get_primers_for_region(region)
            if primer_data:
                click.echo(f"\n🧬 {region}:")
                click.echo(f"  Forward:  {primer_data.get('forward', 'N/A')}")
                click.echo(f"  Reverse:  {primer_data.get('reverse', 'N/A')}")
                click.echo(f"  Length:   F={len(primer_data.get('forward', ''))} / R={len(primer_data.get('reverse', ''))}")
                if 'description' in primer_data:
                    click.echo(f"  Description: {primer_data['description']}")
        
        # Show matching parameters
        params = config.matching_parameters
        click.echo("\n⚙️  Matching Parameters:")
        click.echo(f"  Similarity threshold: {params.get('similarity_threshold', 'N/A')}")
        click.echo(f"  aDNA threshold: {params.get('adna_similarity_threshold', 'N/A')}")
        click.echo(f"  Max mismatches: {params.get('max_mismatches', 'N/A')}")
        click.echo(f"  Search window: {params.get('search_window', 'N/A')}")
    
    click.echo(f"\n📊 Total regions configured: {len(config.get_regions())}")


@click.command()
@click.argument("output_file", type=click.Path())
@click.option(
    "--template-type", 
    type=click.Choice(['basic', 'comprehensive']), 
    default='basic',
    help="Template type: basic (common regions) or comprehensive (with examples)"
)
def generate_primer_config(output_file, template_type):
    """Generate a primer configuration template file."""
    from ...core.primer_config import PrimerConfig
    
    output_path = Path(output_file)
    
    if output_path.exists():
        if not click.confirm(f"File {output_file} already exists. Overwrite?"):
            click.echo("Operation cancelled.")
            return
    
    try:
        if template_type == 'basic':
            # Create basic template with default primers
            config = PrimerConfig()
            config.save_to_yaml(output_path)
            click.echo(f"✅ Basic primer configuration template saved to: {output_file}")
        else:
            # Create comprehensive template with examples
            import shutil
            template_source = Path(__file__).parent.parent.parent.parent.parent / 'config' / 'primers.yaml'
            shutil.copy(template_source, output_path)
            click.echo(f"✅ Comprehensive primer configuration template saved to: {output_file}")
        
        click.echo("\n📝 Next steps:")
        click.echo("1. Edit the primer sequences to match your laboratory primers")
        click.echo("2. Validate the configuration with: validate-primers --primer-config <file>")
        click.echo("3. Use in conversion with: convert-ab1-enhanced --primer-config <file> <input> <output>")
        
    except Exception as e:
        click.echo(f"❌ Error generating template: {e}", err=True)
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
def convert_to_hsd(consensus_dir, output, reference):
    """Convert HVS consensus sequences to HSD format using BWA-MEM alignment."""
    from ...scripts.bwa_aligned_hsd_converter import BWAAlignedHSDConverter
    from pathlib import Path

    click.echo("🧬 Converting consensus sequences to HSD format using BWA-MEM")
    click.echo(f"Input directory: {consensus_dir}")
    click.echo(f"Output file: {output}")
    click.echo(f"Reference: {reference}")

    try:
        # Use BWA-MEM for proper alignment mapping
        converter = BWAAlignedHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output)

        click.echo(f"✅ HSD conversion completed: {output}")
        click.echo("💡 Upload the HSD file to HaploGrep for haplogroup analysis:")
        click.echo("   https://haplogrep.i-med.ac.at/")

    except Exception as e:
        click.echo(f"❌ Error during HSD conversion: {e}", err=True)
        raise click.Abort()
