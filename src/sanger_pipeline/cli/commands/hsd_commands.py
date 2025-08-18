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
@click.option('--consensus-dir', '-i', required=True, help='Directory containing consensus FASTA files')
@click.option('--output', '-o', required=True, help='Output HSD file')
@click.option('--method', '-m', default='direct', type=click.Choice(['direct', 'aligned']), 
              help='Conversion method (default: direct)')
def hsd_enhanced(consensus_dir, output, method):
    """Enhanced HSD converter with alignment options."""
    try:
        from ...scripts.enhanced_hsd_converter import EnhancedHSDConverter
        
        click.echo(f"🧬 Converting consensus files to HSD format using enhanced converter ({method} method)...")
        
        converter = EnhancedHSDConverter()
        converter.convert_to_hsd(consensus_dir, output, method)
        
        click.echo(f"✅ Enhanced HSD conversion completed: {output}")
        
    except Exception as e:
        click.echo(f"❌ Error during enhanced HSD conversion: {e}", err=True)
        raise click.Abort()


@hsd.command("regional")
@click.option('--consensus-dir', '-i', required=True, help='Directory containing consensus FASTA files')
@click.option('--output', '-o', required=True, help='Output HSD file')
def hsd_regional(consensus_dir, output):
    """Regional HSD converter (align each HVS region separately)."""
    try:
        from ...scripts.regional_hsd_converter import RegionalHSDConverter
        
        click.echo("🧬 Converting consensus files to HSD format using regional converter...")
        
        converter = RegionalHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output)
        
        click.echo(f"✅ Regional HSD conversion completed: {output}")
        
    except Exception as e:
        click.echo(f"❌ Error during regional HSD conversion: {e}", err=True)
        raise click.Abort()


@hsd.command("hybrid")
@click.option('--consensus-dir', '-i', required=True, help='Directory containing consensus FASTA files')
@click.option('--output', '-o', required=True, help='Output HSD file')
def hsd_hybrid(consensus_dir, output):
    """Hybrid HSD converter (direct comparison without alignment)."""
    try:
        from ...scripts.hybrid_regional_hsd_converter import HybridRegionalHSDConverter
        
        click.echo("🧬 Converting consensus files to HSD format using hybrid converter...")
        
        converter = HybridRegionalHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output)
        
        click.echo(f"✅ Hybrid HSD conversion completed: {output}")
        
    except Exception as e:
        click.echo(f"❌ Error during hybrid HSD conversion: {e}", err=True)
        raise click.Abort()


@hsd.command("pipeline")
@click.option('--input-dir', '-i', required=True, help='Pipeline output directory containing FASTA files')
@click.option('--output', '-o', required=True, help='Output HSD file')
@click.option('--reference', '-r', default='ref/rCRS.fasta', help='Reference sequence file')
def hsd_pipeline(input_dir, output, reference):
    """Convert pipeline output to HSD format."""
    try:
        from ...scripts.convert_pipeline_to_hsd import main as pipeline_main
        
        click.echo("🧬 Converting pipeline output to HSD format...")
        
        # Temporarily replace sys.argv for the script
        original_argv = sys.argv
        sys.argv = ['convert_pipeline_to_hsd.py', input_dir, output, '--reference', reference]
        
        pipeline_main()
        
        sys.argv = original_argv
        
        click.echo(f"✅ Pipeline HSD conversion completed: {output}")
        
    except Exception as e:
        click.echo(f"❌ Error during pipeline HSD conversion: {e}", err=True)
        raise click.Abort()
