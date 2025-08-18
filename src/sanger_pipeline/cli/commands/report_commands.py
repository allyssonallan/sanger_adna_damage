"""
Report command module - Report generation utilities.

This module contains CLI commands for generating comprehensive QC reports.
"""

import click
import webbrowser
from pathlib import Path


@click.command()
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
        from ...utils.report_generator import QCReportGenerator
        
        report_generator = QCReportGenerator(Path(output_dir))
        report_file = report_generator.generate_report()
        
        click.echo("✅ QC report generated successfully!")
        click.echo(f"📄 Report file: {report_file}")
        
        if open_browser:
            webbrowser.open(f"file://{report_file.absolute()}")
            click.echo("🌐 Report opened in browser")
        else:
            click.echo(f"🌐 Open in browser: file://{report_file.absolute()}")
            
    except Exception as e:
        click.echo(f"❌ Error generating report: {e}", err=True)
        raise click.ClickException(f"Report generation failed: {e}")
