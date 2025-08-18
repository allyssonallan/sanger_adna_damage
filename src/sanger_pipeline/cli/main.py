"""
Refactored command-line interface for the Sanger pipeline.

This module provides a clean, modular CLI structure using separate command modules.
"""

import click

from .commands import (
    run_pipeline,
    status,
    convert_ab1,
    convert_to_hsd,
    analyze_damage,
    damage_summary,
    hsd,
    generate_report,
)
from ..utils.helpers import setup_logging


@click.group()
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def cli(verbose):
    """Sanger DNA damage analysis pipeline."""
    level = "DEBUG" if verbose else "INFO"
    setup_logging(level)


# Add individual commands to the main CLI group
cli.add_command(run_pipeline)
cli.add_command(status)
cli.add_command(convert_ab1)
cli.add_command(convert_to_hsd)
cli.add_command(analyze_damage)
cli.add_command(damage_summary)
cli.add_command(generate_report)

# Add command groups
cli.add_command(hsd)


if __name__ == "__main__":
    cli()
