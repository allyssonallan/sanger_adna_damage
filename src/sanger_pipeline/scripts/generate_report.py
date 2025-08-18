#!/usr/bin/env python3
"""
Generate a beautiful QC report for the Sanger pipeline.

Usage:
    python -m sanger_pipeline.scripts.generate_report [output_directory]
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sanger_pipeline.utils.report_generator import QCReportGenerator


def main():
    """Generate QC report."""

    # Get output directory from command line or use default
    if len(sys.argv) > 1:
        output_dir = Path(sys.argv[1])
    else:
        output_dir = Path("output")

    if not output_dir.exists():
        print(f"Error: Output directory {output_dir} does not exist")
        print("Please run the pipeline first or specify a valid output directory")
        return 1

    print(f"Generating QC report for: {output_dir}")

    try:
        # Generate report
        report_generator = QCReportGenerator(output_dir)
        report_file = report_generator.generate_report()

        print("✅ QC report generated successfully!")
        print(f"📄 Report file: {report_file}")
        print(f"🌐 Open in browser: file://{report_file.absolute()}")

        return 0

    except Exception as e:
        print(f"❌ Error generating report: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
