#!/usr/bin/env python3
"""
Test script to verify all QC report features are working correctly.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from sanger_pipeline.utils.report_generator import QCReportGenerator


def test_report_generation():
    """Test QC report generation with sample data."""
    output_dir = Path("output")
    
    if not output_dir.exists():
        print("❌ Output directory 'output' not found. Please run the pipeline first.")
        return False
    
    print("🧪 Testing QC Report Generation...")
    
    try:
        # Test report generator
        report_generator = QCReportGenerator(output_dir)
        
        # Collect statistics
        print("📊 Collecting pipeline statistics...")
        stats = report_generator.collect_pipeline_statistics()
        
        # Check key statistics
        total_samples = len(stats['samples'])
        damage_files = stats['damage_analysis'].get('files_analyzed', 0)
        hvs_combinations = len(stats['hvs_combinations'].get('combinations', {}))
        final_files = stats['hvs_combinations'].get('total_merged_files', 0)
        
        print(f"  ✓ Found {total_samples} samples")
        print(f"  ✓ Found {damage_files} damage analysis files")
        print(f"  ✓ Found {hvs_combinations} HVS combinations")
        print(f"  ✓ Found {final_files} final merged files")
        
        # Generate report
        print("📝 Generating HTML report...")
        report_file = report_generator.generate_report()
        
        # Verify report file
        if report_file.exists():
            file_size = report_file.stat().st_size / 1024  # KB
            print(f"  ✓ Report generated: {report_file.name}")
            print(f"  ✓ File size: {file_size:.1f} KB")
            
            if file_size > 10:  # Should be substantial
                print("✅ QC Report generation test PASSED!")
                print(f"🌐 Open report: file://{report_file.absolute()}")
                return True
            else:
                print("❌ Report file too small - may be incomplete")
                return False
        else:
            print("❌ Report file not created")
            return False
            
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cli_commands():
    """Test CLI commands are available."""
    print("\n🔧 Testing CLI Commands...")
    
    try:
        from sanger_pipeline.cli.main import cli
        
        # Check that generate-report command exists
        command_names = [cmd.name for cmd in cli.commands.values() if cmd.name]
        
        if 'generate-report' in command_names:
            print("  ✓ generate-report command available")
            return True
        else:
            print("  ❌ generate-report command not found")
            print(f"  Available commands: {', '.join(command_names)}")
            return False
            
    except Exception as e:
        print(f"  ❌ CLI test failed: {e}")
        return False


def main():
    """Run all tests."""
    print("🧪 Testing Sanger Pipeline QC Report Features")
    print("=" * 50)
    
    # Test report generation
    report_test = test_report_generation()
    
    # Test CLI commands
    cli_test = test_cli_commands()
    
    # Final results
    print("\n" + "=" * 50)
    if report_test and cli_test:
        print("🎉 All tests PASSED! QC Report features are working correctly.")
        print("\n📋 Quick Usage Guide:")
        print("  Generate report via CLI:")
        print("    python -m src.sanger_pipeline.cli.main generate-report -o output --open-browser")
        print("  Generate report via script:")
        print("    python generate_report.py output")
        return 0
    else:
        print("❌ Some tests FAILED. Please check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
