#!/usr/bin/env python3
"""
Simple test runner for the Sanger aDNA pipeline.
This script provides a quick way to run tests and check basic functionality.
"""

import sys
import subprocess
from pathlib import Path


def run_pytest():
    """Run pytest with appropriate arguments."""
    print("🧪 Running Sanger aDNA Pipeline Tests with pytest")
    print("=" * 60)

    # Run pytest with coverage and verbose output
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "tests/",
        "-v",
        "--tb=short",
        "--cov=src/sanger_pipeline",
        "--cov-report=term-missing",
    ]

    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except Exception as e:
        print(f"❌ Failed to run pytest: {e}")
        return 1


def run_smoke_tests():
    """Run only the smoke/integration tests."""
    print("🚀 Running Smoke Tests")
    print("=" * 30)

    # Check if pytest is available
    try:
        subprocess.run(
            [sys.executable, "-m", "pytest", "--version"],
            check=True,
            capture_output=True,
        )
        pytest_available = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        pytest_available = False

    if pytest_available:
        cmd = [
            sys.executable,
            "-m",
            "pytest",
            "tests/test_integration_smoke.py",
            "-v",
            "--no-cov",  # Disable coverage for smoke tests
        ]
        try:
            result = subprocess.run(cmd, check=False)
            return result.returncode
        except Exception as e:
            print(f"❌ Failed to run smoke tests with pytest: {e}")
            return 1
    else:
        print("⚠️  pytest not available, running basic smoke tests...")
        return run_basic_smoke_tests()


def run_basic_smoke_tests():
    """Run basic smoke tests without pytest dependency."""
    print("Running basic import and functionality checks...")

    # Add src to path for imports
    sys.path.insert(0, str(Path(__file__).parent / "src"))

    tests_passed = 0
    total_tests = 0

    # Test 1: Import check
    total_tests += 1
    print("1. Testing imports...")
    try:
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        from sanger_pipeline.core.pipeline import SangerPipeline
        from sanger_pipeline.cli.main import cli

        print("   ✅ All imports successful")
        tests_passed += 1

        # Test 2: Basic initialization
        total_tests += 1
        print("2. Testing component initialization...")
        analyzer = ADNADamageAnalyzer()
        assert hasattr(analyzer, "analyze_sequence_damage")
        print("   ✅ Component initialization successful")
        tests_passed += 1

        # Test 3: CLI commands check
        total_tests += 1
        print("3. Testing CLI commands...")
        command_names = [cmd.name for cmd in cli.commands.values()]
        assert "run" in command_names
        assert "analyze-damage" in command_names
        print("   ✅ CLI commands available")
        tests_passed += 1

    except Exception as e:
        print(f"   ❌ Test failed: {e}")

    print(f"\nResults: {tests_passed}/{total_tests} basic smoke tests passed")
    return 0 if tests_passed == total_tests else 1


def check_imports():
    """Quick import check without pytest."""
    print("📦 Checking imports...")

    try:
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        from sanger_pipeline.core.pipeline import SangerPipeline
        from sanger_pipeline.cli.main import cli

        print("✅ All imports successful")
        return True
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Run Sanger aDNA pipeline tests")
    parser.add_argument(
        "--mode",
        choices=["full", "smoke", "imports"],
        default="full",
        help="Test mode: full (all tests), smoke (integration only), imports (import check only)",
    )

    args = parser.parse_args()

    if args.mode == "imports":
        success = check_imports()
        return 0 if success else 1
    elif args.mode == "smoke":
        return run_smoke_tests()
    else:
        return run_pytest()


if __name__ == "__main__":
    sys.exit(main())
