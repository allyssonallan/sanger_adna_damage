#!/usr/bin/env python3
"""
Simple test runner for the Sanger aDNA pipeline.
This script provides a quick way to run tests and check basic functionality.
"""

import subprocess
import sys
from pathlib import Path

UNICODE_FALLBACKS = {
    "🧪": "[TEST]",
    "🚀": "[SMOKE]",
    "⚠️": "[WARN]",
    "✅": "[OK]",
    "❌": "[ERROR]",
    "📦": "[IMPORTS]",
}


def console_print(text=""):
    """Print safely on consoles without Unicode support."""
    try:
        print(text)
    except UnicodeEncodeError:
        # Keep logs readable on Windows terminals that cannot render emoji.
        fallback_text = text
        for original, replacement in UNICODE_FALLBACKS.items():
            fallback_text = fallback_text.replace(original, replacement)
        print(fallback_text)


def run_pytest():
    """Run pytest with appropriate arguments."""
    console_print("🧪 Running Sanger aDNA Pipeline Tests with pytest")
    console_print("=" * 60)

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
        console_print(f"❌ Failed to run pytest: {e}")
        return 1


def run_smoke_tests():
    """Run only the smoke/integration tests."""
    console_print("🚀 Running Smoke Tests")
    console_print("=" * 30)

    try:
        # Probe pytest first so we can fall back to import-level smoke checks.
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
            "--no-cov",
        ]
        try:
            result = subprocess.run(cmd, check=False)
            return result.returncode
        except Exception as e:
            console_print(f"❌ Failed to run smoke tests with pytest: {e}")
            return 1

    console_print("⚠️  pytest not available, running basic smoke tests...")
    return run_basic_smoke_tests()


def run_basic_smoke_tests():
    """Run basic smoke tests without pytest dependency."""
    console_print("Running basic import and functionality checks...")

    # Allow running this script directly from repo root without editable install.
    sys.path.insert(0, str(Path(__file__).parent / "src"))

    tests_passed = 0
    total_tests = 0

    total_tests += 1
    console_print("1. Testing imports...")
    try:
        from sanger_pipeline.cli.main import cli
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        from sanger_pipeline.core.pipeline import SangerPipeline

        console_print("   ✅ All imports successful")
        tests_passed += 1

        total_tests += 1
        console_print("2. Testing component initialization...")
        analyzer = ADNADamageAnalyzer()
        assert hasattr(analyzer, "analyze_sequence_damage")
        console_print("   ✅ Component initialization successful")
        tests_passed += 1

        total_tests += 1
        console_print("3. Testing CLI commands...")
        command_names = [cmd.name for cmd in cli.commands.values()]
        assert "run" in command_names
        assert "analyze-damage" in command_names
        console_print("   ✅ CLI commands available")
        tests_passed += 1

    except Exception as e:
        console_print(f"   ❌ Test failed: {e}")

    console_print(f"\nResults: {tests_passed}/{total_tests} basic smoke tests passed")
    return 0 if tests_passed == total_tests else 1


def check_imports():
    """Quick import check without pytest."""
    console_print("📦 Checking imports...")

    try:
        from sanger_pipeline.cli.main import cli
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        from sanger_pipeline.core.pipeline import SangerPipeline

        console_print("✅ All imports successful")
        return True
    except Exception as e:
        console_print(f"❌ Import failed: {e}")
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
    if args.mode == "smoke":
        return run_smoke_tests()
    return run_pytest()


if __name__ == "__main__":
    sys.exit(main())
