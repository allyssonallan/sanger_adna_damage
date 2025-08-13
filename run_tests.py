#!/usr/bin/env python3
"""
Simple test runner for the Sanger aDNA pipeline without external dependencies.
"""

import sys
import tempfile
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    
    try:
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        from sanger_pipeline.core.pipeline import SangerPipeline
        from sanger_pipeline.core.ab1_converter import AB1Converter
        from sanger_pipeline.core.consensus_builder import ConsensusBuilder
        from sanger_pipeline.cli.main import cli
        print("✅ All imports successful")
        return True
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False

def test_adna_analyzer():
    """Test ADNADamageAnalyzer initialization."""
    print("Testing aDNA analyzer...")
    
    try:
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        analyzer = ADNADamageAnalyzer()
        
        # Check required methods exist
        assert hasattr(analyzer, 'analyze_sequence_damage')
        assert hasattr(analyzer, 'bootstrap_damage_analysis')
        assert hasattr(analyzer, 'assess_authenticity')
        assert hasattr(analyzer, 'generate_damage_plots')
        
        print("✅ ADNADamageAnalyzer tests passed")
        return True
    except Exception as e:
        print(f"❌ aDNA analyzer test failed: {e}")
        return False

def test_pipeline_integration():
    """Test pipeline integration with aDNA analyzer."""
    print("Testing pipeline integration...")
    
    try:
        from sanger_pipeline.core.pipeline import SangerPipeline
        from sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
        
        with tempfile.TemporaryDirectory() as temp_dir:
            input_dir = Path(temp_dir) / 'input'
            output_dir = Path(temp_dir) / 'output'
            input_dir.mkdir()
            
            pipeline = SangerPipeline(input_dir, output_dir)
            
            # Check pipeline has damage analyzer
            assert hasattr(pipeline, 'damage_analyzer')
            assert hasattr(pipeline, '_step_4_adna_damage_analysis')
            
            # Check analyzer is correct type
            assert isinstance(pipeline.damage_analyzer, ADNADamageAnalyzer)
            
        print("✅ Pipeline integration tests passed")
        return True
    except Exception as e:
        print(f"❌ Pipeline integration test failed: {e}")
        return False

def test_cli_commands():
    """Test CLI commands are available."""
    print("Testing CLI commands...")
    
    try:
        from sanger_pipeline.cli.main import cli
        
        command_names = [cmd.name for cmd in cli.commands.values()]
        
        # Check required commands exist
        assert 'run-pipeline' in command_names
        assert 'analyze-damage' in command_names
        assert 'damage-summary' in command_names
        assert 'status' in command_names
        
        print("✅ CLI command tests passed")
        return True
    except Exception as e:
        print(f"❌ CLI command test failed: {e}")
        return False

def test_configuration():
    """Test configuration loading."""
    print("Testing configuration...")
    
    try:
        from sanger_pipeline.utils.helpers import load_config
        
        config_file = Path(__file__).parent / "config" / "default_config.yaml"
        if config_file.exists():
            config = load_config(config_file)
            
            # Check aDNA-specific settings
            assert 'bootstrap' in config
            assert 'damage' in config
            assert 'quality' in config
            
            print("✅ Configuration tests passed")
            return True
        else:
            print("⚠️  Configuration file not found, skipping test")
            return True
    except Exception as e:
        print(f"❌ Configuration test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Running Sanger aDNA Pipeline Tests")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_adna_analyzer,
        test_pipeline_integration,
        test_cli_commands,
        test_configuration
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 50)
    print(f"Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed!")
        return 0
    else:
        print(f"❌ {total - passed} tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
