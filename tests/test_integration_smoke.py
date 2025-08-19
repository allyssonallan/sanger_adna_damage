"""
Integration tests that mirror the functionality from run_tests.py
These tests verify the basic functionality of the pipeline components.
"""

import tempfile
from pathlib import Path
import pytest

from src.sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer
from src.sanger_pipeline.core.pipeline import SangerPipeline
from src.sanger_pipeline.core.enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter
from src.sanger_pipeline.core.consensus_builder import ConsensusBuilder
from src.sanger_pipeline.cli.main import cli


@pytest.mark.integration
class TestPipelineIntegration:
    """Integration tests for the pipeline components."""

    @pytest.mark.unit
    def test_imports(self):
        """Test that all modules can be imported."""
        # All imports done at module level, so if we get here, imports worked
        assert ADNADamageAnalyzer is not None
        assert SangerPipeline is not None
        assert AB1Converter is not None
        assert ConsensusBuilder is not None
        assert cli is not None

    @pytest.mark.unit
    def test_adna_analyzer_initialization(self):
        """Test ADNADamageAnalyzer initialization and basic interface."""
        analyzer = ADNADamageAnalyzer()
        
        # Check required methods exist
        assert hasattr(analyzer, 'analyze_sequence_damage')
        assert hasattr(analyzer, 'bootstrap_damage_analysis')
        assert hasattr(analyzer, 'assess_authenticity')
        assert hasattr(analyzer, 'generate_damage_plots')

    @pytest.mark.integration
    def test_pipeline_integration_with_damage_analyzer(self):
        """Test pipeline integration with aDNA analyzer."""
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

    @pytest.mark.unit
    def test_cli_commands_availability(self):
        """Test CLI commands are available."""
        command_names = [cmd.name for cmd in cli.commands.values()]
        
        # Check required commands exist
        assert 'run' in command_names
        assert 'analyze-damage' in command_names
        assert 'damage-summary' in command_names
        assert 'status' in command_names

    @pytest.mark.unit
    def test_configuration_loading(self):
        """Test configuration loading."""
        from src.sanger_pipeline.utils.helpers import load_config
        
        config_file = Path(__file__).parent.parent.parent / "config" / "default_config.yaml"
        if config_file.exists():
            config = load_config(config_file)
            
            # Check aDNA-specific settings
            assert 'bootstrap' in config
            assert 'damage' in config
            assert 'quality' in config
        else:
            pytest.skip("Configuration file not found")
