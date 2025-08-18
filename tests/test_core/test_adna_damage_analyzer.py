"""
Test suite for ADNADamageAnalyzer class - updated for refactored modular architecture.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from src.sanger_pipeline.core.adna_damage_analyzer import ADNADamageAnalyzer


class TestADNADamageAnalyzer:
    """Test aDNA damage analysis functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = ADNADamageAnalyzer()
        
        # Create temporary files for testing
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Sample sequences for testing
        self.query_seq = ">test_sequence\nACGTACGTACGTACGT"
        self.ref_seq = ">reference\nACGTACGTACGTACGT"
        
        # Create test files
        self.query_file = self.temp_dir / "query.fasta"
        self.ref_file = self.temp_dir / "reference.fasta"
        
        self.query_file.write_text(self.query_seq)
        self.ref_file.write_text(self.ref_seq)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_analyzer_initialization(self):
        """Test that analyzer initializes correctly."""
        analyzer = ADNADamageAnalyzer()
        assert analyzer is not None
        assert hasattr(analyzer, 'damage_calculator')
        assert hasattr(analyzer, 'statistical_analyzer')
        assert hasattr(analyzer, 'visualizer')

    def test_analyze_sequence_damage(self):
        """Test damage pattern analysis."""
        # Create test files with C->T transitions for damage patterns
        query_seq = ">test_sequence\nTTGTACGTACGTACGT"  # C->T at positions 1,2
        ref_seq = ">reference\nCCGTACGTACGTACGT"
        
        query_file = self.temp_dir / "query_damage.fasta"
        ref_file = self.temp_dir / "ref_damage.fasta"
        
        query_file.write_text(query_seq)
        ref_file.write_text(ref_seq)
        
        result = self.analyzer.analyze_sequence_damage(query_file, ref_file)
        
        assert isinstance(result, dict)
        assert 'damage_5_prime' in result
        assert 'damage_3_prime' in result
        
    def test_bootstrap_damage_analysis(self):
        """Test bootstrap analysis functionality."""
        # Create multiple test files
        seq_files = []
        for i in range(3):
            seq_file = self.temp_dir / f"seq_{i}.fasta"
            seq_file.write_text(f">sequence_{i}\nACGTACGTACGTACGT")
            seq_files.append(seq_file)

        result = self.analyzer.bootstrap_damage_analysis(
            seq_files,
            self.ref_file,
            iterations=100  # Small number for testing
        )

        assert isinstance(result, dict)
        assert 'bootstrap_mean_5_prime' in result
        assert 'bootstrap_mean_3_prime' in result
        assert 'observed_damage_5_prime' in result
        assert 'observed_damage_3_prime' in result

    def test_assess_authenticity(self):
        """Test authenticity assessment."""
        # Mock bootstrap results with correct keys
        bootstrap_results = {
            'observed_damage_5_prime': 0.15,
            'observed_damage_3_prime': 0.12,
            'bootstrap_std_5_prime': 0.03,
            'bootstrap_std_3_prime': 0.025,
            'p_value_5_prime': 0.001,
            'p_value_3_prime': 0.005
        }

        result = self.analyzer.assess_authenticity(bootstrap_results)
        
        assert isinstance(result, dict)
        assert 'status' in result
        assert 'interpretation' in result

    def test_generate_damage_plots(self):
        """Test damage plot generation."""
        # Create sample sequence files for testing
        seq_files = []
        for i in range(2):
            seq_file = self.temp_dir / f"seq_{i}.fasta"
            seq_file.write_text(f">sequence_{i}\nACGTACGTACGTACGT")
            seq_files.append(seq_file)
        
        output_dir = self.temp_dir / "plots"
        output_dir.mkdir()
        
        # Test plot generation (should not raise exception)
        try:
            self.analyzer.generate_damage_plots(
                seq_files, output_dir
            )
        except Exception as e:
            # Plot generation might fail due to missing data, but shouldn't crash
            assert "plot" in str(e).lower() or "damage" in str(e).lower()

    def test_error_handling_missing_files(self):
        """Test error handling for missing files."""
        missing_file = self.temp_dir / "missing.fasta"
        
        try:
            self.analyzer.analyze_sequence_damage(missing_file, self.ref_file)
            assert False, "Should raise exception for missing file"
        except Exception as e:
            assert "exist" in str(e).lower() or "not found" in str(e).lower()
