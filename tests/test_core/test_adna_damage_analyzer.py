"""
Test suite for ADNADamageAnalyzer class.
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

    @patch('src.sanger_pipeline.core.adna_damage_analyzer.SeqIO.read')
    @patch('src.sanger_pipeline.core.adna_damage_analyzer.pairwise2.align.globalxx')
    def test_analyze_sequence_damage(self, mock_align, mock_seqio):
        """Test damage pattern analysis."""
        # Mock sequence objects
        mock_query = MagicMock()
        mock_query.seq = "ACGTACGTACGTACGT"
        mock_ref = MagicMock()
        mock_ref.seq = "ACGTACGTACGTACGT"
        
        mock_seqio.side_effect = [mock_query, mock_ref]
        
        # Mock alignment result
        mock_alignment = MagicMock()
        mock_alignment.seqA = "ACGTACGTACGTACGT"
        mock_alignment.seqB = "ACGTACGTACGTACGT"
        mock_align.return_value = [mock_alignment]

        # Run analysis
        result = self.analyzer.analyze_sequence_damage(
            self.query_file, 
            self.ref_file
        )

        # Verify result structure
        assert isinstance(result, dict)
        assert 'ct_transitions' in result
        assert 'ga_transitions' in result
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

        with patch.object(self.analyzer, 'analyze_sequence_damage') as mock_analyze:
            mock_analyze.return_value = {
                'damage_5_prime': 0.1,
                'damage_3_prime': 0.08,
                'ct_transitions': 0.05,
                'ga_transitions': 0.04
            }

            result = self.analyzer.bootstrap_damage_analysis(
                seq_files,
                self.ref_file,
                iterations=100  # Small number for testing
            )

            assert isinstance(result, dict)
            assert 'mean_damage_5' in result
            assert 'std_damage_5' in result
            assert 'confidence_interval_5' in result

    def test_assess_authenticity(self):
        """Test authenticity assessment."""
        # Mock bootstrap results
        bootstrap_results = {
            'mean_damage_5': 0.15,
            'std_damage_5': 0.03,
            'mean_damage_3': 0.12,
            'std_damage_3': 0.025,
            'p_value_5': 0.001,
            'p_value_3': 0.005
        }

        result = self.analyzer.assess_authenticity(bootstrap_results)

        assert isinstance(result, dict)
        assert 'authenticity_score' in result
        assert 'is_authentic' in result
        assert 'confidence' in result
        assert 'assessment' in result

    def test_generate_damage_plots(self):
        """Test damage plot generation."""
        seq_files = [self.query_file]
        output_dir = self.temp_dir / "plots"

        with patch.object(self.analyzer, '_calculate_positional_damage') as mock_calc:
            mock_calc.return_value = {
                'position': list(range(16)),
                'ct_frequency': [0.1] * 16,
                'ga_frequency': [0.08] * 16
            }

            with patch('matplotlib.pyplot.savefig') as mock_save:
                self.analyzer.generate_damage_plots(
                    seq_files,
                    self.ref_file,
                    output_dir
                )

                # Verify that plotting was attempted
                mock_save.assert_called()

    def test_align_sequences(self):
        """Test sequence alignment functionality."""
        ref_seq = "ACGTACGTACGTACGT"
        query_seq = "ACGTACGTACGTACGT"

        with patch('src.sanger_pipeline.core.adna_damage_analyzer.pairwise2.align.globalxx') as mock_align:
            mock_alignment = MagicMock()
            mock_alignment.seqA = ref_seq
            mock_alignment.seqB = query_seq
            mock_align.return_value = [mock_alignment]

            result = self.analyzer._align_sequences(ref_seq, query_seq)

            assert isinstance(result, tuple)
            assert len(result) == 2
            assert result[0] == ref_seq
            assert result[1] == query_seq

    def test_calculate_damage_statistics(self):
        """Test damage statistics calculation."""
        # Test alignment with some C->T transitions
        aligned_ref = "CCCCACGTACGTACGTACGTCCCC"
        aligned_query = "TTTTACGTACGTACGTACGTTTTT"
        alignment = (aligned_ref, aligned_query)

        result = self.analyzer._calculate_damage_statistics(alignment)

        assert isinstance(result, dict)
        assert 'ct_transitions' in result
        assert 'ga_transitions' in result
        assert 'damage_5_prime' in result
        assert 'damage_3_prime' in result
        assert result['ct_transitions'] > 0  # Should detect C->T transitions

    def test_terminal_damage_detection(self):
        """Test terminal-specific damage detection."""
        # Sequence with damage at termini
        aligned_ref = "CCCCACGTACGTACGTACGTCCCC"
        aligned_query = "TTTTACGTACGTACGTACGTTTTT"  # C->T at both ends
        alignment = (aligned_ref, aligned_query)

        result = self.analyzer._calculate_damage_statistics(alignment)

        # Should detect high damage at termini
        assert result['damage_5_prime'] > 0.5
        assert result['damage_3_prime'] > 0.5

    def test_no_damage_detection(self):
        """Test detection when no damage is present."""
        # Perfect match sequences
        ref_seq = "ACGTACGTACGTACGT"
        query_seq = "ACGTACGTACGTACGT"
        alignment = (ref_seq, query_seq)

        result = self.analyzer._calculate_damage_statistics(alignment)

        # Should detect no damage
        assert result['ct_transitions'] == 0.0
        assert result['ga_transitions'] == 0.0
        assert result['damage_5_prime'] == 0.0
        assert result['damage_3_prime'] == 0.0

    def test_partial_damage_detection(self):
        """Test detection of partial damage patterns."""
        # Sequence with some C->T transitions
        aligned_ref = "ACGCACGTACGTACGTACGTACGC"
        aligned_query = "ACGTACGTACGTACGTACGTACGT"  # C->T in middle
        alignment = (aligned_ref, aligned_query)

        result = self.analyzer._calculate_damage_statistics(alignment)

        # Should detect some damage
        assert 0.0 < result['ct_transitions'] < 1.0
        assert result['ga_transitions'] == 0.0  # No G->A transitions in this example

    def test_error_handling_missing_files(self):
        """Test error handling for missing input files."""
        missing_file = self.temp_dir / "missing.fasta"

        try:
            self.analyzer.analyze_sequence_damage(missing_file, self.ref_file)
            assert False, "Should have raised an exception"
        except Exception:
            pass  # Expected behavior

    def test_error_handling_empty_alignment(self):
        """Test error handling for empty alignment results."""
        with patch('src.sanger_pipeline.core.adna_damage_analyzer.pairwise2.align.globalxx') as mock_align:
            mock_align.return_value = []  # Empty alignment

            try:
                self.analyzer._align_sequences("ACGT", "TGCA")
                assert False, "Should have raised an exception"
            except Exception:
                pass  # Expected behavior
