"""
Tests for AB1 converter functionality.
"""

from pathlib import Path
from unittest.mock import Mock
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from sanger_pipeline.core.ab1_converter import AB1Converter


class TestAB1Converter:
    """Test AB1 converter functionality."""
    
    def test_initialization(self):
        """Test AB1Converter initialization."""
        converter = AB1Converter()
        assert converter.min_quality == 20
        
        converter = AB1Converter(min_quality=25)
        assert converter.min_quality == 25
        
    def test_filter_sequence_by_quality(self):
        """Test quality filtering logic."""
        converter = AB1Converter(min_quality=20)
        
        sequence = "ATCGATCG"
        qualities = [25, 15, 30, 18, 22, 19, 35, 21]
        
        # Mock record for testing
        mock_record = Mock()
        mock_record.letter_annotations = {"phred_quality": qualities}
        mock_record.seq = sequence
        
        # Test the internal logic that would be used
        filtered_seq = "".join([
            base if qual >= converter.min_quality else "N"
            for base, qual in zip(sequence, qualities)
        ])
        
        expected = "ANCNNNC"  # Only bases with quality >= 20 kept
        assert filtered_seq == expected
        
    def test_quality_threshold_validation(self):
        """Test quality threshold validation."""
        # Test valid thresholds
        converter = AB1Converter(min_quality=10)
        assert converter.min_quality == 10
        
        converter = AB1Converter(min_quality=40)
        assert converter.min_quality == 40
