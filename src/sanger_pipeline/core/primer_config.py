"""
Primer management and configuration system for the Sanger Pipeline.

This module provides functionality for loading, validating, and managing primer
configurations from YAML files, CLI arguments, and default settings.
"""

import logging
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
import re

logger = logging.getLogger(__name__)


class PrimerConfig:
    """
    Manages primer configurations for the Sanger pipeline.
    
    Supports loading primers from:
    - YAML configuration files
    - CLI arguments
    - Default hardcoded values
    """
    
    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize primer configuration.
        
        Args:
            config_path: Optional path to primer configuration YAML file
        """
        self.primers = {}
        self.matching_parameters = self._get_default_matching_parameters()
        self.quality_control = self._get_default_quality_control()
        
        # Convert string to Path if needed
        if isinstance(config_path, str):
            config_path = Path(config_path)
        
        # Load configuration
        if config_path and config_path.exists():
            self.load_from_yaml(config_path)
        else:
            self.load_default_primers()
            if config_path:
                logger.warning(f"Primer config file not found: {config_path}. Using default primers.")
    
    def load_from_yaml(self, config_path: Path) -> None:
        """
        Load primer configuration from YAML file.
        
        Args:
            config_path: Path to YAML configuration file
        """
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            # Load primers
            for region, primer_data in config.items():
                if region in ['matching_parameters', 'quality_control']:
                    continue
                    
                if isinstance(primer_data, dict) and 'forward' in primer_data and 'reverse' in primer_data:
                    self.primers[region] = self._validate_primer_pair(region, primer_data)
            
            # Load matching parameters if present
            if 'matching_parameters' in config:
                self.matching_parameters.update(config['matching_parameters'])
            
            # Load quality control settings if present  
            if 'quality_control' in config:
                self.quality_control.update(config['quality_control'])
                
            logger.info(f"Loaded {len(self.primers)} primer pairs from {config_path}")
            
        except Exception as e:
            logger.error(f"Failed to load primer config from {config_path}: {e}")
            logger.info("Falling back to default primers")
            self.load_default_primers()
    
    def load_default_primers(self) -> None:
        """Load default hardcoded primer sequences."""
        # Updated with more realistic laboratory primers
        default_primers = {
            "HVS1": {
                "forward": "CACCATTAGCACCCAAAGCT",      # L15996 
                "reverse": "TGATTTCACGGAGGATGGTG",       # H16401
                "description": "HVS1 region (16024-16365)"
            },
            "HVS2": {
                "forward": "GGTCTATCACCCTATTAACCAC",    # L48
                "reverse": "CTGTTAAAAGTGCATACCGCCA",     # H408  
                "description": "HVS2 region (57-372)"
            },
            "HVS3": {
                "forward": "CCGCTTCTGGCCACAGCACT",      # L16209
                "reverse": "GGTGATGTGAGCCCGTCTAAAC",     # H599
                "description": "HVS3 region (438-574)"
            }
        }
        
        for region, primer_data in default_primers.items():
            self.primers[region] = self._validate_primer_pair(region, primer_data)
            
        logger.info(f"Loaded {len(self.primers)} default primer pairs")
    
    def add_custom_primers(self, custom_forward: Optional[Dict[str, str]] = None, 
                          custom_reverse: Optional[Dict[str, str]] = None) -> None:
        """
        Add or update primers with custom sequences.
        
        Args:
            custom_forward: Dictionary of {region: forward_sequence}
            custom_reverse: Dictionary of {region: reverse_sequence}
        """
        if custom_forward:
            for region, forward_seq in custom_forward.items():
                if region not in self.primers:
                    self.primers[region] = {}
                self.primers[region]["forward"] = forward_seq.upper().strip()
                logger.info(f"Updated forward primer for {region}")
        
        if custom_reverse:
            for region, reverse_seq in custom_reverse.items():
                if region not in self.primers:
                    self.primers[region] = {}
                self.primers[region]["reverse"] = reverse_seq.upper().strip()
                logger.info(f"Updated reverse primer for {region}")
        
        # Validate updated primers
        for region in self.primers:
            if "forward" in self.primers[region] and "reverse" in self.primers[region]:
                self.primers[region] = self._validate_primer_pair(region, self.primers[region])
    
    def get_primers_for_region(self, region: str) -> Optional[Dict[str, str]]:
        """
        Get primer pair for specific region.
        
        Args:
            region: HVS region name
            
        Returns:
            Dictionary with primer information or None if not found
        """
        return self.primers.get(region)
    
    def get_all_primers(self) -> Dict[str, Dict[str, str]]:
        """Get all loaded primer pairs."""
        return self.primers.copy()
    
    def get_regions(self) -> List[str]:
        """Get list of available region names."""
        return list(self.primers.keys())
    
    def validate_primers(self) -> Tuple[bool, List[str]]:
        """
        Validate all loaded primers.
        
        Returns:
            Tuple of (all_valid, list_of_issues)
        """
        issues = []
        
        for region, primer_data in self.primers.items():
            # Check required fields
            if "forward" not in primer_data or "reverse" not in primer_data:
                issues.append(f"{region}: Missing forward or reverse primer")
                continue
            
            # Validate sequences
            forward_issues = self._validate_sequence(primer_data["forward"], f"{region} forward")
            reverse_issues = self._validate_sequence(primer_data["reverse"], f"{region} reverse")
            
            issues.extend(forward_issues)
            issues.extend(reverse_issues)
        
        return len(issues) == 0, issues
    
    def _validate_primer_pair(self, region: str, primer_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate and clean a primer pair.
        
        Args:
            region: Region name
            primer_data: Raw primer data dictionary
            
        Returns:
            Validated primer data dictionary
        """
        validated = primer_data.copy()
        
        # Clean and validate sequences
        if "forward" in validated:
            validated["forward"] = self._clean_sequence(validated["forward"])
        if "reverse" in validated:  
            validated["reverse"] = self._clean_sequence(validated["reverse"])
        
        # Add reverse complement for matching
        if "reverse" in validated:
            validated["reverse_complement"] = self._reverse_complement(validated["reverse"])
        
        # Add metadata
        validated.setdefault("region", region)
        validated.setdefault("description", f"Primers for {region} region")
        
        return validated
    
    def _validate_sequence(self, sequence: str, name: str) -> List[str]:
        """
        Validate a primer sequence.
        
        Args:
            sequence: DNA sequence to validate
            name: Name for error reporting
            
        Returns:
            List of validation issues
        """
        issues = []
        
        if not sequence:
            issues.append(f"{name}: Empty sequence")
            return issues
        
        # Check length
        min_len = self.quality_control.get("min_primer_length", 15)
        max_len = self.quality_control.get("max_primer_length", 35)
        
        if len(sequence) < min_len:
            issues.append(f"{name}: Sequence too short ({len(sequence)} < {min_len})")
        elif len(sequence) > max_len:
            issues.append(f"{name}: Sequence too long ({len(sequence)} > {max_len})")
        
        # Check for valid nucleotides
        valid_bases = set("ATCGNRYSWKMBDHV")  # Include IUPAC codes
        invalid_bases = set(sequence.upper()) - valid_bases
        if invalid_bases:
            issues.append(f"{name}: Invalid bases found: {', '.join(invalid_bases)}")
        
        # Check for ambiguous bases if warning enabled
        if self.quality_control.get("warn_ambiguous_bases", True):
            ambiguous_bases = set(sequence.upper()) & set("NRYSWKMBDHV")
            if ambiguous_bases:
                issues.append(f"{name}: Contains ambiguous bases: {', '.join(ambiguous_bases)} (may affect primer specificity)")
        
        return issues
    
    def _clean_sequence(self, sequence: str) -> str:
        """Clean and normalize a DNA sequence."""
        if not sequence:
            return ""
        
        # Remove whitespace and convert to uppercase
        cleaned = sequence.strip().upper()
        
        # Remove any non-nucleotide characters except IUPAC codes
        cleaned = re.sub(r'[^ATCGNRYSWKMBDHV]', '', cleaned)
        
        return cleaned
    
    def _reverse_complement(self, sequence: str) -> str:
        """Generate reverse complement of DNA sequence."""
        complement_map = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W',
            'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'
        }
        
        complement = ''.join(complement_map.get(base, base) for base in sequence.upper())
        return complement[::-1]
    
    def _get_default_matching_parameters(self) -> Dict[str, Any]:
        """Get default primer matching parameters."""
        return {
            "similarity_threshold": 0.7,
            "adna_similarity_threshold": 0.4,
            "max_mismatches": 3,
            "search_window": 50
        }
    
    def _get_default_quality_control(self) -> Dict[str, Any]:
        """Get default quality control parameters."""
        return {
            "min_primer_length": 15,
            "max_primer_length": 35,
            "warn_ambiguous_bases": True,
            "check_primer_dimers": False
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Export configuration as dictionary for serialization."""
        return {
            "primers": self.primers,
            "matching_parameters": self.matching_parameters,
            "quality_control": self.quality_control
        }
    
    def save_to_yaml(self, output_path: Path) -> None:
        """
        Save current configuration to YAML file.
        
        Args:
            output_path: Path where to save the configuration
        """
        config_data = self.to_dict()
        
        with open(output_path, 'w') as f:
            yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Saved primer configuration to {output_path}")


def parse_cli_primers(primer_string: str) -> Dict[str, str]:
    """
    Parse primer sequences from CLI string format.
    
    Expected format: "region1:sequence1,region2:sequence2"
    
    Args:
        primer_string: Comma-separated primer definitions
        
    Returns:
        Dictionary of {region: sequence}
    """
    primers = {}
    
    if not primer_string:
        return primers
    
    for primer_def in primer_string.split(','):
        primer_def = primer_def.strip()
        if ':' in primer_def:
            region, sequence = primer_def.split(':', 1)
            primers[region.strip()] = sequence.strip()
        else:
            logger.warning(f"Invalid primer format: {primer_def}. Expected 'region:sequence'")
    
    return primers


def validate_primer_file(config_path: Path) -> Tuple[bool, List[str]]:
    """
    Validate a primer configuration file without loading it.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    try:
        config = PrimerConfig(config_path)
        valid, validation_issues = config.validate_primers()
        issues.extend(validation_issues)
        
        if len(config.primers) == 0:
            issues.append("No valid primer pairs found in configuration")
        
        return valid and len(config.primers) > 0, issues
        
    except Exception as e:
        issues.append(f"Failed to load configuration: {e}")
        return False, issues
