#!/usr/bin/env python3
"""
Diagnose primer detection in AB1 files.

Usage:
    python -m sanger_pipeline.scripts.diagnose_primers
"""

from pathlib import Path
from ..core.enhanced_ab1_converter import EnhancedAB1Converter
from Bio import SeqIO

def diagnose_primer_detection():
    """Diagnose why primer detection isn't working."""
    
    # Find AB1 files
    input_dir = Path("input")
    ab1_files = list(input_dir.glob("*HVS*-F.ab1"))[:3]  # Test 3 forward files
    
    converter = EnhancedAB1Converter()
    
    print("🔍 Primer Detection Diagnosis\n")
    
    for ab1_file in ab1_files:
        print(f"📁 File: {ab1_file.name}")
        
        # Read AB1 file
        try:
            record = SeqIO.read(ab1_file, "abi")
            sequence = str(record.seq)
            
            print(f"   Length: {len(sequence)} bp")
            print(f"   First 60bp: {sequence[:60]}")
            print(f"   Last 60bp:  {sequence[-60:]}")
            
            # Check each HVS region
            for hvs_region, primers in converter.primers.items():
                forward_primer = primers['forward']
                reverse_primer = primers['reverse']
                
                print(f"\n   🧬 Checking {hvs_region}:")
                print(f"      Forward primer: {forward_primer}")
                print(f"      Reverse primer: {reverse_primer}")
                
                # Check forward primer
                if forward_primer in sequence.upper():
                    print(f"      ✅ Forward primer found (exact match)")
                else:
                    # Check for partial matches
                    best_match = 0
                    best_pos = -1
                    for i in range(min(len(forward_primer) + 10, len(sequence))):
                        window = sequence[i:i + len(forward_primer)].upper()
                        if len(window) == len(forward_primer):
                            matches = sum(1 for a, b in zip(window, forward_primer) if a == b)
                            if matches > best_match:
                                best_match = matches
                                best_pos = i
                    
                    match_percent = (best_match / len(forward_primer)) * 100
                    print(f"      🔍 Best forward match: {best_match}/{len(forward_primer)} ({match_percent:.1f}%) at position {best_pos}")
                    if best_pos >= 0:
                        window = sequence[best_pos:best_pos + len(forward_primer)].upper()
                        print(f"      Sequence: {window}")
                        print(f"      Expected: {forward_primer}")
                
                # Check reverse primer
                if reverse_primer in sequence.upper():
                    print(f"      ✅ Reverse primer found (exact match)")
                else:
                    # Check for partial matches from the end
                    best_match = 0
                    best_pos = -1
                    for i in range(min(len(reverse_primer) + 10, len(sequence))):
                        start_pos = len(sequence) - len(reverse_primer) - i
                        if start_pos >= 0:
                            window = sequence[start_pos:start_pos + len(reverse_primer)].upper()
                            if len(window) == len(reverse_primer):
                                matches = sum(1 for a, b in zip(window, reverse_primer) if a == b)
                                if matches > best_match:
                                    best_match = matches
                                    best_pos = start_pos
                    
                    match_percent = (best_match / len(reverse_primer)) * 100
                    print(f"      🔍 Best reverse match: {best_match}/{len(reverse_primer)} ({match_percent:.1f}%) at position {best_pos}")
                    if best_pos >= 0:
                        window = sequence[best_pos:best_pos + len(reverse_primer)].upper()
                        print(f"      Sequence: {window}")
                        print(f"      Expected: {reverse_primer}")
            
            # Test detection function
            detected = converter.detect_hvs_region(sequence)
            print(f"\n   📍 Auto-detected region: {detected}")
            
            print("\n" + "="*80 + "\n")
            
        except Exception as e:
            print(f"   ❌ Error: {e}\n")

if __name__ == "__main__":
    diagnose_primer_detection()
