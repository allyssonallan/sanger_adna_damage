"""
Individual plot generators for damage visualization.

This module contains specific plot generation functions for different
types of damage analysis visualizations.
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import base64
from io import BytesIO

logger = logging.getLogger(__name__)

# Set matplotlib style for publication-quality plots
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3


class PlotGenerator:
    """Generates individual damage visualization plots."""
    
    def __init__(self, plots_dir: Path):
        """
        Initialize plot generator.
        
        Args:
            plots_dir: Directory for saving plots
        """
        self.plots_dir = Path(plots_dir)
        self.plots_dir.mkdir(exist_ok=True)
        
        # Color scheme for plots
        self.colors = {
            'damage': '#e74c3c',      # Red for damage
            'control': '#3498db',     # Blue for control
            'significant': '#f39c12', # Orange for significant
            'background': '#95a5a6',  # Gray for background
            'high_quality': '#27ae60', # Green for high quality
            'low_quality': '#e67e22'  # Orange for low quality
        }
    
    def create_damage_distribution_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """
        Create damage distribution plot.
        
        Args:
            damage_data: List of damage data dictionaries
            
        Returns:
            Base64 encoded plot string or None
        """
        if not damage_data:
            return None
        
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            # Extract damage values
            damage_5_prime = [d['damage_5_prime'] for d in damage_data]
            damage_3_prime = [d['damage_3_prime'] for d in damage_data]
            overall_damage = [d['overall_damage_rate'] for d in damage_data]
            
            # 5' vs 3' damage scatter plot
            scatter = ax1.scatter(damage_5_prime, damage_3_prime, 
                                c=overall_damage, cmap='Reds', alpha=0.7, s=60)
            ax1.set_xlabel("5' Damage Rate (%)")
            ax1.set_ylabel("3' Damage Rate (%)")
            ax1.set_title("5' vs 3' Damage Distribution")
            
            # Add diagonal line for reference
            max_val = max(max(damage_5_prime), max(damage_3_prime))
            ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Equal damage')
            ax1.legend()
            
            # Colorbar for overall damage
            cbar = plt.colorbar(scatter, ax=ax1)
            cbar.set_label('Overall Damage Rate (%)')
            
            # Overall damage histogram
            ax2.hist(overall_damage, bins=15, color=self.colors['damage'], alpha=0.7, edgecolor='black')
            ax2.set_xlabel('Overall Damage Rate (%)')
            ax2.set_ylabel('Number of Samples')
            ax2.set_title('Overall Damage Distribution')
            ax2.axvline(np.mean(overall_damage), color='blue', linestyle='--', 
                       label=f'Mean: {np.mean(overall_damage):.1f}%')
            ax2.legend()
            
            plt.tight_layout()
            
            # Convert to base64
            result = self._plot_to_base64(fig)
            plt.close(fig)
            
            return result
            
        except Exception as e:
            logger.error(f"Error creating damage distribution plot: {e}")
            plt.close('all')
            return None
    
    def create_damage_correlation_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """
        Create damage-quality correlation plot.
        
        Args:
            damage_data: List of damage data dictionaries
            
        Returns:
            Base64 encoded plot string or None
        """
        if not damage_data:
            return None
        
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # Extract values
            overall_damage = [d['overall_damage_rate'] for d in damage_data]
            valid_percentage = [d['valid_percentage'] for d in damage_data]
            sequence_length = [d['sequence_length'] for d in damage_data]
            damage_5_prime = [d['damage_5_prime'] for d in damage_data]
            
            # Damage vs Quality
            ax1.scatter(overall_damage, valid_percentage, alpha=0.7, 
                       color=self.colors['damage'], s=60)
            ax1.set_xlabel('Overall Damage Rate (%)')
            ax1.set_ylabel('Valid Sequence (%)')
            ax1.set_title('Damage vs Quality Correlation')
            
            # Add trend line
            if len(overall_damage) > 1:
                z = np.polyfit(overall_damage, valid_percentage, 1)
                p = np.poly1d(z)
                ax1.plot(sorted(overall_damage), p(sorted(overall_damage)), 
                        "r--", alpha=0.8, label=f'Trend (r={np.corrcoef(overall_damage, valid_percentage)[0,1]:.2f})')
                ax1.legend()
            
            # Damage vs Length
            ax2.scatter(overall_damage, sequence_length, alpha=0.7, 
                       color=self.colors['control'], s=60)
            ax2.set_xlabel('Overall Damage Rate (%)')
            ax2.set_ylabel('Sequence Length (bp)')
            ax2.set_title('Damage vs Sequence Length')
            
            # 5' Damage vs Quality
            ax3.scatter(damage_5_prime, valid_percentage, alpha=0.7, 
                       color=self.colors['significant'], s=60)
            ax3.set_xlabel("5' Damage Rate (%)")
            ax3.set_ylabel('Valid Sequence (%)')
            ax3.set_title("5' Damage vs Quality")
            
            # Quality distribution by damage tier
            high_damage_quality = [d['valid_percentage'] for d in damage_data 
                                 if d['overall_damage_rate'] > 15]
            low_damage_quality = [d['valid_percentage'] for d in damage_data 
                                if d['overall_damage_rate'] <= 15]
            
            ax4.hist([low_damage_quality, high_damage_quality], 
                    bins=15, alpha=0.7, label=['Low Damage (≤15%)', 'High Damage (>15%)'],
                    color=[self.colors['high_quality'], self.colors['low_quality']])
            ax4.set_xlabel('Valid Sequence (%)')
            ax4.set_ylabel('Number of Samples')
            ax4.set_title('Quality Distribution by Damage Level')
            ax4.legend()
            
            plt.tight_layout()
            
            # Convert to base64
            result = self._plot_to_base64(fig)
            plt.close(fig)
            
            return result
            
        except Exception as e:
            logger.error(f"Error creating correlation plot: {e}")
            plt.close('all')
            return None
    
    def create_status_summary_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """
        Create status summary plot.
        
        Args:
            damage_data: List of damage data dictionaries
            
        Returns:
            Base64 encoded plot string or None
        """
        if not damage_data:
            return None
        
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # Quality tiers
            quality_tiers = {}
            for d in damage_data:
                tier = self._get_quality_tier(d['valid_percentage'])
                quality_tiers[tier] = quality_tiers.get(tier, 0) + 1
            
            if quality_tiers:
                colors_q = [self.colors['high_quality'], self.colors['significant'], 
                           self.colors['damage'], self.colors['background']]
                ax1.pie(quality_tiers.values(), labels=quality_tiers.keys(), 
                       autopct='%1.1f%%', colors=colors_q[:len(quality_tiers)])
                ax1.set_title('Quality Distribution')
            
            # Damage tiers
            damage_tiers = {}
            for d in damage_data:
                tier = self._get_damage_tier(d['overall_damage_rate'])
                damage_tiers[tier] = damage_tiers.get(tier, 0) + 1
            
            if damage_tiers:
                colors_d = [self.colors['high_quality'], self.colors['control'], 
                           self.colors['significant'], self.colors['damage'], self.colors['background']]
                ax2.pie(damage_tiers.values(), labels=damage_tiers.keys(), 
                       autopct='%1.1f%%', colors=colors_d[:len(damage_tiers)])
                ax2.set_title('Damage Level Distribution')
            
            # Confidence levels
            confidence_levels = {}
            for d in damage_data:
                level = self._get_confidence_level(d.get('p_value_5_prime', 1.0), 
                                                 d.get('p_value_3_prime', 1.0))
                confidence_levels[level] = confidence_levels.get(level, 0) + 1
            
            if confidence_levels:
                ax3.bar(confidence_levels.keys(), confidence_levels.values(), 
                       color=[self.colors['high_quality'], self.colors['significant'], 
                             self.colors['damage'], self.colors['background']][:len(confidence_levels)])
                ax3.set_title('Statistical Confidence Levels')
                ax3.set_ylabel('Number of Samples')
                plt.setp(ax3.get_xticklabels(), rotation=45)
            
            # Sample count summary
            summary_data = {
                'Total Samples': len(damage_data),
                'High Quality': len([d for d in damage_data if d['valid_percentage'] > 75]),
                'Significant Damage': len([d for d in damage_data if d['overall_damage_rate'] > 10]),
                'High Confidence': len([d for d in damage_data 
                                      if min(d.get('p_value_5_prime', 1.0), 
                                           d.get('p_value_3_prime', 1.0)) < 0.05])
            }
            
            ax4.bar(summary_data.keys(), summary_data.values(), 
                   color=[self.colors['control'], self.colors['high_quality'], 
                         self.colors['damage'], self.colors['significant']])
            ax4.set_title('Analysis Summary')
            ax4.set_ylabel('Count')
            plt.setp(ax4.get_xticklabels(), rotation=45)
            
            plt.tight_layout()
            
            # Convert to base64
            result = self._plot_to_base64(fig)
            plt.close(fig)
            
            return result
            
        except Exception as e:
            logger.error(f"Error creating status summary plot: {e}")
            plt.close('all')
            return None
    
    def create_quality_damage_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """
        Create quality vs damage analysis plot.
        
        Args:
            damage_data: List of damage data dictionaries
            
        Returns:
            Base64 encoded plot string or None
        """
        if not damage_data:
            return None
        
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            
            # Quality vs Overall Damage Scatter with categories
            overall_damage = [d['overall_damage_rate'] for d in damage_data]
            valid_percentage = [d['valid_percentage'] for d in damage_data]
            
            # Color points by quality tier
            colors = []
            for d in damage_data:
                tier = self._get_quality_tier(d['valid_percentage'])
                if tier == 'excellent':
                    colors.append(self.colors['high_quality'])
                elif tier == 'good':
                    colors.append(self.colors['control'])
                elif tier == 'fair':
                    colors.append(self.colors['significant'])
                else:
                    colors.append(self.colors['damage'])
            
            ax1.scatter(overall_damage, valid_percentage, c=colors, alpha=0.7, s=80)
            ax1.set_xlabel('Overall Damage Rate (%)')
            ax1.set_ylabel('Valid Sequence (%)')
            ax1.set_title('Quality vs Damage Analysis')
            
            # Add quality zones
            ax1.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='Excellent (>90%)')
            ax1.axhline(y=75, color='blue', linestyle='--', alpha=0.5, label='Good (>75%)')
            ax1.axhline(y=50, color='orange', linestyle='--', alpha=0.5, label='Fair (>50%)')
            ax1.legend()
            
            # Sequence length distribution
            sequence_lengths = [d['sequence_length'] for d in damage_data]
            
            ax2.hist(sequence_lengths, bins=20, color=self.colors['control'], alpha=0.7, edgecolor='black')
            ax2.set_xlabel('Sequence Length (bp)')
            ax2.set_ylabel('Number of Samples')
            ax2.set_title('Sequence Length Distribution')
            ax2.axvline(np.mean(sequence_lengths), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(sequence_lengths):.0f} bp')
            ax2.axvline(np.median(sequence_lengths), color='blue', linestyle='--', 
                       label=f'Median: {np.median(sequence_lengths):.0f} bp')
            ax2.legend()
            
            plt.tight_layout()
            
            # Convert to base64
            result = self._plot_to_base64(fig)
            plt.close(fig)
            
            return result
            
        except Exception as e:
            logger.error(f"Error creating quality damage plot: {e}")
            plt.close('all')
            return None
    
    def _plot_to_base64(self, fig) -> str:
        """
        Convert matplotlib figure to base64 string.
        
        Args:
            fig: Matplotlib figure
            
        Returns:
            Base64 encoded image string
        """
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight', facecolor='white')
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        
        graphic = base64.b64encode(image_png)
        return graphic.decode('utf-8')
    
    def _get_quality_tier(self, valid_percentage: float) -> str:
        """Classify sequence quality into tiers."""
        if valid_percentage >= 90:
            return 'excellent'
        elif valid_percentage >= 75:
            return 'good'
        elif valid_percentage >= 50:
            return 'fair'
        else:
            return 'poor'
    
    def _get_damage_tier(self, damage_rate: float) -> str:
        """Classify damage rate into tiers."""
        if damage_rate >= 30:
            return 'extreme'
        elif damage_rate >= 20:
            return 'high'
        elif damage_rate >= 10:
            return 'moderate'
        elif damage_rate >= 5:
            return 'low'
        else:
            return 'minimal'
    
    def _get_confidence_level(self, p_value_5: float, p_value_3: float) -> str:
        """Determine statistical confidence level."""
        min_p = min(p_value_5, p_value_3)
        
        if min_p < 0.001:
            return 'very_high'
        elif min_p < 0.01:
            return 'high'
        elif min_p < 0.05:
            return 'moderate'
        else:
            return 'low'
