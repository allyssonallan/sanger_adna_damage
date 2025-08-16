"""
Damage Visualization Utilities for Ancient DNA Analysis.

This module provides comprehensive plotting functions for visualizing
ancient DNA damage patterns, including position-specific plots,
summary charts, and dashboard visualizations.
"""

import json
import logging
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
import base64
from io import BytesIO

logger = logging.getLogger(__name__)

# Set matplotlib style for publication-quality plots
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3


class DamagePlotGenerator:
    """Generate comprehensive damage visualization plots."""
    
    def __init__(self, output_dir: Path):
        """
        Initialize damage plot generator.
        
        Args:
            output_dir: Pipeline output directory
        """
        self.output_dir = Path(output_dir)
        self.damage_dir = self.output_dir / "damage_analysis"
        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True)
        
        # Color scheme for plots
        self.colors = {
            'damage': '#e74c3c',      # Red for damage
            'control': '#3498db',     # Blue for control
            'significant': '#f39c12', # Orange for significant
            'background': '#95a5a6',  # Gray for background
            'high_damage': '#c0392b', # Dark red
            'low_damage': '#f1c40f'   # Yellow
        }
    
    def generate_comprehensive_damage_plots(self) -> List[str]:
        """
        Generate all damage plots for the report.
        
        Returns:
            List of base64-encoded plot images for embedding in HTML
        """
        plot_images = []
        
        # Collect all damage data
        damage_data = self._collect_damage_data()
        
        if not damage_data:
            logger.warning("No damage data found for plotting")
            return []
        
        # 1. Damage Score Distribution
        dist_plot = self._create_damage_distribution_plot(damage_data)
        if dist_plot:
            plot_images.append(dist_plot)
        
        # 2. 5' vs 3' Damage Correlation
        correlation_plot = self._create_damage_correlation_plot(damage_data)
        if correlation_plot:
            plot_images.append(correlation_plot)
        
        # 3. Sample Status Summary
        status_plot = self._create_status_summary_plot(damage_data)
        if status_plot:
            plot_images.append(status_plot)
        
        # 4. Quality vs Damage Relationship
        quality_plot = self._create_quality_damage_plot(damage_data)
        if quality_plot:
            plot_images.append(quality_plot)
        
        return plot_images

    def get_dashboard_data(self, base_output_dir):
        """Extract comprehensive dashboard data from damage analysis results"""
        try:
            damage_dir = os.path.join(base_output_dir, 'damage_analysis')
            if not os.path.exists(damage_dir):
                return {'has_data': False, 'summary': {}, 'samples': []}
            
            samples = []
            for json_file in glob.glob(os.path.join(damage_dir, '*_damage_results.json')):
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    
                    sample_id = os.path.basename(json_file).replace('_damage_results.json', '')
                    
                    # Extract comprehensive sample information with ALL JSON data
                    damage_patterns = data.get('damage_patterns', {})
                    bootstrap_analysis = data.get('bootstrap_analysis', {})
                    damage_assessment = data.get('damage_assessment', {})
                    
                    sample_info = {
                        'sample_id': sample_id,
                        'total_bases': damage_patterns.get('total_bases', 0),
                        'valid_bases': damage_patterns.get('valid_bases', 0),
                        'n_content': damage_patterns.get('n_content', 0),
                        'ambiguous_content': damage_patterns.get('ambiguous_content', 0),
                        'damage_5_prime': damage_patterns.get('damage_5_prime', 0),
                        'damage_3_prime': damage_patterns.get('damage_3_prime', 0),
                        'overall_damage_rate': damage_patterns.get('overall_damage_rate', 0),
                        'total_ct_transitions': damage_patterns.get('total_ct_transitions', 0),
                        'total_ga_transitions': damage_patterns.get('total_ga_transitions', 0),
                        'sequence_quality': damage_patterns.get('sequence_quality', {}),
                        'p_value_5_prime': bootstrap_analysis.get('p_value_5_prime', 0),
                        'p_value_3_prime': bootstrap_analysis.get('p_value_3_prime', 0),
                        'bootstrap_mean_5_prime': bootstrap_analysis.get('bootstrap_mean_5_prime', 0),
                        'bootstrap_mean_3_prime': bootstrap_analysis.get('bootstrap_mean_3_prime', 0),
                        'bootstrap_std_5_prime': bootstrap_analysis.get('bootstrap_std_5_prime', 0),
                        'bootstrap_std_3_prime': bootstrap_analysis.get('bootstrap_std_3_prime', 0),
                        'damage_status': damage_assessment.get('status', 'UNKNOWN'),
                        'interpretation': damage_assessment.get('interpretation', 'No interpretation available'),
                        'confidence': damage_assessment.get('confidence', 0.0),
                        '5_prime_damage_indicated': damage_assessment.get('5_prime_damage_indicated', False),
                        '3_prime_damage_indicated': damage_assessment.get('3_prime_damage_indicated', False),
                        'damage_patterns': damage_patterns,
                        'bootstrap_analysis': bootstrap_analysis
                    }
                    
                    samples.append(sample_info)
                    
                except (json.JSONDecodeError, FileNotFoundError) as e:
                    print(f"Warning: Could not process {json_file}: {e}")
                    continue
            
            if not samples:
                return {'has_data': False, 'summary': {}, 'samples': []}
            
            # Calculate summary statistics
            total_samples = len(samples)
            samples_with_damage = sum(1 for s in samples if s['damage_status'] == 'ANCIENT_DNA_CONFIRMED')
            samples_partial_damage = sum(1 for s in samples if s['damage_status'] == 'PARTIAL_DAMAGE_SIGNATURE')
            samples_no_damage = sum(1 for s in samples if s['damage_status'] == 'NO_DAMAGE_DETECTED')
            samples_insufficient = sum(1 for s in samples if s['damage_status'] == 'INSUFFICIENT_DATA')
            
            summary = {
                'total_samples': total_samples,
                'samples_with_damage': samples_with_damage,
                'samples_partial_damage': samples_partial_damage,
                'samples_no_damage': samples_no_damage,
                'samples_insufficient': samples_insufficient,
                'damage_indication_rate': ((samples_with_damage + samples_partial_damage) / total_samples * 100) if total_samples > 0 else 0,
                'high_quality_samples': sum(1 for s in samples if s['total_bases'] > 0 and (s['valid_bases'] / s['total_bases'] * 100) >= 80),
                'total_bases': sum(s['total_bases'] for s in samples),
                'total_valid_bases': sum(s['valid_bases'] for s in samples),
                'total_n_content': sum(s['n_content'] for s in samples)
            }
            
            return {
                'has_data': True,
                'summary': summary,
                'samples': samples
            }
            
        except Exception as e:
            print(f"Error processing dashboard data: {e}")
            return {'has_data': False, 'summary': {}, 'samples': []}

    def generate_individual_sample_plots(self, base_output_dir):
        """Generate individual damage plots for each sample"""
        try:
            damage_dir = os.path.join(base_output_dir, 'damage_analysis')
            if not os.path.exists(damage_dir):
                return {}
            
            sample_plots = {}
            for json_file in glob.glob(os.path.join(damage_dir, '*_damage_results.json')):
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    
                    sample_id = os.path.basename(json_file).replace('_damage_results.json', '')
                    
                    # Create individual sample plot
                    fig, ax = plt.subplots(figsize=(8, 6))
                    
                    # Get damage data
                    damage_patterns = data.get('damage_patterns', {})
                    damage_5_prime = damage_patterns.get('damage_5_prime', 0)
                    damage_3_prime = damage_patterns.get('damage_3_prime', 0)
                    total_bases = damage_patterns.get('total_bases', 0)
                    valid_bases = damage_patterns.get('valid_bases', 0)
                    
                    # Create bar plot for this sample
                    categories = ['5\' Damage', '3\' Damage', 'Valid %', 'N Content %']
                    n_content = damage_patterns.get('n_content', 0)
                    valid_pct = (valid_bases / total_bases * 100) if total_bases > 0 else 0
                    n_pct = (n_content / total_bases * 100) if total_bases > 0 else 0
                    
                    values = [damage_5_prime * 100, damage_3_prime * 100, valid_pct, n_pct]
                    colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#f9ca24']
                    
                    bars = ax.bar(categories, values, color=colors, alpha=0.8)
                    
                    # Add value labels on bars
                    for bar, value in zip(bars, values):
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                               f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')
                    
                    ax.set_title(f'Sample Analysis: {sample_id}', fontsize=14, fontweight='bold', pad=20)
                    ax.set_ylabel('Percentage (%)', fontsize=12)
                    ax.set_ylim(0, max(100, max(values) * 1.2))
                    ax.grid(True, alpha=0.3)
                    
                    # Add sample info as text
                    info_text = f'Total Bases: {total_bases:,}\\nValid Bases: {valid_bases:,}\\nStatus: {data.get("damage_assessment", {}).get("status", "UNKNOWN")}'
                    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
                           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
                    
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    
                    # Convert to base64
                    buffer = BytesIO()
                    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                    buffer.seek(0)
                    plot_data = base64.b64encode(buffer.getvalue()).decode()
                    buffer.close()
                    plt.close()
                    
                    sample_plots[sample_id] = plot_data
                    
                except (json.JSONDecodeError, FileNotFoundError) as e:
                    logger.warning(f"Could not process {json_file}: {e}")
                    continue
            
            return sample_plots
            
        except Exception as e:
            logger.error(f"Error generating individual sample plots: {e}")
            return {}
    
    def _collect_damage_data(self) -> List[Dict[str, Any]]:
        """Collect damage data from all JSON files."""
        damage_data = []
        
        if not self.damage_dir.exists():
            return damage_data
        
        for json_file in self.damage_dir.glob("*_damage_results.json"):
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                # Extract sample information
                sample_name = json_file.stem.replace('_damage_results', '')
                damage_patterns = data.get('damage_patterns', {})
                bootstrap_analysis = data.get('bootstrap_analysis', {})
                damage_assessment = data.get('damage_assessment', {})
                
                sample_data = {
                    'sample_name': sample_name,
                    'damage_5_prime': damage_patterns.get('damage_5_prime', 0),
                    'damage_3_prime': damage_patterns.get('damage_3_prime', 0),
                    'overall_damage_rate': damage_patterns.get('overall_damage_rate', 0),
                    'total_bases': damage_patterns.get('total_bases', 0),
                    'valid_bases': damage_patterns.get('valid_bases', 0),
                    'n_content': damage_patterns.get('n_content', 0),
                    'n_percentage': damage_patterns.get('sequence_quality', {}).get('n_percentage', 0),
                    'valid_percentage': damage_patterns.get('sequence_quality', {}).get('valid_percentage', 100),
                    'p_value_5_prime': bootstrap_analysis.get('p_value_5_prime', 1.0),
                    'p_value_3_prime': bootstrap_analysis.get('p_value_3_prime', 1.0),
                    'status': damage_assessment.get('status', 'UNKNOWN'),
                    'interpretation': damage_assessment.get('interpretation', 'No assessment'),
                    'ct_transitions': damage_patterns.get('total_ct_transitions', 0),
                    'ga_transitions': damage_patterns.get('total_ga_transitions', 0)
                }
                
                damage_data.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error reading damage file {json_file}: {e}")
        
        return damage_data
    
    def _create_damage_distribution_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """Create damage score distribution plot."""
        if not damage_data:
            return None
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Extract damage rates
        damage_5_prime = [d['damage_5_prime'] for d in damage_data]
        damage_3_prime = [d['damage_3_prime'] for d in damage_data]
        
        # 5' damage distribution
        ax1.hist(damage_5_prime, bins=15, alpha=0.7, color=self.colors['damage'], 
                label="5' End Damage", edgecolor='black')
        ax1.axvline(np.mean(damage_5_prime), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(damage_5_prime):.3f}')
        ax1.set_xlabel('Damage Rate')
        ax1.set_ylabel('Frequency')
        ax1.set_title("5' End Damage Distribution")
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 3' damage distribution
        ax2.hist(damage_3_prime, bins=15, alpha=0.7, color=self.colors['control'], 
                label="3' End Damage", edgecolor='black')
        ax2.axvline(np.mean(damage_3_prime), color='blue', linestyle='--', 
                   label=f'Mean: {np.mean(damage_3_prime):.3f}')
        ax2.set_xlabel('Damage Rate')
        ax2.set_ylabel('Frequency')
        ax2.set_title("3' End Damage Distribution")
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle('Ancient DNA Damage Rate Distributions', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        return self._plot_to_base64(fig)
    
    def _create_damage_correlation_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """Create 5' vs 3' damage correlation plot."""
        if not damage_data:
            return None
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        # Extract damage rates and statuses
        damage_5_prime = [d['damage_5_prime'] for d in damage_data]
        damage_3_prime = [d['damage_3_prime'] for d in damage_data]
        statuses = [d['status'] for d in damage_data]
        sample_names = [d['sample_name'] for d in damage_data]
        
        # Color mapping for statuses
        color_map = {
            'DAMAGE_INDICATED': self.colors['high_damage'],
            'PARTIAL_DAMAGE_SIGNATURE': self.colors['significant'],
            'NO_DAMAGE_SIGNATURE': self.colors['background'],
            'UNKNOWN': self.colors['control']
        }
        
        colors = [color_map.get(status, self.colors['control']) for status in statuses]
        
        # Create scatter plot
        ax.scatter(damage_5_prime, damage_3_prime, c=colors, 
                  alpha=0.7, s=100, edgecolors='black', linewidth=0.5)
        
        # Add sample names as annotations (for first few samples to avoid clutter)
        for i, (x, y, name) in enumerate(zip(damage_5_prime[:10], damage_3_prime[:10], sample_names[:10])):
            ax.annotate(name, (x, y), xytext=(5, 5), textcoords='offset points', 
                       fontsize=8, alpha=0.8)
        
        # Add correlation line with error handling
        if len(damage_5_prime) > 1:
            try:
                correlation = np.corrcoef(damage_5_prime, damage_3_prime)[0, 1]
                if not np.isnan(correlation) and len(set(damage_5_prime)) > 1 and len(set(damage_3_prime)) > 1:
                    z = np.polyfit(damage_5_prime, damage_3_prime, 1)
                    p = np.poly1d(z)
                    ax.plot(damage_5_prime, p(damage_5_prime), "r--", alpha=0.8, 
                           label=f'Correlation: {correlation:.3f}')
                else:
                    ax.text(0.5, 0.95, f'Correlation: {correlation:.3f}' if not np.isnan(correlation) else 'Correlation: N/A', 
                           transform=ax.transAxes, ha='center', va='top', fontsize=10, alpha=0.7)
            except (np.linalg.LinAlgError, ValueError, Warning) as e:
                logger.warning(f"Could not compute correlation: {e}")
                ax.text(0.5, 0.95, 'Correlation: N/A', transform=ax.transAxes, 
                       ha='center', va='top', fontsize=10, alpha=0.7)
        
        # Create legend
        legend_elements = [mpatches.Patch(color=color, label=status.replace('_', ' ').title()) 
                          for status, color in color_map.items() if status in statuses]
        ax.legend(handles=legend_elements, loc='upper right')
        
        ax.set_xlabel("5' End Damage Rate")
        ax.set_ylabel("3' End Damage Rate")
        ax.set_title("5' vs 3' End Damage Correlation")
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return self._plot_to_base64(fig)
    
    def _create_status_summary_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """Create damage status summary pie chart."""
        if not damage_data:
            return None
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Status distribution
        status_counts = {}
        for d in damage_data:
            status = d['status'].replace('_', ' ').title()
            status_counts[status] = status_counts.get(status, 0) + 1
        
        if status_counts:
            colors = ['#e74c3c', '#f39c12', '#95a5a6', '#3498db'][:len(status_counts)]
            wedges, texts, autotexts = ax1.pie(status_counts.values(), labels=status_counts.keys(), 
                                              autopct='%1.1f%%', colors=colors, startangle=90)
            ax1.set_title('Damage Assessment Status Distribution')
        
        # Quality vs Damage bar chart
        quality_bins = ['<25%', '25-50%', '50-75%', '>75%']
        damage_status_by_quality = {status: [0, 0, 0, 0] for status in ['High', 'Partial', 'None']}
        
        for d in damage_data:
            valid_pct = d['valid_percentage']
            if valid_pct < 25:
                bin_idx = 0
            elif valid_pct < 50:
                bin_idx = 1
            elif valid_pct < 75:
                bin_idx = 2
            else:
                bin_idx = 3
            
            status = d['status']
            if status == 'DAMAGE_INDICATED':
                damage_status_by_quality['High'][bin_idx] += 1
            elif status == 'PARTIAL_DAMAGE_SIGNATURE':
                damage_status_by_quality['Partial'][bin_idx] += 1
            else:
                damage_status_by_quality['None'][bin_idx] += 1
        
        x = np.arange(len(quality_bins))
        width = 0.25
        
        ax2.bar(x - width, damage_status_by_quality['High'], width, label='High Damage', 
               color='#e74c3c', alpha=0.8)
        ax2.bar(x, damage_status_by_quality['Partial'], width, label='Partial Damage', 
               color='#f39c12', alpha=0.8)
        ax2.bar(x + width, damage_status_by_quality['None'], width, label='No Damage', 
               color='#95a5a6', alpha=0.8)
        
        ax2.set_xlabel('Valid Sequence Percentage')
        ax2.set_ylabel('Number of Samples')
        ax2.set_title('Damage Status by Sequence Quality')
        ax2.set_xticks(x)
        ax2.set_xticklabels(quality_bins)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return self._plot_to_base64(fig)
    
    def _create_quality_damage_plot(self, damage_data: List[Dict]) -> Optional[str]:
        """Create quality vs damage relationship plot."""
        if not damage_data:
            return None
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Extract data
        valid_percentages = [d['valid_percentage'] for d in damage_data]
        n_percentages = [d['n_percentage'] for d in damage_data]
        total_bases = [d['total_bases'] for d in damage_data]
        overall_damage = [d['overall_damage_rate'] for d in damage_data]
        
        # 1. Valid percentage vs Overall damage
        ax1.scatter(valid_percentages, overall_damage, alpha=0.7, color=self.colors['damage'])
        ax1.set_xlabel('Valid Sequence Percentage (%)')
        ax1.set_ylabel('Overall Damage Rate')
        ax1.set_title('Sequence Quality vs Damage Rate')
        ax1.grid(True, alpha=0.3)
        
        # 2. N content distribution
        ax2.hist(n_percentages, bins=15, alpha=0.7, color=self.colors['control'], 
                edgecolor='black')
        ax2.axvline(np.mean(n_percentages), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(n_percentages):.1f}%')
        ax2.set_xlabel('N Content Percentage (%)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('N Content Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Total bases vs Damage
        ax3.scatter(total_bases, overall_damage, alpha=0.7, color=self.colors['significant'])
        ax3.set_xlabel('Total Bases')
        ax3.set_ylabel('Overall Damage Rate')
        ax3.set_title('Sequence Length vs Damage Rate')
        ax3.grid(True, alpha=0.3)
        
        # 4. Sample summary statistics
        stats_data = {
            'Total Samples': len(damage_data),
            'High Damage': sum(1 for d in damage_data if d['status'] == 'DAMAGE_INDICATED'),
            'Partial Damage': sum(1 for d in damage_data if d['status'] == 'PARTIAL_DAMAGE_SIGNATURE'),
            'No Damage': sum(1 for d in damage_data if d['status'] == 'NO_DAMAGE_SIGNATURE'),
            'High Quality\n(>50% valid)': sum(1 for d in damage_data if d['valid_percentage'] > 50)
        }
        
        bars = ax4.bar(range(len(stats_data)), list(stats_data.values()), 
                      color=[self.colors['control'], self.colors['high_damage'], 
                            self.colors['significant'], self.colors['background'], 
                            self.colors['damage']][:len(stats_data)])
        ax4.set_xticks(range(len(stats_data)))
        ax4.set_xticklabels(list(stats_data.keys()), rotation=45, ha='right')
        ax4.set_ylabel('Count')
        ax4.set_title('Summary Statistics')
        ax4.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, stats_data.values()):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{value}', ha='center', va='bottom')
        
        plt.suptitle('Comprehensive Quality and Damage Analysis', fontsize=14, fontweight='bold')
        plt.tight_layout()
        return self._plot_to_base64(fig)
    
    def _plot_to_base64(self, fig) -> str:
        """Convert matplotlib figure to base64 string for HTML embedding."""
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        buffer.seek(0)
        plot_data = buffer.getvalue()
        buffer.close()
        plt.close(fig)
        
        return base64.b64encode(plot_data).decode('utf-8')
    
    def generate_dashboard_data(self) -> Dict[str, Any]:
        """Generate comprehensive dashboard data for enhanced reporting."""
        damage_data = self._collect_damage_data()
        
        if not damage_data:
            return {'samples': [], 'summary': {}, 'has_data': False}
        
        # Process sample data for dashboard
        dashboard_samples = []
        for d in damage_data:
            sample_info = {
                'sample_name': d['sample_name'],
                'total_bases': d['total_bases'],
                'valid_bases': d['valid_bases'],
                'n_content': d['n_content'],
                'n_percentage': round(d['n_percentage'], 2),
                'valid_percentage': round(d['valid_percentage'], 2),
                'damage_5_prime': round(d['damage_5_prime'], 4),
                'damage_3_prime': round(d['damage_3_prime'], 4),
                'overall_damage_rate': round(d['overall_damage_rate'], 4),
                'status': d['status'],
                'status_display': d['status'].replace('_', ' ').title(),
                'interpretation': d['interpretation'],
                'ct_transitions': d['ct_transitions'],
                'ga_transitions': d['ga_transitions'],
                'p_value_5_prime': round(d['p_value_5_prime'], 4),
                'p_value_3_prime': round(d['p_value_3_prime'], 4),
                'quality_tier': self._get_quality_tier(d['valid_percentage']),
                'damage_tier': self._get_damage_tier(d['overall_damage_rate']),
                'confidence_level': self._get_confidence_level(d['p_value_5_prime'], d['p_value_3_prime'])
            }
            dashboard_samples.append(sample_info)
        
        # Generate summary statistics
        summary = {
            'total_samples': len(damage_data),
            'samples_with_damage': sum(1 for d in damage_data if d['status'] == 'DAMAGE_INDICATED'),
            'samples_partial_damage': sum(1 for d in damage_data if d['status'] == 'PARTIAL_DAMAGE_SIGNATURE'),
            'samples_no_damage': sum(1 for d in damage_data if d['status'] == 'NO_DAMAGE_SIGNATURE'),
            'high_quality_samples': sum(1 for d in damage_data if d['valid_percentage'] > 50),
            'mean_valid_percentage': round(np.mean([d['valid_percentage'] for d in damage_data]), 2),
            'mean_n_percentage': round(np.mean([d['n_percentage'] for d in damage_data]), 2),
            'mean_damage_rate': round(np.mean([d['overall_damage_rate'] for d in damage_data]), 4),
            'total_bases_analyzed': sum(d['total_bases'] for d in damage_data),
            'total_valid_bases': sum(d['valid_bases'] for d in damage_data),
            'damage_indication_rate': round((sum(1 for d in damage_data if d['status'] == 'DAMAGE_INDICATED') / len(damage_data)) * 100, 1)
        }
        
        return {
            'samples': dashboard_samples,
            'summary': summary,
            'has_data': True
        }
    
    def _get_quality_tier(self, valid_percentage: float) -> str:
        """Categorize sequence quality."""
        if valid_percentage >= 75:
            return 'excellent'
        elif valid_percentage >= 50:
            return 'good'
        elif valid_percentage >= 25:
            return 'fair'
        else:
            return 'poor'
    
    def _get_damage_tier(self, damage_rate: float) -> str:
        """Categorize damage level."""
        if damage_rate >= 0.2:
            return 'high'
        elif damage_rate >= 0.1:
            return 'moderate'
        elif damage_rate >= 0.05:
            return 'low'
        else:
            return 'minimal'
    
    def _get_confidence_level(self, p_value_5: float, p_value_3: float) -> str:
        """Determine confidence level based on p-values."""
        min_p = min(p_value_5, p_value_3)
        if min_p < 0.001:
            return 'very_high'
        elif min_p < 0.01:
            return 'high'
        elif min_p < 0.05:
            return 'moderate'
        else:
            return 'low'
