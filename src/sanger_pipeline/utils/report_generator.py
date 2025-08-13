"""
Comprehensive QC Report Generator for Sanger Pipeline.

This module generates beautiful HTML reports with analysis summaries,
interactive tabs, and damage analysis integration.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any
import pandas as pd

logger = logging.getLogger(__name__)


class QCReportGenerator:
    """Generate comprehensive QC reports for Sanger pipeline analysis."""
    
    def __init__(self, output_dir: Path):
        """
        Initialize report generator.
        
        Args:
            output_dir: Pipeline output directory
        """
        self.output_dir = Path(output_dir)
        self.report_dir = self.output_dir / "reports"
        self.report_dir.mkdir(exist_ok=True)
        
        # Define color palette
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72', 
            'success': '#F18F01',
            'warning': '#C73E1D',
            'info': '#6C757D',
            'light': '#F8F9FA',
            'dark': '#343A40'
        }
    
    def collect_pipeline_statistics(self) -> Dict[str, Any]:
        """Collect comprehensive statistics from all pipeline outputs."""
        stats = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'directories': {},
            'samples': {},
            'damage_analysis': {},
            'hvs_combinations': {},
            'quality_metrics': {}
        }
        
        # Analyze each output directory
        directories = ['fasta', 'filtered', 'consensus', 'aligned', 'final', 'damage_analysis', 'plots']
        
        for dir_name in directories:
            dir_path = self.output_dir / dir_name
            if dir_path.exists():
                stats['directories'][dir_name] = self._analyze_directory(dir_path)
            else:
                stats['directories'][dir_name] = {'exists': False, 'file_count': 0}
        
        # Analyze damage analysis results
        if (self.output_dir / "damage_analysis").exists():
            stats['damage_analysis'] = self._analyze_damage_results()
        
        # Analyze HVS combinations
        if (self.output_dir / "final").exists():
            stats['hvs_combinations'] = self._analyze_hvs_combinations()
        
        # Collect sample-level statistics
        stats['samples'] = self._collect_sample_statistics()
        
        return stats
    
    def _analyze_directory(self, dir_path: Path) -> Dict[str, Any]:
        """Analyze contents of a directory."""
        files = list(dir_path.glob('*'))
        file_types = {}
        total_size = 0
        
        for file_path in files:
            if file_path.is_file():
                suffix = file_path.suffix.lower()
                file_types[suffix] = file_types.get(suffix, 0) + 1
                total_size += file_path.stat().st_size
        
        return {
            'exists': True,
            'file_count': len([f for f in files if f.is_file()]),
            'file_types': file_types,
            'total_size_mb': round(total_size / 1024 / 1024, 2),
            'files': [f.name for f in files if f.is_file()][:20]  # Limit for display
        }
    
    def _analyze_damage_results(self) -> Dict[str, Any]:
        """Analyze aDNA damage analysis results."""
        damage_dir = self.output_dir / "damage_analysis"
        json_files = list(damage_dir.glob("*.json"))
        
        if not json_files:
            return {'files_analyzed': 0, 'summary': {}}
        
        damage_data = []
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    damage_patterns = data.get('damage_patterns', {})
                    
                    # Extract key metrics
                    sample_name = json_file.stem.replace('_damage_results', '')
                    damage_data.append({
                        'sample': sample_name,
                        'damage_5_prime': damage_patterns.get('damage_5_prime', 0),
                        'damage_3_prime': damage_patterns.get('damage_3_prime', 0),
                        'overall_damage_rate': damage_patterns.get('overall_damage_rate', 0),
                        'valid_percentage': damage_patterns.get('sequence_quality', {}).get('valid_percentage', 0),
                        'n_percentage': damage_patterns.get('sequence_quality', {}).get('n_percentage', 0)
                    })
            except Exception as e:
                logger.warning(f"Could not parse damage file {json_file}: {e}")
        
        # Calculate summary statistics
        if damage_data:
            df = pd.DataFrame(damage_data)
            summary = {
                'mean_damage_5_prime': df['damage_5_prime'].mean(),
                'mean_damage_3_prime': df['damage_3_prime'].mean(),
                'mean_overall_damage': df['overall_damage_rate'].mean(),
                'mean_valid_percentage': df['valid_percentage'].mean(),
                'samples_with_damage': len(df[df['overall_damage_rate'] > 0.1]),
                'high_quality_samples': len(df[df['valid_percentage'] > 50])
            }
        else:
            summary = {}
        
        return {
            'files_analyzed': len(json_files),
            'summary': summary,
            'individual_results': damage_data[:10]  # Limit for report
        }
    
    def _analyze_hvs_combinations(self) -> Dict[str, Any]:
        """Analyze HVS region combinations in final merged files."""
        final_dir = self.output_dir / "final"
        final_files = list(final_dir.glob("*.fasta"))
        
        combinations = {}
        for file_path in final_files:
            # Extract HVS combination from filename
            name = file_path.stem
            if '_merged' in name:
                # Extract HVS regions from filename
                hvs_regions = []
                if 'HVS1' in name:
                    hvs_regions.append('HVS1')
                if 'HVS2' in name:
                    hvs_regions.append('HVS2')
                if 'HVS3' in name:
                    hvs_regions.append('HVS3')
                
                combo_key = '_'.join(sorted(hvs_regions)) if hvs_regions else 'Unknown'
                combinations[combo_key] = combinations.get(combo_key, 0) + 1
        
        return {
            'total_merged_files': len(final_files),
            'combinations': combinations,
            'files': [f.name for f in final_files]
        }
    
    def _collect_sample_statistics(self) -> Dict[str, Any]:
        """Collect per-sample processing statistics."""
        samples = {}
        
        # Get all samples from consensus directory
        consensus_dir = self.output_dir / "consensus"
        if consensus_dir.exists():
            for consensus_file in consensus_dir.glob("*_consensus.fasta"):
                try:
                    sample_base = consensus_file.stem.replace('_consensus', '')
                    # Remove HVS region suffix if present
                    for hvs in ['_HVS1', '_HVS2', '_HVS3']:
                        if sample_base.endswith(hvs):
                            sample_base = sample_base[:-len(hvs)]
                            break
                    
                    if sample_base not in samples:
                        samples[sample_base] = {
                            'consensus_files': [],
                            'final_files': [],
                            'damage_files': [],
                            'hvs_regions': set()
                        }
                    
                    samples[sample_base]['consensus_files'].append(consensus_file.name)
                    
                    # Determine HVS region
                    for hvs in ['HVS1', 'HVS2', 'HVS3']:
                        if hvs in consensus_file.name:
                            samples[sample_base]['hvs_regions'].add(hvs)
                            break
                    
                except Exception as e:
                    logger.warning(f"Error processing consensus file {consensus_file}: {e}")
        
        # Check for final merged files
        final_dir = self.output_dir / "final"
        if final_dir.exists():
            for final_file in final_dir.glob("*.fasta"):
                # Extract base sample name
                sample_base = final_file.stem.split('_HVS')[0].replace('_merged', '')
                if sample_base in samples:
                    samples[sample_base]['final_files'].append(final_file.name)
        
        # Check for damage analysis files
        damage_dir = self.output_dir / "damage_analysis"
        if damage_dir.exists():
            for damage_file in damage_dir.glob("*.json"):
                sample_base = damage_file.stem.split('_HVS')[0].replace('_damage_results', '')
                if sample_base in samples:
                    samples[sample_base]['damage_files'].append(damage_file.name)
        
        # Convert sets to lists for JSON serialization
        for sample_data in samples.values():
            sample_data['hvs_regions'] = list(sample_data['hvs_regions'])
        
        return samples
    
    def generate_html_report(self, stats: Dict[str, Any]) -> str:
        """Generate beautiful HTML report with tabs."""
        html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sanger Pipeline QC Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {{
            --primary-color: {self.colors['primary']};
            --secondary-color: {self.colors['secondary']};
            --success-color: {self.colors['success']};
            --warning-color: {self.colors['warning']};
        }}
        
        body {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }}
        
        .main-container {{
            background: white;
            margin: 20px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white;
            padding: 2rem;
            text-align: center;
        }}
        
        .header h1 {{
            margin: 0;
            font-size: 2.5rem;
            font-weight: 300;
        }}
        
        .header .timestamp {{
            opacity: 0.9;
            font-size: 1.1rem;
            margin-top: 0.5rem;
        }}
        
        .nav-tabs {{
            border-bottom: 3px solid var(--primary-color);
            background: #f8f9fa;
        }}
        
        .nav-tabs .nav-link {{
            border: none;
            color: var(--primary-color);
            font-weight: 500;
            padding: 1rem 1.5rem;
            margin-right: 0.25rem;
            border-radius: 0;
        }}
        
        .nav-tabs .nav-link.active {{
            background: var(--primary-color);
            color: white;
            border-bottom: 3px solid var(--secondary-color);
        }}
        
        .tab-content {{
            padding: 2rem;
        }}
        
        .stat-card {{
            background: white;
            border-radius: 10px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            transition: transform 0.2s;
        }}
        
        .stat-card:hover {{
            transform: translateY(-2px);
        }}
        
        .stat-number {{
            font-size: 2.5rem;
            font-weight: bold;
            color: var(--primary-color);
        }}
        
        .directory-card {{
            border-left: 4px solid var(--primary-color);
            margin-bottom: 1rem;
        }}
        
        .hvs-badge {{
            background: var(--success-color);
            color: white;
            padding: 0.25rem 0.5rem;
            border-radius: 20px;
            font-size: 0.8rem;
            margin: 0.1rem;
        }}
        
        .progress-ring {{
            transform: rotate(-90deg);
        }}
        
        .chart-container {{
            background: white;
            border-radius: 10px;
            padding: 1rem;
            margin: 1rem 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="header">
            <h1><i class="fas fa-dna"></i> Sanger Pipeline QC Report</h1>
            <div class="timestamp">
                <i class="fas fa-clock"></i> Generated on {stats['timestamp']}
            </div>
        </div>
        
        <ul class="nav nav-tabs" id="reportTabs" role="tablist">
            <li class="nav-item" role="presentation">
                <button class="nav-link active" id="overview-tab" data-bs-toggle="tab" data-bs-target="#overview" type="button" role="tab">
                    <i class="fas fa-chart-pie"></i> Overview
                </button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="directories-tab" data-bs-toggle="tab" data-bs-target="#directories" type="button" role="tab">
                    <i class="fas fa-folder"></i> Directories
                </button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="samples-tab" data-bs-toggle="tab" data-bs-target="#samples" type="button" role="tab">
                    <i class="fas fa-vials"></i> Samples
                </button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="damage-tab" data-bs-toggle="tab" data-bs-target="#damage" type="button" role="tab">
                    <i class="fas fa-radiation"></i> Damage Analysis
                </button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="hvs-tab" data-bs-toggle="tab" data-bs-target="#hvs" type="button" role="tab">
                    <i class="fas fa-project-diagram"></i> HVS Regions
                </button>
            </li>
        </ul>
        
        <div class="tab-content">
            {self._generate_overview_tab(stats)}
            {self._generate_directories_tab(stats)}
            {self._generate_samples_tab(stats)}
            {self._generate_damage_tab(stats)}
            {self._generate_hvs_tab(stats)}
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Initialize charts and interactive elements
        document.addEventListener('DOMContentLoaded', function() {{
            {self._generate_charts_javascript(stats)}
        }});
    </script>
</body>
</html>
        """
        
        return html_template
    
    def _generate_overview_tab(self, stats: Dict[str, Any]) -> str:
        """Generate the overview tab content."""
        total_samples = len(stats['samples'])
        total_final_files = stats['hvs_combinations'].get('total_merged_files', 0)
        damage_files = stats['damage_analysis'].get('files_analyzed', 0)
        
        return f"""
        <div class="tab-pane fade show active" id="overview" role="tabpanel">
            <div class="row">
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{total_samples}</div>
                        <div class="text-muted">Total Samples</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{total_final_files}</div>
                        <div class="text-muted">Merged Files</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{damage_files}</div>
                        <div class="text-muted">Damage Analyses</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{len(stats['hvs_combinations'].get('combinations', {}))}</div>
                        <div class="text-muted">HVS Combinations</div>
                    </div>
                </div>
            </div>
            
            <div class="row mt-4">
                <div class="col-md-6">
                    <div class="chart-container">
                        <h5><i class="fas fa-chart-donut-alt"></i> HVS Combinations Distribution</h5>
                        <canvas id="hvsChart" width="400" height="200"></canvas>
                    </div>
                </div>
                <div class="col-md-6">
                    <div class="chart-container">
                        <h5><i class="fas fa-folder-open"></i> Directory File Counts</h5>
                        <canvas id="directoriesChart" width="400" height="200"></canvas>
                    </div>
                </div>
            </div>
            
            {self._generate_damage_summary_overview(stats)}
        </div>
        """
    
    def _generate_directories_tab(self, stats: Dict[str, Any]) -> str:
        """Generate the directories tab content."""
        directories_html = ""
        
        for dir_name, dir_info in stats['directories'].items():
            if dir_info['exists']:
                status_color = "success" if dir_info['file_count'] > 0 else "warning"
                file_types_html = ""
                
                for file_type, count in dir_info.get('file_types', {}).items():
                    file_types_html += f'<span class="badge bg-secondary me-1">{file_type}: {count}</span>'
                
                directories_html += f"""
                <div class="col-md-6 mb-3">
                    <div class="card directory-card">
                        <div class="card-body">
                            <h5 class="card-title">
                                <i class="fas fa-folder text-{status_color}"></i> {dir_name.title()}
                                <span class="badge bg-{status_color} ms-2">{dir_info['file_count']} files</span>
                            </h5>
                            <p class="card-text">
                                <strong>Size:</strong> {dir_info['total_size_mb']} MB<br>
                                <strong>File Types:</strong> {file_types_html}
                            </p>
                        </div>
                    </div>
                </div>
                """
            else:
                directories_html += f"""
                <div class="col-md-6 mb-3">
                    <div class="card directory-card">
                        <div class="card-body">
                            <h5 class="card-title">
                                <i class="fas fa-folder text-muted"></i> {dir_name.title()}
                                <span class="badge bg-secondary ms-2">Not Found</span>
                            </h5>
                        </div>
                    </div>
                </div>
                """
        
        return f"""
        <div class="tab-pane fade" id="directories" role="tabpanel">
            <h4><i class="fas fa-folder-tree"></i> Pipeline Output Directories</h4>
            <div class="row">
                {directories_html}
            </div>
        </div>
        """
    
    def _generate_samples_tab(self, stats: Dict[str, Any]) -> str:
        """Generate the samples tab content."""
        samples_html = ""
        
        for sample_name, sample_data in stats['samples'].items():
            hvs_badges = ''.join([f'<span class="hvs-badge">{hvs}</span>' for hvs in sample_data['hvs_regions']])
            
            final_status = "✓ Merged" if sample_data['final_files'] else "⚠ Not Merged"
            damage_status = "✓ Analyzed" if sample_data['damage_files'] else "⚠ No Analysis"
            
            samples_html += f"""
            <tr>
                <td><strong>{sample_name}</strong></td>
                <td>{hvs_badges}</td>
                <td>{len(sample_data['consensus_files'])}</td>
                <td>{len(sample_data['final_files'])}</td>
                <td>{len(sample_data['damage_files'])}</td>
                <td>{final_status}</td>
                <td>{damage_status}</td>
            </tr>
            """
        
        return f"""
        <div class="tab-pane fade" id="samples" role="tabpanel">
            <h4><i class="fas fa-vials"></i> Sample Processing Summary</h4>
            <div class="table-responsive">
                <table class="table table-striped table-hover">
                    <thead class="table-primary">
                        <tr>
                            <th>Sample Name</th>
                            <th>HVS Regions</th>
                            <th>Consensus Files</th>
                            <th>Final Files</th>
                            <th>Damage Files</th>
                            <th>Merge Status</th>
                            <th>Damage Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {samples_html}
                    </tbody>
                </table>
            </div>
        </div>
        """
    
    def _generate_damage_tab(self, stats: Dict[str, Any]) -> str:
        """Generate the damage analysis tab content."""
        damage_data = stats['damage_analysis']
        
        if not damage_data.get('files_analyzed', 0):
            return """
            <div class="tab-pane fade" id="damage" role="tabpanel">
                <h4><i class="fas fa-radiation"></i> aDNA Damage Analysis</h4>
                <div class="alert alert-info">
                    <i class="fas fa-info-circle"></i> No damage analysis results found.
                </div>
            </div>
            """
        
        summary = damage_data.get('summary', {})
        individual_results = damage_data.get('individual_results', [])
        
        individual_html = ""
        for result in individual_results[:10]:  # Show top 10
            damage_color = "success" if result['overall_damage_rate'] > 0.1 else "warning"
            individual_html += f"""
            <tr>
                <td>{result['sample']}</td>
                <td>{result['damage_5_prime']:.3f}</td>
                <td>{result['damage_3_prime']:.3f}</td>
                <td><span class="badge bg-{damage_color}">{result['overall_damage_rate']:.3f}</span></td>
                <td>{result['valid_percentage']:.1f}%</td>
            </tr>
            """
        
        return f"""
        <div class="tab-pane fade" id="damage" role="tabpanel">
            <h4><i class="fas fa-radiation"></i> aDNA Damage Analysis</h4>
            
            <div class="row mb-4">
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{summary.get('mean_damage_5_prime', 0):.3f}</div>
                        <div class="text-muted">Mean 5' Damage</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{summary.get('mean_damage_3_prime', 0):.3f}</div>
                        <div class="text-muted">Mean 3' Damage</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{summary.get('samples_with_damage', 0)}</div>
                        <div class="text-muted">Samples with Damage</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{summary.get('high_quality_samples', 0)}</div>
                        <div class="text-muted">High Quality Samples</div>
                    </div>
                </div>
            </div>
            
            <div class="table-responsive">
                <table class="table table-striped">
                    <thead class="table-primary">
                        <tr>
                            <th>Sample</th>
                            <th>5' Damage</th>
                            <th>3' Damage</th>
                            <th>Overall Damage</th>
                            <th>Valid %</th>
                        </tr>
                    </thead>
                    <tbody>
                        {individual_html}
                    </tbody>
                </table>
            </div>
        </div>
        """
    
    def _generate_hvs_tab(self, stats: Dict[str, Any]) -> str:
        """Generate the HVS combinations tab content."""
        combinations = stats['hvs_combinations'].get('combinations', {})
        files = stats['hvs_combinations'].get('files', [])
        
        combinations_html = ""
        for combo, count in sorted(combinations.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(files) * 100) if files else 0
            combinations_html += f"""
            <tr>
                <td><strong>{combo}</strong></td>
                <td>{count}</td>
                <td>
                    <div class="progress">
                        <div class="progress-bar" style="width: {percentage}%">{percentage:.1f}%</div>
                    </div>
                </td>
            </tr>
            """
        
        files_html = ""
        for file_name in files[:20]:  # Show first 20 files
            files_html += f'<li class="list-group-item">{file_name}</li>'
        
        return f"""
        <div class="tab-pane fade" id="hvs" role="tabpanel">
            <h4><i class="fas fa-project-diagram"></i> HVS Region Analysis</h4>
            
            <div class="row">
                <div class="col-md-6">
                    <h5>Region Combinations</h5>
                    <table class="table table-striped">
                        <thead class="table-primary">
                            <tr>
                                <th>HVS Combination</th>
                                <th>Count</th>
                                <th>Percentage</th>
                            </tr>
                        </thead>
                        <tbody>
                            {combinations_html}
                        </tbody>
                    </table>
                </div>
                <div class="col-md-6">
                    <h5>Final Merged Files</h5>
                    <ul class="list-group" style="max-height: 400px; overflow-y: auto;">
                        {files_html}
                    </ul>
                </div>
            </div>
        </div>
        """
    
    def _generate_damage_summary_overview(self, stats: Dict[str, Any]) -> str:
        """Generate damage analysis summary for overview tab."""
        damage_data = stats['damage_analysis']
        
        if not damage_data.get('files_analyzed', 0):
            return ""
        
        summary = damage_data.get('summary', {})
        
        return f"""
        <div class="row mt-4">
            <div class="col-12">
                <div class="chart-container">
                    <h5><i class="fas fa-radiation"></i> Damage Analysis Summary</h5>
                    <div class="row">
                        <div class="col-md-4 text-center">
                            <div class="stat-number text-warning">{summary.get('mean_overall_damage', 0):.3f}</div>
                            <div class="text-muted">Mean Damage Rate</div>
                        </div>
                        <div class="col-md-4 text-center">
                            <div class="stat-number text-success">{summary.get('mean_valid_percentage', 0):.1f}%</div>
                            <div class="text-muted">Mean Valid Sequence</div>
                        </div>
                        <div class="col-md-4 text-center">
                            <div class="stat-number text-info">{summary.get('samples_with_damage', 0)}</div>
                            <div class="text-muted">Samples with Damage</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_charts_javascript(self, stats: Dict[str, Any]) -> str:
        """Generate JavaScript for interactive charts."""
        # HVS combinations chart data
        combinations = stats['hvs_combinations'].get('combinations', {})
        hvs_labels = list(combinations.keys())
        hvs_data = list(combinations.values())
        
        # Directories chart data
        dir_labels = []
        dir_data = []
        for dir_name, dir_info in stats['directories'].items():
            if dir_info['exists']:
                dir_labels.append(dir_name.title())
                dir_data.append(dir_info['file_count'])
        
        return f"""
        // HVS Combinations Pie Chart
        const hvsCtx = document.getElementById('hvsChart');
        if (hvsCtx) {{
            new Chart(hvsCtx, {{
                type: 'doughnut',
                data: {{
                    labels: {hvs_labels},
                    datasets: [{{
                        data: {hvs_data},
                        backgroundColor: ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6C757D'],
                        borderWidth: 2
                    }}]
                }},
                options: {{
                    responsive: true,
                    plugins: {{
                        legend: {{
                            position: 'bottom'
                        }}
                    }}
                }}
            }});
        }}
        
        // Directories Bar Chart
        const dirCtx = document.getElementById('directoriesChart');
        if (dirCtx) {{
            new Chart(dirCtx, {{
                type: 'bar',
                data: {{
                    labels: {dir_labels},
                    datasets: [{{
                        label: 'File Count',
                        data: {dir_data},
                        backgroundColor: '#2E86AB',
                        borderColor: '#1e5f7a',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    scales: {{
                        y: {{
                            beginAtZero: true
                        }}
                    }}
                }}
            }});
        }}
        """
    
    def generate_report(self) -> Path:
        """
        Generate comprehensive QC report.
        
        Returns:
            Path to generated HTML report
        """
        logger.info("Generating comprehensive QC report...")
        
        # Collect all statistics
        stats = self.collect_pipeline_statistics()
        
        # Generate HTML report
        html_content = self.generate_html_report(stats)
        
        # Save report
        report_file = self.report_dir / f"sanger_qc_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logger.info(f"QC report generated: {report_file}")
        
        return report_file
