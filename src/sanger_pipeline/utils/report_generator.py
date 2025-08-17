"""
Comprehensive QC Report Generator for Sanger Pipeline.

This module generates beautiful HTML reports with analysis summaries,
interactive tabs, damage analysis integration, and comprehensive
damage visualization plots.
"""

import json
import logging
import base64
from datetime import datetime
from pathlib import Path
from typing import Dict, Any
import pandas as pd

# Import our new damage plotting utilities
from .damage_plots import DamagePlotGenerator

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
        
        # Initialize damage plot generator
        self.damage_plotter = DamagePlotGenerator(output_dir)
        
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
    
    def _encode_logos(self) -> Dict[str, str]:
        """Encode institutional logos as base64 strings."""
        logos = {}
        config_dir = Path(__file__).parent.parent.parent.parent / "config" / "logos"
        
        logo_files = {
            'funcap': 'funcap.png',
            'labbat': 'labbat.png', 
            'npdm': 'npdm.png',
            'ufc': 'ufc.png'
        }
        
        for logo_key, filename in logo_files.items():
            logo_path = config_dir / filename
            if logo_path.exists():
                try:
                    with open(logo_path, 'rb') as f:
                        logo_data = f.read()
                    logos[logo_key] = base64.b64encode(logo_data).decode('utf-8')
                except Exception as e:
                    logger.warning(f"Could not encode logo {filename}: {e}")
            else:
                logger.warning(f"Logo file not found: {logo_path}")
        
        return logos
    
    def collect_pipeline_statistics(self) -> Dict[str, Any]:
        """Collect comprehensive statistics from all pipeline outputs."""
        stats = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'directories': {},
            'samples': {},
            'damage_analysis': {},
            'hvs_combinations': {},
            'quality_metrics': {},
            'dashboard_data': {},
            'damage_plots': []
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
        
        # Generate comprehensive dashboard data
        # Generate dashboard data
        stats['dashboard_data'] = self.damage_plotter.get_dashboard_data(str(self.output_dir))
        
        # Generate damage data summary
        dashboard_data = stats['dashboard_data']
        if dashboard_data.get('has_data'):
            samples = dashboard_data.get('samples', [])
            stats['damage_data'] = {
                'files_analyzed': len(samples),
                'summary': dashboard_data.get('summary', {})
            }
        else:
            stats['damage_data'] = {'files_analyzed': 0, 'summary': {}}
        
        # Generate damage plots for embedding in report
        stats['damage_plots'] = self.damage_plotter.generate_comprehensive_damage_plots()
        
        # Generate individual sample plots
        stats['individual_sample_plots'] = self.damage_plotter.generate_individual_sample_plots(str(self.output_dir))
        
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
        
        # Encode logos
        logos = self._encode_logos()
        
        # Generate logo HTML
        logo_html_top = ""
        logo_html_bottom = ""
        
        if logos:
            # Top row: UFC and FUNCAP
            if 'ufc' in logos:
                logo_html_top += f'<a href="https://ufc.br" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["ufc"]}" alt="UFC Logo" title="Universidade Federal do Ceará - https://ufc.br"></a>'
            if 'funcap' in logos:
                logo_html_top += f'<a href="https://funcap.ce.gov.br" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["funcap"]}" alt="FUNCAP Logo" title="Fundação Cearense de Apoio ao Desenvolvimento Científico e Tecnológico - https://www.funcap.ce.gov.br"></a>'
            
            # Bottom row: LABBAT and NPDM  
            if 'labbat' in logos:
                logo_html_bottom += f'<a href="https://instagram.com/labbat.npdm.ufc" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["labbat"]}" alt="LABBAT Logo" title="Laboratório de Bioarqueologia Translacional - https://instagram.com/labbat.npdm.ufc"></a>'
            if 'npdm' in logos:
                logo_html_bottom += f'<a href="https://npdm.ufc.br" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["npdm"]}" alt="NPDM Logo" title="Núcleo de Pesquisa e Desenvolvimento de Medicamentos - https://npdm.ufc.br"></a>'
        
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
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .header-content {{
            text-align: center;
            flex: 1;
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
        
        .header-logos {{
            display: flex;
            flex-direction: column;
            gap: 10px;
            align-items: center;
        }}
        
        .logo-row {{
            display: flex;
            gap: 15px;
            align-items: center;
        }}
        
        .header-logos img {{
            height: 40px;
            max-width: 80px;
            object-fit: contain;
            background: rgba(255, 255, 255, 0.9);
            padding: 5px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        
        .header-logos a {{
            display: inline-block;
            text-decoration: none;
        }}
        
        .header-logos a:hover img {{
            transform: scale(1.05);
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        }}
        
        @media (max-width: 768px) {{
            .header {{
                flex-direction: column;
                gap: 1rem;
            }}
            
            .header-logos img {{
                height: 30px;
                max-width: 60px;
            }}
            
            .logo-row {{
                gap: 10px;
            }}
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
            <div class="header-content">
                <h1><i class="fas fa-dna"></i> Sanger Pipeline QC Report</h1>
                <div class="timestamp">
                    <i class="fas fa-clock"></i> Generated on {stats['timestamp']}
                </div>
            </div>
            <div class="header-logos">
                <div class="logo-row">
                    {logo_html_top}
                </div>
                <div class="logo-row">
                    {logo_html_bottom}
                </div>
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
                <button class="nav-link" id="sample-details-tab" data-bs-toggle="tab" data-bs-target="#sample-details" type="button" role="tab">
                    <i class="fas fa-table"></i> Sample Details
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
            {self._generate_sample_details_tab(stats)}
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
        """Generate the enhanced damage analysis tab content with plots and comprehensive data."""
        damage_data = stats['damage_analysis']
        dashboard_data = stats['dashboard_data']
        damage_plots = stats['damage_plots']
        
        # Check if we have damage data and dashboard data
        damage_data = stats.get('damage_data', {})
        dashboard_data = stats.get('dashboard_data', {'has_data': False, 'summary': {}, 'samples': []})
        
        if not damage_data.get('files_analyzed', 0) and not dashboard_data.get('has_data', False):
            return """
            <div class="tab-pane fade" id="damage" role="tabpanel">
                <h4><i class="fas fa-radiation"></i> aDNA Damage Analysis</h4>
                <div class="alert alert-info">
                    <i class="fas fa-info-circle"></i> No damage analysis results found.
                </div>
            </div>
            """
        
        # Get pre-computed dashboard data from stats
        dashboard_data = stats.get('dashboard_data', {'has_data': False, 'summary': {}, 'samples': []})
        dashboard_samples = dashboard_data.get('samples', [])
        
        # Calculate summary statistics from samples
        total_samples = len(dashboard_samples)
        samples_with_damage = sum(1 for s in dashboard_samples if s.get('damage_status') == 'ANCIENT_DNA_CONFIRMED')
        samples_partial_damage = sum(1 for s in dashboard_samples if s.get('damage_status') == 'PARTIAL_DAMAGE_SIGNATURE')
        samples_no_damage = sum(1 for s in dashboard_samples if s.get('damage_status') == 'NO_DAMAGE_DETECTED')
        damage_indication_rate = ((samples_with_damage + samples_partial_damage) / total_samples * 100) if total_samples > 0 else 0
        high_quality_samples = sum(1 for s in dashboard_samples if (s.get('valid_bases', 0) / s.get('total_bases', 1) * 100) >= 80)
        
        # Generate summary statistics section
        summary_html = f"""
        <div class="row mb-4">
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-primary">{total_samples}</div>
                    <div class="text-muted">Total Samples</div>
                </div>
            </div>
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-success">{samples_with_damage}</div>
                    <div class="text-muted">High Damage</div>
                </div>
            </div>
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-warning">{samples_partial_damage}</div>
                    <div class="text-muted">Partial Damage</div>
                </div>
            </div>
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-secondary">{samples_no_damage}</div>
                    <div class="text-muted">No Damage</div>
                </div>
            </div>
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-info">{damage_indication_rate:.1f}%</div>
                    <div class="text-muted">Damage Rate</div>
                </div>
            </div>
            <div class="col-md-2">
                <div class="stat-card text-center">
                    <div class="stat-number text-dark">{high_quality_samples}</div>
                    <div class="text-muted">High Quality</div>
                </div>
            </div>
        </div>
        """
        
        # Generate comprehensive sample table
        sample_table_html = ""
        if dashboard_samples:
            sample_table_html = """
            <div class="table-responsive mb-4">
                <table class="table table-striped table-hover" id="damageTable">
                    <thead class="table-primary">
                        <tr>
                            <th>Sample</th>
                            <th>Total Bases</th>
                            <th>Valid Bases</th>
                            <th>N Content</th>
                            <th>Valid %</th>
                            <th>5' Damage</th>
                            <th>3' Damage</th>
                            <th>Overall Damage</th>
                            <th>Status</th>
                            <th>Interpretation</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for sample in dashboard_samples:
                # Color coding for status
                status_colors = {
                    'ANCIENT_DNA_CONFIRMED': 'success',
                    'PARTIAL_DAMAGE_SIGNATURE': 'warning',
                    'NO_DAMAGE_DETECTED': 'secondary',
                    'INSUFFICIENT_DATA': 'danger'
                }
                # Use correct field name from dashboard data
                sample_status = sample.get('damage_status', 'UNKNOWN')
                status_color = status_colors.get(sample_status, 'info')
                
                # Calculate quality tier based on valid percentage
                total_bases = sample.get('total_bases', 0)
                valid_bases = sample.get('valid_bases', 0)
                valid_percentage = (valid_bases / total_bases * 100) if total_bases > 0 else 0
                
                if valid_percentage >= 80:
                    sample_quality = 'excellent'
                elif valid_percentage >= 60:
                    sample_quality = 'good'
                elif valid_percentage >= 40:
                    sample_quality = 'fair'
                else:
                    sample_quality = 'poor'
                
                # Quality tier badge
                quality_colors = {
                    'excellent': 'success',
                    'good': 'info',
                    'fair': 'warning',
                    'poor': 'danger'
                }
                quality_color = quality_colors.get(sample_quality, 'secondary')
                
                # Format status for display
                status_display = sample_status.replace('_', ' ').title()
                
                sample_table_html += f"""
                <tr>
                    <td><strong>{sample['sample_id']}</strong></td>
                    <td>{sample['total_bases']:,}</td>
                    <td>{sample['valid_bases']:,}</td>
                    <td>{sample['n_content']:,}</td>
                    <td><span class="badge bg-{quality_color}">{valid_percentage:.1f}%</span></td>
                    <td>{sample['damage_5_prime']:.4f}</td>
                    <td>{sample['damage_3_prime']:.4f}</td>
                    <td>{sample['overall_damage_rate']:.4f}</td>
                    <td><span class="badge bg-{status_color}">{status_display}</span></td>
                    <td><small>{sample['interpretation'][:100]}...</small></td>
                </tr>
                """
            
            sample_table_html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Generate damage plots section
        plots_html = ""
        if damage_plots:
            plots_html = """
            <div class="row mb-4">
                <div class="col-12">
                    <h5><i class="fas fa-chart-line"></i> Damage Analysis Visualizations</h5>
                    <div class="row">
            """
            
            plot_titles = [
                "Damage Rate Distributions",
                "5' vs 3' Damage Correlation", 
                "Status Summary & Quality Analysis",
                "Comprehensive Quality-Damage Analysis"
            ]
            
            for i, (plot_data, title) in enumerate(zip(damage_plots, plot_titles)):
                plots_html += f"""
                    <div class="col-md-6 mb-3">
                        <div class="card">
                            <div class="card-header">
                                <h6 class="mb-0">{title}</h6>
                            </div>
                            <div class="card-body text-center">
                                <img src="data:image/png;base64,{plot_data}" 
                                     class="img-fluid" alt="{title}" style="max-height: 400px;">
                            </div>
                        </div>
                    </div>
                """
                
                # Break into rows of 2
                if (i + 1) % 2 == 0 and i < len(damage_plots) - 1:
                    plots_html += """
                    </div>
                    <div class="row">
                    """
            
            plots_html += """
                    </div>
                </div>
            </div>
            """
        
        return f"""
        <div class="tab-pane fade" id="damage" role="tabpanel">
            <h4><i class="fas fa-radiation"></i> aDNA Damage Analysis</h4>
            
            <!-- Important Disclaimer -->
            <div class="alert alert-warning">
                <i class="fas fa-exclamation-triangle"></i>
                <strong>Important:</strong> This tool provides damage pattern screening for ancient DNA research. 
                Results should be validated with additional methods for definitive authentication.
            </div>
            
            {summary_html}
            
            {plots_html}
            
            <h5><i class="fas fa-table"></i> Detailed Sample Analysis</h5>
            {sample_table_html}
        </div>
        """
    
    def _generate_sample_details_tab(self, stats):
        """Generate comprehensive sample details tab with enhanced sample information"""
        # Get dashboard data from pre-computed stats
        dashboard_data = stats.get('dashboard_data', {'has_data': False, 'summary': {}, 'samples': []})
        
        if not dashboard_data.get('has_data', False):
            return """
            <div class="tab-pane fade" id="sample-details" role="tabpanel">
                <h4><i class="fas fa-table"></i> Sample Details</h4>
                <div class="alert alert-info">
                    <i class="fas fa-info-circle"></i> No sample data available for detailed analysis.
                </div>
            </div>
            """
        
        # Get dashboard sample data
        dashboard_samples = dashboard_data.get('samples', [])
        individual_plots = stats.get('individual_sample_plots', {})
        
        # Create summary statistics cards
        total_samples = len(dashboard_samples)
        total_bases = sum(sample.get('total_bases', 0) for sample in dashboard_samples)
        total_valid_bases = sum(sample.get('valid_bases', 0) for sample in dashboard_samples)
        total_n_content = sum(sample.get('n_content', 0) for sample in dashboard_samples)
        
        # Count samples by status
        status_counts = {}
        for sample in dashboard_samples:
            status = sample.get('damage_status', 'UNKNOWN')
            status_counts[status] = status_counts.get(status, 0) + 1
        
        # Generate summary cards HTML
        summary_cards_html = f"""
        <div class="row mb-4">
            <div class="col-md-3">
                <div class="card text-white bg-primary">
                    <div class="card-body">
                        <h5 class="card-title"><i class="fas fa-flask"></i> Total Samples</h5>
                        <h2 class="card-text">{total_samples:,}</h2>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-white bg-info">
                    <div class="card-body">
                        <h5 class="card-title"><i class="fas fa-dna"></i> Total Bases</h5>
                        <h2 class="card-text">{total_bases:,}</h2>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-white bg-success">
                    <div class="card-body">
                        <h5 class="card-title"><i class="fas fa-check-circle"></i> Valid Bases</h5>
                        <h2 class="card-text">{total_valid_bases:,}</h2>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-white bg-warning">
                    <div class="card-body">
                        <h5 class="card-title"><i class="fas fa-question-circle"></i> N Content</h5>
                        <h2 class="card-text">{total_n_content:,}</h2>
                    </div>
                </div>
            </div>
        </div>
        """
        
        # Generate status distribution cards
        status_cards_html = ""
        if status_counts:
            status_cards_html = '<div class="row mb-4">'
            for status, count in status_counts.items():
                badge_class = {
                    'ANCIENT_DNA_CONFIRMED': 'bg-success',
                    'PARTIAL_DAMAGE_SIGNATURE': 'bg-warning', 
                    'NO_DAMAGE_DETECTED': 'bg-danger',
                    'INSUFFICIENT_DATA': 'bg-secondary',
                    'UNKNOWN': 'bg-secondary'
                }.get(status, 'bg-secondary')
                
                percentage = (count / total_samples * 100) if total_samples > 0 else 0
                
                status_cards_html += f"""
                <div class="col-md-6 col-lg-3">
                    <div class="card {badge_class} text-white">
                        <div class="card-body">
                            <h6 class="card-title">{status.replace('_', ' ').title()}</h6>
                            <h3 class="card-text">{count} <small>({percentage:.1f}%)</small></h3>
                        </div>
                    </div>
                </div>
                """
            status_cards_html += '</div>'
        
        # Generate comprehensive sample table with ALL JSON data
        sample_table_html = """
        <div class="card">
            <div class="card-header">
                <h5><i class="fas fa-table"></i> Comprehensive Sample Analysis (All Samples)</h5>
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-striped table-hover table-sm">
                        <thead class="table-dark">
                            <tr>
                                <th>Sample ID</th>
                                <th>Total Bases</th>
                                <th>Valid Bases</th>
                                <th>N Content</th>
                                <th>Valid %</th>
                                <th>N %</th>
                                <th>5' Damage</th>
                                <th>3' Damage</th>
                                <th>Overall Damage Rate</th>
                                <th>CT Transitions</th>
                                <th>GA Transitions</th>
                                <th>P-value 5'</th>
                                <th>P-value 3'</th>
                                <th>Bootstrap Mean 5'</th>
                                <th>Bootstrap Mean 3'</th>
                                <th>5' Indicated</th>
                                <th>3' Indicated</th>
                                <th>Status</th>
                                <th>Quality</th>
                                <th>Interpretation</th>
                            </tr>
                        </thead>
                        <tbody>
        """
        
        for sample in dashboard_samples:
            sample_id = sample.get('sample_id', 'Unknown')
            total_bases = sample.get('total_bases', 0)
            valid_bases = sample.get('valid_bases', 0)
            n_content = sample.get('n_content', 0)
            damage_5_prime = sample.get('damage_5_prime', 0)
            damage_3_prime = sample.get('damage_3_prime', 0)
            overall_damage_rate = sample.get('overall_damage_rate', 0)
            ct_transitions = sample.get('total_ct_transitions', 0)
            ga_transitions = sample.get('total_ga_transitions', 0)
            p_value_5 = sample.get('p_value_5_prime', 0)
            p_value_3 = sample.get('p_value_3_prime', 0)
            bootstrap_mean_5 = sample.get('bootstrap_mean_5_prime', 0)
            bootstrap_mean_3 = sample.get('bootstrap_mean_3_prime', 0)
            damage_5_indicated = sample.get('5_prime_damage_indicated', False)
            damage_3_indicated = sample.get('3_prime_damage_indicated', False)
            damage_status = sample.get('damage_status', 'UNKNOWN')
            interpretation = sample.get('interpretation', 'No interpretation available')
            
            # Calculate percentages
            valid_percent = (valid_bases / total_bases * 100) if total_bases > 0 else 0
            n_percent = (n_content / total_bases * 100) if total_bases > 0 else 0
            
            # Determine quality tier
            if valid_percent >= 80:
                quality_tier = "HIGH"
                quality_class = "badge bg-success"
            elif valid_percent >= 50:
                quality_tier = "MEDIUM"
                quality_class = "badge bg-warning"
            else:
                quality_tier = "LOW"
                quality_class = "badge bg-danger"
            
            # Status badge styling
            status_class = {
                'ANCIENT_DNA_CONFIRMED': 'badge bg-success',
                'PARTIAL_DAMAGE_SIGNATURE': 'badge bg-warning',
                'NO_DAMAGE_DETECTED': 'badge bg-danger',
                'INSUFFICIENT_DATA': 'badge bg-secondary',
                'UNKNOWN': 'badge bg-secondary'
            }.get(damage_status, 'badge bg-secondary')
            
            # Damage indication badges
            damage_5_badge = '<span class="badge bg-success">Yes</span>' if damage_5_indicated else '<span class="badge bg-secondary">No</span>'
            damage_3_badge = '<span class="badge bg-success">Yes</span>' if damage_3_indicated else '<span class="badge bg-secondary">No</span>'
            
            sample_table_html += f"""
            <tr>
                <td><strong>{sample_id}</strong></td>
                <td>{total_bases:,}</td>
                <td>{valid_bases:,}</td>
                <td>{n_content:,}</td>
                <td>{valid_percent:.1f}%</td>
                <td>{n_percent:.1f}%</td>
                <td>{damage_5_prime:.3f}</td>
                <td>{damage_3_prime:.3f}</td>
                <td>{overall_damage_rate:.3f}</td>
                <td>{ct_transitions}</td>
                <td>{ga_transitions}</td>
                <td>{p_value_5:.3f}</td>
                <td>{p_value_3:.3f}</td>
                <td>{bootstrap_mean_5:.3f}</td>
                <td>{bootstrap_mean_3:.3f}</td>
                <td>{damage_5_badge}</td>
                <td>{damage_3_badge}</td>
                <td><span class="{status_class}">{damage_status.replace('_', ' ')}</span></td>
                <td><span class="{quality_class}">{quality_tier}</span></td>
                <td><small class="text-muted">{interpretation[:150]}{'...' if len(interpretation) > 150 else ''}</small></td>
            </tr>
            """
        
        sample_table_html += """
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        """
        
        # Generate individual sample plots section
        individual_plots_html = ""
        if individual_plots:
            individual_plots_html = """
            <div class="card mt-4">
                <div class="card-header">
                    <h5><i class="fas fa-chart-bar"></i> Individual Sample Analysis Plots</h5>
                    <p class="text-muted mb-0">Click on any sample to view its detailed damage analysis plot</p>
                </div>
                <div class="card-body">
                    <div class="accordion" id="samplePlotsAccordion">
            """
            
            for i, sample in enumerate(dashboard_samples):  # Show ALL samples
                sample_id = sample.get('sample_id', 'Unknown')
                if sample_id in individual_plots:
                    plot_data = individual_plots[sample_id]
                    
                    # Create status badge
                    status = sample.get('damage_status', 'UNKNOWN')
                    status_class = {
                        'ANCIENT_DNA_CONFIRMED': 'badge bg-success',
                        'PARTIAL_DAMAGE_SIGNATURE': 'badge bg-warning',
                        'NO_DAMAGE_DETECTED': 'badge bg-danger',
                        'INSUFFICIENT_DATA': 'badge bg-secondary',
                        'UNKNOWN': 'badge bg-secondary'
                    }.get(status, 'badge bg-secondary')
                    
                    # Additional data for display
                    damage_5 = sample.get('damage_5_prime', 0)
                    damage_3 = sample.get('damage_3_prime', 0)
                    ct_transitions = sample.get('total_ct_transitions', 0)
                    ga_transitions = sample.get('total_ga_transitions', 0)
                    
                    individual_plots_html += f"""
                        <div class="accordion-item">
                            <h2 class="accordion-header" id="heading{i}">
                                <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" 
                                        data-bs-target="#collapse{i}" aria-expanded="false" aria-controls="collapse{i}">
                                    <strong>{sample_id}</strong> 
                                    <span class="{status_class} ms-2">{status.replace('_', ' ')}</span>
                                    <small class="text-muted ms-2">
                                        ({sample.get('total_bases', 0):,} bases, {sample.get('valid_bases', 0):,} valid)
                                        - 5':{"%.3f" % damage_5}, 3':{"%.3f" % damage_3} - CT:{ct_transitions}, GA:{ga_transitions}
                                    </small>
                                </button>
                            </h2>
                            <div id="collapse{i}" class="accordion-collapse collapse" 
                                 aria-labelledby="heading{i}" data-bs-parent="#samplePlotsAccordion">
                                <div class="accordion-body">
                                    <div class="row">
                                        <div class="col-md-8">
                                            <div class="text-center">
                                                <img src="data:image/png;base64,{plot_data}" class="img-fluid" 
                                                     alt="Sample {sample_id} Analysis" style="max-width: 100%; height: auto;">
                                            </div>
                                        </div>
                                        <div class="col-md-4">
                                            <h6>Detailed Metrics:</h6>
                                            <ul class="list-unstyled small">
                                                <li><strong>5' Damage:</strong> {damage_5:.3f}</li>
                                                <li><strong>3' Damage:</strong> {damage_3:.3f}</li>
                                                <li><strong>Overall Damage Rate:</strong> {sample.get('overall_damage_rate', 0):.3f}</li>
                                                <li><strong>CT Transitions:</strong> {ct_transitions}</li>
                                                <li><strong>GA Transitions:</strong> {ga_transitions}</li>
                                                <li><strong>P-value 5':</strong> {sample.get('p_value_5_prime', 0):.3f}</li>
                                                <li><strong>P-value 3':</strong> {sample.get('p_value_3_prime', 0):.3f}</li>
                                                <li><strong>5' Indicated:</strong> {'Yes' if sample.get('5_prime_damage_indicated') else 'No'}</li>
                                                <li><strong>3' Indicated:</strong> {'Yes' if sample.get('3_prime_damage_indicated') else 'No'}</li>
                                            </ul>
                                        </div>
                                    </div>
                                    <div class="mt-3">
                                        <small class="text-muted">
                                            <strong>Interpretation:</strong> {sample.get('interpretation', 'No interpretation available')}
                                        </small>
                                    </div>
                                </div>
                            </div>
                        </div>
                    """
            
            individual_plots_html += """
                    </div>
                </div>
            </div>
            """
        
        return f"""
        <div class="tab-pane fade" id="sample-details" role="tabpanel">
            <h4><i class="fas fa-table"></i> Sample Details</h4>
            <p class="lead">Comprehensive overview of all analyzed samples with detailed metrics and interpretations.</p>
            
            {summary_cards_html}
            
            {status_cards_html}
            
            {sample_table_html}
            
            {individual_plots_html}
        </div>
        """

    def _get_dashboard_data(self, base_output_dir):
        """Extract comprehensive dashboard data from damage analysis results"""
        if not self.damage_plotter:
            return {'has_data': False, 'summary': {}, 'samples': []}
        
        return self.damage_plotter.get_dashboard_data(base_output_dir)
    
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
