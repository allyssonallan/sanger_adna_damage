"""
HTML template generation for Sanger Pipeline QC reports.

This module handles the generation of HTML content, templates,
and styling for comprehensive QC reports.
"""

import base64
import logging
from pathlib import Path
from typing import Dict, Any

logger = logging.getLogger(__name__)


class HTMLTemplateGenerator:
    """Generates HTML templates and content for QC reports."""
    
    def __init__(self, output_dir: Path):
        """
        Initialize HTML template generator.
        
        Args:
            output_dir: Pipeline output directory
        """
        self.output_dir = Path(output_dir)
        
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
        """
        Encode logo images to base64 for embedding in HTML.
        
        Returns:
            Dictionary mapping logo names to base64 encoded strings
        """
        logos = {}
        logo_dir = Path(__file__).parent.parent.parent.parent / "config" / "logos"
        
        if logo_dir.exists():
            for logo_file in logo_dir.glob("*.png"):
                try:
                    with open(logo_file, "rb") as f:
                        encoded = base64.b64encode(f.read()).decode('utf-8')
                        logos[logo_file.stem.lower()] = encoded
                except Exception as e:
                    logger.warning(f"Could not encode logo {logo_file}: {e}")
        
        return logos
    
    def generate_html_report(self, stats: Dict[str, Any]) -> str:
        """
        Generate complete HTML report.
        
        Args:
            stats: Statistics dictionary from StatisticsCollector
            
        Returns:
            Complete HTML report as string
        """
        # Encode logos
        logos = self._encode_logos()
        
        # Generate logo HTML
        logo_html_top, logo_html_bottom = self._generate_logo_html(logos)
        
        # Generate main HTML template
        html_template = self._generate_main_template(logo_html_top, logo_html_bottom)
        
        # Generate tab content
        tab_content = self._generate_all_tabs(stats)
        
        # Generate JavaScript for charts
        charts_js = self._generate_charts_javascript(stats)
        
        # Combine everything
        return html_template.format(
            tab_content=tab_content,
            charts_javascript=charts_js
        )
    
    def _generate_logo_html(self, logos: Dict[str, str]) -> tuple[str, str]:
        """
        Generate HTML for logos.
        
        Args:
            logos: Dictionary of base64 encoded logos
            
        Returns:
            Tuple of (top_row_html, bottom_row_html)
        """
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
        
        return logo_html_top, logo_html_bottom
    
    def _generate_main_template(self, logo_html_top: str, logo_html_bottom: str) -> str:
        """
        Generate main HTML template structure.
        
        Args:
            logo_html_top: HTML for top row logos
            logo_html_bottom: HTML for bottom row logos
            
        Returns:
            HTML template string with placeholders
        """
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sanger Pipeline QC Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    {self._generate_css_styles()}
</head>
<body>
    <div class="main-container">
        {self._generate_header_html(logo_html_top, logo_html_bottom)}
        
        <!-- Main Content -->
        <div class="content-area">
            <!-- Navigation Tabs -->
            <ul class="nav nav-tabs custom-tabs" role="tablist">
                <li class="nav-item">
                    <a class="nav-link active" data-bs-toggle="tab" href="#overview" role="tab">
                        <i class="fas fa-chart-line"></i> Overview
                    </a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-bs-toggle="tab" href="#directories" role="tab">
                        <i class="fas fa-folder-open"></i> Directories
                    </a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-bs-toggle="tab" href="#samples" role="tab">
                        <i class="fas fa-dna"></i> Samples
                    </a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-bs-toggle="tab" href="#damage" role="tab">
                        <i class="fas fa-exclamation-triangle"></i> Damage Analysis
                    </a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-bs-toggle="tab" href="#hvs" role="tab">
                        <i class="fas fa-layer-group"></i> HVS Regions
                    </a>
                </li>
                <li class="nav-item">
                    <a class="nav-link" data-bs-toggle="tab" href="#sample-details" role="tab">
                        <i class="fas fa-microscope"></i> Sample Details
                    </a>
                </li>
            </ul>
            
            <!-- Tab Content -->
            <div class="tab-content mt-4">
                {{tab_content}}
            </div>
        </div>
        
        {self._generate_footer_html()}
    </div>
    
    <!-- Bootstrap JS -->
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    
    <!-- Chart JavaScript -->
    <script>
        {{charts_javascript}}
    </script>
</body>
</html>
"""
    
    def _generate_css_styles(self) -> str:
        """Generate CSS styles for the report."""
        return f"""
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
        
        .header h1 {{
            margin: 0;
            font-size: 2.5rem;
            font-weight: 300;
        }}
        
        .logo-container {{
            display: flex;
            flex-direction: column;
            gap: 15px;
        }}
        
        .logo-row {{
            display: flex;
            gap: 15px;
            align-items: center;
        }}
        
        .logo-container img {{
            height: 40px;
            max-width: 120px;
            object-fit: contain;
            transition: transform 0.3s ease;
        }}
        
        .logo-container img:hover {{
            transform: scale(1.1);
        }}
        
        .custom-tabs {{
            border-bottom: 3px solid var(--primary-color);
            background: #f8f9fa;
            border-radius: 10px 10px 0 0;
            padding: 0 1rem;
        }}
        
        .custom-tabs .nav-link {{
            border: none;
            color: #666;
            font-weight: 500;
            padding: 1rem 1.5rem;
            transition: all 0.3s ease;
        }}
        
        .custom-tabs .nav-link:hover {{
            color: var(--primary-color);
            background: rgba(46, 134, 171, 0.1);
        }}
        
        .custom-tabs .nav-link.active {{
            color: var(--primary-color);
            background: white;
            border-bottom: 3px solid var(--primary-color);
        }}
        
        .content-area {{
            padding: 2rem;
        }}
        
        .stat-card {{
            background: white;
            border-radius: 10px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        
        .stat-card:hover {{
            transform: translateY(-2px);
        }}
        
        .stat-number {{
            font-size: 2rem;
            font-weight: bold;
            color: var(--primary-color);
        }}
        
        .stat-label {{
            color: #666;
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .table {{
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .table thead {{
            background: var(--primary-color);
            color: white;
        }}
        
        .badge-custom {{
            background: var(--success-color);
            color: white;
            padding: 0.5rem 1rem;
            border-radius: 20px;
            font-size: 0.8rem;
        }}
        
        .footer {{
            background: #f8f9fa;
            padding: 1rem 2rem;
            text-align: center;
            color: #666;
            border-top: 1px solid #dee2e6;
        }}
        
        .chart-container {{
            position: relative;
            height: 400px;
            margin: 1rem 0;
        }}
        
        .alert-custom {{
            border-radius: 10px;
            border: none;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        @media (max-width: 768px) {{
            .header {{
                flex-direction: column;
                text-align: center;
                gap: 1rem;
            }}
            
            .header h1 {{
                font-size: 2rem;
            }}
            
            .logo-container {{
                align-items: center;
            }}
            
            .main-container {{
                margin: 10px;
            }}
        }}
    </style>
"""
    
    def _generate_header_html(self, logo_html_top: str, logo_html_bottom: str) -> str:
        """Generate header HTML section."""
        return f"""
        <div class="header">
            <div>
                <h1><i class="fas fa-dna"></i> Sanger Pipeline QC Report</h1>
                <p class="mb-0">Comprehensive quality control and damage analysis</p>
            </div>
            <div class="logo-container">
                <div class="logo-row">
                    {logo_html_top}
                </div>
                <div class="logo-row">
                    {logo_html_bottom}
                </div>
            </div>
        </div>
"""
    
    def _generate_footer_html(self) -> str:
        """Generate footer HTML section."""
        return """
        <div class="footer">
            <p class="mb-0">
                Generated by Sanger aDNA Pipeline v2.0.0 | 
                <i class="fas fa-calendar-alt"></i> Report Date: <span id="report-date"></span>
            </p>
        </div>
        <script>
            document.getElementById('report-date').textContent = new Date().toLocaleDateString();
        </script>
"""
    
    def _generate_all_tabs(self, stats: Dict[str, Any]) -> str:
        """
        Generate content for all tabs.
        
        Args:
            stats: Statistics dictionary
            
        Returns:
            HTML content for all tabs
        """
        return f"""
            {self._generate_overview_tab(stats)}
            {self._generate_directories_tab(stats)}
            {self._generate_samples_tab(stats)}
            {self._generate_damage_tab(stats)}
            {self._generate_hvs_tab(stats)}
            {self._generate_sample_details_tab(stats)}
        """
    
    def _generate_overview_tab(self, stats: Dict[str, Any]) -> str:
        """Generate overview tab content."""
        total_files = sum(stats.get(dir_name, {}).get('file_count', 0) 
                         for dir_name in ['input', 'output', 'consensus', 'final'])
        
        return f"""
        <div class="tab-pane fade show active" id="overview" role="tabpanel">
            <div class="row">
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{total_files}</div>
                        <div class="stat-label">Total Files</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{len(stats.get('samples', {}))}</div>
                        <div class="stat-label">Samples Processed</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{stats.get('damage_data', {}).get('files_analyzed', 0)}</div>
                        <div class="stat-label">Damage Files</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{stats.get('hvs_combinations', {}).get('total_merged_files', 0)}</div>
                        <div class="stat-label">Merged Files</div>
                    </div>
                </div>
            </div>
            
            <div class="row mt-4">
                <div class="col-12">
                    <div class="stat-card">
                        <h3><i class="fas fa-info-circle"></i> Pipeline Information</h3>
                        <table class="table">
                            <tr>
                                <td><strong>Report Generated:</strong></td>
                                <td>{stats.get('report_generated', 'Unknown')}</td>
                            </tr>
                            <tr>
                                <td><strong>Pipeline Version:</strong></td>
                                <td>{stats.get('pipeline_version', 'Unknown')}</td>
                            </tr>
                            <tr>
                                <td><strong>Output Directory:</strong></td>
                                <td><code>{stats.get('output_directory', 'Unknown')}</code></td>
                            </tr>
                        </table>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_directories_tab(self, stats: Dict[str, Any]) -> str:
        """Generate directories tab content."""
        directories_html = ""
        
        for dir_name in ['input', 'output', 'consensus', 'final']:
            dir_data = stats.get(dir_name, {})
            if dir_data.get('exists', False):
                file_types_html = ""
                for ext, count in dir_data.get('file_types', {}).items():
                    file_types_html += f'<span class="badge badge-custom me-1">{ext}: {count}</span>'
                
                directories_html += f"""
                <div class="col-md-6 mb-4">
                    <div class="stat-card">
                        <h4><i class="fas fa-folder"></i> {dir_name.title()}</h4>
                        <p><strong>Files:</strong> {dir_data.get('file_count', 0)}</p>
                        <p><strong>Size:</strong> {dir_data.get('total_size_mb', 0)} MB</p>
                        <p><strong>Types:</strong><br>{file_types_html}</p>
                    </div>
                </div>
                """
            else:
                directories_html += f"""
                <div class="col-md-6 mb-4">
                    <div class="stat-card">
                        <h4><i class="fas fa-folder"></i> {dir_name.title()}</h4>
                        <div class="alert alert-warning">Directory does not exist</div>
                    </div>
                </div>
                """
        
        return f"""
        <div class="tab-pane fade" id="directories" role="tabpanel">
            <div class="row">
                {directories_html}
            </div>
        </div>
        """
    
    def _generate_samples_tab(self, stats: Dict[str, Any]) -> str:
        """Generate samples tab content."""
        samples_data = stats.get('samples', {})
        
        if not samples_data:
            return """
            <div class="tab-pane fade" id="samples" role="tabpanel">
                <div class="alert alert-info">No sample data available</div>
            </div>
            """
        
        samples_rows = ""
        for sample_name, sample_info in samples_data.items():
            hvs_regions = ', '.join(sample_info.get('hvs_regions', []))
            samples_rows += f"""
            <tr>
                <td>{sample_name}</td>
                <td>{len(sample_info.get('consensus_files', []))}</td>
                <td>{len(sample_info.get('final_files', []))}</td>
                <td>{len(sample_info.get('damage_files', []))}</td>
                <td>{hvs_regions}</td>
            </tr>
            """
        
        return f"""
        <div class="tab-pane fade" id="samples" role="tabpanel">
            <div class="stat-card">
                <h3><i class="fas fa-dna"></i> Sample Processing Summary</h3>
                <div class="table-responsive">
                    <table class="table table-striped">
                        <thead>
                            <tr>
                                <th>Sample Name</th>
                                <th>Consensus Files</th>
                                <th>Final Files</th>
                                <th>Damage Files</th>
                                <th>HVS Regions</th>
                            </tr>
                        </thead>
                        <tbody>
                            {samples_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        """
    
    def _generate_damage_tab(self, stats: Dict[str, Any]) -> str:
        """Generate damage analysis tab content."""
        damage_data = stats.get('damage_data', {})
        
        if damage_data.get('files_analyzed', 0) == 0:
            return """
            <div class="tab-pane fade" id="damage" role="tabpanel">
                <div class="alert alert-info">No damage analysis data available</div>
            </div>
            """
        
        return f"""
        <div class="tab-pane fade" id="damage" role="tabpanel">
            <div class="stat-card">
                <h3><i class="fas fa-exclamation-triangle"></i> Damage Analysis Summary</h3>
                <p><strong>Files Analyzed:</strong> {damage_data.get('files_analyzed', 0)}</p>
                <!-- Damage analysis charts and details would go here -->
                <div id="damage-chart" class="chart-container"></div>
            </div>
        </div>
        """
    
    def _generate_hvs_tab(self, stats: Dict[str, Any]) -> str:
        """Generate HVS regions tab content."""
        hvs_data = stats.get('hvs_combinations', {})
        
        if not hvs_data.get('combinations'):
            return """
            <div class="tab-pane fade" id="hvs" role="tabpanel">
                <div class="alert alert-info">No HVS combination data available</div>
            </div>
            """
        
        combinations_html = ""
        for combo, count in hvs_data.get('combinations', {}).items():
            combinations_html += f"""
            <tr>
                <td>{combo}</td>
                <td>{count}</td>
            </tr>
            """
        
        return f"""
        <div class="tab-pane fade" id="hvs" role="tabpanel">
            <div class="stat-card">
                <h3><i class="fas fa-layer-group"></i> HVS Region Combinations</h3>
                <div class="table-responsive">
                    <table class="table table-striped">
                        <thead>
                            <tr>
                                <th>HVS Combination</th>
                                <th>Count</th>
                            </tr>
                        </thead>
                        <tbody>
                            {combinations_html}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        """
    
    def _generate_sample_details_tab(self, stats: Dict[str, Any]) -> str:
        """Generate sample details tab content."""
        return """
        <div class="tab-pane fade" id="sample-details" role="tabpanel">
            <div class="stat-card">
                <h3><i class="fas fa-microscope"></i> Detailed Sample Analysis</h3>
                <p>Detailed sample analysis and visualizations would be displayed here.</p>
                <!-- Individual sample plots and detailed analysis -->
            </div>
        </div>
        """
    
    def _generate_charts_javascript(self, stats: Dict[str, Any]) -> str:
        """Generate JavaScript for charts and interactive elements."""
        return """
        // Initialize charts when tab is shown
        document.addEventListener('DOMContentLoaded', function() {
            // Add chart initialization code here
            console.log('Charts initialized');
        });
        """
