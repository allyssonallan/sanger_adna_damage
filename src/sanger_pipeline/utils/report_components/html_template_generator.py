"""
HTML template generation for Sanger Pipeline QC reports.

This module handles the generation of HTML content, templates,
and styling for comprehensive QC reports.
"""

import base64
import html
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

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
            "primary": "#2E86AB",
            "secondary": "#A23B72",
            "success": "#F18F01",
            "warning": "#C73E1D",
            "info": "#6C757D",
            "light": "#F8F9FA",
            "dark": "#343A40",
        }

    def _encode_logos(self) -> Dict[str, str]:
        """
        Encode logo images to base64 for embedding in HTML.

        Returns:
            Dictionary mapping logo names to base64 encoded strings
        """
        logos: Dict[str, str] = {}
        module_path = Path(__file__).resolve()
        project_root = module_path.parents[4]
        candidate_dirs = [
            project_root / "config" / "logos",
            self.output_dir / "config" / "logos",
        ]

        seen_files = set()

        for logo_dir in candidate_dirs:
            if not logo_dir.exists():
                continue

            for logo_file in logo_dir.glob("*.png"):
                lowercase_name = logo_file.stem.lower()
                if lowercase_name in seen_files:
                    continue

                try:
                    with open(logo_file, "rb") as f:
                        encoded = base64.b64encode(f.read()).decode("utf-8")
                        logos[lowercase_name] = encoded
                        seen_files.add(lowercase_name)
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
        safe_template = (
            html_template.replace("{", "{{")
            .replace("}", "}}")
            .replace("{{tab_content}}", "{tab_content}")
            .replace("{{charts_javascript}}", "{charts_javascript}")
        )

        return safe_template.format(
            tab_content=tab_content, charts_javascript=charts_js
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
            if "ufc" in logos:
                logo_html_top += f'<a href="https://ufc.br" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["ufc"]}" alt="UFC Logo" title="Universidade Federal do Ceará - https://ufc.br"></a>'
            if "funcap" in logos:
                logo_html_top += f'<a href="https://funcap.ce.gov.br" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["funcap"]}" alt="FUNCAP Logo" title="Fundação Cearense de Apoio ao Desenvolvimento Científico e Tecnológico - https://www.funcap.ce.gov.br"></a>'

            # Bottom row: LABBAT and NPDM
            if "labbat" in logos:
                logo_html_bottom += f'<a href="https://instagram.com/labbat.npdm.ufc" target="_blank" rel="noopener noreferrer"><img src="data:image/png;base64,{logos["labbat"]}" alt="LABBAT Logo" title="Laboratório de Bioarqueologia Translacional - https://instagram.com/labbat.npdm.ufc"></a>'
            if "npdm" in logos:
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
        css_styles = self._generate_css_styles()
        header_html = self._generate_header_html(logo_html_top, logo_html_bottom)
        footer_html = self._generate_footer_html()

        return (
            """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sanger Pipeline QC Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>"""
            + css_styles
            + """
</head>
<body>
    <div class="main-container">
"""
            + header_html
            + """
        
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
                        <i class="fas fa-exclamation-triangle"></i> N Content
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
                {tab_content}
            </div>
        </div>
        
"""
            + footer_html
            + """
    </div>
    
    <!-- Bootstrap JS -->
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    
    <!-- Chart JavaScript -->
    <script>
        {charts_javascript}
    </script>
</body>
</html>
"""
        )

    def _generate_css_styles(self) -> str:
        """Generate CSS styles for the report."""
        return """
    <style>
        :root {
            --card-radius: 16px;
            --border-soft: #eef2f7;
            --text-muted: #6c7a89;
            --brand: #2e86ab;
            --surface: #ffffff;
            --surface-alt: #f7f9fc;
        }

        body {
            background-color: #f4f6fb;
            color: #1f232a;
            font-family: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
            line-height: 1.6;
        }

        a {
            color: inherit;
        }

        .main-container {
            max-width: 1220px;
            margin: 0 auto;
            padding: 2.5rem 1.5rem 3rem;
        }

        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 2.5rem;
            padding-bottom: 1.5rem;
            border-bottom: 1px solid var(--border-soft);
        }

        .header h1 {
            margin: 0;
            font-weight: 600;
            font-size: 2.15rem;
            letter-spacing: -0.015em;
        }

        .header p {
            color: var(--text-muted);
            font-size: 0.95rem;
        }

        .logo-container {
            display: flex;
            flex-direction: column;
            gap: 0.75rem;
            align-items: flex-end;
        }

        .logo-row {
            display: flex;
            gap: 1rem;
            align-items: center;
        }

        .logo-row img {
            max-height: 56px;
            max-width: 150px;
            width: auto;
            object-fit: contain;
            border-radius: 12px;
            background: rgba(255,255,255,0.85);
            padding: 0.35rem 0.65rem;
            box-shadow: 0 6px 18px rgba(31, 45, 61, 0.08);
        }

        .content-area {
            margin-top: 2.5rem;
        }

        .custom-tabs .nav-link {
            border: none;
            border-bottom: 2px solid transparent;
            color: var(--text-muted);
            font-weight: 500;
            font-size: 0.95rem;
            padding: 0.75rem 1.35rem;
            transition: color 0.2s ease, border-color 0.2s ease;
        }

        .custom-tabs .nav-link:hover,
        .custom-tabs .nav-link:focus {
            color: #1f232a;
        }

        .custom-tabs .nav-link.active {
            color: #1f232a;
            border-bottom-color: var(--brand);
            background: none;
        }

        .tab-content {
            margin-top: 1.75rem;
            background: var(--surface);
            border-radius: var(--card-radius);
            box-shadow: 0 24px 40px rgba(31, 45, 61, 0.08);
            padding: 2.25rem;
        }

        .stat-card {
            background: var(--surface);
            border-radius: var(--card-radius);
            border: 1px solid var(--border-soft);
            padding: 1.75rem;
            box-shadow: none;
        }

        .stat-card + .stat-card {
            margin-top: 1.5rem;
        }

        .stat-number {
            font-size: 2.4rem;
            font-weight: 600;
            color: var(--brand);
            margin-bottom: 0.35rem;
        }

        .stat-label {
            color: var(--text-muted);
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.12em;
        }

        .badge-custom {
            display: inline-flex;
            align-items: center;
            gap: 0.35rem;
            background: #e3f2fd;
            color: #1f4b7a;
            border-radius: 999px;
            padding: 0.32rem 0.75rem;
            font-size: 0.75rem;
        }

        .table {
            margin-bottom: 0;
        }

        .table thead th {
            border-bottom-width: 1px;
            color: #1f232a;
            font-weight: 600;
        }

        .table td {
            color: #2f3b49;
            vertical-align: middle;
        }

        ul {
            padding-left: 1.1rem;
        }

        .footer {
            text-align: center;
            padding-top: 1.75rem;
            color: #94a3b8;
            font-size: 0.85rem;
        }

        .alert {
            border-radius: 12px;
            border: none;
        }

        .config-section {
            background: var(--surface-alt);
            border-radius: 14px;
        }

        .sample-stack {
            display: flex;
            flex-direction: column;
            gap: 1.5rem;
        }

        .sample-card {
            border-radius: var(--card-radius);
            border: 1px solid var(--border-soft);
            background: var(--surface);
            padding: 1.75rem;
            box-shadow: none;
        }

        .sample-card h4 {
            margin-bottom: 0.25rem;
            font-size: 1.15rem;
            font-weight: 600;
        }

        .meta-group {
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
            margin-bottom: 1.25rem;
        }

        .meta-chip {
            display: inline-flex;
            align-items: center;
            gap: 0.3rem;
            background: #eef6ff;
            color: #29649b;
            border-radius: 999px;
            padding: 0.3rem 0.75rem;
            font-size: 0.78rem;
            letter-spacing: 0.02em;
        }

        .meta-chip.muted {
            background: #f1f3f5;
            color: #616e7c;
        }

        .list-section + .list-section {
            margin-top: 1.25rem;
        }

        .list-section h5 {
            font-size: 0.95rem;
            font-weight: 600;
            margin-bottom: 0.5rem;
            color: #1f232a;
        }

        @media (max-width: 992px) {
            .header {
                flex-direction: column;
                align-items: flex-start;
            }

            .logo-container {
                align-items: flex-start;
            }

            .tab-content {
                padding: 1.75rem;
            }
        }

        @media (max-width: 576px) {
            .custom-tabs .nav-link {
                padding: 0.5rem 0.75rem;
                font-size: 0.85rem;
            }

            .stat-number {
                font-size: 1.85rem;
            }

            .logo-row img {
                max-height: 46px;
            }
        }
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

    def _format_report_date(self, value: Any) -> str:
        """Format ISO datetime strings into a readable date/time."""
        if not value:
            return "Unknown"

        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M")

        if isinstance(value, str):
            try:
                parsed = datetime.fromisoformat(value)
                return parsed.strftime("%Y-%m-%d %H:%M")
            except ValueError:
                return html.escape(value)

        return html.escape(str(value))

    def _format_ratio_percentage(self, value: Any, precision: int = 1) -> str:
        """Format a 0-1 ratio as a percentage string."""
        try:
            return f"{float(value) * 100:.{precision}f}%"
        except (TypeError, ValueError):
            return "N/A"

    def _format_plain_percentage(self, value: Any, precision: int = 1) -> str:
        """Format an absolute percentage value (0-100)."""
        try:
            return f"{float(value):.{precision}f}%"
        except (TypeError, ValueError):
            return "N/A"

    def _format_number(self, value: Any, precision: int = 2) -> str:
        """Format a numeric value with the given precision."""
        try:
            return f"{float(value):.{precision}f}"
        except (TypeError, ValueError):
            return "N/A"

    def _render_file_list(self, files: List[str], max_items: int = 5) -> str:
        """Render an HTML unordered list for a collection of file names."""
        if not files:
            return "<em>None</em>"

        items = [
            f"<li>{html.escape(file_name)}</li>" for file_name in files[:max_items]
        ]
        remaining = len(files) - max_items
        if remaining > 0:
            items.append(f"<li>… and {remaining} more</li>")

        return '<ul class="mb-0">' + "".join(items) + "</ul>"

    def _render_configuration_sections(
        self, config: Dict[str, Any], source: Optional[str]
    ) -> str:
        """Render configuration parameters for the overview tab."""

        if not config:
            return """
            <div class="row mt-4">
                <div class="col-12">
                    <div class="stat-card">
                        <div class="alert alert-info mb-0">Configuration details unavailable for this run.</div>
                    </div>
                </div>
            </div>
            """

        sections_html = "".join(
            self._render_config_section(section, values)
            for section, values in sorted(config.items())
        )

        source_html = (
            f'<p class="text-muted small mb-3">Source: <code>{html.escape(source)}</code></p>'
            if source
            else ""
        )

        return f"""
        <div class="row mt-4">
            <div class="col-12">
                <div class="stat-card">
                    <h3><i class="fas fa-cogs"></i> Configuration Parameters</h3>
                    {source_html}
                    <div class="row mt-2">
                        {sections_html}
                    </div>
                </div>
            </div>
        </div>
        """

    def _render_config_section(self, title: str, content: Any) -> str:
        """Render an individual configuration section."""

        safe_title = html.escape(title.replace("_", " ").title())
        body_html = self._render_config_content(content)

        return f"""
            <div class="col-md-6 mb-3">
                <div class="config-section border rounded p-3 bg-light h-100">
                    <h5 class="mb-3"><i class="fas fa-sliders-h"></i> {safe_title}</h5>
                    {body_html}
                </div>
            </div>
        """

    def _render_config_content(self, content: Any) -> str:
        """Render configuration content, supporting nested structures."""

        if isinstance(content, dict):
            rows = "".join(
                f"<tr><td>{html.escape(str(key))}</td><td>{self._format_config_value(value)}</td></tr>"
                for key, value in sorted(content.items())
            )
            return f"""
                <div class="table-responsive">
                    <table class="table table-sm table-borderless mb-0">
                        <tbody>{rows}</tbody>
                    </table>
                </div>
            """

        return f'<p class="mb-0">{self._format_config_value(content)}</p>'

    def _format_config_value(self, value: Any) -> str:
        """Format configuration values into HTML-safe strings."""

        if isinstance(value, dict):
            items = "".join(
                f"<li><strong>{html.escape(str(key))}:</strong> {self._format_config_value(val)}</li>"
                for key, val in sorted(value.items())
            )
            return '<ul class="mb-0">' + items + "</ul>"

        if isinstance(value, (list, tuple, set)):
            if not value:
                return "<em>None</em>"
            return ", ".join(html.escape(str(item)) for item in value)

        if value is None:
            return "<em>None</em>"

        return html.escape(str(value))

    def _generate_all_tabs(self, stats: Dict[str, Any]) -> str:
        """
        Generate content for all tabs.

        Args:
            stats: Statistics dictionary

        Returns:
            HTML content for all tabs
        """
        sections = [
            self._generate_overview_tab(stats),
            self._generate_directories_tab(stats),
            self._generate_samples_tab(stats),
            self._generate_damage_tab(stats),
            self._generate_hvs_tab(stats),
            self._generate_sample_details_tab(stats),
        ]

        return "\n".join(sections)

    def _generate_overview_tab(self, stats: Dict[str, Any]) -> str:
        """Generate overview tab content."""
        directories = stats.get("directories", {})
        total_files = sum(
            directory_stats.get("file_count", 0)
            for directory_stats in directories.values()
            if isinstance(directory_stats, dict)
        )

        config_html = self._render_configuration_sections(
            stats.get("config_parameters", {}), stats.get("config_source")
        )

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
                        <div class="stat-number">{len(stats.get("samples", {}))}</div>
                        <div class="stat-label">Samples Processed</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{stats.get("damage_data", {}).get("files_analyzed", 0)}</div>
                        <div class="stat-label">N Result Files</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stat-card text-center">
                        <div class="stat-number">{stats.get("hvs_combinations", {}).get("total_merged_files", 0)}</div>
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
                                <td>{self._format_report_date(stats.get("report_generated"))}</td>
                            </tr>
                            <tr>
                                <td><strong>Pipeline Version:</strong></td>
                                <td>{stats.get("pipeline_version", "Unknown")}</td>
                            </tr>
                            <tr>
                                <td><strong>Output Directory:</strong></td>
                                <td><code>{stats.get("output_directory", "Unknown")}</code></td>
                            </tr>
                            <tr>
                                <td><strong>Minimum Phred Score:</strong></td>
                                <td>{stats.get("quality_settings", {}).get("min_phred_score", "Unknown")}</td>
                            </tr>
                            <tr>
                                <td><strong>Minimum Sequence Length:</strong></td>
                                <td>{stats.get("quality_settings", {}).get("min_sequence_length", "Unknown")}</td>
                            </tr>
                        </table>
                    </div>
                </div>
            </div>
            {config_html}
        </div>
        """

    def _generate_directories_tab(self, stats: Dict[str, Any]) -> str:
        """Generate directories tab content."""
        directories = stats.get("directories", {})
        directory_order = [
            ("fasta", "FASTA"),
            ("filtered", "Filtered"),
            ("consensus", "Consensus"),
            ("final", "Final"),
            ("aligned", "Aligned"),
            ("damage_analysis", "N Content Results"),
            ("plots", "Plots"),
        ]

        directories_html = ""
        for key, label in directory_order:
            dir_data = directories.get(key, {})
            safe_dir_title = html.escape(label)

            if dir_data.get("exists", False):
                file_types_html = ""
                for ext, count in sorted(dir_data.get("file_types", {}).items()):
                    file_types_html += (
                        f'<span class="badge badge-custom me-1">{html.escape(ext)}: '
                        f"{count}</span>"
                    )

                directories_html += f"""
                <div class="col-md-6 mb-4">
                    <div class="stat-card">
                        <h4><i class="fas fa-folder"></i> {safe_dir_title}</h4>
                        <p><strong>Files:</strong> {dir_data.get("file_count", 0)}</p>
                        <p><strong>Size:</strong> {dir_data.get("total_size_mb", 0)} MB</p>
                        <div class="mb-2"><strong>Types:</strong><br>{file_types_html}</div>
                        <div><strong>Sample files:</strong> {self._render_file_list(dir_data.get("files", []))}</div>
                    </div>
                </div>
                """
            else:
                directories_html += f"""
                <div class="col-md-6 mb-4">
                    <div class="stat-card">
                        <h4><i class="fas fa-folder"></i> {safe_dir_title}</h4>
                        <div class="alert alert-warning mb-0">Directory does not exist</div>
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
        samples_data = stats.get("samples", {})

        if not samples_data:
            return """
            <div class="tab-pane fade" id="samples" role="tabpanel">
                <div class="alert alert-info">No sample data available</div>
            </div>
            """

        samples_rows = ""
        for sample_name, sample_info in samples_data.items():
            safe_name = html.escape(sample_name)
            hvs_regions = ", ".join(
                html.escape(region) for region in sample_info.get("hvs_regions", [])
            )
            samples_rows += f"""
            <tr>
                <td>{safe_name}</td>
                <td>{len(sample_info.get("consensus_files", []))}</td>
                <td>{len(sample_info.get("final_files", []))}</td>
                <td>{len(sample_info.get("damage_files", []))}</td>
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
        """Generate N content tab content."""
        damage_data = stats.get("damage_data", {})

        if damage_data.get("files_analyzed", 0) == 0:
            return """
            <div class="tab-pane fade" id="damage" role="tabpanel">
                <div class="alert alert-info">No N content data available</div>
            </div>
            """

        files_analyzed = damage_data.get("files_analyzed", 0)
        summary = damage_data.get("summary", {})
        individual = damage_data.get("individual_results", [])

        summary_metrics = [
            ("Total N Bases", self._format_number(summary.get("total_n_bases"), 0)),
            (
                "Average N Percentage",
                self._format_plain_percentage(summary.get("average_n_percentage")),
            ),
            (
                "Median N Percentage",
                self._format_plain_percentage(summary.get("median_n_percentage")),
            ),
            (
                "Maximum N Percentage",
                self._format_plain_percentage(summary.get("max_n_percentage")),
            ),
            (
                "Total Valid Bases",
                self._format_number(summary.get("total_valid_bases"), 0),
            ),
            (
                "Average Valid Percentage",
                self._format_plain_percentage(summary.get("average_valid_percentage")),
            ),
        ]

        summary_rows = "".join(f"""
                <tr>
                    <td>{html.escape(label)}</td>
                    <td>{value if value not in (None, "N/A") else "N/A"}</td>
                </tr>
            """ for label, value in summary_metrics if value not in (None, "N/A"))

        if not summary_rows:
            summary_table = (
                '<div class="alert alert-info">No summary data available yet.</div>'
            )
        else:
            summary_table = f"""
                <table class="table table-striped">
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                        </tr>
                    </thead>
                    <tbody>
                        {summary_rows}
                    </tbody>
                </table>
            """

        individual_rows = ""
        for item in individual:
            sample = html.escape(str(item.get("sample", "Unknown")))
            n_bases = self._format_number(item.get("n_content"), 0)
            n_pct = self._format_plain_percentage(item.get("n_percentage"))
            valid_bases = self._format_number(item.get("valid_bases"), 0)
            valid_pct = self._format_plain_percentage(item.get("valid_percentage"))
            ambiguous = self._format_number(item.get("ambiguous_content"), 0)
            total_bases = self._format_number(item.get("total_bases"), 0)

            individual_rows += f"""
                <tr>
                    <td>{sample}</td>
                    <td>{n_bases}</td>
                    <td>{n_pct}</td>
                    <td>{valid_bases}</td>
                    <td>{valid_pct}</td>
                    <td>{ambiguous}</td>
                    <td>{total_bases}</td>
                </tr>
            """

        if not individual_rows:
            individual_table = '<div class="alert alert-info">No sample-level N content available.</div>'
        else:
            individual_table = f"""
                <div class="table-responsive">
                    <table class="table table-sm table-hover">
                        <thead>
                            <tr>
                                <th>Sample</th>
                                <th>N Bases</th>
                                <th>N Percentage</th>
                                <th>Valid Bases</th>
                                <th>Valid Percentage</th>
                                <th>Ambiguous Bases</th>
                                <th>Total Bases</th>
                            </tr>
                        </thead>
                        <tbody>
                            {individual_rows}
                        </tbody>
                    </table>
                </div>
            """

        return f"""
        <div class="tab-pane fade" id="damage" role="tabpanel">
            <div class="stat-card">
                <h3><i class="fas fa-exclamation-triangle"></i> N Content Overview</h3>
                <p><strong>Files Analyzed:</strong> {files_analyzed}</p>
                <h5>N Content Summary</h5>
                {summary_table}
                <h5 class="mt-4">Sample-Level N Content</h5>
                {individual_table}
            </div>
        </div>
        """

    def _generate_hvs_tab(self, stats: Dict[str, Any]) -> str:
        """Generate HVS regions tab content."""
        hvs_data = stats.get("hvs_combinations", {})

        if not hvs_data.get("combinations"):
            return """
            <div class="tab-pane fade" id="hvs" role="tabpanel">
                <div class="alert alert-info">No HVS combination data available</div>
            </div>
            """

        combinations_html = ""
        for combo, count in hvs_data.get("combinations", {}).items():
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
        samples = stats.get("samples", {})

        if not samples:
            return """
            <div class=\"tab-pane fade\" id=\"sample-details\" role=\"tabpanel\">
                <div class=\"stat-card\">
                    <div class=\"alert alert-info mb-0\">No per-sample details collected.</div>
                </div>
            </div>
            """

        sample_panels = ""
        for sample_name in sorted(samples):
            info = samples[sample_name]
            safe_name = html.escape(sample_name)
            consensus_list = info.get("consensus_files", [])
            final_list = info.get("final_files", [])
            damage_list = info.get("damage_files", [])
            hvs_regions = sorted(info.get("hvs_regions", []))

            consensus_files = self._render_file_list(consensus_list)
            final_files = self._render_file_list(final_list)
            damage_files = self._render_file_list(damage_list)

            if hvs_regions:
                hvs_html = "".join(
                    f'<span class="meta-chip"><i class="fas fa-hashtag"></i> {html.escape(region)}</span>'
                    for region in hvs_regions
                )
            else:
                hvs_html = '<span class="meta-chip muted">No HVS regions</span>'

            meta_chips = "".join(
                [
                    f'<span class="meta-chip"><i class="fas fa-layer-group"></i> {len(consensus_list)} consensus</span>',
                    f'<span class="meta-chip"><i class="fas fa-dna"></i> {len(final_list)} final</span>',
                    f'<span class="meta-chip"><i class="fas fa-chart-area"></i> {len(damage_list)} reports</span>',
                ]
            )

            sample_panels += f"""
                <article class=\"sample-card\">
                    <div class=\"d-flex justify-content-between align-items-start flex-wrap gap-3\">
                        <div>
                            <h4><i class=\"fas fa-vial\"></i> {safe_name}</h4>
                            <div class=\"meta-group\">{hvs_html}</div>
                        </div>
                        <div class=\"meta-group\">{meta_chips}</div>
                    </div>
                    <div class=\"list-section\">
                        <h5><i class=\"fas fa-layer-group\"></i> Consensus Files</h5>
                        {consensus_files}
                    </div>
                    <div class=\"list-section\">
                        <h5><i class=\"fas fa-dna\"></i> Final Assemblies</h5>
                        {final_files}
                    </div>
                    <div class=\"list-section\">
                        <h5><i class=\"fas fa-chart-area\"></i> N Content / Damage Reports</h5>
                        {damage_files}
                    </div>
                </article>
            """

        return f"""
        <div class=\"tab-pane fade\" id=\"sample-details\" role=\"tabpanel\">
            <div class=\"sample-stack\">
                {sample_panels}
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
