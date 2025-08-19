# Sanger aDNA Damage Analysis - Modular Browser Application

## 📁 Project Structure

``` bash
browser/
├── index.html              # Main HTML file (modular version)
├── browser_page.html       # Original monolithic version (backup)
├── browser_page_checkpoint.html  # Checkpoint backup
├── css/
│   └── main.css            # Extracted styles
└── js/
    ├── main.js             # Application initialization & coordination
    ├── core.js             # Base utilities, theme, constants, data store
    ├── abif-parser.js      # ABIF file format parser
    ├── sequence-processor.js # Quality filtering, alignment, damage analysis
    ├── chromatogram.js     # Canvas rendering and visualization
    ├── pipeline.js         # Main processing workflow
    ├── ui-manager.js       # User interface components & management
    └── export-manager.js   # ZIP creation and downloads
```

## 🏗️ Module Architecture

### **Core Module** (`core.js`)

- **Theme management**: Light/dark mode toggle
- **Constants**: Quality thresholds, trace colors, damage parameters
- **Data store**: Centralized storage for all pipeline results
- **Utilities**: File parsing, counters, formatting functions
- **Logo management**: Auto-injection of institutional logos

### **ABIF Parser** (`abif-parser.js`)

- **File format parsing**: Reads Sanger .ab1 files
- **Data extraction**: Bases, quality scores, chromatogram traces
- **Error handling**: Validates file format and structure

### **Sequence Processor** (`sequence-processor.js`)

- **Quality filtering**: Q30+ base calling with N replacement
- **Reverse complement**: Proper DNA strand handling
- **Sequence alignment**: Overlap consensus for F+R reads
- **Damage analysis**: C→T and G→A transition detection with bootstrap confidence intervals

### **Chromatogram Renderer** (`chromatogram.js`)

- **Canvas visualization**: Real-time trace plotting
- **Multi-channel display**: A, C, G, T trace rendering with color coding
- **Export functionality**: PNG snapshot generation
- **Error handling**: Graceful fallback for rendering issues

### **Pipeline Controller** (`pipeline.js`)

- **File processing**: Orchestrates the complete analysis workflow
- **Stage management**: Individual and batch processing modes
- **Data flow**: Coordinates between parser, processor, and UI
- **Progress tracking**: Updates counters and statistics

### **UI Manager** (`ui-manager.js`)

- **Interface control**: Manages all user interactions
- **Statistics display**: Dynamic table updates and empty states
- **Modal management**: File summary and detailed views
- **Event handling**: Button clicks, file uploads, selections

### **Export Manager** (`export-manager.js`)

- **ZIP creation**: In-browser archive generation
- **File downloads**: Blob URL handling for multiple formats
- **Format support**: FASTA, CSV, JSON, PNG, HSD
- **Error recovery**: Fallback to individual downloads

### **Main Application** (`main.js`)

- **Initialization**: Coordinates all module startup
- **Dependency injection**: Connects modules with proper dependencies
- **Error handling**: Global error management and user feedback
- **Public API**: Provides external access points

## 🎯 Key Features

### **Modular Design Benefits**

- ✅ **Separation of concerns**: Each module has a single responsibility
- ✅ **Maintainability**: Easy to modify individual components
- ✅ **Testability**: Modules can be unit tested independently
- ✅ **Reusability**: Components can be reused in other projects
- ✅ **Debugging**: Isolated error tracking and resolution

### **Performance Optimizations**

- ✅ **Lazy loading**: Modules load only when needed
- ✅ **Memory management**: Efficient data structures and cleanup
- ✅ **Canvas optimization**: Device pixel ratio handling
- ✅ **File processing**: Chunked and asynchronous operations

### **User Experience**

- ✅ **Responsive design**: Works on desktop and mobile
- ✅ **Theme support**: Light and dark modes
- ✅ **Progress feedback**: Real-time status updates
- ✅ **Error recovery**: Graceful handling of failures

## 🚀 Usage

### **Development**

```bash
# Serve the application locally
python -m http.server 8000
# Navigate to: http://localhost:8000/browser/
```

### **Files**

- Use `index.html` for the modular version
- `browser_page.html` contains the original monolithic version
- All modules are ES6 compatible and use modern JavaScript features

### **Browser Support**

- Chrome/Edge 90+
- Firefox 88+
- Safari 14+
- Requires ES6 module support

## 🔧 Customization

### **Adding New Pipeline Stages**

1. Add stage logic to `pipeline.js`
2. Update UI in `ui-manager.js`
3. Add export support in `export-manager.js`
4. Update constants in `core.js`

### **Modifying Analysis Parameters**

- Edit `CONSTANTS` object in `core.js`
- Quality thresholds, damage windows, bootstrap samples

### **Styling Changes**

- Modify `css/main.css`
- CSS custom properties for easy theming
- Component-specific selectors

## 📊 Data Flow

``` bash
AB1 Files → ABIF Parser → Sequence Processor → Pipeline Controller
    ↓
UI Manager ← Export Manager ← Chromatogram Renderer ← Data Store
```

1. **File Upload**: User selects .ab1 files
2. **Parsing**: ABIF parser extracts sequences and traces
3. **Processing**: Quality filtering, alignment, damage analysis
4. **Visualization**: Chromatogram rendering and statistics
5. **Export**: ZIP archives with results

This modular architecture provides a solid foundation for extending the application with additional analysis methods, visualization options, or integration with external tools.
