/**
 * Main application initialization and coordination
 */

import { Theme, LogoManager, store } from './core.js';
import ChromatogramRenderer from './chromatogram.js';
import UIManager from './ui-manager.js';
import Pipeline from './pipeline.js';
import ExportManager from './export-manager.js';

export class SangerApp {
  constructor() {
    this.components = {};
    this.initialized = false;
  }

  async init() {
    if (this.initialized) return;

    try {
      console.log('🚀 Initializing Sanger aDNA Damage Application...');

      // Initialize theme
      Theme.init();
      console.log('✅ Theme initialized');

      // Initialize core components
      this.components.chromatogram = new ChromatogramRenderer('chrom');
      console.log('✅ Chromatogram renderer initialized');
      
      this.components.uiManager = new UIManager();
      console.log('✅ UI Manager initialized');
      
      this.components.exportManager = new ExportManager();
      console.log('✅ Export Manager initialized');
      
      // Initialize pipeline with dependencies
      this.components.pipeline = new Pipeline(
        this.components.chromatogram,
        this.components.uiManager
      );
      console.log('✅ Pipeline initialized');

      // Setup UI connections
      this.setupUIConnections();
      console.log('✅ UI connections established');

      // Load logos
      await LogoManager.inject();
      console.log('✅ Logos loaded');

      // Initialize stats display
      this.components.uiManager.toggleStatsEmpty();
      
      // Initialize stage counters
      this.components.uiManager.updateStageCounters();

      // Setup paste support for files
      this.setupPasteSupport();

      this.initialized = true;
      console.log('✅ Application initialized successfully');

      // Set global reference for modal interactions
      window.app = this;

    } catch (error) {
      console.error('Failed to initialize application:', error);
      this.showError('Failed to initialize application: ' + error.message);
    }
  }

  // Helper methods for modal interactions
  showChromatogram(key) {
    console.log('Showing chromatogram for:', key);
    const traceSelect = document.getElementById('traceSelect');
    if (traceSelect) {
      traceSelect.value = key;
      traceSelect.dispatchEvent(new Event('change'));
    }
    
    // Close modal
    if (this.components.uiManager) {
      this.components.uiManager.hideFileModal();
    }
  }

  testAllTriggers() {
    console.log('🧪 Testing all pipeline triggers...');
    const stages = ['fasta', 'filter', 'align', 'consensus', 'damage'];
    
    stages.forEach(stage => {
      const button = document.querySelector(`[data-stage="${stage}"]`);
      if (button) {
        console.log(`✓ Found trigger for ${stage}:`, button);
      } else {
        console.warn(`✗ Missing trigger for ${stage}`);
      }
    });
    
    const exportButtons = document.querySelectorAll('[data-exp]');
    console.log(`✓ Found ${exportButtons.length} export buttons:`, 
      Array.from(exportButtons).map(b => b.getAttribute('data-exp')));
      
    alert(`Trigger Test Complete!\nCheck console for details.\n\nFound:\n- ${stages.filter(s => document.querySelector(`[data-stage="${s}"]`)).length}/${stages.length} stage triggers\n- ${exportButtons.length} export buttons`);
  }

  validatePipeline() {
    console.log('🔍 Validating pipeline state...');
    
    const stages = ['fasta', 'filtered', 'aligned', 'consensus', 'damage_analysis'];
    const results = stages.map(stage => ({
      stage,
      count: store[stage]?.size || 0,
      hasData: (store[stage]?.size || 0) > 0
    }));
    
    console.table(results);
    
    const summary = results.map(r => `${r.stage}: ${r.count} items`).join('\n');
    alert(`Pipeline Validation:\n\n${summary}\n\nSee console for detailed table.`);
  }

  setupUIConnections() {
    const { uiManager, chromatogram, pipeline } = this.components;

    // Setup file input
    uiManager.setupFileInput('fileInput', pipeline);

    // Setup pipeline control buttons
    uiManager.setupPipelineButtons(pipeline);

    // Export buttons are automatically initialized in ExportManager constructor

    // Setup chromatogram PNG download
    chromatogram.setupPNGDownload('dlPNG');

    console.log('UI connections established');
  }

  showError(message) {
    // Create or update error display
    let errorDiv = document.getElementById('app-error');
    if (!errorDiv) {
      errorDiv = document.createElement('div');
      errorDiv.id = 'app-error';
      errorDiv.style.cssText = `
        position: fixed;
        top: 20px;
        right: 20px;
        background: #fee;
        color: #c33;
        padding: 16px;
        border-radius: 8px;
        border: 1px solid #fcc;
        max-width: 400px;
        z-index: 10000;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
      `;
      document.body.appendChild(errorDiv);
    }

    errorDiv.innerHTML = `
      <strong>Error:</strong> ${message}
      <button onclick="this.parentElement.remove()" style="
        float: right;
        background: none;
        border: none;
        color: #c33;
        cursor: pointer;
        font-size: 18px;
        padding: 0;
        margin-left: 10px;
      ">&times;</button>
    `;

    // Auto-remove after 10 seconds
    setTimeout(() => {
      if (errorDiv.parentElement) {
        errorDiv.remove();
      }
    }, 10000);
  }

  // Public API methods for external access
  getComponent(name) {
    return this.components[name];
  }

  async processFiles(files) {
    if (!this.initialized) {
      await this.init();
    }
    return this.components.pipeline.processFiles(files);
  }

  async runStage(stage) {
    if (!this.initialized) {
      await this.init();
    }
    return this.components.pipeline.runStage(stage);
  }

  // Test the file loading system
  async testFileLoading() {
    console.log('🧪 Testing file loading system...');
    
    // Test if file input exists and is properly connected
    const fileInput = document.getElementById('fileInput');
    if (!fileInput) {
      console.error('❌ File input element not found');
      return;
    }
    console.log('✅ File input element found');
    
    // Test if pipeline is initialized
    if (!this.components.pipeline) {
      console.error('❌ Pipeline not initialized');
      return;
    }
    console.log('✅ Pipeline initialized');
    
    // Test ABIF parser
    try {
      const { ABIFParser } = await import('./abif-parser.js');
      console.log('✅ ABIF parser accessible');
    } catch (error) {
      console.error('❌ ABIF parser import failed:', error);
    }
    
    console.log('📂 File loading test complete - check results above');
  }

  // Direct method to load AB1 files from input directory
  async loadSampleFiles() {
    console.log('🔄 Attempting to load sample AB1 files...');
    
    // Create a manual file selection interface
    const modal = document.createElement('div');
    modal.style.cssText = `
      position: fixed; top: 0; left: 0; width: 100%; height: 100%;
      background: rgba(0,0,0,0.8); display: flex; align-items: center;
      justify-content: center; z-index: 2000;
    `;
    
    const sampleFiles = [
      '10_CON_TA_PB_B2_B1_HVS1-F.ab1',
      '1_SJT_P1_PB_A1_A1_HVS1-F.ab1', 
      '2_SJC_FO_PB_B_B1_HVS1-F.ab1',
      'C_A8_A1_HVS1-F.ab1'
    ];
    
    modal.innerHTML = `
      <div style="
        background: var(--bg); border-radius: 12px; padding: 24px;
        max-width: 600px; max-height: 80%; overflow-y: auto;
      ">
        <h3 style="margin: 0 0 16px 0; color: var(--ink);">Load Sample AB1 Files</h3>
        <p style="color: var(--muted); margin-bottom: 16px;">
          Select sample AB1 files to load for testing:
        </p>
        <div id="sampleFilesList"></div>
        <div style="margin-top: 16px; text-align: right;">
          <button id="cancelLoad" style="margin-right: 8px; padding: 8px 16px; background: var(--muted); color: white; border: none; border-radius: 6px; cursor: pointer;">Cancel</button>
          <button id="loadSelected" style="padding: 8px 16px; background: var(--accent); color: white; border: none; border-radius: 6px; cursor: pointer;">Load Selected</button>
        </div>
      </div>
    `;
    
    // Add file checkboxes
    const filesList = modal.querySelector('#sampleFilesList');
    sampleFiles.forEach(fileName => {
      const div = document.createElement('div');
      div.style.cssText = 'margin-bottom: 8px; padding: 8px; border: 1px solid var(--panel); border-radius: 6px;';
      div.innerHTML = `
        <label style="display: flex; align-items: center; cursor: pointer;">
          <input type="checkbox" value="${fileName}" style="margin-right: 8px;">
          <span style="color: var(--ink);">${fileName}</span>
        </label>
      `;
      filesList.appendChild(div);
    });
    
    document.body.appendChild(modal);
    
    // Handle close
    modal.querySelector('#cancelLoad').onclick = () => modal.remove();
    
    // Handle load
    modal.querySelector('#loadSelected').onclick = async () => {
      const selected = Array.from(modal.querySelectorAll('input[type="checkbox"]:checked'))
        .map(cb => cb.value);
      
      modal.remove();
      
      if (selected.length > 0) {
        console.log('📂 Loading selected files:', selected);
        await this.loadFilesFromPaths(selected);
      }
    };
  }

  async loadFilesFromPaths(filePaths) {
    console.log('🔄 Attempting to load files from paths:', filePaths);
    
    for (const filePath of filePaths) {
      try {
        const response = await fetch(`../input/${filePath}`);
        if (response.ok) {
          const arrayBuffer = await response.arrayBuffer();
          const file = new File([arrayBuffer], filePath, { type: 'application/octet-stream' });
          
          console.log(`✅ Loaded ${filePath}:`, file.size, 'bytes');
          await this.components.pipeline.processSingleFile(file);
        } else {
          console.error(`❌ Failed to fetch ${filePath}:`, response.status);
        }
      } catch (error) {
        console.error(`❌ Error loading ${filePath}:`, error);
      }
    }
    
    // Update UI after loading
    if (this.components.uiManager) {
      this.components.uiManager.updateStageCounters();
      this.components.uiManager.toggleStatsEmpty();
    }
  }

  // Alternative file loading using File System Access API (Chrome/Edge)
  async openFileSystemDialog() {
    console.log('🗂️ Opening File System Access API dialog...');
    
    try {
      // Check if File System Access API is supported
      if (!window.showOpenFilePicker) {
        this.showError('File System Access API not supported in this browser. Try Chrome or Edge.');
        return;
      }

      const fileHandles = await window.showOpenFilePicker({
        types: [{
          description: 'AB1 Files',
          accept: {
            'application/octet-stream': ['.ab1']
          }
        }],
        multiple: true
      });

      console.log(`Selected ${fileHandles.length} files via File System API`);

      for (const fileHandle of fileHandles) {
        const file = await fileHandle.getFile();
        console.log(`Processing ${file.name} (${file.size} bytes)`);
        await this.components.pipeline.processSingleFile(file);
      }

      // Update UI
      if (this.components.uiManager) {
        this.components.uiManager.updateStageCounters();
        this.components.uiManager.toggleStatsEmpty();
      }

    } catch (error) {
      if (error.name !== 'AbortError') {
        console.error('File System API error:', error);
        this.showError('Failed to open files: ' + error.message);
      }
    }
  }

  // Load demo data with synthetic sequences
  async loadDemoData() {
    console.log('🎯 Loading demo data...');
    
    const demoSequences = [
      {
        name: 'DEMO_Sample1_HVS1-F.ab1',
        sequence: 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC',
        qualities: Array(60).fill().map(() => Math.floor(Math.random() * 40) + 10)
      },
      {
        name: 'DEMO_Sample1_HVS1-R.ab1', 
        sequence: 'GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT',
        qualities: Array(60).fill().map(() => Math.floor(Math.random() * 40) + 10)
      },
      {
        name: 'DEMO_Sample2_HVS2-F.ab1',
        sequence: 'TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGG',
        qualities: Array(60).fill().map(() => Math.floor(Math.random() * 40) + 10)
      }
    ];

    for (const demo of demoSequences) {
      try {
        // Create a synthetic ABIF-like structure
        const demoFile = {
          name: demo.name,
          size: demo.sequence.length,
          arrayBuffer: () => Promise.resolve(new ArrayBuffer(100)) // Dummy buffer
        };

        // Mock the ABIF parsing for demo data
        const originalParseSingle = this.components.pipeline.processSingleFile;
        this.components.pipeline.processSingleFile = async function(file) {
          if (file.name.startsWith('DEMO_')) {
            console.log('Processing demo file:', file.name);
            
            const info = Utils.parseFileName(file.name);
            if (!info) {
              console.warn('Invalid demo filename format:', file.name);
              return;
            }

            const key = `${info.sample}|${info.region}|${info.dir}`;
            const demoData = demoSequences.find(d => d.name === file.name);

            // Store sequences directly
            store.fasta.set(key, demoData.sequence);
            store.filtered.set(key, demoData.sequence);

            // Update UI
            this.updateTraceSelector(key, info, { 
              bases: demoData.sequence,
              quals: demoData.qualities,
              traces: null 
            }, demoData.sequence);
            
            this.addStatsRow(info, demoData.sequence, demoData.qualities);
            
            if (this.uiManager) {
              this.uiManager.updateStageCounters();
            }

            console.log('Demo file processed:', file.name, 'Key:', key);
          } else {
            return originalParseSingle.call(this, file);
          }
        };

        await this.components.pipeline.processSingleFile(demoFile);

      } catch (error) {
        console.error('Error loading demo data:', error);
      }
    }

    // Restore original method
    setTimeout(() => {
      delete this.components.pipeline.processSingleFile.isDemo;
    }, 100);

    console.log('✅ Demo data loaded successfully');
  }

  // Paste functionality for file data
  async setupPasteSupport() {
    document.addEventListener('paste', async (e) => {
      console.log('📋 Paste event detected');
      
      const items = e.clipboardData?.items;
      if (!items) return;

      for (let item of items) {
        if (item.kind === 'file') {
          const file = item.getAsFile();
          if (file && file.name.toLowerCase().endsWith('.ab1')) {
            console.log('📋 Processing pasted AB1 file:', file.name);
            await this.components.pipeline.processSingleFile(file);
            
            // Update UI
            if (this.components.uiManager) {
              this.components.uiManager.updateStageCounters();
              this.components.uiManager.toggleStatsEmpty();
            }
          }
        }
      }
    });
    
    console.log('📋 Paste support initialized');
  }

  reset() {
    if (this.components.uiManager) {
      this.components.uiManager.resetApplication();
    }
  }
}

// Global instance
let appInstance = null;

export function getApp() {
  if (!appInstance) {
    appInstance = new SangerApp();
  }
  return appInstance;
}

// Auto-initialize when DOM is ready
if (typeof document !== 'undefined') {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => {
      getApp().init();
    });
  } else {
    // DOM already loaded
    getApp().init();
  }
}

export default SangerApp;
