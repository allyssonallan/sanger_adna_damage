/**
 * UI components and interface management
 */

import { store, Utils, COUNTERS } from './core.js';

export class UIManager {
  constructor() {
    this.statsTable = document.getElementById('statsTable');
    this.statsBody = this.statsTable?.querySelector('tbody');
    this.statsEmpty = document.getElementById('statsEmpty');
    this.modal = document.getElementById('fileSummaryModal');
    
    this.initializeEventHandlers();
  }

  initializeEventHandlers() {
    // View All Files button
    const viewAllFilesBtn = document.getElementById('viewAllFiles');
    if (viewAllFilesBtn) {
      viewAllFilesBtn.onclick = () => this.showFileModal();
    }

    // Modal close functionality
    if (this.modal) {
      const closeBtn = this.modal.querySelector('.close');
      if (closeBtn) {
        closeBtn.onclick = () => this.hideFileModal();
      }
      
      // Close on outside click
      this.modal.onclick = (e) => {
        if (e.target === this.modal) {
          this.hideFileModal();
        }
      };
    }

    // Reset button
    const resetBtn = document.getElementById('reset');
    if (resetBtn) {
      resetBtn.onclick = () => this.resetApplication();
    }

    // Trace selector
    const traceSelect = document.getElementById('traceSelect');
    if (traceSelect) {
      traceSelect.addEventListener('change', () => {
        const key = traceSelect.value;
        console.log('Trace selector changed to:', key);
        
        // Update sequence display
        const filtered = store.filtered.get(key);
        if (filtered) {
          const seqText = document.getElementById('seqText');
          if (seqText) {
            seqText.textContent = filtered;  // filtered is now just the sequence
          }
        }
        
        // Update chromatogram if trace data is available
        this.updateChromatogramForKey(key);
      });
    }

    // List uploaded files button
    const listUploadedBtn = document.getElementById('listUploaded');
    if (listUploadedBtn) {
      listUploadedBtn.onclick = () => this.showUploadedFilesModal();
    }

    // List temp files button
    const listTempFilesBtn = document.getElementById('listTempFiles');
    if (listTempFilesBtn) {
      listTempFilesBtn.onclick = () => this.showTempFilesModal();
    }

    // Test file loading button
    const testFileLoadingBtn = document.getElementById('testFileLoading');
    if (testFileLoadingBtn) {
      testFileLoadingBtn.onclick = () => {
        if (window.app) {
          window.app.testFileLoading();
        }
      };
    }

    // Load sample files button
    const loadSampleFilesBtn = document.getElementById('loadSampleFiles');
    if (loadSampleFilesBtn) {
      loadSampleFilesBtn.onclick = () => {
        if (window.app) {
          window.app.loadSampleFiles();
        }
      };
    }

    // Load demo data button
    const loadDemoDataBtn = document.getElementById('loadDemoData');
    if (loadDemoDataBtn) {
      loadDemoDataBtn.onclick = () => {
        if (window.app) {
          window.app.loadDemoData();
        }
      };
    }

    // File System API button
    const openFileSystemBtn = document.getElementById('openFileSystemDialog');
    if (openFileSystemBtn) {
      openFileSystemBtn.onclick = () => {
        if (window.app) {
          window.app.openFileSystemDialog();
        }
      };
    }
  }

  updateChromatogramForKey(key) {
    // Check if we have stored trace data for this key
    if (window.chromatogramData && window.chromatogramData[key]) {
      console.log('Rendering stored chromatogram for:', key);
      const chromatogramRenderer = document.querySelector('canvas#chrom')?._renderer;
      if (chromatogramRenderer) {
        chromatogramRenderer.draw(window.chromatogramData[key]);
      }
    } else {
      console.warn('No chromatogram data found for key:', key);
    }
  }

  clearStats() {
    if (this.statsBody) {
      this.statsBody.innerHTML = '';
    }
  }

  updateStageCounters() {
    // Update all stage counters based on store data
    const stages = [
      { id: 'fastaCount', key: 'fasta', suffix: 'files' },
      { id: 'filteredCount', key: 'filtered', suffix: 'files' },
      { id: 'alignedCount', key: 'aligned', suffix: '' },
      { id: 'consensusCount', key: 'consensus', suffix: '' },
      { id: 'finalCount', key: 'final', suffix: 'files' },
      { id: 'plotsCount', key: 'plots', suffix: '' },
      { id: 'damageCount', key: 'damage_analysis', suffix: '' }
    ];
    
    stages.forEach(stage => {
      const element = document.getElementById(stage.id);
      const count = store[stage.key]?.size || 0;
      if (element) {
        element.textContent = `${count} ${stage.suffix}`;
        
        // Add visual status indicator
        const parentStage = element.closest('.stage');
        if (parentStage) {
          parentStage.classList.toggle('has-data', count > 0);
          parentStage.classList.toggle('no-data', count === 0);
        }
      }
    });
    
    console.log('Stage counters updated:', {
      fasta: store.fasta?.size || 0,
      filtered: store.filtered?.size || 0,
      aligned: store.aligned?.size || 0,
      consensus: store.consensus?.size || 0,
      damage: store.damage_analysis?.size || 0
    });
  }

  addStatsRow(data) {
    if (!this.statsBody) return;

    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${data.sample}</td>
      <td>${data.region}</td>
      <td>${data.dir}</td>
      <td>${data.len}</td>
      <td>${data.Ns}</td>
      <td>${data.q30}${typeof data.q30 === 'number' ? '%' : ''}</td>
      <td>${data.rate5}</td>
      <td>${data.rate3}</td>
    `;
    this.statsBody.appendChild(tr);
  }

  updateDamageRates(sample, region, damage) {
    if (!this.statsBody) return;

    const rows = [...this.statsBody.querySelectorAll('tr')];
    const consensusRow = rows.find(row => {
      const cells = row.querySelectorAll('td');
      return cells[0]?.textContent === sample && 
             cells[1]?.textContent === region && 
             cells[2]?.textContent === 'C';
    });

    if (consensusRow) {
      const cells = consensusRow.querySelectorAll('td');
      cells[6].textContent = Utils.formatPercentage(damage.rate5, damage.ci5);
      cells[7].textContent = Utils.formatPercentage(damage.rate3, damage.ci3);
    }
  }

  toggleStatsEmpty() {
    if (!this.statsBody || !this.statsEmpty || !this.statsTable) return;

    const hasData = this.statsBody.children.length > 0;
    
    this.statsTable.style.display = hasData ? '' : 'none';
    this.statsEmpty.style.display = hasData ? 'none' : 'block';
  }

  showFileModal() {
    if (!this.modal) return;
    
    this.updateFilesList();
    this.modal.style.display = 'block';
  }

  hideFileModal() {
    if (!this.modal) return;
    this.modal.style.display = 'none';
  }

  updateFilesList() {
    const filesList = document.getElementById('loadedFilesList');
    if (!filesList) return;

    filesList.innerHTML = '';

    console.log('DEBUG: updateFilesList called');
    console.log('DEBUG: store.fasta size:', store.fasta.size);
    console.log('DEBUG: store.fasta contents:', [...store.fasta.entries()]);

    // Collect all loaded files with validation
    const validFiles = [];
    const invalidFiles = [];
    
    for (const [key, sequence] of store.fasta) {
      console.log('DEBUG: Processing key:', key, 'sequence length:', sequence?.length || 0);
      const [sample, region, dir] = key.split('|');
      const length = sequence ? sequence.length : 0;
      const isValid = length > 0;
      
      const fileInfo = { 
        sample, 
        region, 
        dir, 
        key, 
        length,
        sequence,
        hasFiltered: store.filtered.has(key),
        hasAligned: store.aligned.has(key),
        hasDamage: store.damage_analysis.has(key)
      };
      
      console.log('DEBUG: fileInfo:', fileInfo);
      
      if (isValid) {
        validFiles.push(fileInfo);
      } else {
        invalidFiles.push(fileInfo);
      }
    }

    if (validFiles.length === 0 && invalidFiles.length === 0) {
      filesList.innerHTML = '<p style="text-align: center; color: var(--muted); padding: 20px;">No files loaded yet.</p>';
      return;
    }

    // Add header with summary
    const summaryDiv = document.createElement('div');
    summaryDiv.style.cssText = 'margin-bottom: 16px; padding: 12px; background: var(--panel); border-radius: 8px; border: 1px solid #e5e7eb;';
    summaryDiv.innerHTML = `
      <h4 style="margin: 0 0 8px 0; color: var(--ink);">Files Summary</h4>
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr)); gap: 8px; font-size: 13px;">
        <div>✅ Valid: <strong>${validFiles.length}</strong></div>
        <div>❌ Invalid: <strong>${invalidFiles.length}</strong></div>
        <div>🔄 Filtered: <strong>${validFiles.filter(f => f.hasFiltered).length}</strong></div>
        <div>📊 Damage: <strong>${validFiles.filter(f => f.hasDamage).length}</strong></div>
      </div>
      <div style="margin-top: 10px;">
        <button class="secondary" style="font-size: 12px; margin-right: 8px;" onclick="window.app?.testAllTriggers?.()">🧪 Test All Triggers</button>
        <button class="secondary" style="font-size: 12px;" onclick="window.app?.validatePipeline?.()">🔍 Validate Pipeline</button>
      </div>
    `;
    filesList.appendChild(summaryDiv);

    // Valid files section
    if (validFiles.length > 0) {
      const validHeader = document.createElement('h4');
      validHeader.textContent = `✅ Valid Files (${validFiles.length})`;
      validHeader.style.cssText = 'color: var(--good); margin: 16px 0 8px 0;';
      filesList.appendChild(validHeader);

      validFiles.forEach(file => {
        const fileDiv = document.createElement('div');
        fileDiv.style.cssText = `
          padding: 12px;
          border: 1px solid #d1fae5;
          border-radius: 6px;
          margin-bottom: 8px;
          background: #f0fdf4;
        `;
        
        const statusIcons = [
          file.hasFiltered ? '🔄' : '⚪',
          file.hasAligned ? '🔗' : '⚪', 
          file.hasDamage ? '📊' : '⚪'
        ].join(' ');
        
        fileDiv.innerHTML = `
          <div style="display: flex; justify-content: space-between; align-items: center;">
            <div>
              <div style="font-weight: 600; color: var(--ink);">${file.sample}_${file.region}-${file.dir}</div>
              <div style="font-size: 12px; color: var(--muted); margin-top: 2px;">
                Length: ${file.length} bp | Status: ${statusIcons}
              </div>
            </div>
            <button class="ghost" style="font-size: 11px; padding: 4px 8px;" onclick="window.app?.showChromatogram?.('${file.key}')">View</button>
          </div>
        `;
        
        filesList.appendChild(fileDiv);
      });
    }

    // Invalid files section  
    if (invalidFiles.length > 0) {
      const invalidHeader = document.createElement('h4');
      invalidHeader.textContent = `❌ Invalid Files (${invalidFiles.length})`;
      invalidHeader.style.cssText = 'color: var(--bad); margin: 16px 0 8px 0;';
      filesList.appendChild(invalidHeader);

      invalidFiles.forEach(file => {
        const fileDiv = document.createElement('div');
        fileDiv.style.cssText = `
          padding: 12px;
          border: 1px solid #fecaca;
          border-radius: 6px;
          margin-bottom: 8px;
          background: #fef2f2;
        `;
        
        fileDiv.innerHTML = `
          <div style="font-weight: 600; color: var(--bad);">${file.sample}_${file.region}-${file.dir}</div>
          <div style="font-size: 12px; color: var(--muted); margin-top: 2px;">
            Issue: ${file.length === 0 ? 'Zero length sequence' : 'Unknown parsing error'}
          </div>
        `;
        
        filesList.appendChild(fileDiv);
      });
    }

    // Update counters
    const fileCount = document.getElementById('fileCount');
    const validCount = document.getElementById('validCount'); 
    const processedCount = document.getElementById('processedCount');
    
    if (fileCount) fileCount.textContent = validFiles.length + invalidFiles.length;
    if (validCount) validCount.textContent = validFiles.length;
    if (processedCount) processedCount.textContent = validFiles.filter(f => f.hasFiltered).length;
  }

  resetApplication() {
    // Clear all stores
    for (const key of Object.keys(store)) {
      if (store[key] instanceof Map) {
        store[key].clear();
      } else if (Array.isArray(store[key])) {
        store[key].length = 0;
      } else {
        store[key] = typeof store[key] === 'string' ? '' : null;
      }
    }

    // Reset counters
    Object.keys(COUNTERS).forEach(key => {
      COUNTERS[key] = 0;
    });

    // Clear UI elements
    this.clearStats();
    
    const traceSelect = document.getElementById('traceSelect');
    if (traceSelect) traceSelect.innerHTML = '';
    
    const seqText = document.getElementById('seqText');
    if (seqText) seqText.textContent = '';

    // Clear chromatogram data
    if (window.chromatogramData) {
      window.chromatogramData = {};
    }
    
    // Clear canvas
    const canvas = document.getElementById('chrom');
    if (canvas) {
      const ctx = canvas.getContext('2d');
      if (ctx) {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
      }
    }

    // Update displays
    Utils.updateCounts();
    this.toggleStatsEmpty();
    this.updateStageCounters();
    
    console.log('Application reset completed');
  }

  showUploadedFilesModal() {
    const uploadedFiles = Array.from(store.fasta.keys());
    
    const modal = this.createModal('Uploaded AB1 Files', uploadedFiles.length === 0 
      ? '<p style="color: var(--muted); text-align: center; padding: 20px;">No files uploaded yet</p>'
      : this.generateUploadedFilesList(uploadedFiles)
    );
    
    document.body.appendChild(modal);
  }

  showTempFilesModal() {
    const tempFiles = this.collectTempFiles();
    
    const modal = this.createModal('Temporary Processing Files', tempFiles.length === 0 
      ? '<p style="color: var(--muted); text-align: center; padding: 20px;">No temporary files generated yet</p>'
      : this.generateTempFilesList(tempFiles)
    );
    
    document.body.appendChild(modal);
  }

  createModal(title, content) {
    const modal = document.createElement('div');
    modal.className = 'modal';
    modal.style.cssText = `
      position: fixed; top: 0; left: 0; width: 100%; height: 100%; 
      background: rgba(0,0,0,0.5); display: flex; align-items: center; 
      justify-content: center; z-index: 1000;
    `;
    
    modal.innerHTML = `
      <div class="modal-content" style="
        background: var(--bg); border-radius: 12px; padding: 24px; 
        max-width: 80%; max-height: 80%; overflow-y: auto;
        box-shadow: 0 20px 40px rgba(0,0,0,0.15);
        min-width: 500px;
      ">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px;">
          <h3 style="margin: 0; color: var(--ink);">${title}</h3>
          <button class="close" style="
            background: none; border: none; font-size: 24px; cursor: pointer;
            color: var(--muted); padding: 0; width: 32px; height: 32px;
            display: flex; align-items: center; justify-content: center;
            border-radius: 6px;
          ">×</button>
        </div>
        <div>${content}</div>
      </div>
    `;
    
    // Close handlers
    const closeBtn = modal.querySelector('.close');
    closeBtn.onclick = () => modal.remove();
    modal.onclick = (e) => {
      if (e.target === modal) modal.remove();
    };
    
    return modal;
  }

  generateUploadedFilesList(files) {
    let html = `<div style="margin-bottom: 12px; color: var(--muted);">
      <strong>${files.length}</strong> uploaded AB1 files:
    </div>`;
    
    files.forEach(file => {
      const parsed = Utils.parseFileName(file);
      const seqLength = store.fasta.get(file)?.length || 0;
      
      html += `
        <div style="
          padding: 12px; border: 1px solid var(--panel); border-radius: 8px; 
          margin-bottom: 8px; background: var(--panel);
        ">
          <div style="font-weight: 600; color: var(--ink);">${file}</div>
          ${parsed ? `
            <div style="font-size: 13px; color: var(--muted); margin-top: 4px;">
              Sample: ${parsed.sample} | Region: ${parsed.region} | Direction: ${parsed.dir}
            </div>
          ` : ''}
          <div style="font-size: 13px; color: var(--muted); margin-top: 4px;">
            Sequence length: ${seqLength} bp
          </div>
        </div>
      `;
    });
    
    return html;
  }

  generateTempFilesList(tempFiles) {
    let html = `<div style="margin-bottom: 12px; color: var(--muted);">
      <strong>${tempFiles.length}</strong> temporary processing files:
    </div>`;
    
    tempFiles.forEach(file => {
      html += `
        <div style="
          padding: 12px; border: 1px solid var(--panel); border-radius: 8px; 
          margin-bottom: 8px; background: var(--panel);
        ">
          <div style="font-weight: 600; color: var(--ink);">${file.name}</div>
          <div style="font-size: 13px; color: var(--muted); margin-top: 4px;">
            Stage: ${file.stage} | Type: ${file.type} | Size: ${file.size}
          </div>
        </div>
      `;
    });
    
    return html;
  }

  collectTempFiles() {
    const tempFiles = [];
    
    // Collect from all processing stages
    const stages = [
      { stage: 'Filtered', store: store.filtered, type: 'FASTA sequence' },
      { stage: 'Aligned', store: store.aligned, type: 'Alignment result' },
      { stage: 'Consensus', store: store.consensus, type: 'Consensus sequence' },
      { stage: 'Final', store: store.final, type: 'Final sequence' },
      { stage: 'Damage Analysis', store: store.damage_analysis, type: 'Analysis data' },
      { stage: 'Plots', store: store.plots, type: 'Plot data' }
    ];
    
    stages.forEach(({ stage, store: stageStore, type }) => {
      if (stageStore instanceof Map) {
        for (const [key, value] of stageStore.entries()) {
          tempFiles.push({
            name: key,
            stage: stage,
            type: type,
            size: typeof value === 'string' ? `${value.length} chars` : 'Data object'
          });
        }
      }
    });
    
    return tempFiles;
  }

  setupFileInput(inputId, pipeline) {
    const fileInput = document.getElementById(inputId);
    if (!fileInput) {
      console.error('File input not found:', inputId);
      return;
    }

    console.log('Setting up file input:', inputId);

    fileInput.addEventListener('change', async (e) => {
      console.log('File input change event triggered');
      const files = [...e.target.files];
      console.log('Selected files:', files.length, files.map(f => f.name));
      
      if (files.length > 0 && pipeline) {
        console.log('Processing files with pipeline...');
        await pipeline.processFiles(files);
      } else {
        if (files.length === 0) {
          console.log('No files selected');
        }
        if (!pipeline) {
          console.error('Pipeline not available');
        }
      }
    });
    
    console.log('File input event listener attached');
    
    // Setup drag and drop functionality
    this.setupDragAndDrop(pipeline);
  }

  setupDragAndDrop(pipeline) {
    const dragDropZone = document.getElementById('dragDropZone');
    if (!dragDropZone) return;

    console.log('Setting up drag and drop functionality');

    // Prevent default drag behaviors
    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
      dragDropZone.addEventListener(eventName, (e) => {
        e.preventDefault();
        e.stopPropagation();
      });
    });

    // Highlight drop zone when item is dragged over it
    ['dragenter', 'dragover'].forEach(eventName => {
      dragDropZone.addEventListener(eventName, () => {
        dragDropZone.classList.add('drag-over');
      });
    });

    ['dragleave', 'drop'].forEach(eventName => {
      dragDropZone.addEventListener(eventName, () => {
        dragDropZone.classList.remove('drag-over');
      });
    });

    // Handle dropped files
    dragDropZone.addEventListener('drop', async (e) => {
      console.log('Files dropped');
      const files = [...e.dataTransfer.files];
      console.log('Dropped files:', files.length, files.map(f => f.name));
      
      // Filter for AB1 files
      const ab1Files = files.filter(file => 
        file.name.toLowerCase().endsWith('.ab1') || 
        file.type === 'application/octet-stream'
      );
      
      if (ab1Files.length === 0) {
        this.showMessage('No AB1 files found in dropped files', 'warning');
        return;
      }
      
      if (ab1Files.length !== files.length) {
        this.showMessage(`Found ${ab1Files.length} AB1 files out of ${files.length} dropped files`, 'info');
      }
      
      // Process the AB1 files
      if (pipeline) {
        await pipeline.processFiles(ab1Files);
      } else {
        console.error('Pipeline not available for drag and drop');
      }
    });

    // Click to open file dialog
    dragDropZone.addEventListener('click', () => {
      const fileInput = document.getElementById('fileInput');
      if (fileInput) {
        fileInput.click();
      }
    });
  }

  showMessage(message, type = 'info') {
    const messageDiv = document.createElement('div');
    messageDiv.style.cssText = `
      position: fixed;
      top: 20px;
      right: 20px;
      padding: 12px 16px;
      border-radius: 8px;
      color: white;
      font-weight: 600;
      z-index: 1000;
      max-width: 300px;
      background: ${type === 'warning' ? '#f59e0b' : type === 'error' ? '#ef4444' : '#3b82f6'};
    `;
    messageDiv.textContent = message;
    
    document.body.appendChild(messageDiv);
    
    setTimeout(() => {
      messageDiv.remove();
    }, 4000);
  }

  setupPipelineButtons(pipeline) {
    // Individual stage buttons
    const stageButtons = document.querySelectorAll('[data-stage]');
    stageButtons.forEach(button => {
      button.addEventListener('click', async () => {
        const stage = button.getAttribute('data-stage');
        if (pipeline) {
          // Show loading state
          const originalText = button.textContent;
          button.textContent = '⏳ Processing...';
          button.disabled = true;
          
          try {
            await pipeline.runStage(stage);
            
            // Show success state
            button.textContent = '✅ Complete';
            setTimeout(() => {
              button.textContent = originalText;
              button.disabled = false;
            }, 2000);
            
          } catch (error) {
            console.error('Stage error:', error);
            button.textContent = '❌ Error';
            setTimeout(() => {
              button.textContent = originalText;
              button.disabled = false;
            }, 3000);
          }
        }
      });
    });

    // Run all button
    const runAllBtn = document.getElementById('runAll');
    if (runAllBtn && pipeline) {
      runAllBtn.addEventListener('click', async () => {
        const originalText = runAllBtn.textContent;
        runAllBtn.textContent = '⏳ Running All...';
        runAllBtn.disabled = true;
        
        try {
          await pipeline.runStage('align');
          await pipeline.runStage('consensus');
          await pipeline.runStage('damage');
          
          runAllBtn.textContent = '✅ All Complete';
          setTimeout(() => {
            runAllBtn.textContent = originalText;
            runAllBtn.disabled = false;
          }, 2000);
          
        } catch (error) {
          console.error('Run all error:', error);
          runAllBtn.textContent = '❌ Error';
          setTimeout(() => {
            runAllBtn.textContent = originalText;
            runAllBtn.disabled = false;
          }, 3000);
        }
      });
    }
  }
}

export default UIManager;
