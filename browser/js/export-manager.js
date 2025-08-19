/**
 * Export functionality for downloading results
 */

import { store } from './core.js';

export class ExportManager {
  constructor() {
    this.initializeExportButtons();
  }

  initializeExportButtons() {
    // Generic export buttons
    document.querySelectorAll('button[data-exp]').forEach(btn => {
      btn.addEventListener('click', () => {
        try {
          console.log('Export button clicked for:', btn.getAttribute('data-exp'));
          
          const key = btn.getAttribute('data-exp');
          const data = store[key];
          
          console.log('Store data for key:', key, data);
          
          if (!(data instanceof Map) || data.size === 0) {
            console.log('No data to export for key:', key);
            alert("Nothing to export for " + key);
            return;
          }
          
          this.exportDataAsZip(key, data);
        } catch (error) {
          console.error('Error in export:', error);
          alert('Export failed: ' + error.message);
        }
      });
    });

    // Summary CSV export
    const exportSummaryBtn = document.getElementById('exportSummary');
    if (exportSummaryBtn) {
      exportSummaryBtn.onclick = () => this.exportSummaryCSV();
    }

    // HSD export
    const exportHSDBtn = document.getElementById('exportHSD');
    if (exportHSDBtn) {
      exportHSDBtn.onclick = () => this.exportHSD();
    }
  }

  async exportDataAsZip(key, dataMap) {
    const parts = [];
    
    for (const [k, v] of dataMap) {
      const fileName = this.generateFileName(key, k);
      console.log('Adding file to export:', fileName, 'size:', v?.length || (v instanceof Blob ? v.size : 'unknown'));
      
      let fileContent;
      
      // Convert sequence data to FASTA format for export
      if ((key === 'fasta' || key === 'filtered') && typeof v === 'string' && !v.startsWith('>')) {
        const [sample, region, dir] = k.split('|');
        const fastaContent = `>${sample}_${region}-${dir}\n${v}\n`;
        fileContent = new Blob([fastaContent]);
      } else {
        fileContent = v instanceof Blob ? v : new Blob([v]);
      }
      
      parts.push(new File([fileContent], fileName));
    }
    
    console.log('Exporting', parts.length, 'files as ZIP');
    await this.downloadZip(parts, `${key}.zip`);
  }

  generateFileName(exportType, key) {
    switch (exportType) {
      case 'plots':
        return key; // Already has extension
      case 'damage_analysis':
        return key.replaceAll("|", "_") + "_damage_results.json";
      case 'consensus':
        return key.replaceAll("|", "_") + "_consensus.fasta";
      case 'aligned':
        return key.replaceAll("|", "_") + "_aligned.txt";
      default:
        return key.replaceAll("|", "_") + ".fasta";
    }
  }

  exportSummaryCSV() {
    const statsBody = document.querySelector('#statsTable tbody');
    if (!statsBody || !statsBody.children.length) {
      alert("No stats yet.");
      return;
    }

    let csv = "sample,region,dir,length,Ns,Q30,rate5,rate3\n";
    
    for (const tr of statsBody.children) {
      const cells = [...tr.children].map(td => td.textContent);
      csv += cells.join(",") + "\n";
    }
    
    const file = new File([csv], "summary.csv");
    this.downloadBlob(file);
  }

  exportHSD() {
    const txt = "# HSD preview built from consensus names only\n" +
      Array.from(store.consensus.keys())
        .map(k => k.replaceAll("|", "_"))
        .join("\n");
    
    this.downloadBlob(new File([txt], "haplogroup_analysis.hsd"));
  }

  downloadBlob(file) {
    try {
      console.log('Downloading file:', file.name, 'size:', file.size);
      
      const url = URL.createObjectURL(file);
      const a = document.createElement('a');
      a.href = url;
      a.download = file.name;
      
      // Add to DOM temporarily for some browsers
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      
      setTimeout(() => {
        URL.revokeObjectURL(url);
        console.log('Download cleanup completed for:', file.name);
      }, 1000);
    } catch (error) {
      console.error('Error downloading file:', error);
      alert('Download failed: ' + error.message);
    }
  }

  async downloadZip(files, name) {
    // For many files, use individual downloads to avoid memory issues
    if (files.length > 120) {
      for (const file of files) {
        this.downloadBlob(file);
      }
      return;
    }
    
    try {
      const zipBlob = await this.createZip(files);
      this.downloadBlob(new File([zipBlob], name));
    } catch (error) {
      console.error('Error creating ZIP:', error);
      // Fallback to individual downloads
      for (const file of files) {
        this.downloadBlob(file);
      }
    }
  }

  async createZip(files) {
    // Simple ZIP creation (STORE method only)
    const encoder = new TextEncoder();
    const chunks = [];
    const centralDir = [];
    let offset = 0;

    for (const file of files) {
      const content = await file.arrayBuffer();
      const fileName = encoder.encode(file.name);
      
      // Local file header
      const localHeader = new ArrayBuffer(30 + fileName.length);
      const localView = new DataView(localHeader);
      
      localView.setUint32(0, 0x04034b50, true); // Local file header signature
      localView.setUint16(4, 10, true); // Version needed to extract
      localView.setUint16(6, 0, true); // General purpose bit flag
      localView.setUint16(8, 0, true); // Compression method (stored)
      localView.setUint16(10, 0, true); // Last mod file time
      localView.setUint16(12, 0, true); // Last mod file date
      localView.setUint32(14, this.crc32(content), true); // CRC-32
      localView.setUint32(18, content.byteLength, true); // Compressed size
      localView.setUint32(22, content.byteLength, true); // Uncompressed size
      localView.setUint16(26, fileName.length, true); // File name length
      localView.setUint16(28, 0, true); // Extra field length
      
      const fileHeader = new Uint8Array(localHeader);
      fileHeader.set(fileName, 30);
      
      chunks.push(fileHeader);
      chunks.push(new Uint8Array(content));
      
      // Central directory entry
      centralDir.push({
        fileName,
        crc: this.crc32(content),
        size: content.byteLength,
        offset
      });
      
      offset += fileHeader.length + content.byteLength;
    }

    // Central directory
    const centralDirStart = offset;
    for (const entry of centralDir) {
      const centralHeader = new ArrayBuffer(46 + entry.fileName.length);
      const centralView = new DataView(centralHeader);
      
      centralView.setUint32(0, 0x02014b50, true); // Central directory signature
      centralView.setUint16(4, 20, true); // Version made by
      centralView.setUint16(6, 10, true); // Version needed to extract
      centralView.setUint16(8, 0, true); // General purpose bit flag
      centralView.setUint16(10, 0, true); // Compression method
      centralView.setUint16(12, 0, true); // Last mod file time
      centralView.setUint16(14, 0, true); // Last mod file date
      centralView.setUint32(16, entry.crc, true); // CRC-32
      centralView.setUint32(20, entry.size, true); // Compressed size
      centralView.setUint32(24, entry.size, true); // Uncompressed size
      centralView.setUint16(28, entry.fileName.length, true); // File name length
      centralView.setUint16(30, 0, true); // Extra field length
      centralView.setUint16(32, 0, true); // File comment length
      centralView.setUint16(34, 0, true); // Disk number start
      centralView.setUint16(36, 0, true); // Internal file attributes
      centralView.setUint32(38, 0, true); // External file attributes
      centralView.setUint32(42, entry.offset, true); // Relative offset
      
      const centralHeaderBytes = new Uint8Array(centralHeader);
      centralHeaderBytes.set(entry.fileName, 46);
      chunks.push(centralHeaderBytes);
    }

    const centralDirSize = offset - centralDirStart;

    // End of central directory record
    const endRecord = new ArrayBuffer(22);
    const endView = new DataView(endRecord);
    
    endView.setUint32(0, 0x06054b50, true); // End of central dir signature
    endView.setUint16(4, 0, true); // Number of this disk
    endView.setUint16(6, 0, true); // Disk where central directory starts
    endView.setUint16(8, centralDir.length, true); // Number of central directory records on this disk
    endView.setUint16(10, centralDir.length, true); // Total number of central directory records
    endView.setUint32(12, centralDirSize, true); // Size of central directory
    endView.setUint32(16, centralDirStart, true); // Offset of start of central directory
    endView.setUint16(20, 0, true); // ZIP file comment length
    
    chunks.push(new Uint8Array(endRecord));

    return new Blob(chunks, { type: 'application/zip' });
  }

  crc32(data) {
    // Simple CRC32 implementation
    const table = [];
    for (let i = 0; i < 256; i++) {
      let c = i;
      for (let j = 0; j < 8; j++) {
        c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      }
      table[i] = c;
    }

    let crc = 0 ^ (-1);
    const bytes = new Uint8Array(data);
    
    for (let i = 0; i < bytes.length; i++) {
      crc = (crc >>> 8) ^ table[(crc ^ bytes[i]) & 0xFF];
    }
    
    return (crc ^ (-1)) >>> 0;
  }
}

export default ExportManager;
