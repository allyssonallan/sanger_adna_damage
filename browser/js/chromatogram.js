/**
 * Chromatogram rendering and visualization
 */

import { CONSTANTS } from './core.js';

export class ChromatogramRenderer {
  constructor(canvasId) {
    this.canvas = document.getElementById(canvasId);
    this.ctx = this.canvas?.getContext('2d');
    
    if (!this.canvas) {
      console.error(`Canvas element with id "${canvasId}" not found!`);
      return;
    }
    
    if (!this.ctx) {
      console.error('Failed to get canvas 2D context');
      return;
    }

    console.log('Chromatogram renderer initialized');
  }

  draw(traces) {
    try {
      console.log('Drawing chromatogram with traces:', traces);
      
      if (!traces || typeof traces !== 'object') {
        console.error('Invalid traces data:', traces);
        return;
      }

      if (!this.canvas || !this.ctx) {
        console.error('Canvas or context not available');
        return;
      }

      const dpr = window.devicePixelRatio || 1;
      this.canvas.width = this.canvas.clientWidth * dpr;
      this.canvas.height = this.canvas.clientHeight * dpr;
      this.ctx.scale(dpr, dpr);
      this.ctx.clearRect(0, 0, this.canvas.clientWidth, this.canvas.clientHeight);

      const W = this.canvas.clientWidth;
      const H = this.canvas.clientHeight;
      const pad = 16;

      console.log('Canvas dimensions:', { W, H, clientWidth: this.canvas.clientWidth, clientHeight: this.canvas.clientHeight });

      const keys = ["A", "C", "G", "T"];
      
      // Find maximum value across all traces
      let max = 1;
      keys.forEach(key => {
        if (traces[key] && traces[key].length > 0) {
          const traceMax = Math.max(...traces[key]);
          max = Math.max(max, traceMax);
          console.log(`Trace ${key}: length=${traces[key].length}, max=${traceMax}`);
        }
      });

      console.log('Max value found:', max);

      // Draw each trace
      keys.forEach(key => {
        const arr = traces[key];
        if (!arr || arr.length === 0) {
          console.log(`No data for trace ${key}`);
          return;
        }

        console.log(`Drawing trace ${key} with ${arr.length} points`);
        
        this.ctx.beginPath();
        this.ctx.strokeStyle = CONSTANTS.TRACE_COLORS[key];
        this.ctx.lineWidth = 1.2;
        
        for (let x = 0; x < arr.length; x++) {
          const y = H - pad - (arr[x] / max) * (H - 2 * pad);
          if (x === 0) {
            this.ctx.moveTo(pad, y);
          } else {
            this.ctx.lineTo(pad + x * (W - 2 * pad) / arr.length, y);
          }
        }
        
        this.ctx.stroke();
      });

      console.log('Chromatogram drawing completed');
    } catch (error) {
      console.error('Error drawing chromatogram:', error);
    }
  }

  createSnapshot() {
    if (!this.canvas) return null;
    
    try {
      // Create a smaller snapshot canvas
      const snapshotCanvas = document.createElement('canvas');
      snapshotCanvas.width = 200;
      snapshotCanvas.height = 100;
      const snapshotCtx = snapshotCanvas.getContext('2d');
      
      // Draw scaled-down version
      snapshotCtx.drawImage(this.canvas, 0, 0, 200, 100);
      
      return snapshotCanvas.toDataURL('image/png');
    } catch (error) {
      console.error('Error creating chromatogram snapshot:', error);
      return null;
    }
  }

  setupPNGDownload(buttonId) {
    const button = document.getElementById(buttonId);
    if (!button) return;

    button.onclick = () => {
      try {
        if (!this.canvas) return;
        
        const url = this.canvas.toDataURL('image/png');
        const a = document.createElement('a');
        a.href = url;
        a.download = 'chromatogram.png';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      } catch (error) {
        console.error('Error downloading PNG:', error);
        alert('Failed to download PNG: ' + error.message);
      }
    };
  }
}

export default ChromatogramRenderer;
