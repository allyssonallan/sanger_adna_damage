/**
 * Core utilities, theme management, and constants
 */

// Theme management
export const Theme = {
  init() {
    const root = document.documentElement;
    const themeToggle = document.getElementById('themeToggle');
    
    themeToggle.addEventListener('change', () => 
      root.classList.toggle('dark', themeToggle.checked)
    );
    
    if (window.matchMedia('(prefers-color-scheme: dark)').matches) {
      themeToggle.checked = true;
      root.classList.add('dark');
    }
  }
};

// Constants
export const CONSTANTS = {
  QUALITY_THRESHOLD: 30,
  MAX_BOOTSTRAP_SAMPLES: 1000,
  DAMAGE_WINDOW_SIZE: 20,
  ALIGNMENT_BAND: 30,
  RC_MAP: {A:"T", C:"G", G:"C", T:"A", N:"N"},
  TRACE_COLORS: {
    "A": "#2b8a3e",
    "C": "#2563eb", 
    "G": "#f59e0b",
    "T": "#ef4444"
  }
};

// Data store with counters
export const store = {
  fasta: new Map(),
  filtered: new Map(),
  aligned: new Map(),
  consensus: new Map(),
  final: new Map(),
  plots: new Map(),
  damage_analysis: new Map(),
  hsd: "",
  stats: []
};

export const COUNTERS = {
  totalFiles: 0,
  validFiles: 0,
  processedFiles: 0
};

// Utility functions
export const Utils = {
  async fileToBase64(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result);
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  },

  updateCounts() {
    document.getElementById('fileCount').textContent = COUNTERS.totalFiles;
    document.getElementById('validCount').textContent = COUNTERS.validFiles;
    document.getElementById('processedCount').textContent = COUNTERS.processedFiles;
  },

  parseFileName(fileName) {
    // Updated regex to match actual filename format: 10_CON_TA_PB_B2_B1_HVS1-F.ab1
    const regex = /^(.+)_([A-Z0-9]+)-([FR])\.ab1$/i;
    const match = fileName.match(regex);
    console.log('Parsing filename:', fileName, 'Match result:', match);
    if (!match) return null;
    
    // Extract the last underscore-separated part as region (HVS1, HVS2, etc.)
    const fullSample = match[1]; // Everything before the region
    const region = match[2]; // The region part (HVS1, HVS2, etc.)
    const direction = match[3].toUpperCase(); // F or R
    
    return {
      sample: fullSample,
      region: region, 
      dir: direction,
      file: null
    };
  },

  formatPercentage(p, ci) {
    return `${(100*p).toFixed(1)}% [${(100*ci[0]).toFixed(1)}–${(100*ci[1]).toFixed(1)}]`;
  }
};

// Logo injection
export const LogoManager = {
  async inject() {
    console.log('LogoManager.inject() called');
    const logoRowTop = document.getElementById('logoRowTop');
    const logoRowBottom = document.getElementById('logoRowBottom');
    
    console.log('logoRowTop:', logoRowTop);
    console.log('logoRowBottom:', logoRowBottom);
    
    if (!logoRowTop || !logoRowBottom) {
      console.warn('Logo containers not found');
      return;
    }

    console.log('Starting logo loading...');
    const logos = [
      {url:'https://raw.githubusercontent.com/allyssonallan/sanger_adna_damage/main/config/logos/ufc.png',    href:'https://ufc.br', title:'Universidade Federal do Ceará - https://ufc.br', h:28},
      {url:'https://raw.githubusercontent.com/allyssonallan/sanger_adna_damage/main/config/logos/funcap.png', href:'https://funcap.ce.gov.br', title:'FUNCAP - https://www.funcap.ce.gov.br', h:28},
      {url:'https://raw.githubusercontent.com/allyssonallan/sanger_adna_damage/main/config/logos/labbat.png', href:'https://instagram.com/labbat.npdm.ufc', title:'LABBAT - https://instagram.com/labbat.npdm.ufc', h:22},
      {url:'https://raw.githubusercontent.com/allyssonallan/sanger_adna_damage/main/config/logos/npdm.png',   href:'https://npdm.ufc.br', title:'NPDM - https://npdm.ufc.br', h:22},
    ];

    const toDataURL = async (u) => {
      const r = await fetch(u, {mode:'cors'}); 
      const b = await r.blob();
      return new Promise(res => {
        const fr = new FileReader(); 
        fr.onload = () => res(fr.result); 
        fr.readAsDataURL(b);
      });
    };

    try {
      console.log('Fetching logos...');
      const [ufc, funcap, labbat, npdm] = await Promise.all(logos.map(l => toDataURL(l.url)));
      
      console.log('Logos loaded, injecting HTML...');
      logoRowTop.innerHTML =
        `<a href="${logos[0].href}" target="_blank" rel="noopener"><img src="${ufc}" alt="UFC" title="${logos[0].title}" style="height:${logos[0].h}px; background: white; padding: 4px; border-radius: 4px; margin: 2px;"></a>` +
        `<a href="${logos[1].href}" target="_blank" rel="noopener"><img src="${funcap}" alt="FUNCAP" title="${logos[1].title}" style="height:${logos[1].h}px; background: white; padding: 4px; border-radius: 4px; margin: 2px;"></a>`;
        
      logoRowBottom.innerHTML =
        `<a href="${logos[2].href}" target="_blank" rel="noopener"><img src="${labbat}" alt="LABBAT" title="${logos[2].title}" style="height:${logos[2].h}px; background: white; padding: 4px; border-radius: 4px; margin: 2px;"></a>` +
        `<a href="${logos[3].href}" target="_blank" rel="noopener"><img src="${npdm}" alt="NPDM" title="${logos[3].title}" style="height:${logos[3].h}px; background: white; padding: 4px; border-radius: 4px; margin: 2px;"></a>`;
      
      console.log('Logos injected successfully');
    } catch (error) {
      console.warn('Error loading logos:', error);
    }
  }
};
