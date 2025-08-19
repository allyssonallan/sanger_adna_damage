/**
 * Main processing pipeline for Sanger aDNA damage analysis
 */

import { store, COUNTERS, Utils } from './core.js';
import { ABIFParser } from './abif-parser.js';
import { SequenceProcessor } from './sequence-processor.js';

export class Pipeline {
  constructor(chromatogramRenderer, uiManager) {
    this.chromatogramRenderer = chromatogramRenderer;
    this.uiManager = uiManager;
  }

  async processFiles(files) {
    console.log('Processing', files.length, 'files');
    
    for (const file of files) {
      await this.processSingleFile(file);
    }
    
    Utils.updateCounts();
    this.uiManager.toggleStatsEmpty();
  }

  async processSingleFile(file) {
    try {
      console.log('Processing file:', file.name);
      
      COUNTERS.totalFiles++;
      const info = Utils.parseFileName(file.name);
      
      console.log('DEBUG: parseFileName result for', file.name, ':', info);
      
      if (!info) {
        console.warn('Invalid filename format:', file.name);
        return;
      }
      
      info.file = file;
      COUNTERS.validFiles++;

      // Read and parse ABIF file
      const arrayBuffer = await file.arrayBuffer();
      const abifData = ABIFParser.readABIF(arrayBuffer);
      
      console.log('ABIF data parsed:', {
        basesLength: abifData.bases?.length || 0,
        basesType: typeof abifData.bases,
        basesContent: abifData.bases?.substring(0, 50) + '...',
        qualsLength: abifData.quals?.length || 0,
        hasTraces: !!abifData.traces
      });
      
      if (!abifData.bases || abifData.bases.length === 0) {
        console.error('No sequence data found in ABIF file:', file.name);
        throw new Error('No sequence data found in ABIF file');
      }
      
      const seqRaw = abifData.bases.toUpperCase().replace(/[^ACGT]/g, "N");
      const seqFiltered = SequenceProcessor.qualityFilter(seqRaw, abifData.quals);
      
      console.log('Processed sequences:', {
        rawLength: seqRaw.length,
        filteredLength: seqFiltered.length,
        rawSample: seqRaw.substring(0, 20),
        filteredSample: seqFiltered.substring(0, 20)
      });
      
      const key = `${info.sample}|${info.region}|${info.dir}`;

      // Store FASTA and filtered sequences
      store.fasta.set(key, seqRaw);  // Store just the sequence, not FASTA format
      store.filtered.set(key, seqFiltered);
      
      console.log('DEBUG: Stored in store.fasta:', key, '- length:', seqRaw.length);
      console.log('DEBUG: store.fasta size now:', store.fasta.size);
      
      // Create chromatogram snapshot
      if (this.chromatogramRenderer) {
        const snapshot = this.chromatogramRenderer.createSnapshot();
        if (snapshot) {
          store.plots.set(`${info.sample}_${info.region}-${info.dir}_quality.png`, snapshot);
        }
      }

      // Update UI
      this.updateTraceSelector(key, info, abifData, seqFiltered);
      this.addStatsRow(info, seqFiltered, abifData.quals);
      
      // Update stage counters
      if (this.uiManager) {
        this.uiManager.updateStageCounters();
      }
      
      COUNTERS.processedFiles++;
      
      console.log('File processed successfully:', file.name, 'Key:', key);
    } catch (error) {
      console.error('Error processing file:', file.name, error);
    }
  }

  updateTraceSelector(key, info, abifData, sequence) {
    const traceSelect = document.getElementById('traceSelect');
    const seqText = document.getElementById('seqText');
    
    // Store trace data globally for later access
    if (!window.chromatogramData) {
      window.chromatogramData = {};
    }
    window.chromatogramData[key] = abifData.traces;
    
    if (traceSelect) {
      traceSelect.add(new Option(`${info.sample} ${info.region}-${info.dir}`, key));
      traceSelect.value = key;
    }
    
    if (this.chromatogramRenderer) {
      console.log('About to draw chromatogram for key:', key);
      console.log('Traces data:', abifData.traces);
      
      // Associate renderer with canvas for later access
      const canvas = document.getElementById('chrom');
      if (canvas) {
        canvas._renderer = this.chromatogramRenderer;
      }
      
      this.chromatogramRenderer.draw(abifData.traces);
    }
    
    if (seqText) {
      seqText.textContent = sequence;
    }
  }

  addStatsRow(info, sequence, qualities) {
    console.log('Adding stats row:', {
      sample: info.sample,
      region: info.region,
      dir: info.dir,
      sequenceLength: sequence?.length || 0,
      qualitiesLength: qualities?.length || 0,
      sequenceSample: sequence?.substring(0, 20)
    });
    
    const Ns = (sequence.match(/N/g) || []).length;
    const q30 = qualities.length ? 
      Math.round(100 * qualities.filter(q => q >= 30).length / qualities.length) : 0;
    
    this.uiManager.addStatsRow({
      sample: info.sample,
      region: info.region,
      dir: info.dir,
      len: sequence.length,
      Ns,
      q30,
      rate5: "–",
      rate3: "–"
    });
  }

  async runStage(stage, file = null, fileInfo = null) {
    console.log('Running pipeline stage:', stage);
    
    try {
      this.uiManager.clearStats();
      
      switch (stage) {
        case 'fasta':
          // Already processed during file upload
          break;
          
        case 'filter':
          // Already done during FASTA processing
          break;
          
        case 'align':
          await this.runAlignmentStage();
          break;
          
        case 'consensus':
          await this.runConsensusStage();
          break;
          
        case 'damage':
          await this.runDamageStage();
          break;
          
        default:
          console.warn('Unknown pipeline stage:', stage);
      }
      
      this.uiManager.toggleStatsEmpty();
      Utils.updateCounts();
      
      // Update stage counters after any pipeline operation
      if (this.uiManager) {
        this.uiManager.updateStageCounters();
      }
      
    } catch (error) {
      console.error('Error in pipeline stage:', stage, error);
      throw error; // Re-throw for button error handling
    }
  }

  async runAlignmentStage() {
    console.log('Running alignment stage');
    
    // Group by sample and region for alignment
    const groups = {};
    
    for (const [key, fasta] of store.filtered) {
      const [sample, region] = key.split('|');
      const groupKey = `${sample}|${region}`;
      
      if (!groups[groupKey]) {
        groups[groupKey] = {};
      }
      
      const sequence = fasta.split('\n')[1] || '';
      groups[groupKey][key] = sequence;
    }

    // Align forward and reverse sequences
    for (const [groupKey, sequences] of Object.entries(groups)) {
      const keys = Object.keys(sequences);
      
      if (keys.length >= 2) {
        const forward = sequences[keys.find(k => k.endsWith('|F'))];
        const reverse = sequences[keys.find(k => k.endsWith('|R'))];
        
        if (forward && reverse) {
          const reverseComp = SequenceProcessor.reverseComplement(reverse);
          const aligned = SequenceProcessor.alignOverlap(forward, reverseComp);
          
          store.aligned.set(groupKey, `>${groupKey.replace('|', '_')}_aligned\n${aligned}\n`);
        }
      }
    }
  }

  async runConsensusStage() {
    console.log('Running consensus stage');
    
    // Build consensus from aligned sequences
    for (const [key, alignedFasta] of store.aligned) {
      const consensus = alignedFasta.split('\n')[1] || '';
      store.consensus.set(key, `>${key.replace('|', '_')}_consensus\n${consensus}\n`);
      
      // Add consensus row to stats
      const [sample, region] = key.split('|');
      this.uiManager.addStatsRow({
        sample,
        region,
        dir: 'C',
        len: consensus.length,
        Ns: (consensus.match(/N/g) || []).length,
        q30: '–',
        rate5: '–',
        rate3: '–'
      });
    }
    
    this.buildFinalSequences();
  }

  async runDamageStage() {
    console.log('Running damage analysis stage');
    
    // Analyze damage patterns
    for (const [regionKey, consensusText] of store.consensus) {
      const [sample, region] = regionKey.split('|');
      const consensus = consensusText.split('\n')[1] || '';
      const damage = SequenceProcessor.analyzeDamage(consensus);
      
      store.damage_analysis.set(regionKey, JSON.stringify(damage, null, 2));
      
      // Update consensus row with damage rates
      this.uiManager.updateDamageRates(sample, region, damage);
    }
  }

  buildFinalSequences() {
    const bySample = {};
    
    for (const [key, fasta] of store.consensus) {
      const [sample, region] = key.split("|");
      bySample[sample] = bySample[sample] || {};
      bySample[sample][region] = fasta.split("\n")[1];
    }
    
    for (const [sample, regions] of Object.entries(bySample)) {
      const merged = [regions.HVS1, regions.HVS2, regions.HVS3]
        .filter(Boolean)
        .join("NNNN");
      
      if (merged) {
        store.final.set(sample, `>${sample}_merged\n${merged}\n`);
      }
    }
  }
}

export default Pipeline;
