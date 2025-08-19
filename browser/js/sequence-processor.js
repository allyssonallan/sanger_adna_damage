/**
 * Sequence processing utilities: quality filtering, alignment, damage analysis
 */

import { CONSTANTS } from './core.js';

export class SequenceProcessor {
  static reverseComplement(sequence) {
    return sequence
      .split("")
      .reverse()
      .map(base => CONSTANTS.RC_MAP[base] || "N")
      .join("");
  }

  static qualityFilter(sequence, qualities) {
    if (!qualities || qualities.length === 0) return sequence;
    
    return sequence
      .split("")
      .map((base, i) => {
        const q = qualities[i] || 0;
        return q >= CONSTANTS.QUALITY_THRESHOLD ? base : "N";
      })
      .join("");
  }

  static alignOverlap(seqA, seqB, band = CONSTANTS.ALIGNMENT_BAND) {
    const m = seqA.length;
    const n = seqB.length;
    const maxShift = Math.min(m, n, band);
    
    let bestScore = -Infinity;
    let bestOverlap = "";
    let bestI = 0;

    // Try different overlap positions
    for (let shift = -maxShift; shift <= maxShift; shift++) {
      let score = 0;
      let matches = 0;
      let overlapStart, overlapEnd;

      if (shift >= 0) {
        overlapStart = shift;
        overlapEnd = Math.min(m, n + shift);
      } else {
        overlapStart = 0;
        overlapEnd = Math.min(m + shift, n);
      }

      const overlapLen = overlapEnd - overlapStart;
      if (overlapLen <= 0) continue;

      // Score the overlap
      for (let i = 0; i < overlapLen; i++) {
        const posA = shift >= 0 ? shift + i : i;
        const posB = shift >= 0 ? i : -shift + i;
        
        if (posA < m && posB < n) {
          if (seqA[posA] === seqB[posB] && seqA[posA] !== 'N') {
            score += 2;
            matches++;
          } else if (seqA[posA] === 'N' || seqB[posB] === 'N') {
            score += 0;
          } else {
            score -= 1;
          }
        }
      }

      if (score > bestScore && matches > overlapLen * 0.6) {
        bestScore = score;
        bestI = shift;
        
        // Build consensus
        let consensus = "";
        const totalLen = Math.max(m, n + shift, m - shift);
        
        for (let pos = 0; pos < totalLen; pos++) {
          const aPos = shift >= 0 ? pos - shift : pos;
          const bPos = shift >= 0 ? pos : pos + shift;
          
          const aChar = (aPos >= 0 && aPos < m) ? seqA[aPos] : "";
          const bChar = (bPos >= 0 && bPos < n) ? seqB[bPos] : "";
          
          if (aChar && bChar) {
            consensus += (aChar !== 'N') ? aChar : bChar;
          } else {
            consensus += aChar || bChar;
          }
        }
        
        bestOverlap = consensus;
      }
    }

    return bestOverlap || seqA; // fallback to first sequence
  }

  static analyzeDamage(consensus) {
    const windowSize = CONSTANTS.DAMAGE_WINDOW_SIZE;
    const len = consensus.length;
    
    if (len < windowSize * 2) {
      return { rate5: 0, rate3: 0, ci5: [0, 0], ci3: [0, 0] };
    }

    // 5' end: count C->T transitions
    const start5 = consensus.slice(0, windowSize);
    let ct5 = 0, totalC5 = 0;
    
    for (let i = 0; i < start5.length; i++) {
      if (start5[i] === 'C' || start5[i] === 'T') {
        totalC5++;
        if (start5[i] === 'T') ct5++;
      }
    }

    // 3' end: count G->A transitions  
    const start3 = consensus.slice(-windowSize);
    let ga3 = 0, totalG3 = 0;
    
    for (let i = 0; i < start3.length; i++) {
      if (start3[i] === 'G' || start3[i] === 'A') {
        totalG3++;
        if (start3[i] === 'A') ga3++;
      }
    }

    const rate5 = totalC5 > 0 ? ct5 / totalC5 : 0;
    const rate3 = totalG3 > 0 ? ga3 / totalG3 : 0;

    // Bootstrap confidence intervals
    const bootstraps5 = [];
    const bootstraps3 = [];
    
    for (let b = 0; b < CONSTANTS.MAX_BOOTSTRAP_SAMPLES; b++) {
      let bootCt5 = 0, bootTotalC5 = 0;
      let bootGa3 = 0, bootTotalG3 = 0;
      
      // Bootstrap sample 5' end
      for (let i = 0; i < windowSize; i++) {
        const idx = Math.floor(Math.random() * start5.length);
        if (start5[idx] === 'C' || start5[idx] === 'T') {
          bootTotalC5++;
          if (start5[idx] === 'T') bootCt5++;
        }
      }
      
      // Bootstrap sample 3' end
      for (let i = 0; i < windowSize; i++) {
        const idx = Math.floor(Math.random() * start3.length);
        if (start3[idx] === 'G' || start3[idx] === 'A') {
          bootTotalG3++;
          if (start3[idx] === 'A') bootGa3++;
        }
      }
      
      bootstraps5.push(bootTotalC5 > 0 ? bootCt5 / bootTotalC5 : 0);
      bootstraps3.push(bootTotalG3 > 0 ? bootGa3 / bootTotalG3 : 0);
    }

    // Calculate confidence intervals
    function calculateCI(arr) {
      arr.sort((a, b) => a - b);
      const lo = arr[Math.floor(0.025 * arr.length)] || 0;
      const hi = arr[Math.floor(0.975 * arr.length)] || 0;
      return [lo, hi];
    }

    return {
      rate5,
      rate3,
      ci5: calculateCI(bootstraps5),
      ci3: calculateCI(bootstraps3)
    };
  }
}

export default SequenceProcessor;
