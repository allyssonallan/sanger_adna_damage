/**
 * ABIF file format parser for Sanger sequencing data
 */

export class ABIFParser {
  static bufStr(buf, off, len) {
    return new TextDecoder("ascii").decode(buf.slice(off, off + len));
  }

  static getUint32(v, i) {
    return (new DataView(v.buffer || v, v.byteOffset || 0)).getUint32(i, false);
  }

  static readABIF(arrayBuffer) {
    const buf = new Uint8Array(arrayBuffer);
    
    if (ABIFParser.bufStr(buf, 0, 4) !== "ABIF") {
      throw new Error("Not a valid ABIF file");
    }

    const dirOff = ABIFParser.getUint32(buf, 26);
    const numEntries = ABIFParser.getUint32(buf, 18);

    function readDir(off) {
      return {
        tag: ABIFParser.bufStr(buf, off, 4),
        id: ABIFParser.getUint32(buf, off + 4),
        type: ABIFParser.getUint32(buf, off + 8),
        size: ABIFParser.getUint32(buf, off + 12),
        numItems: ABIFParser.getUint32(buf, off + 16),
        dataSize: ABIFParser.getUint32(buf, off + 20),
        dataOff: ABIFParser.getUint32(buf, off + 24)
      };
    }

    const dirs = [];
    for (let i = 0; i < numEntries; i++) {
      dirs.push(readDir(dirOff + i * 28));
    }

    function readTag(tag, id) {
      const d = dirs.find(x => x.tag === tag && x.id === id);
      if (!d) return null;
      
      if (d.dataSize <= 4) {
        // Data stored in offset field
        const dv = new DataView(arrayBuffer);
        return d.type === 2 ? String.fromCharCode(dv.getUint8(dirOff + dirs.indexOf(d) * 28 + 24)) :
               d.type === 4 ? dv.getUint16(dirOff + dirs.indexOf(d) * 28 + 24, false) :
               d.type === 5 ? dv.getUint32(dirOff + dirs.indexOf(d) * 28 + 24, false) : null;
      }
      
      const start = d.dataOff;
      if (d.type === 18) { // pString
        return ABIFParser.bufStr(buf, start + 1, buf[start]);
      } else if (d.type === 19) { // cString
        return ABIFParser.bufStr(buf, start, d.dataSize - 1);
      } else if (d.type === 4) { // short array
        const arr = new Uint16Array(arrayBuffer, start, d.dataSize / 2);
        return Array.from(arr);
      }
      return null;
    }

    function readTrace(tagId) {
      const d = dirs.find(x => x.tag === "DATA" && x.id === tagId);
      if (!d) return null;
      
      const start = d.dataOff;
      const len = d.dataSize / 2;
      const arr = new Uint16Array(arrayBuffer, start, len);
      return Array.from(arr);
    }

    // Extract data
    const bases = readTag("PBAS", 1) || "";
    const quals = readTag("PCON", 2) || [];
    
    const traces = {
      A: readTrace(9),
      C: readTrace(10), 
      G: readTrace(11),
      T: readTrace(12)
    };

    return { bases, quals, traces };
  }
}

export default ABIFParser;
