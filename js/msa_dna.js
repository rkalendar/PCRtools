// Multiple Sequence Alignment for DNA sequences

class MultipleSequenceAlignment {
  constructor(sequences, matchScore = 1, mismatchPenalty = -1, gapPenalty = -2) {
    this.sequences = sequences;
    this.matchScore = matchScore;
    this.mismatchPenalty = mismatchPenalty;
    this.gapPenalty = gapPenalty;
  }

  // Needleman-Wunsch pairwise alignment
  pairwiseAlign(seq1, seq2) {
    const m = seq1.length;
    const n = seq2.length;
    const dp = Array(m + 1).fill(null).map(() => Array(n + 1).fill(0));
    
    // Initialize first row and column
    for (let i = 0; i <= m; i++) dp[i][0] = i * this.gapPenalty;
    for (let j = 0; j <= n; j++) dp[0][j] = j * this.gapPenalty;
    
    // Fill the matrix
    for (let i = 1; i <= m; i++) {
      for (let j = 1; j <= n; j++) {
        const match = dp[i-1][j-1] + (seq1[i-1] === seq2[j-1] ? this.matchScore : this.mismatchPenalty);
        const del = dp[i-1][j] + this.gapPenalty;
        const ins = dp[i][j-1] + this.gapPenalty;
        dp[i][j] = Math.max(match, del, ins);
      }
    }
    
    // Traceback
    let aligned1 = '', aligned2 = '';
    let i = m, j = n;
    
    while (i > 0 || j > 0) {
      if (i > 0 && j > 0 && dp[i][j] === dp[i-1][j-1] + (seq1[i-1] === seq2[j-1] ? this.matchScore : this.mismatchPenalty)) {
        aligned1 = seq1[i-1] + aligned1;
        aligned2 = seq2[j-1] + aligned2;
        i--;
        j--;
      } else if (i > 0 && dp[i][j] === dp[i-1][j] + this.gapPenalty) {
        aligned1 = seq1[i-1] + aligned1;
        aligned2 = '-' + aligned2;
        i--;
      } else {
        aligned1 = '-' + aligned1;
        aligned2 = seq2[j-1] + aligned2;
        j--;
      }
    }
    
    return { seq1: aligned1, seq2: aligned2, score: dp[m][n] };
  }

  // Calculate similarity score between two sequences
  calculateSimilarity(seq1, seq2) {
    const alignment = this.pairwiseAlign(seq1, seq2);
    return alignment.score;
  }

  // Build distance matrix
  buildDistanceMatrix() {
    const n = this.sequences.length;
    const matrix = Array(n).fill(null).map(() => Array(n).fill(0));
    
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const score = this.calculateSimilarity(this.sequences[i], this.sequences[j]);
        matrix[i][j] = matrix[j][i] = -score; // Convert to distance
      }
    }
    
    return matrix;
  }

  // Progressive alignment using guide tree
  align() {
    if (this.sequences.length === 0) return [];
    if (this.sequences.length === 1) return [this.sequences[0]];
    if (this.sequences.length === 2) {
      const result = this.pairwiseAlign(this.sequences[0], this.sequences[1]);
      return [result.seq1, result.seq2];
    }

    // Build distance matrix
    const distMatrix = this.buildDistanceMatrix();
    
    // Find most similar pair
    let minDist = Infinity, pair = [0, 1];
    for (let i = 0; i < this.sequences.length; i++) {
      for (let j = i + 1; j < this.sequences.length; j++) {
        if (distMatrix[i][j] < minDist) {
          minDist = distMatrix[i][j];
          pair = [i, j];
        }
      }
    }
    
    // Align the most similar pair
    const aligned = this.pairwiseAlign(this.sequences[pair[0]], this.sequences[pair[1]]);
    
    // Create new sequence list with aligned pair as consensus
    const newSequences = [];
    for (let i = 0; i < this.sequences.length; i++) {
      if (i !== pair[0] && i !== pair[1]) {
        newSequences.push(this.sequences[i]);
      }
    }
    
    // Add aligned sequences
    const alignedPair = [aligned.seq1, aligned.seq2];
    
    // If there are more sequences, align them progressively
    if (newSequences.length > 0) {
      const consensus = this.getConsensus([aligned.seq1, aligned.seq2]);
      newSequences.push(consensus);
      
      const msa = new MultipleSequenceAlignment(newSequences, this.matchScore, this.mismatchPenalty, this.gapPenalty);
      const result = msa.align();
      
      // Replace consensus with original aligned pair
      const consensusIdx = result.length - 1;
      result.splice(consensusIdx, 1, ...this.alignToProfile(alignedPair, result[consensusIdx]));
      
      return result;
    }
    
    return alignedPair;
  }

  // Get consensus sequence
  getConsensus(alignedSeqs) {
    if (alignedSeqs.length === 0) return '';
    
    const len = alignedSeqs[0].length;
    let consensus = '';
    
    for (let i = 0; i < len; i++) {
      const counts = { A: 0, T: 0, G: 0, C: 0, '-': 0 };
      
      for (const seq of alignedSeqs) {
        if (seq[i] in counts) counts[seq[i]]++;
      }
      
      delete counts['-'];
      const max = Math.max(...Object.values(counts));
      const base = Object.keys(counts).find(k => counts[k] === max) || '-';
      consensus += base;
    }
    
    return consensus;
  }

  // Align sequences to a profile
  alignToProfile(seqs, profile) {
    return seqs.map(seq => {
      const result = this.pairwiseAlign(seq, profile);
      return result.seq1;
    });
  }

  // Format alignment output
  formatAlignment(aligned) {
    const maxLen = Math.max(...aligned.map(s => s.length));
    let output = '';
    
    aligned.forEach((seq, idx) => {
      output += `Seq${idx + 1}: ${seq}\n`;
    });
    
    // Add conservation line
    output += 'Cons: ';
    for (let i = 0; i < maxLen; i++) {
      const bases = aligned.map(s => s[i]).filter(b => b !== '-');
      const allSame = bases.every(b => b === bases[0]);
      output += allSame && bases.length > 1 ? '*' : ' ';
    }
    
    return output;
  }
}

// Example usage
const sequences = [
  'ATGCTAGCT',
  'ATGCTAGC',
  'ATGCTGCT',
  'ATGCTAGCTA'
];

const msa = new MultipleSequenceAlignment(sequences);
const aligned = msa.align();

console.log('Multiple Sequence Alignment Result:');
console.log(msa.formatAlignment(aligned));
console.log('\nAligned sequences:');
aligned.forEach((seq, i) => console.log(`Sequence ${i + 1}: ${seq}`));