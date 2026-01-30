/**
 * MUSCLE-JS: JavaScript Implementation of Sequence Alignment Algorithms
 * Based on MUSCLE (MUltiple Sequence Comparison by Log-Expectation)
 * 
 * This is a simplified implementation focusing on core algorithms:
 * - Needleman-Wunsch (global alignment)
 * - Smith-Waterman (local alignment)
 * - Progressive multiple sequence alignment
 * 
 * Original MUSCLE by Robert C. Edgar
 * JS implementation for educational and lightweight use cases
 */

// ============================================================================
// SUBSTITUTION MATRICES
// ============================================================================

/**
 * BLOSUM62 substitution matrix for protein sequences
 * Standard matrix for detecting distant evolutionary relationships
 */
const BLOSUM62 = {
  'A': { 'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4 },
  'R': { 'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4 },
  'N': { 'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4 },
  'D': { 'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4 },
  'C': { 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2, '*': -4 },
  'Q': { 'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4 },
  'E': { 'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4 },
  'G': { 'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4 },
  'H': { 'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4 },
  'I': { 'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4 },
  'L': { 'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4 },
  'K': { 'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4 },
  'M': { 'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4 },
  'F': { 'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4 },
  'P': { 'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -2, '*': -4 },
  'S': { 'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0, '*': -4 },
  'T': { 'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0, '*': -4 },
  'W': { 'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2, '*': -4 },
  'Y': { 'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4 },
  'V': { 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4 },
  'B': { 'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4 },
  'Z': { 'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4 },
  'X': { 'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4 },
  '*': { 'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4, '*': 1 }
};

/**
 * Simple nucleotide scoring matrix
 * Match: +2, Mismatch: -1
 */
const DNA_MATRIX = (() => {
  const bases = ['A', 'T', 'G', 'C', 'N'];
  const matrix = {};
  for (const b1 of bases) {
    matrix[b1] = {};
    for (const b2 of bases) {
      if (b1 === 'N' || b2 === 'N') {
        matrix[b1][b2] = 0;
      } else {
        matrix[b1][b2] = (b1 === b2) ? 2 : -1;
      }
    }
  }
  return matrix;
})();

// ============================================================================
// FASTA PARSER
// ============================================================================

/**
 * Parse FASTA format string into array of sequences
 * @param {string} fastaText - FASTA formatted text
 * @returns {Array<{id: string, description: string, sequence: string}>}
 */
function parseFasta(fastaText) {
  const sequences = [];
  const lines = fastaText.trim().split('\n');
  let currentSeq = null;

  for (const line of lines) {
    const trimmedLine = line.trim();
    if (trimmedLine.startsWith('>')) {
      if (currentSeq) {
        sequences.push(currentSeq);
      }
      const header = trimmedLine.substring(1);
      const spaceIndex = header.indexOf(' ');
      currentSeq = {
        id: spaceIndex > 0 ? header.substring(0, spaceIndex) : header,
        description: spaceIndex > 0 ? header.substring(spaceIndex + 1) : '',
        sequence: ''
      };
    } else if (currentSeq && trimmedLine) {
      currentSeq.sequence += trimmedLine.toUpperCase().replace(/\s/g, '');
    }
  }

  if (currentSeq) {
    sequences.push(currentSeq);
  }

  return sequences;
}

/**
 * Convert aligned sequences to FASTA format
 * @param {Array<{id: string, sequence: string}>} sequences
 * @param {number} lineWidth - Characters per line (default 60)
 * @returns {string}
 */
function toFasta(sequences, lineWidth = 60) {
  let output = '';
  for (const seq of sequences) {
    output += `>${seq.id}`;
    if (seq.description) {
      output += ` ${seq.description}`;
    }
    output += '\n';
    
    for (let i = 0; i < seq.sequence.length; i += lineWidth) {
      output += seq.sequence.substring(i, i + lineWidth) + '\n';
    }
  }
  return output;
}

// ============================================================================
// SCORING FUNCTIONS
// ============================================================================

/**
 * Get substitution score for two residues
 * @param {string} a - First residue
 * @param {string} b - Second residue
 * @param {Object} matrix - Substitution matrix
 * @returns {number}
 */
function getScore(a, b, matrix) {
  if (matrix[a] && matrix[a][b] !== undefined) {
    return matrix[a][b];
  }
  // Handle gaps and unknown characters
  if (a === '-' || b === '-') return 0;
  return a === b ? 1 : -1;
}

/**
 * Detect sequence type (protein or nucleotide)
 * @param {string} sequence
 * @returns {'protein' | 'dna'}
 */
function detectSequenceType(sequence) {
  const dnaChars = new Set(['A', 'T', 'G', 'C', 'N', 'U']);
  const seq = sequence.toUpperCase().replace(/-/g, '');
  let dnaCount = 0;
  
  for (const char of seq) {
    if (dnaChars.has(char)) dnaCount++;
  }
  
  return (dnaCount / seq.length) > 0.9 ? 'dna' : 'protein';
}

// ============================================================================
// NEEDLEMAN-WUNSCH ALGORITHM (Global Alignment)
// ============================================================================

/**
 * Needleman-Wunsch global alignment algorithm
 * @param {string} seq1 - First sequence
 * @param {string} seq2 - Second sequence
 * @param {Object} options - Alignment options
 * @returns {{alignedSeq1: string, alignedSeq2: string, score: number, identity: number}}
 */
function needlemanWunsch(seq1, seq2, options = {}) {
  const {
    gapOpen = -10,
    gapExtend = -1,
    matrix = null
  } = options;

  // Auto-detect matrix if not provided
  const scoringMatrix = matrix || 
    (detectSequenceType(seq1) === 'dna' ? DNA_MATRIX : BLOSUM62);

  const m = seq1.length;
  const n = seq2.length;

  // Initialize scoring matrices for affine gap penalties
  // M[i][j] = best score ending with match/mismatch
  // X[i][j] = best score ending with gap in seq1
  // Y[i][j] = best score ending with gap in seq2
  const M = Array(m + 1).fill(null).map(() => Array(n + 1).fill(-Infinity));
  const X = Array(m + 1).fill(null).map(() => Array(n + 1).fill(-Infinity));
  const Y = Array(m + 1).fill(null).map(() => Array(n + 1).fill(-Infinity));
  
  // Traceback matrices
  const traceM = Array(m + 1).fill(null).map(() => Array(n + 1).fill(null));
  const traceX = Array(m + 1).fill(null).map(() => Array(n + 1).fill(null));
  const traceY = Array(m + 1).fill(null).map(() => Array(n + 1).fill(null));

  // Initialize
  M[0][0] = 0;
  
  for (let i = 1; i <= m; i++) {
    X[i][0] = gapOpen + (i - 1) * gapExtend;
    traceX[i][0] = 'X';
  }
  
  for (let j = 1; j <= n; j++) {
    Y[0][j] = gapOpen + (j - 1) * gapExtend;
    traceY[0][j] = 'Y';
  }

  // Fill matrices
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const matchScore = getScore(seq1[i - 1], seq2[j - 1], scoringMatrix);

      // Update M (match/mismatch)
      const mFromM = M[i - 1][j - 1] + matchScore;
      const mFromX = X[i - 1][j - 1] + matchScore;
      const mFromY = Y[i - 1][j - 1] + matchScore;
      
      if (mFromM >= mFromX && mFromM >= mFromY) {
        M[i][j] = mFromM;
        traceM[i][j] = 'M';
      } else if (mFromX >= mFromY) {
        M[i][j] = mFromX;
        traceM[i][j] = 'X';
      } else {
        M[i][j] = mFromY;
        traceM[i][j] = 'Y';
      }

      // Update X (gap in seq2)
      const xFromM = M[i - 1][j] + gapOpen;
      const xFromX = X[i - 1][j] + gapExtend;
      
      if (xFromM >= xFromX) {
        X[i][j] = xFromM;
        traceX[i][j] = 'M';
      } else {
        X[i][j] = xFromX;
        traceX[i][j] = 'X';
      }

      // Update Y (gap in seq1)
      const yFromM = M[i][j - 1] + gapOpen;
      const yFromY = Y[i][j - 1] + gapExtend;
      
      if (yFromM >= yFromY) {
        Y[i][j] = yFromM;
        traceY[i][j] = 'M';
      } else {
        Y[i][j] = yFromY;
        traceY[i][j] = 'Y';
      }
    }
  }

  // Find best ending score
  const finalScores = [M[m][n], X[m][n], Y[m][n]];
  const maxScore = Math.max(...finalScores);
  let currentMatrix = finalScores.indexOf(maxScore) === 0 ? 'M' : 
                      finalScores.indexOf(maxScore) === 1 ? 'X' : 'Y';

  // Traceback
  let alignedSeq1 = '';
  let alignedSeq2 = '';
  let i = m, j = n;

  while (i > 0 || j > 0) {
    if (currentMatrix === 'M' && i > 0 && j > 0) {
      alignedSeq1 = seq1[i - 1] + alignedSeq1;
      alignedSeq2 = seq2[j - 1] + alignedSeq2;
      const trace = traceM[i][j];
      i--; j--;
      currentMatrix = trace;
    } else if (currentMatrix === 'X' && i > 0) {
      alignedSeq1 = seq1[i - 1] + alignedSeq1;
      alignedSeq2 = '-' + alignedSeq2;
      const trace = traceX[i][j];
      i--;
      currentMatrix = trace;
    } else if (currentMatrix === 'Y' && j > 0) {
      alignedSeq1 = '-' + alignedSeq1;
      alignedSeq2 = seq2[j - 1] + alignedSeq2;
      const trace = traceY[i][j];
      j--;
      currentMatrix = trace;
    } else if (i > 0) {
      alignedSeq1 = seq1[i - 1] + alignedSeq1;
      alignedSeq2 = '-' + alignedSeq2;
      i--;
    } else {
      alignedSeq1 = '-' + alignedSeq1;
      alignedSeq2 = seq2[j - 1] + alignedSeq2;
      j--;
    }
  }

  // Calculate identity
  let matches = 0;
  let alignLength = alignedSeq1.length;
  for (let k = 0; k < alignLength; k++) {
    if (alignedSeq1[k] === alignedSeq2[k] && alignedSeq1[k] !== '-') {
      matches++;
    }
  }
  const identity = matches / alignLength;

  return {
    alignedSeq1,
    alignedSeq2,
    score: maxScore,
    identity,
    matches,
    alignLength
  };
}

// ============================================================================
// SMITH-WATERMAN ALGORITHM (Local Alignment)
// ============================================================================

/**
 * Smith-Waterman local alignment algorithm
 * @param {string} seq1 - First sequence
 * @param {string} seq2 - Second sequence
 * @param {Object} options - Alignment options
 * @returns {{alignedSeq1: string, alignedSeq2: string, score: number, start1: number, start2: number, end1: number, end2: number}}
 */
function smithWaterman(seq1, seq2, options = {}) {
  const {
    gapOpen = -10,
    gapExtend = -1,
    matrix = null
  } = options;

  const scoringMatrix = matrix || 
    (detectSequenceType(seq1) === 'dna' ? DNA_MATRIX : BLOSUM62);

  const m = seq1.length;
  const n = seq2.length;

  // Initialize scoring matrix
  const H = Array(m + 1).fill(null).map(() => Array(n + 1).fill(0));
  const trace = Array(m + 1).fill(null).map(() => Array(n + 1).fill(null));

  let maxScore = 0;
  let maxI = 0, maxJ = 0;

  // Fill matrix
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const matchScore = getScore(seq1[i - 1], seq2[j - 1], scoringMatrix);
      
      const scores = [
        0,  // Start new alignment
        H[i - 1][j - 1] + matchScore,  // Match/mismatch
        H[i - 1][j] + gapOpen,  // Gap in seq2
        H[i][j - 1] + gapOpen   // Gap in seq1
      ];

      const maxIdx = scores.indexOf(Math.max(...scores));
      H[i][j] = scores[maxIdx];
      
      if (maxIdx === 1) trace[i][j] = 'D';  // Diagonal
      else if (maxIdx === 2) trace[i][j] = 'U';  // Up
      else if (maxIdx === 3) trace[i][j] = 'L';  // Left
      else trace[i][j] = null;  // Stop

      if (H[i][j] > maxScore) {
        maxScore = H[i][j];
        maxI = i;
        maxJ = j;
      }
    }
  }

  // Traceback from max score position
  let alignedSeq1 = '';
  let alignedSeq2 = '';
  let i = maxI, j = maxJ;
  const end1 = maxI - 1;
  const end2 = maxJ - 1;

  while (i > 0 && j > 0 && H[i][j] > 0) {
    if (trace[i][j] === 'D') {
      alignedSeq1 = seq1[i - 1] + alignedSeq1;
      alignedSeq2 = seq2[j - 1] + alignedSeq2;
      i--; j--;
    } else if (trace[i][j] === 'U') {
      alignedSeq1 = seq1[i - 1] + alignedSeq1;
      alignedSeq2 = '-' + alignedSeq2;
      i--;
    } else if (trace[i][j] === 'L') {
      alignedSeq1 = '-' + alignedSeq1;
      alignedSeq2 = seq2[j - 1] + alignedSeq2;
      j--;
    } else {
      break;
    }
  }

  const start1 = i;
  const start2 = j;

  // Calculate identity
  let matches = 0;
  for (let k = 0; k < alignedSeq1.length; k++) {
    if (alignedSeq1[k] === alignedSeq2[k] && alignedSeq1[k] !== '-') {
      matches++;
    }
  }

  return {
    alignedSeq1,
    alignedSeq2,
    score: maxScore,
    identity: matches / alignedSeq1.length,
    start1,
    start2,
    end1,
    end2,
    matches
  };
}

// ============================================================================
// DISTANCE CALCULATION
// ============================================================================

/**
 * Calculate pairwise distance using k-mer method (fast approximation)
 * @param {string} seq1 
 * @param {string} seq2 
 * @param {number} k - k-mer size
 * @returns {number} - Distance (0-1)
 */
function kmerDistance(seq1, seq2, k = 3) {
  const getKmers = (seq) => {
    const kmers = new Set();
    for (let i = 0; i <= seq.length - k; i++) {
      kmers.add(seq.substring(i, i + k));
    }
    return kmers;
  };

  const kmers1 = getKmers(seq1);
  const kmers2 = getKmers(seq2);

  let intersection = 0;
  for (const kmer of kmers1) {
    if (kmers2.has(kmer)) intersection++;
  }

  const union = kmers1.size + kmers2.size - intersection;
  const similarity = union > 0 ? intersection / union : 0;
  
  return 1 - similarity;
}

/**
 * Calculate distance matrix for sequences
 * @param {Array<{id: string, sequence: string}>} sequences 
 * @returns {Array<Array<number>>}
 */
function calculateDistanceMatrix(sequences) {
  const n = sequences.length;
  const matrix = Array(n).fill(null).map(() => Array(n).fill(0));

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dist = kmerDistance(sequences[i].sequence, sequences[j].sequence);
      matrix[i][j] = dist;
      matrix[j][i] = dist;
    }
  }

  return matrix;
}

// ============================================================================
// UPGMA TREE BUILDING
// ============================================================================

/**
 * Build guide tree using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
 * @param {Array<Array<number>>} distMatrix 
 * @param {Array<string>} labels 
 * @returns {Object} - Tree structure
 */
function buildUPGMATree(distMatrix, labels) {
  const n = labels.length;
  
  // Clone distance matrix
  const D = distMatrix.map(row => [...row]);
  
  // Initialize clusters
  const clusters = labels.map((label, i) => ({
    id: i,
    label,
    size: 1,
    height: 0,
    left: null,
    right: null
  }));
  
  const active = new Set(clusters.map((_, i) => i));
  let nextId = n;

  while (active.size > 1) {
    // Find minimum distance
    let minDist = Infinity;
    let minI = -1, minJ = -1;
    
    const activeArray = [...active];
    for (let ii = 0; ii < activeArray.length; ii++) {
      for (let jj = ii + 1; jj < activeArray.length; jj++) {
        const i = activeArray[ii];
        const j = activeArray[jj];
        if (D[i][j] < minDist) {
          minDist = D[i][j];
          minI = i;
          minJ = j;
        }
      }
    }

    // Create new cluster
    const newCluster = {
      id: nextId,
      label: `Node${nextId}`,
      size: clusters[minI].size + clusters[minJ].size,
      height: minDist / 2,
      left: clusters[minI],
      right: clusters[minJ]
    };
    
    clusters.push(newCluster);

    // Update distance matrix
    for (const k of active) {
      if (k !== minI && k !== minJ) {
        const newDist = (D[minI][k] * clusters[minI].size + 
                        D[minJ][k] * clusters[minJ].size) / 
                       (clusters[minI].size + clusters[minJ].size);
        D[nextId] = D[nextId] || [];
        D[nextId][k] = newDist;
        D[k] = D[k] || [];
        D[k][nextId] = newDist;
      }
    }
    D[nextId] = D[nextId] || [];
    D[nextId][nextId] = 0;

    active.delete(minI);
    active.delete(minJ);
    active.add(nextId);
    nextId++;
  }

  return clusters[clusters.length - 1];
}

/**
 * Get leaf order from tree (for progressive alignment)
 * @param {Object} tree 
 * @returns {Array<number>}
 */
function getTreeOrder(tree) {
  const order = [];
  
  function traverse(node) {
    if (node.left === null && node.right === null) {
      order.push(node.id);
    } else {
      if (node.left) traverse(node.left);
      if (node.right) traverse(node.right);
    }
  }
  
  traverse(tree);
  return order;
}

// ============================================================================
// PROFILE ALIGNMENT
// ============================================================================

/**
 * Create sequence profile from aligned sequences
 * @param {Array<string>} alignedSeqs 
 * @returns {Array<Object>} - Profile (frequency at each position)
 */
function createProfile(alignedSeqs) {
  if (alignedSeqs.length === 0) return [];
  
  const length = alignedSeqs[0].length;
  const profile = [];

  for (let i = 0; i < length; i++) {
    const counts = { '-': 0 };
    for (const seq of alignedSeqs) {
      const char = seq[i] || '-';
      counts[char] = (counts[char] || 0) + 1;
    }
    
    // Convert to frequencies
    const freq = {};
    for (const [char, count] of Object.entries(counts)) {
      freq[char] = count / alignedSeqs.length;
    }
    profile.push(freq);
  }

  return profile;
}

/**
 * Score two profile positions
 * @param {Object} prof1 - Profile position 1
 * @param {Object} prof2 - Profile position 2
 * @param {Object} matrix - Substitution matrix
 * @returns {number}
 */
function profileScore(prof1, prof2, matrix) {
  let score = 0;
  
  for (const [char1, freq1] of Object.entries(prof1)) {
    for (const [char2, freq2] of Object.entries(prof2)) {
      if (char1 !== '-' && char2 !== '-') {
        score += freq1 * freq2 * getScore(char1, char2, matrix);
      }
    }
  }
  
  return score;
}

/**
 * Align two profiles (groups of aligned sequences)
 * @param {Array<string>} group1 
 * @param {Array<string>} group2 
 * @param {Object} options 
 * @returns {{aligned1: Array<string>, aligned2: Array<string>}}
 */
function alignProfiles(group1, group2, options = {}) {
  const {
    gapOpen = -10,
    gapExtend = -1,
    matrix = BLOSUM62
  } = options;

  const profile1 = createProfile(group1);
  const profile2 = createProfile(group2);
  
  const m = profile1.length;
  const n = profile2.length;

  // Simple gap scoring for profile alignment
  const gapScore = gapOpen;

  // Initialize DP matrix
  const H = Array(m + 1).fill(null).map(() => Array(n + 1).fill(0));
  const trace = Array(m + 1).fill(null).map(() => Array(n + 1).fill(null));

  for (let i = 1; i <= m; i++) {
    H[i][0] = i * gapScore;
    trace[i][0] = 'U';
  }
  for (let j = 1; j <= n; j++) {
    H[0][j] = j * gapScore;
    trace[0][j] = 'L';
  }

  // Fill matrix
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const matchScore = profileScore(profile1[i - 1], profile2[j - 1], matrix);
      
      const scores = [
        H[i - 1][j - 1] + matchScore,
        H[i - 1][j] + gapScore,
        H[i][j - 1] + gapScore
      ];

      const maxIdx = scores.indexOf(Math.max(...scores));
      H[i][j] = scores[maxIdx];
      
      if (maxIdx === 0) trace[i][j] = 'D';
      else if (maxIdx === 1) trace[i][j] = 'U';
      else trace[i][j] = 'L';
    }
  }

  // Traceback
  const alignment1 = [];
  const alignment2 = [];
  let i = m, j = n;

  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && trace[i][j] === 'D') {
      alignment1.unshift(i - 1);
      alignment2.unshift(j - 1);
      i--; j--;
    } else if (i > 0 && (j === 0 || trace[i][j] === 'U')) {
      alignment1.unshift(i - 1);
      alignment2.unshift(-1);  // Gap
      i--;
    } else {
      alignment1.unshift(-1);  // Gap
      alignment2.unshift(j - 1);
      j--;
    }
  }

  // Apply alignment to sequences
  const aligned1 = group1.map(seq => {
    let result = '';
    for (const idx of alignment1) {
      result += idx >= 0 ? seq[idx] : '-';
    }
    return result;
  });

  const aligned2 = group2.map(seq => {
    let result = '';
    for (const idx of alignment2) {
      result += idx >= 0 ? seq[idx] : '-';
    }
    return result;
  });

  return { aligned1, aligned2 };
}

// ============================================================================
// PROGRESSIVE MULTIPLE SEQUENCE ALIGNMENT
// ============================================================================

/**
 * Perform progressive multiple sequence alignment
 * @param {Array<{id: string, sequence: string}>} sequences 
 * @param {Object} options 
 * @returns {Array<{id: string, sequence: string}>}
 */
function progressiveAlignment(sequences, options = {}) {
  if (sequences.length === 0) return [];
  if (sequences.length === 1) return sequences;

  const {
    gapOpen = -10,
    gapExtend = -1,
    matrix = null
  } = options;

  // Auto-detect matrix
  const scoringMatrix = matrix || 
    (detectSequenceType(sequences[0].sequence) === 'dna' ? DNA_MATRIX : BLOSUM62);

  // Step 1: Calculate distance matrix
  const distMatrix = calculateDistanceMatrix(sequences);

  // Step 2: Build guide tree using UPGMA
  const labels = sequences.map(s => s.id);
  const tree = buildUPGMATree(distMatrix, labels);

  // Step 3: Progressive alignment following tree
  const alignedGroups = {};  // Maps node id to aligned sequences

  function alignNode(node) {
    if (node.left === null && node.right === null) {
      // Leaf node
      alignedGroups[node.id] = [sequences[node.id].sequence];
      return;
    }

    // Recursively align children
    alignNode(node.left);
    alignNode(node.right);

    // Align the two groups
    const group1 = alignedGroups[node.left.id];
    const group2 = alignedGroups[node.right.id];

    const { aligned1, aligned2 } = alignProfiles(group1, group2, {
      gapOpen,
      gapExtend,
      matrix: scoringMatrix
    });

    alignedGroups[node.id] = [...aligned1, ...aligned2];
  }

  alignNode(tree);

  // Get final alignment in original order
  const treeOrder = getTreeOrder(tree);
  const orderMap = {};
  treeOrder.forEach((origIdx, alignIdx) => {
    orderMap[origIdx] = alignIdx;
  });

  const finalAlignment = alignedGroups[tree.id];
  
  // Map back to original sequence order
  const result = sequences.map((seq, i) => ({
    id: seq.id,
    description: seq.description || '',
    sequence: finalAlignment[orderMap[i]]
  }));

  return result;
}

// ============================================================================
// MAIN ALIGNMENT FUNCTION (MUSCLE-like interface)
// ============================================================================

/**
 * Main alignment function with MUSCLE-like interface
 * @param {string|Array} input - FASTA string or array of sequences
 * @param {Object} options - Alignment options
 * @returns {Object} - Alignment result
 */
function align(input, options = {}) {
  // Parse input
  let sequences;
  if (typeof input === 'string') {
    sequences = parseFasta(input);
  } else {
    sequences = input;
  }

  if (sequences.length === 0) {
    throw new Error('No sequences provided');
  }

  const {
    type = 'auto',  // 'global', 'local', or 'auto'
    gapOpen = -10,
    gapExtend = -1,
    matrix = null,
    outputFormat = 'object'  // 'object', 'fasta', or 'clustal'
  } = options;

  let result;

  if (sequences.length === 2 && type !== 'auto') {
    // Pairwise alignment
    const seq1 = sequences[0].sequence;
    const seq2 = sequences[1].sequence;
    
    if (type === 'local') {
      const alignment = smithWaterman(seq1, seq2, { gapOpen, gapExtend, matrix });
      result = {
        type: 'local',
        sequences: [
          { id: sequences[0].id, sequence: alignment.alignedSeq1 },
          { id: sequences[1].id, sequence: alignment.alignedSeq2 }
        ],
        score: alignment.score,
        identity: alignment.identity,
        start: [alignment.start1, alignment.start2],
        end: [alignment.end1, alignment.end2]
      };
    } else {
      const alignment = needlemanWunsch(seq1, seq2, { gapOpen, gapExtend, matrix });
      result = {
        type: 'global',
        sequences: [
          { id: sequences[0].id, sequence: alignment.alignedSeq1 },
          { id: sequences[1].id, sequence: alignment.alignedSeq2 }
        ],
        score: alignment.score,
        identity: alignment.identity
      };
    }
  } else {
    // Multiple sequence alignment
    const aligned = progressiveAlignment(sequences, { gapOpen, gapExtend, matrix });
    
    // Calculate overall statistics
    let totalIdentity = 0;
    let comparisons = 0;
    const alignLength = aligned[0]?.sequence.length || 0;
    
    for (let i = 0; i < aligned.length; i++) {
      for (let j = i + 1; j < aligned.length; j++) {
        let matches = 0;
        for (let k = 0; k < alignLength; k++) {
          if (aligned[i].sequence[k] === aligned[j].sequence[k] && 
              aligned[i].sequence[k] !== '-') {
            matches++;
          }
        }
        totalIdentity += matches / alignLength;
        comparisons++;
      }
    }

    result = {
      type: 'multiple',
      sequences: aligned,
      alignmentLength: alignLength,
      averageIdentity: comparisons > 0 ? totalIdentity / comparisons : 0,
      numSequences: aligned.length
    };
  }

  // Format output
  if (outputFormat === 'fasta') {
    return toFasta(result.sequences);
  } else if (outputFormat === 'clustal') {
    return toClustal(result.sequences);
  }

  return result;
}

/**
 * Convert to CLUSTAL format
 * @param {Array<{id: string, sequence: string}>} sequences 
 * @returns {string}
 */
function toClustal(sequences) {
  let output = 'CLUSTAL W (1.83) multiple sequence alignment\n\n';
  
  const maxIdLen = Math.max(...sequences.map(s => s.id.length), 10);
  const blockSize = 60;
  const alignLen = sequences[0]?.sequence.length || 0;

  for (let start = 0; start < alignLen; start += blockSize) {
    for (const seq of sequences) {
      const id = seq.id.padEnd(maxIdLen);
      const block = seq.sequence.substring(start, start + blockSize);
      output += `${id} ${block}\n`;
    }
    
    // Conservation line
    output += ' '.repeat(maxIdLen) + ' ';
    for (let i = start; i < Math.min(start + blockSize, alignLen); i++) {
      const chars = sequences.map(s => s.sequence[i]);
      const allSame = chars.every(c => c === chars[0] && c !== '-');
      output += allSame ? '*' : ' ';
    }
    output += '\n\n';
  }

  return output;
}

// ============================================================================
// EXPORTS
// ============================================================================

// For Node.js
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    // Main functions
    align,
    needlemanWunsch,
    smithWaterman,
    progressiveAlignment,
    
    // Parsing
    parseFasta,
    toFasta,
    toClustal,
    
    // Utilities
    detectSequenceType,
    calculateDistanceMatrix,
    buildUPGMATree,
    kmerDistance,
    
    // Matrices
    BLOSUM62,
    DNA_MATRIX
  };
}

// For browser
if (typeof window !== 'undefined') {
  window.MuscleJS = {
    align,
    needlemanWunsch,
    smithWaterman,
    progressiveAlignment,
    parseFasta,
    toFasta,
    toClustal,
    detectSequenceType,
    calculateDistanceMatrix,
    buildUPGMATree,
    kmerDistance,
    BLOSUM62,
    DNA_MATRIX
  };
}
