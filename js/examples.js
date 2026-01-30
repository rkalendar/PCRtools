/**
 * MUSCLE-JS Examples and Tests
 */

const {
  align,
  needlemanWunsch,
  smithWaterman,
  progressiveAlignment,
  parseFasta,
  toFasta,
  toClustal,
  detectSequenceType,
  BLOSUM62,
  DNA_MATRIX
} = require('./muscle.js');

// ============================================================================
// EXAMPLE 1: Simple Pairwise Global Alignment (Needleman-Wunsch)
// ============================================================================

console.log('=' .repeat(70));
console.log('EXAMPLE 1: Pairwise Global Alignment (Needleman-Wunsch)');
console.log('=' .repeat(70));

const seq1 = 'HEAGAWGHEE';
const seq2 = 'PAWHEAE';

const globalResult = needlemanWunsch(seq1, seq2);

console.log('\nSequence 1:', seq1);
console.log('Sequence 2:', seq2);
console.log('\nAlignment:');
console.log('  ', globalResult.alignedSeq1);
console.log('  ', globalResult.alignedSeq2);
console.log('\nScore:', globalResult.score);
console.log('Identity:', (globalResult.identity * 100).toFixed(1) + '%');

// ============================================================================
// EXAMPLE 2: Local Alignment (Smith-Waterman)
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 2: Local Alignment (Smith-Waterman)');
console.log('=' .repeat(70));

const longSeq1 = 'XXXXHEAGAWGHEEXXXX';
const longSeq2 = 'YYYYPAWHEAEYYY';

const localResult = smithWaterman(longSeq1, longSeq2);

console.log('\nSequence 1:', longSeq1);
console.log('Sequence 2:', longSeq2);
console.log('\nLocal Alignment:');
console.log('  ', localResult.alignedSeq1);
console.log('  ', localResult.alignedSeq2);
console.log('\nScore:', localResult.score);
console.log('Position in Seq1:', localResult.start1 + 1, '-', localResult.end1 + 1);
console.log('Position in Seq2:', localResult.start2 + 1, '-', localResult.end2 + 1);

// ============================================================================
// EXAMPLE 3: DNA Sequence Alignment
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 3: DNA Sequence Alignment');
console.log('=' .repeat(70));

const dna1 = 'ATCGATCGATCG';
const dna2 = 'ATCGTTCGATC';

const dnaResult = needlemanWunsch(dna1, dna2);

console.log('\nDNA Sequence 1:', dna1);
console.log('DNA Sequence 2:', dna2);
console.log('\nAlignment:');
console.log('  ', dnaResult.alignedSeq1);
console.log('  ', dnaResult.alignedSeq2);

// Show match line
let matchLine = '  ';
for (let i = 0; i < dnaResult.alignedSeq1.length; i++) {
  if (dnaResult.alignedSeq1[i] === dnaResult.alignedSeq2[i]) {
    matchLine += '|';
  } else if (dnaResult.alignedSeq1[i] === '-' || dnaResult.alignedSeq2[i] === '-') {
    matchLine += ' ';
  } else {
    matchLine += '.';
  }
}
console.log(matchLine);
console.log('\nIdentity:', (dnaResult.identity * 100).toFixed(1) + '%');

// ============================================================================
// EXAMPLE 4: Multiple Sequence Alignment
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 4: Multiple Sequence Alignment');
console.log('=' .repeat(70));

const fastaInput = `
>Sequence1
MKFLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGVQLTLF
>Sequence2
MKTLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGVQLTLF
>Sequence3
MKFLILLNILCLFPVLADNHGVGPQGAGVDPITFDINSNQTGVQLTLF
>Sequence4
MKFLILLFNILCPVLAADNHGVGPQGASGVDPITFDINNQTGVQLTLF
`;

const sequences = parseFasta(fastaInput);
console.log('\nInput sequences:');
sequences.forEach(seq => {
  console.log(`  ${seq.id}: ${seq.sequence.length} aa`);
});

const msaResult = align(fastaInput);

console.log('\nMultiple Sequence Alignment:');
console.log('-'.repeat(70));
msaResult.sequences.forEach(seq => {
  console.log(`${seq.id.padEnd(12)} ${seq.sequence}`);
});

console.log('\nStatistics:');
console.log('  Alignment length:', msaResult.alignmentLength);
console.log('  Number of sequences:', msaResult.numSequences);
console.log('  Average identity:', (msaResult.averageIdentity * 100).toFixed(1) + '%');

// ============================================================================
// EXAMPLE 5: Output Formats
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 5: Output Formats');
console.log('=' .repeat(70));

// FASTA format
console.log('\n--- FASTA Format ---');
const fastaOutput = align(fastaInput, { outputFormat: 'fasta' });
console.log(fastaOutput);

// CLUSTAL format
console.log('\n--- CLUSTAL Format ---');
const clustalOutput = align(fastaInput, { outputFormat: 'clustal' });
console.log(clustalOutput);

// ============================================================================
// EXAMPLE 6: Custom Gap Penalties
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 6: Custom Gap Penalties');
console.log('=' .repeat(70));

const testSeq1 = 'ACDEFGHIKLMNPQRSTVWY';
const testSeq2 = 'ACDEFIKLMNPQRSTVWY';

console.log('\nSequence 1:', testSeq1);
console.log('Sequence 2:', testSeq2);

// Default gap penalties
const defaultGaps = needlemanWunsch(testSeq1, testSeq2);
console.log('\nDefault gap penalties (open=-10, extend=-1):');
console.log('  ', defaultGaps.alignedSeq1);
console.log('  ', defaultGaps.alignedSeq2);
console.log('  Score:', defaultGaps.score);

// Severe gap penalty
const severeGaps = needlemanWunsch(testSeq1, testSeq2, { gapOpen: -20, gapExtend: -5 });
console.log('\nSevere gap penalties (open=-20, extend=-5):');
console.log('  ', severeGaps.alignedSeq1);
console.log('  ', severeGaps.alignedSeq2);
console.log('  Score:', severeGaps.score);

// Mild gap penalty
const mildGaps = needlemanWunsch(testSeq1, testSeq2, { gapOpen: -5, gapExtend: -0.5 });
console.log('\nMild gap penalties (open=-5, extend=-0.5):');
console.log('  ', mildGaps.alignedSeq1);
console.log('  ', mildGaps.alignedSeq2);
console.log('  Score:', mildGaps.score);

// ============================================================================
// EXAMPLE 7: Hemoglobin Family Alignment
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 7: Hemoglobin Family Alignment');
console.log('=' .repeat(70));

const hemoglobinFasta = `
>HBA_HUMAN Alpha-globin
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLL
SHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
>HBB_HUMAN Beta-globin
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLST
PDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPE
NFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH
>MYG_HUMAN Myoglobin
MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLK
SEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPI
KYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG
`;

const hemoResult = align(hemoglobinFasta);

console.log('\nHemoglobin Family Multiple Sequence Alignment:');
console.log('-'.repeat(70));

// Print alignment in blocks
const blockSize = 60;
const alignLen = hemoResult.sequences[0].sequence.length;

for (let start = 0; start < alignLen; start += blockSize) {
  console.log(`\nPositions ${start + 1}-${Math.min(start + blockSize, alignLen)}:`);
  for (const seq of hemoResult.sequences) {
    const block = seq.sequence.substring(start, start + blockSize);
    console.log(`${seq.id.substring(0, 10).padEnd(12)} ${block}`);
  }
}

console.log('\n\nAlignment Statistics:');
console.log('  Total length:', hemoResult.alignmentLength, 'positions');
console.log('  Average pairwise identity:', (hemoResult.averageIdentity * 100).toFixed(1) + '%');

// ============================================================================
// EXAMPLE 8: Sequence Type Detection
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('EXAMPLE 8: Sequence Type Detection');
console.log('=' .repeat(70));

const proteinSeq = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH';
const dnaSeq = 'ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG';

console.log('\nProtein sequence detected as:', detectSequenceType(proteinSeq));
console.log('DNA sequence detected as:', detectSequenceType(dnaSeq));

// ============================================================================
// PERFORMANCE TEST
// ============================================================================

console.log('\n' + '=' .repeat(70));
console.log('PERFORMANCE TEST');
console.log('=' .repeat(70));

// Generate random sequences
function randomSequence(length, alphabet = 'ACDEFGHIKLMNPQRSTVWY') {
  let seq = '';
  for (let i = 0; i < length; i++) {
    seq += alphabet[Math.floor(Math.random() * alphabet.length)];
  }
  return seq;
}

// Test pairwise alignment performance
const lengths = [50, 100, 200, 500];
console.log('\nPairwise alignment (Needleman-Wunsch):');

for (const len of lengths) {
  const s1 = randomSequence(len);
  const s2 = randomSequence(len);
  
  const start = Date.now();
  needlemanWunsch(s1, s2);
  const elapsed = Date.now() - start;
  
  console.log(`  Length ${len}: ${elapsed}ms`);
}

// Test MSA performance
console.log('\nMultiple sequence alignment:');
const numSeqs = [3, 5, 10];
const seqLen = 100;

for (const n of numSeqs) {
  let fasta = '';
  for (let i = 0; i < n; i++) {
    fasta += `>Seq${i + 1}\n${randomSequence(seqLen)}\n`;
  }
  
  const start = Date.now();
  align(fasta);
  const elapsed = Date.now() - start;
  
  console.log(`  ${n} sequences x ${seqLen} aa: ${elapsed}ms`);
}

console.log('\n' + '=' .repeat(70));
console.log('All examples completed successfully!');
console.log('=' .repeat(70));
