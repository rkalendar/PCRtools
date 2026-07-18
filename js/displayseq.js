'use strict';
// Wrapped in an IIFE so the helpers below (parseFasta/parseSequences/
// parseDelimited) stay private and cannot collide with same-named globals
// from other scripts loaded on the same page (e.g. muscle.js also declares a
// global `parseFasta` with a different signature). Only `Display` is exported.
(function (global) {
function Display(s, mode = "full") {
    if (!s || s.length < 2) return "";

    s = s.toLowerCase();
    const { names, seqs } = parseSequences(s);

    let totalSeq = "";
    let count = 0;
    for (const seq of seqs) {
        const dna = DNA(seq);
        if (dna.length > 0) { totalSeq += dna; count++; }
    }

    if (count === 0) return "";

    const len = totalSeq.length;
    const cg  = CG(totalSeq).toFixed(1);

    const formats = {
        count:  ` ${count} : ${len} nt`,
        length: ` ${len} nt`,
        cg:     ` ${cg}%CG`,
        full:   ` ${count} : ${len} nt (${cg}%CG)`,
    };

    return formats[mode] ?? formats.full;
}

function parseSequences(s) {
    // FASTA: starts with ">"
    if (s.includes(">")) return parseFasta(s);
    // Tab-separated
    if (s.includes("\t")) return parseDelimited(s, "\t");
    // Space-separated
    if (s.includes(" "))  return parseDelimited(s, " ");
    // Plain sequence — no header
    return { names: ["Seq1"], seqs: [s] };
}

function parseFasta(s) {
    const names = [], seqs = [];
    const blocks = s.split(/(?=>)/);         // split before each ">"
    for (const block of blocks) {
        if (!block.startsWith(">")) continue;
        const nl = block.search(/[\r\n ]/);  // end of header line
        if (nl === -1) continue;
        names.push(block.substring(1, nl).trim());
        seqs.push(block.substring(nl + 1));
    }
    return { names, seqs };
}

function parseDelimited(s, delimiter) {
    const names = [], seqs = [];
    const lines = s.split(/\r?\n/);
    for (const line of lines) {
        const idx = line.indexOf(delimiter);
        if (idx === -1) continue;
        names.push(line.substring(0, idx).trim());
        seqs.push(line.substring(idx));
    }
    return { names, seqs };
}

// --- Protein / mixed-alphabet support (used by the MSA page) ---------------
// Classify one sequence as "dna" or "protein" using the same heuristic as
// muscle.js detectSequenceType(): >90% plain nucleotide letters => DNA.
function classifySeq(seq) {
    const residues = seq.toLowerCase().replace(/[^a-z]/g, "");
    if (!residues.length) return null;
    const nt = (residues.match(/[atgcnu]/g) || []).length;
    return nt / residues.length > 0.9 ? "dna" : "protein";
}

// Count amino-acid residues (every letter is a valid IUPAC residue code;
// gaps, stop "*", digits and whitespace are ignored).
function countResidues(seq) {
    return (seq.toLowerCase().match(/[a-z]/g) || []).length;
}

// Like Display(), but reports DNA and protein sequences separately:
//   DNA     ->  " N : LEN nt (CG%CG)"
//   protein ->  " N : LEN aa"
// When both kinds are present, the two summaries are concatenated.
function DisplayBio(s) {
    if (!s || s.length < 2) return "";
    const { seqs } = parseSequences(s.toLowerCase());

    let dnaCount = 0, dnaConcat = "";
    let aaCount = 0, aaLen = 0;
    for (const seq of seqs) {
        const type = classifySeq(seq);
        if (type === "dna") {
            const dna = DNA(seq);
            if (dna.length) { dnaConcat += dna; dnaCount++; }
        } else if (type === "protein") {
            const n = countResidues(seq);
            if (n) { aaLen += n; aaCount++; }
        }
    }

    const parts = [];
    if (dnaCount) parts.push(` ${dnaCount} : ${dnaConcat.length} nt (${CG(dnaConcat).toFixed(1)}%CG)`);
    if (aaCount)  parts.push(` ${aaCount} : ${aaLen} aa`);
    return parts.join(" ");
}

global.Display = Display;
global.DisplayBio = DisplayBio;
})(typeof window !== "undefined" ? window : globalThis);
