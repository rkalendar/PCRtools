// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const T_TAIL = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"; // 30-nt T-tail added to reverse primers
const MIN_DIMER_BASES = 5;

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------
function analysis() {
    let probe    = numToDouble(document.getElementById('probe').value);
    let distance = numToDouble(document.getElementById('distance').value);

    // Clamp probe length to [35, 40]
    probe = Math.min(Math.max(probe, 35), 40);

    // Clamp minimum distance to 5
    if (distance < 5) distance = 5;

    const rawInput = document.getElementById('inputText').value;
    const normalized = rawInput.toLowerCase().trim();

    if (normalized.length < 2) return;

    // Parse input into parallel name / sequence arrays
    const { nameSeqs, seqs } = parseSequences(normalized);

    const result = ["", ""];

    for (let n = 0; n < seqs.length; n++) {
        const header = "Location(ID)\tLinguistic_Complexity(%)\tSequence(5'-3')\n"
                     + nameSeqs[n].toUpperCase() + ":\n";

        const boundsArr  = [];
        const maskedArr  = [];
        const cleanSeq   = parseSeq(normalized === seqs[n] ? seqs[n] : seqs[n], boundsArr, maskedArr);
        const seqLen     = cleanSeq.length;

        // Build repeat-mask and apply manually-masked regions
        // NOTE: RepeatMask and LingComplexity are defined in the external dna.js
        let mask = RepeatMask(cleanSeq);
        for (let y = 0; y < maskedArr.length; y += 2) {
            mask = mask.fill(1, maskedArr[y], maskedArr[y + 1]);
        }

        // Forward-pass bounds: [b0, b1, b2, b3] as-is
        const fwdBounds = boundsArr.slice();

        // Reverse-pass bounds: coordinates must be mirrored onto the complement strand
        // Original logic: x1=lseqs-b[3], x2=lseqs-b[2]+1-probe (≡ end before probe subtraction = lseqs-b[2])
        //                 x3=lseqs-b[1], x4=lseqs-b[0]+1-probe (≡ lseqs-b[0])
        const b = boundsArr;
        const revBounds = [seqLen - b[3], seqLen - b[2], seqLen - b[1], seqLen - b[0]];

        // Forward direction
        result[0] += header + scanPrimerPairs(
            cleanSeq, mask, fwdBounds, probe, distance, seqLen,
            /* isReverse */ false
        );

        // Reverse (complement) direction
        const complement = complementDNA(cleanSeq);
        const reverseMask = [...mask].reverse();

        result[1] += header + scanPrimerPairs(
            complement, reverseMask, revBounds, probe, distance, seqLen,
            /* isReverse */ true
        );
    }

    return result;
}

// ---------------------------------------------------------------------------
// Sequence parsing helpers
// ---------------------------------------------------------------------------

/**
 * Detect input format (FASTA / tab-delimited / space-delimited / plain) and
 * return parallel arrays of sequence names and sequences.
 */
function parseSequences(s) {
    // FASTA format (lines starting with '>')
    if (s.includes('>')) {
        return parseFasta(s);
    }

    // Tab-delimited
    if (s.includes('\t')) {
        return parseDelimited(s, '\t');
    }

    // Space-delimited
    if (s.includes(' ')) {
        return parseDelimited(s, ' ');
    }

    // Plain sequence with no name
    return { nameSeqs: ['SEQ1'], seqs: [s] };
}

function parseFasta(s) {
    const nameSeqs = [];
    const seqs     = [];
    let pos = s.indexOf('>');

    while (pos !== -1) {
        // Find end of header line
        let lineEnd = s.indexOf('\r', pos + 1);
        const lf    = s.indexOf('\n', pos + 1);
        const sp    = s.indexOf(' ',  pos + 1);

        if (lineEnd === -1 && lf !== -1) lineEnd = lf;
        if (lineEnd === -1 && sp  > -1)  lineEnd = sp;

        if (lineEnd === -1) break;

        const name   = s.substring(pos + 1, lineEnd).trim();
        const nextGt = s.indexOf('>', lineEnd + 1);
        const seqStr = nextGt !== -1
            ? s.substring(lineEnd + 1, nextGt).trim()
            : s.substring(lineEnd + 1).trim();

        nameSeqs.push(name);
        seqs.push(seqStr);
        pos = nextGt;
    }

    return { nameSeqs, seqs };
}

function parseDelimited(s, delimiter) {
    const nameSeqs = [];
    const seqs     = [];
    const lines    = s.split(/[\r\n]+/);

    for (const line of lines) {
        const idx = line.indexOf(delimiter);
        if (idx === -1) continue;
        nameSeqs.push(line.substring(0, idx).trim());
        seqs.push(line.substring(idx));
    }

    return { nameSeqs, seqs };
}

// ---------------------------------------------------------------------------
// Core primer-scanning logic (shared by forward & reverse passes)
// ---------------------------------------------------------------------------
function scanPrimerPairs(seq, mask, bounds, probe, distance, seqLen, isReverse) {
    const fwd = collectCandidates(seq, mask, bounds[0], bounds[1] - probe + 1, probe, seqLen, isReverse);
    const rev = collectCandidates(seq, mask, bounds[2], bounds[3] - probe + 1, probe, seqLen, isReverse);

    let output = "";
    const sameRegion = bounds[0] === bounds[2] && bounds[1] === bounds[3];

    for (let i = 0; i < fwd.primers.length; i++) {
        for (let j = 0; j < rev.primers.length; j++) {
            const p1 = fwd, p2 = rev;

            if (sameRegion && p2.x1[j] + 1 <= p1.x2[i] + distance) continue;

            if (dimerLook(p1.primers[i], p2.primers[j], MIN_DIMER_BASES) !== -1) continue;

            if (isReverse) {
                output += p2.names[j] + "\t" + p2.lc[j] + "\t" + p2.primers[j] + "\n";
                output += p1.names[i] + "\t" + p1.lc[i] + "\t" + T_TAIL + p1.primers[i] + T_TAIL + "\t\n\n";
            } else {
                output += p1.names[i] + "\t" + p1.lc[i] + "\t" + p1.primers[i] + "\n";
                output += p2.names[j] + "\t" + p2.lc[j] + "\t" + T_TAIL + p2.primers[j] + T_TAIL + "\t\n\n";
            }
            break;
        }
    }

    return output;
}

/**
 * Slide a probe-length window across [start, end) and collect valid primer candidates.
 */
function collectCandidates(seq, mask, start, end, probe, seqLen, isReverse) {
    const primers = [], names = [], lc = [], x1 = [], x2 = [];

    // Pre-fill sliding window count for masked positions
    let maskedCount = 0;
    for (let y = start; y < start + probe - 1; y++) {
        if (mask[y] > 0) maskedCount++;
    }

    for (let pos = start; pos < end; pos++) {
        if (mask[pos + probe - 1] > 0) maskedCount++;

        if (maskedCount === 0) {
            const candidate = seq.substring(pos, pos + probe);

            if (dimerLook(candidate, candidate, MIN_DIMER_BASES) === -1 && hasNoRepeat(candidate)) {
                const complexity = LingComplexity(candidate);
                const label = isReverse
                    ? "R_" + (seqLen - pos) + "-" + (seqLen - (pos + probe) + 1)
                    : "F_" + (pos + 1)      + "-" + (pos + probe + 1);

                primers.push(candidate);
                names.push(label);
                lc.push(complexity);
                x1.push(pos + 1);
                x2.push(pos + probe + 1);
            }
        }

        if (mask[pos] > 0) maskedCount--;
    }

    return { primers, names, lc, x1, x2 };
}

// ---------------------------------------------------------------------------
// Sequence utilities
// ---------------------------------------------------------------------------

// Valid IUPAC nucleotide characters (lowercase)
const VALID_BASES = new Set(['a','t','c','g','r','y','m','k','w','b','d','v','h','s','n']);

/**
 * Strip non-nucleotide characters and resolve [bracket] / slash annotations.
 * Populates `bounds` (region boundaries) and `masked` (hard-masked ranges).
 */
function parseSeq(str, bounds, masked) {
    let clean = "";
    const bracketPositions = [];

    for (const chr of str) {
        if (VALID_BASES.has(chr)) {
            clean += chr;
        } else if (chr === 'u') {
            clean += 't';
        } else if (chr === 'i') {
            clean += 'g';
        } else if (chr === '[') {
            bracketPositions.push(clean.length);        // positive = open
        } else if (chr === ']') {
            bracketPositions.push(-clean.length);       // negative = close
        } else if (chr === '/') {
            masked.push(clean.length);
        }
    }

    if (bracketPositions.length === 0) {
        bounds.push(0, clean.length, 0, clean.length);
    } else {
        for (let i = 0; i < bracketPositions.length - 1; i++) {
            if (bracketPositions[i] >= 0) {
                bounds.push(bracketPositions[i]);
                for (let j = i + 1; j < bracketPositions.length; j++) {
                    if (bracketPositions[j] < 0) {
                        bounds.push(-bracketPositions[j]);
                        break;
                    }
                }
            }
        }
        if (bounds.length === 2) {
            bounds.push(bounds[0], bounds[1]);
        }
    }

    return clean;
}

function reverseDNA(source) {
    return source.split("").reverse().join("");
}

// IUPAC complement lookup table (char-code indexed)
const COMPLEMENT_TABLE = (() => {
    const d = new Array(128).fill(0);
    d[97]  = 116; // a -> t
    d[98]  = 118; // b -> v
    d[99]  = 103; // c -> g
    d[100] = 104; // d -> h
    d[103] = 99;  // g -> c
    d[104] = 100; // h -> d
    d[105] = 99;  // i -> c (inosine)
    d[107] = 109; // k -> m
    d[109] = 107; // m -> k
    d[110] = 110; // n -> n
    d[114] = 121; // r -> y
    d[115] = 115; // s -> s
    d[116] = 97;  // t -> a
    d[117] = 97;  // u -> a
    d[118] = 98;  // v -> b
    d[119] = 119; // w -> w
    d[121] = 114; // y -> r
    return d;
})();

function antiSenseDNA(source) {
    const result = [];
    for (let i = 0; i < source.length; i++) {
        const mapped = COMPLEMENT_TABLE[source.charCodeAt(i)];
        if (mapped) result.push(String.fromCharCode(mapped));
    }
    return result.join('');
}

function complementDNA(source) {
    return antiSenseDNA(reverseDNA(source));
}

// ---------------------------------------------------------------------------
// Repeat / complexity filters
// ---------------------------------------------------------------------------

const REPEAT_PATTERNS = ['gggg', 'ccccc', 'aaaaa', 'ttttt', 'gcgcgc', 'cgcgcg', 'atatat', 'tatata'];

/**
 * Returns -1 if the sequence contains no problematic repeat motifs, 0 if it does.
 * (Consistent with DimerLook convention: -1 = pass, 0 = fail.)
 */
function hasNoRepeat(str) {
    // BUG FIX: original loop used `i < t.length - 1`, skipping the last pattern 'tatata'
    for (const pattern of REPEAT_PATTERNS) {
        if (str.includes(pattern)) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Dimer detection
// ---------------------------------------------------------------------------

/**
 * Check for primer-dimer formation between r1 and r2.
 * Returns -1 if no dimer found (safe), 0 if a dimer of >= n_base is found.
 */
function dimerLook(r1, r2, n_base) {
    n_base = Math.min(Math.max(n_base, 3), 10);

    let seq1, seq2;
    if (r1.length < r2.length) {
        seq1 = r2;
        seq2 = antiSenseDNA(reverseDNA(r1));
    } else {
        seq1 = r1;
        seq2 = antiSenseDNA(reverseDNA(r2));
    }

    const l1     = seq1.length;
    const l2     = seq2.length;
    const initK  = n_base - 1;

    for (let y = 0; y < l2 - initK; y++) {
        const fragment = seq2.substring(y, y + initK);
        let x = -1;

        while (true) {
            x = seq1.indexOf(fragment, x + 1);
            if (x === -1) break;

            let matchLen  = initK;
            let baseCount = initK;

            // Extend match forward
            while (x + matchLen <= l1 - 1 && y + matchLen <= l2 - 1) {
                if (seq1[x + matchLen] === seq2[y + matchLen]) {
                    matchLen++;
                    baseCount++;
                } else if (
                    x + matchLen + 1 <= l1 - 1 &&
                    y + matchLen + 1 <= l2 - 1 &&
                    seq1[x + matchLen + 1] === seq2[y + matchLen + 1]
                ) {
                    matchLen  += 2;
                    baseCount += 1;
                } else {
                    break;
                }
            }

            // Extend match backward
            let backSteps = 0;
            while (x - backSteps - 1 >= 0 && y - backSteps - 1 >= 0) {
                if (seq1[x - backSteps - 1] === seq2[y - backSteps - 1]) {
                    matchLen++;
                    baseCount++;
                    backSteps++;
                } else if (
                    x - backSteps - 2 >= 0 &&
                    y - backSteps - 2 >= 0 &&
                    seq1[x - backSteps - 2] === seq2[y - backSteps - 2]
                ) {
                    matchLen  += 2;
                    baseCount += 1;
                    backSteps += 2;
                } else {
                    break;
                }
            }

            if (baseCount > n_base) return 0;
        }
    }

    return -1;
}


// ---------------------------------------------------------------------------
// Numeric parser
// ---------------------------------------------------------------------------

/**
 * Parse a string to a float, accepting both '.' and ',' as decimal separators.
 * BUG FIX: original used replace(',', '.') which only replaced the *first* comma.
 */
function numToDouble(val) {
    val = val.replace(/,/g, '.').replace(/[^0-9.-]/g, '');
    const out = parseFloat(val);
    return isNaN(out) ? 0 : out;
}
