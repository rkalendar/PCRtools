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