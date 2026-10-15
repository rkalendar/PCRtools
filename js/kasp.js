'use strict';
function analysis() {
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);
    let snp = parseInt(document.getElementById('snp').value);

    let overlapping = document.getElementById('overlapping').checked;
    let maskrepeats = document.getElementById('repeats').checked;
    let ctconvert = document.getElementById('ctconvert').checked;
    const lowcomplexity = document.getElementById('lowcomplexity').checked;
    // Allele sizing reads the genotype off a gel, so it replaces the fluorophore
    // tails rather than adding to them.
    const sizing = {
        on: document.getElementById('sizingmode')?.checked || false,
        step: parseInt(document.getElementById('sizestep')?.value ?? '5', 10),
        maxlen: parseInt(document.getElementById('sizemaxlen')?.value ?? '45', 10),
        mismatch: document.getElementById('sizemismatch')?.checked || false
    };
    if (!Number.isFinite(sizing.step) || sizing.step < 1) { sizing.step = 1; }
    if (!Number.isFinite(sizing.maxlen) || sizing.maxlen < 20) { sizing.maxlen = 20; }


    let ends3 = document.getElementById('end3com').value.toLowerCase().trim();   // "swh ssw wsh sww www";
    if (DNA(ends3).length < 1) { ends3 = "n"; }
    const raw = document.getElementById('minlc')?.value;
    let minlc = Number.parseInt(raw ?? '70', 10);
    if (!Number.isFinite(minlc)) minlc = 70;
    minlc = Math.min(90, Math.max(10, minlc));

    if (minlen < 12) { minlen = 12; }
    if (maxlen < 12) { maxlen = 12; }
    if (minlen > 100) { minlen = 100; }
    if (maxlen > 100) { maxlen = 100; }
    if (minlen > maxlen) { maxlen = minlen; }
    if (mintm < 37) { mintm = 37; }
    if (mintm > 80) { mintm = 80; }
    if (snp < 1) { snp = 1; }


    let primeroverlap = 12;
    if (overlapping) { primeroverlap = minlen; }

    let ReadResult = ReadingSeq(document.getElementById('inputTails').value);
    let tails = ReadResult.seqs;
    let n_tails = tails.length;
    for (let n = 0; n < n_tails; n++) {
        tails[n] = DNA(tails[n]).toUpperCase();
    }

    ReadResult = ReadingSeq(document.getElementById('inputText').value);
    let name_seq = ReadResult.name_seq;
    let seqs = ReadResult.seqs;
    let n_seq = seqs.length;

    const result = ["", ""];
    let resultarea1 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n";
    let resultarea2 = "PrimerID\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\tFragment_Size(bp)/Tm(°C)\n\n";
    for (let n = 0; n < n_seq; n++) {
        resultarea1 += name_seq[n] + ":\n";
        // The validator promises that a record without a variant bracket is skipped: say so,
        // rather than print empty primer rows for it.
        if (seqs[n].indexOf("[") < 0) {
            const note = "// " + name_seq[n] + ": no variant is marked in square brackets, so no assay can be designed; the record is skipped.\n\n";
            resultarea1 += note;
            resultarea2 += note;
            continue;
        }

        let resl1 = AsPCR(name_seq[n], "Forward", snp, seqs[n], tails, maskrepeats, lowcomplexity, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3, sizing);
        resultarea1 += resl1.resultarea1;
        resultarea2 += resl1.resultarea2;

        //REVERSE DIRECTION
        let resl2 = AsPCR(name_seq[n], "Reverse", snp, ComplementDNA(seqs[n]), tails, maskrepeats, lowcomplexity, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3, sizing);
        resultarea2 += "\n";
        resultarea1 += resl2.resultarea1;
        resultarea2 += resl2.resultarea2;
    }
    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

// Allele-specific primers, one per allele, all ending at the same place: the variant's 3'
// window. Only the 5' end is free, so each primer grows from minlen towards maxlen until
// its Tm enters the window, and it stops at the start of the sequence or at a masked base
// (msk: one mask per allele template; 1 repeat or low complexity, 2 a /.../ exclusion,
// 4 another marker). If an allele gets no primer, or the primers pair with each other, no
// set can be made: fail says why and nothing from it is printed. A primer whose Tm cannot
// be brought into the window is kept; why[j] records where its growth ended. ends names
// the two ends of the sequence as the user sees them: the reverse pass reads it backwards.
function ProbeDesignForward(snp, tails, s1, s2, pol, minlen, maxlen, mintm, maxtm, nm, tpl, msk, ends = ["start", "end"]) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const ntails = tails.length;

    let n = pol.length;
    let minpol = 999;
    let z = 0;
    let k = snp;
    let kmer = minlen;
    let x1 = s1.length;
    let s = new Array(n).fill("");
    let t = new Array(n).fill("");
    let p = new Array(n).fill("");
    let fpr = new Array(n).fill("");
    let fpn = new Array(n).fill("");
    let ftm = new Array(n).fill(0.0);
    let fcg = new Array(n).fill(0.0);
    let flc = new Array(n).fill(0);
    let fpx = new Array(n).fill(0);
    let why = new Array(n).fill(null);
    let fail = "";

    for (let i = 0; i < n; i++) {
        if (minpol > pol[i].length) { minpol = pol[i].length; }
        // Bisulphite mode reads each allele from its own converted template (tpl);
        // pol still names the alleles as they were entered.
        s[i] = tpl ? (tpl[i][1] + tpl[i][2]) : (pol[i] + s2);
        t[i] = tpl ? tpl[i].join("") : (s1 + pol[i] + s2);
    }
    if (kmer > minpol) { kmer = minpol; }
    if (kmer == 0) { kmer = 12; }

    if (minlen <= minpol) { // long insertions - Haplotypes
        let w = 0;
        for (let i = 0; i < minpol - minlen + 1; i++) {
            for (let j = 0; j < n; j++) {
                p[j] = s[j].substring(i, i + minlen);
            }
            let j = countUniqueSNPs(p);
            if (j > w) {
                w = j; z = i;
                k = minlen;
            }
        }
        x1 = x1 + z + k;
    }
    else {
        let w = 0;
        for (let i = 0; i < kmer + snp; i++) {
            for (let j = 0; j < n; j++) {
                p[j] = s[j].substring(i, i + snp);
            }
            let j = countUniqueSNPs(p);
            if (j > w) { w = j; z = i; }
            if (j > 2) { break; }
        }
        x1 = x1 + z + snp;
    }

    for (let j = 0; j < n; j++) {
        const lab = pol[j] === "" ? "-" : pol[j].toUpperCase();
        // The 3' window may run past the variant; it needs that much sequence after it.
        if (x1 > t[j].length) {
            fail = "the variant is too close to the " + ends[1] + " of the sequence: the allele-specific primer for " + lab + " would need " + (x1 - t[j].length) + " more nt past it";
            break;
        }
        let first = 0;  // shortest length above the Tm floor
        let pick = 0;   // shortest length inside the Tm window
        let last = 0;   // longest length tried
        let stop = "";  // what ended the growth before maxlen: "start" or a mask value
        for (let i = minlen; i <= maxlen; i++) {
            const x = x1 - i;
            if (x < 0) { stop = "start"; break; }
            // No base of the primer may be masked, not only its 3' end.
            let m = 0;
            for (let y = x; y < (i === minlen ? x1 : x + 1); y++) { m |= msk[j][y]; }
            if (m > 0) { stop = m; break; }
            const tm = DNA_Tm(t[j].substring(x, x1), salt, Mg_M, p_mkM);
            last = i;
            if (tm > mintm) {
                if (first === 0) { first = i; }
                if (tm <= maxtm) { pick = i; break; }
            }
        }
        if (last === 0) {
            fail = stop === "start"
                ? "the variant is too close to the " + ends[0] + " of the sequence: only " + x1 + " nt are available for the allele-specific primer for " + lab + ", which needs at least " + minlen
                : "the allele-specific primer for " + lab + " cannot avoid " + MaskWhat(stop) + ": even at " + minlen + " nt it covers some";
            break;
        }
        // Outside the window the position still fixes the primer: the one nearest the
        // floor is kept and reported with a warning.
        const len = pick || first || last;
        why[j] = pick ? null : { last, stop };
        fpx[j] = x1 - len;
        fpr[j] = t[j].substring(fpx[j], x1);
        ftm[j] = DNA_Tm(fpr[j], salt, Mg_M, p_mkM);
        flc[j] = LingComplexity2(fpr[j]);
        fcg[j] = CG(fpr[j]);
        fpn[j] = nm + (fpx[j] + 1) + "-" + (fpx[j] + fpr[j].length) + "_ASP" + "_" + pol[j].toUpperCase();
        // The discriminating base(s) are what the primer is for, so they are the
        // one part of the body written in upper case.
        const va = Math.max(fpx[j], s1.length) - fpx[j];
        const vb = Math.min(fpx[j] + fpr[j].length, s1.length + pol[j].length) - fpx[j];
        if (vb > va) { fpr[j] = fpr[j].substring(0, va) + fpr[j].substring(va, vb).toUpperCase() + fpr[j].substring(vb); }
    }

    if (!fail) {
        if (ntails > 0) {
            for (let i = 0; i < Math.min(ntails, n); i++) {
                fpr[i] = tails[i] + fpr[i];
            }
        }
        // Primers that pair are rejected, not shown: the whole set, since one allele
        // cannot be called without the others.
        const d = DimerPair(fpr);
        if (d) {
            fail = DimerWhy(fpn, d) + (ntails > 0 ? " (tails included)" : "") + ", so no set is reported";
        }
    }
    return { fpr, fpn, ftm, flc, fpx, fcg, x1, why, fail };
}

// The first pair of primers (a self-pair when i === j) that forms a dimer, or null.
function DimerPair(pr) {
    const min3 = 6;
    for (let i = 0; i < pr.length; i++) {
        for (let j = i; j < pr.length; j++) {
            if (DimerLook2(pr[i].toLowerCase(), pr[j].toLowerCase(), min3) > 0) { return [i, j]; }
        }
    }
    return null;
}

// Names the primers of a DimerPair() result, for messages.
function DimerWhy(fpn, d) {
    return d[0] === d[1] ? "the allele-specific primer " + fpn[d[0]] + " forms a self-dimer" : "the allele-specific primers " + fpn[d[0]] + " and " + fpn[d[1]] + " pair with each other";
}

// What a mask value stands for, for messages.
function MaskWhat(m) {
    if (m & 4) { return "another marker"; }
    if (m & 2) { return "a /.../ excluded region"; }
    return "repeat or low-complexity sequence";
}

function PrimerDesignReverse(cs, msk, prlist, mplist, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, nm, e5, e3) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 6;


    const le5 = e5.length;
    const s5 = GeneratorVariants(e5).join(' ');

    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;


    const lseqs = cs.length;


    let rn = -1;
    let q = 0;

    for (let y = x1; y < x1 + minlen - 1; y++) {
        if (msk[y] > 0) { q++; }
    }

    for (let x = x1; x < x2; x++) {
        if (msk[x + minlen - 1] > 0) { q++; }
        if (q === 0) {
            let q1 = 1;
            let x5 = x;
            let tm1 = maxtm + 1;
            let lc = 0;

            let s1 = cs.substring(x5 + minlen - le3, x5 + minlen);

            if (s3.indexOf(s1) > -1) {
                for (let d = 0; d < 1 + maxlen - minlen; d++) {
                    x5 = x - d;
                    // The 5' extension may not reach a masked base either.
                    if (x5 >= x1 && !(msk[x5] > 0)) {
                        s1 = cs.substring(x5, x5 + minlen + d);
                        if (s5.indexOf(s1.substring(0, le5)) > -1) {
                            tm1 = DNA_Tm(s1, salt, Mg_M, p_mkM);
                        }
                        if (tm1 >= mintm) {
                            break;
                        }
                    } else { break; }
                }

                if (tm1 >= mintm && tm1 <= maxtm) {
                    if (ssrrepeat(s1) === -1) {
                        lc = LingComplexity2(s1);
                        if (lc >= minlc) {
                            if (DimerLook2(s1, s1, min3) === 0) {
                                q1 = 0;
                            }
                        }
                    }
                }

                if (q1 === 0) {
                    if (rn > -1) {
                        const p2 = rpx[rn] + rpr[rn].length - primeroverlap;
                        if (x5 < p2) {
                            if (rlc[rn] < lc) {
                                rn--; q1 = 1;
                            }
                        } else { q1 = 0; }
                    }

                    if (mplist.length > 0 && q1 === 0) {
                        for (let j = 0; j < mplist.length; j++) {
                            if (DimerLook2(s1, prlist[j], min3) > 0) {
                                q1 = 1; break;
                            }
                        }
                    }

                    if (prlist.length > 0 && q1 === 0) {
                        for (let j = 0; j < prlist.length; j++) {
                            if (DimerLook2(s1, prlist[j].toLowerCase(), min3) > 0) {
                                q1 = 1; break;
                            }
                        }
                    }


                    if (q1 === 0) {
                        rn++;
                        rpr[rn] = s1;
                        rtm[rn] = tm1;
                        rlc[rn] = lc;
                        rcg[rn] = CG(s1);
                        rpx[rn] = x5;
                        rpn[rn] = nm + (lseqs - x5) + "-" + (lseqs - (x5 + s1.length) + 1);
                    }
                }
            }
        }
        if (msk[x] > 0) { q--; }
    }
    return rn;
}

// Reads one record. The first bracket is the marker this pass genotypes: s2 holds its
// alleles as written, ref the first one entered, s1 the flank before it and s3 everything
// after it. Later markers are written into s3 at their first entered allele, so s1 + ref +
// s3 is the reference template, and each is listed in mk (at, len in s3; lo, hi its
// shortest and longest allele) so that primers can be kept off it. A '/' outside brackets
// opens or closes a /.../ exclusion (inside brackets it only separates alleles); bx holds
// those bounds on the reference template. nb counts the brackets. rev: the text is the
// reverse complement of the record, so a list of alleles reads backwards.
function Seqs(str, rev) {
    let s1 = "";
    let s2 = "";
    let s3 = "";
    let ref = "";
    let nb = 0;
    const mk = [];
    const bx = [];
    // Without a complete bracket before anything else there is nothing to genotype.
    const none = { s1: "", s2: "", s3: "", ref: "", mk: [], bx: [], nb: 0 };
    const first = (a) => rev ? ComplementDNA(GeneratorVariants(ComplementDNA(a))[0]) : GeneratorVariants(a)[0];
    const parts = str.match(/\[[^\[\]]*\]|\/|[^\[\]\/]+|[\[\]]/g) || [];
    for (const w of parts) {
        if (w === "[" || w === "]") {
            if (nb === 0) { return none; }
            continue;
        }
        if (w.charAt(0) === "[") {
            // DNASeq, not DNA: U, I and LNA letters are recoded in place rather than
            // dropped, so every coordinate stays that of the sequence as entered.
            const a = w.slice(1, -1).trim().replace(/[uiefjl]/g, c => DNASeq(c));
            const r = first(a);
            if (nb === 0) { s2 = a; ref = r; }
            else {
                const len = GeneratorVariants(a).map(v => v.length);
                mk.push({ at: s3.length, len: r.length, lo: Math.min(...len), hi: Math.max(...len) });
                s3 += r;
            }
            nb++;
        }
        else if (w === "/") { bx.push(s1.length + ref.length + s3.length); }
        else if (nb === 0) { s1 += DNASeq(w); }
        else { s3 += DNASeq(w); }
    }
    return nb === 0 ? none : { s1, s2, s3, ref, mk, bx, nb };
}

// One allele's template after C->T conversion of the input's top strand, as
// [s1, allele, s2] in the sense of the pass. The reverse pass reads the reverse
// complement, so its pieces are converted in top-strand sense and turned back.
function BisulfitePassAllele(s1, allele, s2, rev) {
    if (!rev) { return BisulfiteAllele(s1, allele, s2); }
    const c = BisulfiteAllele(ComplementDNA(s2), ComplementDNA(allele), ComplementDNA(s1));
    return [ComplementDNA(c[2]), ComplementDNA(c[1]), ComplementDNA(c[0])];
}

function AsPCR(name_seq, direction, snp, inseq, tails, maskrepeats, lowcomplexity, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3, sizing) {
    sizing = sizing || { on: false };
    const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    const Kmax = 14;  // Kmax=11 for SSR  Kmax=14 for telomers     
    const kmer2 = 21;  // SSR min length kmer2 = 21 nt 
    let resultarea1 = "";
    let resultarea2 = "";
    // The reverse pass runs on the reverse complement of the input; in bisulphite
    // mode it must still read the input's converted top strand, not the bottom one.
    const rev = direction === "Reverse";
    let rst = Seqs(inseq, rev);
    let s1 = rst.s1;
    let s2 = rst.s3;
    let pol = rst.s2;
    let st = GeneratorVariants(pol);
    // With several markers this pass genotypes one of them: the first going forward,
    // the last going back.
    const several = rst.nb > 1;
    const which = several ? "// " + name_seq + " (" + direction.toLowerCase() + " set): the record has " + rst.nb + " markers; this set genotypes the " +
        (rev ? "last" : "first") + ", [" + (rev ? ComplementDNA(pol) : pol).toUpperCase() + "], and no primer covers the others.\n" : "";
    // A pass that cannot give a usable set says why, in both tabs, in place of the set.
    const refuse = (why) => {
        const note = which + "// " + name_seq + " (" + direction.toLowerCase() + " set): " + why + ".\n\n";
        return { resultarea1: note, resultarea2: note };
    };
    if (rst.nb === 0) { return refuse("no complete [ ] variant bracket"); }
    if (st.every(a => a === "")) { return refuse("the variant bracket holds no allele"); }
    // Messages name the alleles as entered: the reverse pass holds them complemented
    // and, for a list, in reverse order.
    const entered = GeneratorVariants(rev ? ComplementDNA(pol) : pol);
    const asEntered = (i) => rev ? ComplementDNA(st[i]) : st[i];
    const rank = (i) => entered.indexOf(asEntered(i));
    // Positions, product sizes and Ta refer to the reference template: the record with
    // the first entered allele at every marker.
    const ref = rst.ref;
    const iref = Math.max(0, st.indexOf(ref));

    // Bases no primer may cover, on the two flanks: /.../ exclusions (2) and the other
    // markers (4). Repeats and low complexity (1) are added per allele template below.
    const l1 = s1.length;
    const l3 = s2.length;
    const m1 = new Array(l1).fill(0);
    const m3 = new Array(l3).fill(0);
    for (let i = 0; i + 1 < rst.bx.length; i += 2) {
        const a = rst.bx[i];
        const b = rst.bx[i + 1];
        if (a < l1 + ref.length && b > l1) { return refuse("the variant lies inside a /.../ excluded region"); }
        for (let r = a; r < b; r++) {
            if (r < l1) { m1[r] |= 2; }
            else if (r >= l1 + ref.length) { m3[r - l1 - ref.length] |= 2; }
        }
    }
    for (const m of rst.mk) {
        for (let y = m.at; y < m.at + m.len; y++) { m3[y] |= 4; }
        // An insertion point has no letters: masking the base after it keeps a primer
        // from spanning it.
        if (m.len === 0 && m.at < l3) { m3[m.at] |= 4; }
        // After conversion the top-strand base before a marker is C or T by its allele.
        if (ctconvert) {
            if (!rev && m.at > 0) { m3[m.at - 1] |= 4; }
            if (rev && m.at + m.len < l3) { m3[m.at + m.len] |= 4; }
        }
    }
    if (ctconvert && !rev && rst.mk.some(m => m.at === 0) && st.some(a => (s1 + a).endsWith("c"))) {
        return refuse("the next marker follows the variant directly, and after C->T conversion the variant's last base depends on it");
    }

    let tpl = null;
    if (ctconvert) {
        const groups = rev
            ? BisulfiteCollapsedAlleles(ComplementDNA(s2), st.map(a => ComplementDNA(a)), ComplementDNA(s1))
            : BisulfiteCollapsedAlleles(s1, st, s2);
        if (groups.length > 0) {
            // With one marker both passes read it and share the verdict, so the forward pass
            // alone reports it. With several, the reverse pass reads the last marker: it reports
            // its own verdict, with the alleles turned back to the top strand and in the order
            // they were entered.
            const where = several ? (rev ? " (last marker, reverse pass)" : " (first marker, forward pass)") : "";
            const same = groups.map(g => g.slice().sort((a, b) => rank(a) - rank(b))).sort((a, b) => rank(a[0]) - rank(b[0]))
                .map(g => g.map(i => { const a = asEntered(i); return a === "" ? "-" : a.toUpperCase(); }).join("/")).join(", ");
            const note = rev && !several ? "" : "// " + name_seq + where + ": after C->T conversion the alleles " + same +
                " are indistinguishable on the converted top strand (a non-CpG C/T variant); no allele-specific oligonucleotide can be designed on this strand.\n\n";
            return { resultarea1: note, resultarea2: note };
        }
        tpl = st.map(a => BisulfitePassAllele(s1, a, s2, rev));
    }
    // The reference template; each allele's own is s1 + allele + s2.
    const seq = s1 + ref + s2;
    const lseqs = seq.length;

    let sq = "";
    let cs = "";
    const end5 = "n";
    const end3 = ends3;

    let mplist = []; //probe list

    // Repeats and low complexity are looked for in each allele's own template, since the
    // alleles differ in sequence, and a common primer must avoid them in every one. Both
    // passes take them from the top strand as entered: the low-complexity search is not
    // strand-symmetric, and one record should have one mask.
    const msk = st.map(a => {
        const t = rev ? ComplementDNA(s1 + a + s2) : s1 + a + s2;
        let m = new Array(t.length).fill(0);
        if (maskrepeats) {
            m = RepeatMask(t).map(v => (v > 0 ? 1 : 0));
        }
        if (lowcomplexity) {
            const lowsequence = new LowComplexitySequence(t, Kmax, kmer2);
            const bloks = lowsequence.getBlocks();
            const dataL = bloks[0];
            for (let i = 0; i < dataL.length; i += 2) {
                for (let k = dataL[i]; k < dataL[i] + dataL[i + 1]; k++) {
                    m[k] |= 1;
                }
            }
        }
        if (rev) { m.reverse(); }
        for (let i = 0; i < l1; i++) { m[i] |= m1[i]; }
        for (let i = 0; i < l3; i++) { m[l1 + a.length + i] |= m3[i]; }
        return m;
    });

    if (ctconvert) {
        // DNActBisulfiteForward() gives the converted bottom strand, which no longer
        // pairs with the top one. Both passes read the input's converted top strand:
        // cs in the common primer's orientation, sq as its reverse complement.
        cs = rev ? BisulfiteTop(ComplementDNA(seq)) : DNActBisulfiteReverse(seq);
        sq = ComplementDNA(cs);
    }
    else {
        sq = seq;
        cs = ComplementDNA(seq);
    }

    //Forward direction - Probe design
    // In sizing mode the tails are computed from the primers themselves, so the
    // fluorophore tails from the tails tab are left out of this pass.
    const ends = rev ? ["end", "start"] : ["start", "end"];
    let result = ProbeDesignForward(snp, sizing.on ? [] : tails, s1, s2, st, minlen, maxlen, mintm, maxtm, "F_", tpl, msk, ends);
    if (result.fail) { return refuse(result.fail); }
    // The common primer binds wholly downstream of every allele's primer: for a short
    // allele that primer's 3' end reaches further into s2.
    const reach = Math.max(1, ...st.map(a => result.x1 - l1 - a.length));
    let x2 = l3 - minlen - reach + 1;
    if (x2 < 1) {
        return refuse("only " + Math.max(0, l3 - reach) + " nt of sequence lie beyond the allele-specific primers, and a common primer needs at least " + minlen);
    }
    let sizingNotes = sizing.on ? AlleleSizing(result, s1.length, st, snp, sizing) : null;
    if (result.fail) { return refuse(result.fail); }
    let fpr = result.fpr;
    let fpn = result.fpn;
    let ftm = result.ftm;
    let flc = result.flc;
    let fcg = result.fcg;
    let fpx = result.fpx;
    let fn = fpr.length;
    // A primer whose Tm could not be brought into the window is shown with the reason.
    let warn = "";
    for (let d = 0; d < fn; d++) {
        const w = result.why[d];
        if (!w || (ftm[d] > mintm && ftm[d] <= maxtm)) { continue; }
        warn += "// " + fpn[d] + ": Tm " + ftm[d].toFixed(1) + " °C is " + (ftm[d] > maxtm
            ? "above the " + maxtm + " °C ceiling; with its 3' end fixed at the variant, no length from " + minlen + " to " + w.last + " nt falls within " + mintm + "-" + maxtm + " °C"
            : "below the " + mintm + " °C floor; with its 3' end fixed at the variant, it cannot grow past " +
            (w.stop === "" ? "the maximum length of " + maxlen + " nt" : w.stop === "start" ? "the " + ends[0] + " of the sequence" : MaskWhat(w.stop))) + ".\n";
    }
    // The fluorescent readout sees an allele only through its primer's tail.
    const untailed = sizing.on || tails.length === 0 ? [] : fpn.slice(tails.length);
    if (untailed.length > 0) {
        warn += "// " + untailed.join(", ") + ": only " + tails.length + " tails are given for " + fn + " alleles, so " +
            (untailed.length === 1 ? "this primer carries none and a fluorescent readout cannot report its allele" : "these primers carry none and a fluorescent readout cannot report their alleles") + ".\n";
    }
    resultarea1 += which;
    for (let d = 0; d < fn; d++) {
        resultarea1 += fpn[d] + "\t" + fpr[d] + "\t" + fpr[d].length + "\t" + ftm[d].toFixed(1) + "\t" + fcg[d].toFixed(1) + "\t" + flc[d] + "\n";
        mplist.push(fpr[d].toLocaleLowerCase());
    }
    resultarea1 += warn;

    // reverse direction    
    // The common primer reads the reference template; a base masked in any allele's
    // template is masked for it.
    let mskR = new Array(lseqs).fill(1);
    for (let i = 0; i < l3; i++) {
        let v = 0;
        for (let j = 0; j < st.length; j++) { v |= msk[j][l1 + st[j].length + i]; }
        mskR[l1 + ref.length + i] = v;
    }
    mskR.reverse();
    let rpr = [];
    let rpn = [];
    let rtm = [];
    let rlc = [];
    let rpx = [];
    let rcg = [];
    let x1 = 0;
    let rn = PrimerDesignReverse(cs, mskR, fpr, mplist, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, "R_", end5, end3);
    for (let d = rn; d > -1; d--) {
        rpx[d] = lseqs - (rpx[d] + rpr[d].length) + 1;
        resultarea1 += rpn[d] + "\t" + rpr[d] + "\t" + rpr[d].length + "\t" + rtm[d].toFixed(1) + "\t" + rcg[d].toFixed(1) + "\t" + rlc[d] + "\n";
    }
    resultarea1 += "\n\n";
    // combinations
    if (rn > -1) {
        resultarea2 += "[KASP " + direction + " sets] " + name_seq + ":\n" + which;

        for (let f = 0; f < fn; f++) {
            resultarea2 += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
        }
        if (sizingNotes) { for (const t of sizingNotes) { resultarea2 += "// " + t + "\n"; } }
        resultarea2 += warn;
        resultarea2 += "\n";

        // Each allele has a product of its own, from its primer's 5' end to the common
        // primer's on that allele's template. rpx is on the reference template, so a
        // position past the variant moves by the allele's difference in length from ref.
        const tailLen = (j) => (result.tail ? result.tail[j].length : 0);
        const bandOf = (r, j) => rpx[r] + rpr[r].length - 1 - (ref.length - st[j].length) - fpx[j] + tailLen(j);
        // Other markers inside the product lengthen or shorten it with their alleles.
        const spread = (r) => {
            const end = rpx[r] + rpr[r].length - 1 - l1 - ref.length;
            let lo = 0;
            let hi = 0;
            for (const m of rst.mk) { if (m.at < end) { lo += m.lo - m.len; hi += m.hi - m.len; } }
            return [lo, hi];
        };
        const sameLength = st.every(a => a.length === ref.length);

        for (let r = rn; r > -1; r--) {
            let x3r = rpx[r];
            // Per-allele band lengths matter only when the gel is the readout; the
            // fluorescent assay gets one fragment size, that of the reference template, or
            // the range when the alleles differ in length.
            let p;
            if (sizing.on) {
                p = st.map((a, j) => (a === "" ? "del" : a.toUpperCase()) + ":" + bandOf(r, j)).join(" ");
            }
            else {
                const sp = spread(r);
                const all = st.map((a, j) => bandOf(r, j));
                const lo = Math.min(...all) + sp[0];
                const hi = Math.max(...all) + sp[1];
                p = sameLength && sp[0] === 0 && sp[1] === 0 ? String(bandOf(r, iref)) : (lo === hi ? String(lo) : lo + "-" + hi);
            }
            let ta = Tm77(sq.substring(fpx[iref], x3r + rpr[r].length - 1), KMg_M, 0);
            resultarea2 += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n";
        }
        resultarea2 += "\n";
    }
    return { resultarea1, resultarea2 };
}


// ── Allele sizing (gel readout) ──────────────────────────────────────────────
// A KASP assay tells the alleles apart by which fluorophore lights up.  Without a
// reader they have to differ in something a gel can see, and one such difference is
// already there: the allele-specific primers share a 3'-end at the variant and share
// one reverse primer, so each product is exactly as long as its own primer reaches.
// A 5'-tail lengthens that allele's band and nothing else, because the tail has no
// counterpart in the template and takes no part in discrimination.
//
// The shortest product is the baseline and keeps a bare primer; the others are padded
// onto the rungs of a ladder with a fixed step.  Where the alleles already differ by
// more than the step — InDels and haplotypes usually do — the tail stays empty.

// Filler for the tails: AT-leaning so it adds little to the melting temperature, free
// of runs that would slip during synthesis, and screened before use.
const sizingDonors = [
    "acggactattatgccacggcttgttatcaagcgt",
    "atctcttcatccacacgtactaataaggtacaag",
    "acgacgtcatacaagtgaagcgtatataaggcgt",
    "acgctagactaggacgttgatataagtgtgccta",
    "tggaatgtaggacggcaattataacgccggttat",
    "tgaatggccaccatagcaatactgcattacctca"
];

// The fillers for this particular set of primers, quietest first.  Taking the first one
// that happens to pass would fall back to an arbitrary sequence when none is perfectly
// clear; scoring them all means the assay gets the best available either way, and the
// next ones are there if the tails cut from the best still make the primers pair.
function SizingDonor(avoid) {
    const min3 = 6;
    const all = [];
    for (const d of sizingDonors) {
        let s = DimerLook2(d, d, min3) * 4;
        for (const o of avoid) {
            if (!o) { continue; }
            const low = o.toLowerCase();
            s += DimerLook2(d, low, min3) + DimerLook2(low, d, min3);
        }
        all.push({ seq: d, score: s });
    }
    // A stable sort: among equal scores the listed order decides, as before.
    return all.sort((a, b) => a.score - b.score);
}

// spans[j] = template bases the allele-j product carries before any tail.  Returns one
// tail per allele plus the longest tail that could not be cut from the donor.
function SizingLadder(spans, step, donor) {
    const n = spans.length;
    const rank = spans.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const tails = new Array(n).fill("");
    let over = 0;
    let prev = 0;
    for (let k = 0; k < n; k++) {
        const span = rank[k][0];
        const j = rank[k][1];
        const target = (k === 0) ? span : Math.max(prev + step, span);
        const need = target - span;
        if (need > donor.length) { over = Math.max(over, need); }
        tails[j] = need > 0 ? donor.substring(donor.length - Math.min(need, donor.length)) : "";
        prev = target;
    }
    return { tails: tails, over: over };
}

// One mismatch at the 3'-end leaks, and on a gel a leak looks like a real band; a
// second, deliberate mismatch just behind the variant stalls the polymerase on the
// wrong allele.  It goes in only where every allele primer carries the same base, so
// that it cannot start discriminating on its own.
function SizingMismatch(fpr, offset) {
    if (fpr.length < 2) { return null; }
    let from = "";
    for (const p of fpr) {
        const i = p.length - offset;
        if (i < 1) { return null; }
        const c = p.charAt(i).toLowerCase();
        if (from === "") { from = c; }
        else if (from !== c) { return null; }
    }
    const to = AntiSenseDNA(from);
    if (!to || to === from) { return null; }
    for (let j = 0; j < fpr.length; j++) {
        const i = fpr[j].length - offset;
        fpr[j] = fpr[j].substring(0, i) + to + fpr[j].substring(i + 1);
    }
    return { at: offset, from: from.toUpperCase(), to: to.toUpperCase() };
}

// Turns a finished set of allele-specific primers into a length-coded one, in place.
function AlleleSizing(result, s1len, alleles, snp, opt) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const fpr = result.fpr;
    const n = fpr.length;
    const notes = [];

    if (opt.mismatch) {
        // The mismatch sits immediately behind the discriminating window, so it only
        // lands on the classic second position when the variant is on the last base.
        if (snp > 1) {
            notes.push("The variant is " + snp + " bases from the 3'-end, so the second mismatch went to position " + (snp + 1) + ". Allele-specific PCR discriminates best with the variant on the last base — set the ASP offset to 1.");
        }
        const mm = SizingMismatch(fpr, snp + 1);
        notes.push(mm
            ? "Second mismatch " + mm.at + " nt from the 3'-end of every allele primer: " + mm.from + "→" + mm.to + "."
            : "Second mismatch not placed: the base behind the variant is not the same in every allele primer.");
    }

    // Tm, CG and complexity describe the part that anneals in the first cycle, so they
    // are taken from the body — after the mismatch, before the tail.
    for (let j = 0; j < n; j++) {
        const body = fpr[j].toLowerCase();
        result.ftm[j] = DNA_Tm(body, salt, Mg_M, p_mkM);
        result.fcg[j] = CG(body);
        result.flc[j] = LingComplexity2(body);
    }

    const spans = [];
    for (let j = 0; j < n; j++) { spans.push(s1len + alleles[j].length - result.fpx[j]); }

    // The mismatch changes the bodies, so they are screened again before any tail goes on.
    const d0 = DimerPair(fpr);
    if (d0) {
        result.fail = "with the second mismatch in, " + DimerWhy(result.fpn, d0) + ", so no set is reported";
        return notes;
    }

    // Tails are new sequence in the tube, so the tailed primers are screened against each
    // other exactly as the untailed ones were; primers that pair are not reported.
    let donor = null;
    let lad = null;
    for (const d of SizingDonor(fpr)) {
        const l = SizingLadder(spans, opt.step, d.seq);
        if (!DimerPair(fpr.map((b, j) => l.tails[j] + b))) { donor = d; lad = l; break; }
    }
    if (!lad) {
        result.fail = "with every filler tried, the tailed allele-sizing primers pair with each other, so no set is reported";
        return notes;
    }
    for (let j = 0; j < n; j++) { fpr[j] = lad.tails[j] + fpr[j]; }
    result.tail = lad.tails;

    const pad = [];
    let longest = 0;
    for (let j = 0; j < n; j++) {
        pad.push((alleles[j] === "" ? "del" : alleles[j].toUpperCase()) + (lad.tails[j].length ? " +" + lad.tails[j].length : " bare"));
        if (fpr[j].length > longest) { longest = fpr[j].length; }
    }
    notes.unshift("Bands step by at least " + opt.step + " bp; 5'-tails (nt): " + pad.join(", ") + ".");
    if (lad.over > 0) {
        notes.push("A tail of " + lad.over + " nt was needed but the filler is only " + donor.seq.length + " nt — lower the step.");
    }
    if (longest > opt.maxlen) {
        notes.push("Longest primer is " + longest + " nt, over the " + opt.maxlen + " nt limit — lower the step or shorten the primer length range.");
    }
    return notes;
}
