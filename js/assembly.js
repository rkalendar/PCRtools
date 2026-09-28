'use strict';
function analysis() {
    // A design field that is empty or not a number stops the run, named, rather than being
    // replaced by a default the user never sees.
    let badField = "";
    const readInt = (id, def, label) => {
        const raw = String(document.getElementById(id).value ?? "").trim();
        const v = parseInt(raw, 10);
        if (raw === "" || !Number.isFinite(Number(raw)) || !Number.isFinite(v)) {
            if (!badField) { badField = label; }
            return def;
        }
        return v;
    };

    let minlc  = readInt('minlc', 70, "Min. linguistic complexity");
    let mintm  = readInt('mintm', 60, "Min. Tm");
    let maxtm  = readInt('maxtm', 62, "Max. Tm");
    let minlen = readInt('minlen', 18, "Min. length");
    if (badField) {
        const msg = "Run stopped: the field \u201C" + badField + "\u201D is empty or not a number.\n";
        return [msg, msg];
    }

    let ends3 = document.getElementById('end3com').value.toLowerCase().trim();   // e.g. "swh ssw wsh sww www"
    if (DNA(ends3).length < 1) { ends3 = "n"; }

    if (minlc < 10) { minlc = 10; }
    if (minlc > 90) { minlc = 90; }
    if (minlen < 12) { minlen = 12; }
    if (minlen > 100) { minlen = 100; }
    if (mintm < 37) { mintm = 37; }
    if (mintm > 80) { mintm = 80; }
    // maxtm is an enforced ceiling: a primer's Tm must never exceed it. There is no maximum
    // length — the primer is extended as far as needed to reach mintm — but if reaching mintm
    // would push Tm above maxtm, the primer is shortened instead (even below minlen) so it stays
    // under maxtm. minlen is therefore a soft preference, not a hard floor.
    if (maxtm < mintm) { maxtm = mintm; }
    if (maxtm > 90)    { maxtm = 90; }

    // My Primers lists oligos, one per line or named.
    const plistRead = ReadingOligos(document.getElementById('inputPrimerList').value);
    let plist = plistRead.seqs;
    const n_plist = plist.length;
    for (let n = 0; n < n_plist; n++) {
        plist[n] = DNA(plist[n]).toLowerCase();
    }

    const ReadResult = ReadingSeq(document.getElementById('inputText').value);
    const name_seq = ReadResult.name_seq;
    const seqs = ReadResult.seqs;
    const n_seq = seqs.length;

    // A box that cannot be read, or no sequence at all, stops the run with the reason, also
    // on the pages that run without the input check (primerdigital).
    const stop = [ReadingStop(ReadResult), ReadingStop(plistRead, "My Primers")].filter(Boolean).join("\n");
    if (stop) {
        return [stop + "\n", stop + "\n"];
    }

    // Gibson Assembly joins fragments in order and needs at least 3
    // (vector-left, insert(s)…, vector-right). With fewer, the overlap step below would
    // read out-of-range entries (e.g. fpr[-1]); return a clear message instead of throwing.
    if (n_seq < 3) {
        const msg = "Gibson Assembly needs at least 3 fragments in assembly order " +
                    "(e.g. vector-left, insert(s)\u2026, vector-right). Found " + n_seq + ".\n" +
                    "Add the missing fragments and run again.\n";
        return [msg, msg];
    }

    const min3 = 8;
    const result = ["", ""];
    let resultarea1 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(LC%)\tLinguistic_Complexity(YR%)\n";
    let resultarea2 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(LC%)\tLinguistic_Complexity(YR%)\n";

    let fpr = new Array(n_seq).fill("");
    let fpn = new Array(n_seq).fill("");
    let ftm = new Array(n_seq).fill(0.0);
    let flc = new Array(n_seq).fill(0);
    let flci = new Array(n_seq).fill(0);
    let fcg = new Array(n_seq).fill(0.0);

    let rpr = new Array(n_seq).fill("");
    let rpn = new Array(n_seq).fill("");
    let rtm = new Array(n_seq).fill(0.0);
    let rlc = new Array(n_seq).fill(0);
    let rlci = new Array(n_seq).fill(0);
    let rcg = new Array(n_seq).fill(0.0);


    let fn = 0;
    let rn = 0;

    // Auto-relax the linguistic-complexity floor when no good primer is found.
    // Strategy: first hunt for an in-window primer (mintm <= Tm <= maxtm), stepping minlc down
    // to MINLC_FLOOR; only if none is achievable at any complexity do we accept the best
    // primer that stays under maxtm (which may be below mintm — the "shorten to respect maxtm"
    // case). Returns the minlc that succeeded, or -1 if even the floor produced nothing.
    const MINLC_FLOOR = 10;
    const MINLC_STEP = 10;
    // tail/pl: a junction oligo is designed with its 5' overlap tail in place and with the partner
    // oligo added to the primer list, so the self- and cross-dimer filters judge what is ordered.
    const designRelax = (n, sq, name, nm, lseqs, prArr, pnArr, tmArr, lcArr, lciArr, cgArr, tail = "", pl = plist) => {
        // 1) Prefer an in-window primer (mintm <= Tm <= maxtm), keeping complexity as high as
        //    possible: try the user's minlc first, then step down to the floor.
        for (let lc = minlc; lc >= MINLC_FLOOR; lc -= MINLC_STEP) {
            if (PrimerDesign(n, sq, name, pl, min3, lc, minlen, mintm, maxtm,
                             prArr, pnArr, tmArr, lcArr, lciArr, cgArr, nm, ends3, tail, lseqs, true)) {
                return { ok: true, diag: null };
            }
        }
        // 2) No in-window primer is reachable (Tm window unreachable). Take the best compromise
        //    UNDER maxtm at the loosest complexity, so its Tm is as close to mintm as possible.
        if (PrimerDesign(n, sq, name, pl, min3, MINLC_FLOOR, minlen, mintm, maxtm,
                         prArr, pnArr, tmArr, lcArr, lciArr, cgArr, nm, ends3, tail, lseqs, false)) {
            return { ok: true, diag: null };
        }
        // 3) Nothing usable. Run once more to capture WHY (e.g. self-dimer on in-window lengths).
        const diag = {};
        PrimerDesign(n, sq, name, pl, min3, MINLC_FLOOR, minlen, mintm, maxtm,
                     prArr, pnArr, tmArr, lcArr, lciArr, cgArr, nm, ends3, tail, lseqs, false, diag);
        return { ok: false, diag };
    };

    const notes = [];
    // Why a fragment has no forward or no reverse primer, for Report 2 as well as the notes.
    const fwhy = new Array(n_seq).fill(""), rwhy = new Array(n_seq).fill("");
    const failReason = (d) =>
        (d && d.selfDimerInWindow)
            ? "in-window primers exist but self-dimerise (often a palindromic restriction site at this end) \u2014 prime this end manually or move the fragment boundary"
            : "no primer with Tm \u2264 " + maxtm + "\u00B0C; try lowering minlen/min. complexity or widening the Tm range";

    for (let n = 0; n < n_seq; n++) {
        let RestSeq = Seq(seqs[n]);
        const seq = RestSeq.returnString;
        const cs = ComplementDNA(seq);

        const fRes = designRelax(n, seq, name_seq[n], "F_", 0, fpr, fpn, ftm, flc, flci, fcg);
        if (fRes.ok) {
            resultarea1 += fpn[n] + "\t" + fpr[n] + "\t" + fpr[n].length + "\t" + ftm[n].toFixed(1) + "\t" + fcg[n].toFixed(1) + "\t" + flc[n] + "\t" + flci[n] + "\n";
            if (flc[n] < minlc) { notes.push("Forward " + name_seq[n] + ": complexity " + flc[n] + "% is below the requested " + minlc + "% (relaxed to find a primer)."); }
            if (ftm[n] < mintm) { notes.push("Forward " + name_seq[n] + ": Tm " + ftm[n].toFixed(1) + "\u00B0C is below mintm (" + mintm + ") \u2014 shortened to keep Tm under maxtm (" + maxtm + ")."); }
        } else {
            fwhy[n] = failReason(fRes.diag);
            notes.push("Forward " + name_seq[n] + ": NO primer \u2014 " + fwhy[n] + ".");
        }

        // Reverse names count from the end of the template (the bases Seq() read), not of the
        // record text with its spaces and line breaks.
        const rRes = designRelax(n, cs, name_seq[n], "R_", seq.length, rpr, rpn, rtm, rlc, rlci, rcg);
        if (rRes.ok) {
            resultarea1 += rpn[n] + "\t" + rpr[n] + "\t" + rpr[n].length + "\t" + rtm[n].toFixed(1) + "\t" + rcg[n].toFixed(1) + "\t" + rlc[n] + "\t" + rlci[n] + "\n";
            if (rlc[n] < minlc) { notes.push("Reverse " + name_seq[n] + ": complexity " + rlc[n] + "% is below the requested " + minlc + "% (relaxed to find a primer)."); }
            if (rtm[n] < mintm) { notes.push("Reverse " + name_seq[n] + ": Tm " + rtm[n].toFixed(1) + "\u00B0C is below mintm (" + mintm + ") \u2014 shortened to keep Tm under maxtm (" + maxtm + ")."); }
        } else {
            rwhy[n] = failReason(rRes.diag);
            notes.push("Reverse " + name_seq[n] + ": NO primer \u2014 " + rwhy[n] + ".");
        }
    }

    if (notes.length) {
        resultarea1 += "\n# Notes\n# " + notes.join("\n# ") + "\n";
    }


    // combinations
    // One junction: the forward oligo of fragment n carries the end of fragment n-1 as a 5' tail
    // and its reverse oligo the start of fragment n+1, so these composed oligos are what is
    // ordered and what the dimer test must judge. A pair that fails it is not printed with its
    // complexity overwritten by 0 (S5): the annealing parts are designed again with the tails in
    // place and the partner oligo excluded, and if no dimer-free pair exists the junction is
    // reported as having none.
    const pairDimer = (f1, r1) => DimerLook2(f1, f1, min3) > 0 || DimerLook2(r1, r1, min3) > 0 ||
        DimerLook2(f1, r1, min3) > 0;
    // Design one oligo of a junction again, with its overlap tail and with `avoid` treated as an
    // existing primer, so the tail's dimers are caught by the same filters. Tm, CG and complexity
    // describe the annealing part, as in the per-fragment rows; null when nothing passes. The new
    // one must reach Min. Tm, or the Tm of the oligo it replaces when that one is already below
    // (floor): designRelax's last resort, the longest primer under Max. Tm, can lie far below the
    // range and would not anneal beside its partner.
    const redesign = (n, sq, nm, lseqs, tail, avoid, floor) => {
        const pr = [""], pn = [""], tm = [0], lc = [0], lci = [0], cg = [0];
        if (!designRelax(0, sq, name_seq[n], nm, lseqs, pr, pn, tm, lc, lci, cg, tail, plist.concat([avoid])).ok ||
            tm[0] < floor) {
            return null;
        }
        const core = pr[0].substring(tail.length);
        return {
            seq: pr[0], len: core.length, tm: tm[0], cg: CG(core), lc: lc[0], lci: lci[0],
            name: lseqs > 0 ? name_seq[n] + "_" + nm + lseqs + "-" + (lseqs - core.length + 1)
                : name_seq[n] + "_" + nm + "1-" + core.length,
        };
    };
    const row = o => o.name + "\t" + o.seq + "\t" + o.len + "\t" + o.tm.toFixed(1) + "\t" + o.cg.toFixed(1) + "\t" + o.lc + "\t" + o.lci + "\n";
    // A redesigned oligo below Min. Tm (it replaces one that was below too) is noted, as the
    // per-fragment rows are in the Report tab.
    const lowTm = o => o.tm < mintm ? "# " + o.name + ": Tm " + o.tm.toFixed(1) + "\u00B0C is below mintm (" + mintm + "), as for the primer it replaces.\n" : "";

    // The two vector-arm primers meet when the construct is circularised: they are used in one
    // reaction, so they are judged as a junction pair is (S5). A pair that dimerises is designed
    // again, each primer with the other treated as an existing one, and if no dimer-free pair
    // exists the vector arms are reported as having none.
    if (rpr[0].length > 0 && fpr[n_seq - 1].length > 0) {
        let r = { seq: rpr[0], len: rpr[0].length, tm: rtm[0], cg: rcg[0], lc: rlc[0], lci: rlci[0], name: rpn[0] };
        let f = { seq: fpr[n_seq - 1], len: fpr[n_seq - 1].length, tm: ftm[n_seq - 1], cg: fcg[n_seq - 1], lc: flc[n_seq - 1], lci: flci[n_seq - 1], name: fpn[n_seq - 1] };
        if (pairDimer(r.seq, f.seq)) {
            const left = Seq(seqs[0]).returnString, right = Seq(seqs[n_seq - 1]).returnString;
            const nr = redesign(0, ComplementDNA(left), "R_", left.length, "", f.seq, Math.min(mintm, r.tm));
            if (nr) { r = nr; }
            const nf = redesign(n_seq - 1, right, "F_", 0, "", r.seq, Math.min(mintm, f.tm));
            if (nf) { f = nf; }
        }
        if (pairDimer(r.seq, f.seq)) {
            resultarea2 += "No primer pair for the vector arms (" + name_seq[0] + " / " + name_seq[n_seq - 1] +
                "): the reverse primer of " + name_seq[0] + " and the forward primer of " + name_seq[n_seq - 1] +
                " dimerise with each other, and no other annealing part within the Tm range avoids it. They are used in the same " +
                "reaction when the construct is circularised, so move a vector boundary or prime the vector manually.\n";
        } else {
            if (r.name !== rpn[0] || f.name !== fpn[n_seq - 1]) {
                resultarea2 += "# Vector arms (" + name_seq[0] + " / " + name_seq[n_seq - 1] + "): designed again to avoid a primer-dimer, so they differ from the Report (Primer list) tab.\n";
                if (r.name !== rpn[0]) { resultarea2 += lowTm(r); }
                if (f.name !== fpn[n_seq - 1]) { resultarea2 += lowTm(f); }
            }
            resultarea2 += row(r);
            resultarea2 += row(f);
        }
    } else {
        if (rpr[0].length > 0) {
            resultarea2 += rpn[0] + "\t" + rpr[0] + "\t" + rpr[0].length + "\t" + rtm[0].toFixed(1) + "\t" + rcg[0].toFixed(1) + "\t" + rlc[0] + "\t" + rlci[0] + "\n";
        } else {
            resultarea2 += "# Reverse primer of " + name_seq[0] + ": none found (" + rwhy[0] + ").\n";
        }
        if (fpr[n_seq - 1].length > 0) {
            resultarea2 += fpn[n_seq - 1] + "\t" + fpr[n_seq - 1] + "\t" + fpr[n_seq - 1].length + "\t" + ftm[n_seq - 1].toFixed(1) + "\t" + fcg[n_seq - 1].toFixed(1) + "\t" + flc[n_seq - 1] + "\t" + flci[n_seq - 1] + "\n";
        } else {
            resultarea2 += "# Forward primer of " + name_seq[n_seq - 1] + ": none found (" + fwhy[n_seq - 1] + ").\n";
        }
    }

    for (let n = 1; n < n_seq - 1; n++) {
        if (fpr[n].length > 0 && rpr[n].length > 0 && rpr[n - 1].length > 0 && fpr[n + 1].length > 0) {
            const ftail = ComplementDNA(rpr[n - 1]).toUpperCase();
            const rtail = ComplementDNA(fpr[n + 1]).toUpperCase();
            let f = { seq: ftail + fpr[n], len: fpr[n].length, tm: ftm[n], cg: fcg[n], lc: flc[n], lci: flci[n], name: fpn[n] };
            let r = { seq: rtail + rpr[n], len: rpr[n].length, tm: rtm[n], cg: rcg[n], lc: rlc[n], lci: rlci[n], name: rpn[n] };

            if (pairDimer(f.seq, r.seq)) {
                const sq = Seq(seqs[n]).returnString;
                const nf = redesign(n, sq, "F_", 0, ftail, r.seq, Math.min(mintm, f.tm));
                if (nf) { f = nf; }
                const nr = redesign(n, ComplementDNA(sq), "R_", sq.length, rtail, f.seq, Math.min(mintm, r.tm));
                if (nr) { r = nr; }
            }

            if (pairDimer(f.seq, r.seq)) {
                // Nothing to order for this junction; say so instead of printing the pair.
                const why = DimerLook2(f.seq, r.seq, min3) > 0
                    ? "the forward and the reverse oligo dimerise with each other"
                    : DimerLook2(f.seq, f.seq, min3) > 0
                        ? "the forward oligo self-dimerises" : "the reverse oligo self-dimerises";
                resultarea2 += "\nNo primer pair for junction " + n + " (" + name_seq[n] + " / " + name_seq[n - 1] +
                    " overlap): " + why + ", and no other annealing part within the Tm range avoids it. The overlap tails are the " +
                    "neighbouring fragments' primers, so move the fragment boundary or prime this junction manually.\n";
                continue;
            }
            if (f.name !== fpn[n] || r.name !== rpn[n]) {
                resultarea2 += "# Junction " + n + " (" + name_seq[n] + "): the annealing part was designed again with the overlap tail in place to avoid a primer-dimer, so it differs from the Report (Primer list) tab.\n";
                if (f.name !== fpn[n]) { resultarea2 += lowTm(f); }
                if (r.name !== rpn[n]) { resultarea2 += lowTm(r); }
            }
            resultarea2 += row(f);
            resultarea2 += row(r);
        }
        else {
            // The neighbours' primers are this junction's overlap tails, so any one missing leaves
            // it without oligos: name each missing primer and why it is missing.
            const missing = [];
            if (!fpr[n].length) { missing.push(name_seq[n] + " has no forward primer (" + fwhy[n] + ")"); }
            if (!rpr[n].length) { missing.push(name_seq[n] + " has no reverse primer (" + rwhy[n] + ")"); }
            if (!rpr[n - 1].length) { missing.push(name_seq[n - 1] + " has no reverse primer for the forward oligo's overlap tail (" + rwhy[n - 1] + ")"); }
            if (!fpr[n + 1].length) { missing.push(name_seq[n + 1] + " has no forward primer for the reverse oligo's overlap tail (" + fwhy[n + 1] + ")"); }
            resultarea2 += "\nNo primer pair for junction " + n + " (" + name_seq[n] + " / " + name_seq[n - 1] +
                " overlap): " + missing.join("; ") + ".\n";
        }
    }

    // ---- Amplification primers for the assembled construct (vector-anchored pair) ----
    // Forward primer placed ANYWHERE on the left vector arm (fragment 0) and reverse primer
    // ANYWHERE on the right vector arm (last fragment), both pointing inward toward the joined
    // inserts. Together they flank the assembled cassette so it can be PCR-amplified / verified
    // after joining. Unlike the per-fragment primers, the position is free (it slides to find a
    // good primer), so it is not forced onto the palindromic cloning site at the very end.
    {
        const leftVec  = Seq(seqs[0]).returnString;
        const rightVec = Seq(seqs[n_seq - 1]).returnString;
        const fp = VectorPrimerScan(leftVec, plist, min3, minlc, minlen, mintm, maxtm, ends3);
        const rp = VectorPrimerScan(ComplementDNA(rightVec), plist, min3, minlc, minlen, mintm, maxtm, ends3);

        let mid = 0;
        for (let n = 1; n < n_seq - 1; n++) mid += Seq(seqs[n]).returnString.length;

        resultarea2 += "\n# Amplification of the assembled construct \u2014 vector-anchored pair (flanks the joined inserts)\n";
        resultarea2 += "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(\u00B0C)\tCG(%)\tLinguistic_Complexity(LC%)\tLinguistic_Complexity(YR%)\n";

        if (fp.ok) {
            const id = name_seq[0] + "_ampF_" + (fp.start + 1) + "-" + fp.end;
            resultarea2 += id + "\t" + fp.seq + "\t" + fp.len + "\t" + fp.tm.toFixed(1) + "\t" + fp.cg.toFixed(1) + "\t" + fp.lc + "\t" + fp.lci + "\n";
            if (fp.len < minlen) { resultarea2 += "# (forward amplification primer is " + fp.len + " nt, below minlen " + minlen + ", to avoid cloning-site self-dimers)\n"; }
        } else {
            resultarea2 += "# Forward amplification primer: none found on " + name_seq[0] + " (" + failReason(fp.diag) + ").\n";
        }

        if (rp.ok) {
            const Lr = rightVec.length;
            const a = Lr - rp.end + 1;   // 1-based start on the right vector arm
            const b = Lr - rp.start;     // 1-based end on the right vector arm
            const id = name_seq[n_seq - 1] + "_ampR_" + a + "-" + b;
            resultarea2 += id + "\t" + rp.seq + "\t" + rp.len + "\t" + rp.tm.toFixed(1) + "\t" + rp.cg.toFixed(1) + "\t" + rp.lc + "\t" + rp.lci + "\n";
            if (rp.len < minlen) { resultarea2 += "# (reverse amplification primer is " + rp.len + " nt, below minlen " + minlen + ", to avoid cloning-site self-dimers)\n"; }
        } else {
            resultarea2 += "# Reverse amplification primer: none found on " + name_seq[n_seq - 1] + " (" + failReason(rp.diag) + ").\n";
        }

        if (fp.ok && rp.ok) {
            const leftReach  = leftVec.length - fp.start;
            const rightReach = rightVec.length - rp.start;   // rp.start is on the reverse-complement strand
            const amplicon = leftReach + mid + rightReach;
            resultarea2 += "# Approx. amplicon \u2248 " + amplicon + " bp (left arm " + leftReach + " + inserts " + mid + " + right arm " + rightReach + ").\n";
            if (DimerLook2(fp.seq, rp.seq, min3) > 0) {
                resultarea2 += "# Warning: the two amplification primers may form a primer-dimer.\n";
            }
        }
    }

    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

// The page's sequence badge (displayseq.js) counts the bases this engine reads: u and i, but
// not the LNA letters e, f, j and l.
function engineBases(s) { return Seq(s).returnString; }

function Seq(str) {
    let returnString = "";
    let b1 = [];
    let b = [];
    let bx = [];
    for (let i = 0; i < str.length; i++) {
        let chr = str.charAt(i);
        if (chr === 'a' || chr === 't' || chr === 'c' || chr === 'g' || chr === 'r' || chr === 'y' || chr === 'm' || chr === 'k' || chr === 'w' || chr === 'b' || chr === 'd' || chr === 'v' || chr === 'h' || chr === 's' || chr === 'n') {
            returnString += chr;
        }
        if (chr === 'u') {
            returnString += 't';
        }
        if (chr === 'i') {
            returnString += 'g';
        }
        if (chr === '[') {
            b1.push(returnString.length);
        }
        if (chr === ']') {
            b1.push(-returnString.length);
        }
        if (chr === '/') {
            bx.push(returnString.length);
        }
    }
    if (b1.length === 0) {
        b[0] = 0;
        b[1] = returnString.length;
        b[2] = 0;
        b[3] = returnString.length;
    }
    else {
        for (let i = 0; i < b1.length - 1; i++) {
            if (b1[i] >= 0) {
                b.push(b1[i]);
                for (let j = i + 1; j < b1.length; j++) {
                    if (b1[j] < 0) {
                        b.push(-b1[j]);
                        break;
                    }
                }
            }
        }
        if (b.length === 2) { //  [SNP]
            b1[0] = b[0];
            b1[1] = b[1];
            b[0] = 0;
            b[1] = b1[0] - 1;
            b[2] = b1[1] + 1;
            b[3] = returnString.length;
        }
    }
    return { returnString, b, bx };
}

function PrimerDesign(n, sq, seqname, plist, min3, minlc, minlen, mintm, maxtm, fpr, fpn, ftm, flc, flci, fcg, nm, e3, tail, lseqs, inWindowOnly, diag) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;
    const nplist = plist.length;

    // Absolute shortest primer we will ever return. minlen is only a preference: the primer
    // may be shortened below it (down to this floor) so that Tm stays under maxtm.
    const absMin = Math.max(le3, 10);
    const xmax = sq.length;                 // no maximum length — extend as far as needed

    // Candidate stored only when it passes every non-Tm filter AND Tm <= maxtm.
    let best = null;                        // chosen primer
    let bestUnder = null;                   // best fallback with Tm < mintm (highest Tm = longest)
    let sawInWindow = false;                // some 3'-valid length had mintm <= Tm <= maxtm
    let selfDimerInWindow = false;          // such an in-window length was blocked by a self-dimer

    // Tm rises with length, so once Tm > maxtm no longer length can qualify -> stop.
    for (let x = absMin; x <= xmax; x++) {
        const core = sq.substring(0, x);
        const s2 = core.substring(x - le3);
        if (s3.indexOf(s2) === -1) continue;            // 3'-end pattern must match

        const tm1 = DNA_Tm(core, salt, Mg_M, p_mkM);
        if (tm1 > maxtm) break;                          // ceiling reached; nothing longer fits
        const inWin = (tm1 >= mintm);                    // Tm is inside [mintm, maxtm]
        if (inWin) sawInWindow = true;

        if (ssrrepeat(core) !== -1) continue;            // reject simple-sequence repeats
        const lc = LingComplexity2(core);
        const lci = LingComplexityRY(core);
        if (lc < minlc) continue;                        // complexity floor

        const cand = tail + core;
        if (DimerLook2(cand, cand, min3) !== 0) {        // self-dimer
            if (inWin) selfDimerInWindow = true;
            continue;
        }
        if (nplist > 0) {
            let clash = false;
            for (let j = 0; j < nplist; j++) {
                if (DimerLook2(cand, plist[j], min3) > 0) { clash = true; break; }
            }
            if (clash) continue;                          // dimer with an existing primer
        }

        // cand passes all filters and Tm <= maxtm.
        if (inWin) {
            // In-window candidate. Prefer length >= minlen and the shortest such (lowest Tm,
            // i.e. the safest margin below maxtm). While still below minlen, keep moving up
            // toward it; once at/above minlen, freeze the first one found.
            if (best === null || best.x < minlen) {
                best = { x, tm: tm1, lc, lci, cand };
            }
            if (best.x >= minlen) break;                  // shortest in-window >= minlen found
        } else {
            // Below the window but under the ceiling — remember the highest-Tm (longest) one.
            bestUnder = { x, tm: tm1, lc, lci, cand };
        }
    }

    if (diag) { diag.sawInWindow = sawInWindow; diag.selfDimerInWindow = selfDimerInWindow; }

    // In-window primer is always preferred. The under-ceiling fallback is used ONLY when the
    // Tm window was genuinely unreachable (Tm steps over [mintm,maxtm] between two lengths) —
    // i.e. the real "shorten to respect maxtm" case. If in-window lengths existed but were
    // blocked by a filter (e.g. a palindromic self-dimer), we do NOT emit a misleading short
    // primer; we report no primer so the reason can be surfaced.
    const pick = best || ((!inWindowOnly && !sawInWindow) ? bestUnder : null);
    if (!pick) return false;

    fpr[n] = pick.cand;
    ftm[n] = pick.tm;
    flc[n] = pick.lc;
    flci[n] = pick.lci;
    fcg[n] = CG(pick.cand);
    fpn[n] = (lseqs > 0)
        ? seqname + "_" + nm + lseqs + "-" + (lseqs - pick.cand.length + 1)
        : seqname + "_" + nm + "1-" + pick.cand.length;
    return true;
}

// Find a good FORWARD primer ANYWHERE on `strand`, preferring a 3'-end as close as possible to
// the strand's 3' end (used to flank the joined inserts). Position is free, so the primer can
// slide away from an awkward end (e.g. a palindromic cloning site). minlen is a preference: a
// shorter primer (down to an absolute floor) is accepted when that is what fits the Tm window or
// avoids a self-dimer. Same Tm/complexity/dimer rules as PrimerDesign; maxtm is a hard ceiling.
// Returns { ok, seq, tm, lc, cg, start, end, len, diag }  (start/end are 0-based on `strand`).
function VectorPrimerScan(strand, plist, min3, minlc, minlen, mintm, maxtm, e3) {
    const salt = 0.055, Mg_M = 0.001, p_mkM = 0.2;
    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;
    const nplist = plist.length;
    const L = strand.length;
    const hardMin = Math.max(le3, 10);

    let bestUnder = null;                       // highest-Tm primer with Tm < mintm (fallback)
    let sawInWindow = false, selfDimerInWindow = false;

    // One pass at a given minimum length. Scans the 3'-end position high -> low (junction first)
    // and, for each, grows the length until the Tm window is reached. Returns the junction-most,
    // shortest in-window primer, or null. Updates the shared sawInWindow / bestUnder state.
    const scan = (minLenLocal) => {
        const lo = Math.max(le3, minLenLocal);
        for (let e = L; e >= lo; e--) {
            const s2 = strand.substring(e - le3, e);
            if (s3.indexOf(s2) === -1) continue;        // 3'-end pattern (depends only on the 3' end)
            for (let len = minLenLocal; e - len >= 0; len++) {
                const core = strand.substring(e - len, e);
                const tm = DNA_Tm(core, salt, Mg_M, p_mkM);
                if (tm > maxtm) break;                   // hotter than ceiling; stop for this 3' end
                const inWin = (tm >= mintm);
                if (inWin) sawInWindow = true;
                if (ssrrepeat(core) !== -1) continue;
                const lc = LingComplexity2(core);
        const lci = LingComplexityRY(core);
                if (lc < minlc) continue;
                if (DimerLook2(core, core, min3) !== 0) { if (inWin) selfDimerInWindow = true; continue; }
                let clash = false;
                for (let j = 0; j < nplist; j++) { if (DimerLook2(core, plist[j], min3) > 0) { clash = true; break; } }
                if (clash) continue;
                if (inWin) {
                    return { ok: true, seq: core, tm, lc, lci, cg: CG(core), start: e - len, end: e, len };
                } else if (!bestUnder || tm > bestUnder.tm) {
                    bestUnder = { ok: true, seq: core, tm, lc, lci, cg: CG(core), start: e - len, end: e, len };
                }
            }
        }
        return null;
    };

    // Prefer a primer of length >= minlen; if none is in-window, allow a shorter one (down to
    // hardMin) — the vector cloning site often only yields short dimer-free primers.
    const hit = scan(minlen) || scan(hardMin);
    if (hit) { hit.diag = { sawInWindow, selfDimerInWindow }; return hit; }

    // No in-window primer anywhere. Use the under-window fallback only if the Tm window was
    // never reachable (the genuine narrow-window case), not when filters blocked the good lengths.
    if (bestUnder && !sawInWindow) {
        bestUnder.diag = { sawInWindow, selfDimerInWindow };
        return bestUnder;
    }
    return { ok: false, diag: { sawInWindow, selfDimerInWindow } };
}
