let dataStore = {};
function addRecord(category, record) {
    if (!dataStore[category]) {
        dataStore[category] = [];
    }
    dataStore[category].push(record);
}
function analysis() {
    dataStore = {};
    let minlc = parseInt(document.getElementById('minlc').value);
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);

    let multiplex = document.getElementById('multiplex').checked;

    const minpcr = parseInt(document.getElementById('minpcr').value);
    const maxpcr = parseInt(document.getElementById('maxpcr').value);
    const overlapping = document.getElementById('overlapping').checked;
    const maskrepeats = document.getElementById('repeats').checked;
    const ctconvert = document.getElementById('ctconvert').checked;
    const mgbprobe = document.getElementById('mgbprobe').checked;
    const taqmanprobe = document.getElementById('taqmanprobe').checked;
    const invertedPCR = document.getElementById('invertedPCR').checked; 
    const ftail = DNA(document.getElementById('ftail').value.toLowerCase().trim()); // "tttcacacaggaaacagctatgac";
    const rtail = DNA(document.getElementById('rtail').value.toLowerCase().trim()); // "tttcacacaggaaacagctatgac";

    const ends5 = "n";
    let ends3 = document.getElementById('end3com').value.toLowerCase().trim();   // "swh ssw wsh sww www";
    if (DNA(ends3).length < 1) { ends3 = "n"; }

    if (minlc < 10) { minlc = 10; }
    if (minlc > 90) { minlc = 90; }
    if (minlen < 12) { minlen = 12; }
    if (maxlen < 12) { maxlen = 12; }
    if (minlen > 500) { minlen = 500; }
    if (maxlen > 500) { maxlen = 500; }
    if (minlen > maxlen) { maxlen = minlen; }
    if (mintm < 37) { mintm = 37; }
    if (mintm > 80) { mintm = 80; }

    let primeroverlap = 12;
    if (overlapping) { primeroverlap = minlen; }

    let probetm = 0;
    let probedesign = false;
    let probename = "";
    let probeminlen = 12;
    let probemaxlen = 21;
    if (mgbprobe) { probetm = -2; probename = "MGB_"; probedesign = true; probeminlen = 12; probemaxlen = 21; }
    if (taqmanprobe) { probetm = 8; probename = "TaqMan_"; probedesign = true; probeminlen = 18; probemaxlen = 25; }

    let mplist = []; //multiplex list

    ReadResult = ReadingSeq(document.getElementById('inputPrimerList').value);
    let plist = ReadResult.seqs;
    const n_plist = plist.length;
    for (let n = 0; n < n_plist; n++) {
        plist[n] = DNA(plist[n]).toUpperCase();
    }

    ReadResult = ReadingSeq(document.getElementById('inputText').value);
    const name_seq = ReadResult.name_seq;
    const seqs = ReadResult.seqs;
    const n_seq = seqs.length;
    const result = ["", ""];

    if (n_seq === 0) {
        return result;
    }
    if (n_seq === 1) {
        multiplex = false;
    }
    let resultarea1 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n";
    let resultarea2 = "PrimerID\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\tFragment_Size(bp)/Tm(°C)\n\n";
    let startTime = performance.now();// Start time
    /*
        for (let i = 0; i < 10001; i++) {
            const result = qDimer("cgacgttgtaaaacgacggccagtctatggtcatgg", "agtaaagatggttgaggatgtaaggtcttgacccatgatgatacgtctcctatagtgagtcgtattaggatcc");
        }
    */

    for (let n = 0; n < n_seq; n++) {
        resultarea1 += name_seq[n] + ":\n";

        let RestSeq = Seq(seqs[n]);
        const b = RestSeq.b;
        const bx = RestSeq.bx;
        const seq = RestSeq.returnString;
        const lseqs = seq.length;
        let sq = "";
        let cs = "";

        let msk = new Array(lseqs).fill(0);
        if (maskrepeats) {
            msk = RepeatMask(seq);
        }
        if (bx.length > 1) {
            for (let i = 0; i < bx.length; i += 2) {
                // msk = msk.fill(1, bx[i], bx[i + 1]);
                for (let j = bx[i]; j < bx[i + 1]; j++) {
                    msk[j]++;
                }
            }
        }

        if (ctconvert) {
            sq = DNActBisulfiteForward(seq);
            cs = DNActBisulfiteReverse(seq);
        }
        else {
            sq = seq;
            cs = ComplementDNA(seq);
        }

        let fpr = [];
        let fpn = [];
        let ftm = [];
        let flc = [];
        let fcg = [];
        let fpx = [];
        let x1 = b[0];
        let x2 = b[1] - minlen + 1;
        if (invertedPCR) {
            x1 = b[2];
            x2 = b[3] - minlen + 1;
        }
        let fn = PrimerDesignForward(sq, plist, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, fpr, fpn, ftm, flc, fpx, fcg, (n + 1) + "F_", ends5, ends3, ftail, mplist);
        for (let d = 0; d < fn + 1; d++) {
            if (multiplex) {
                mplist.push(fpr[d]);
            }
            resultarea1 += fpn[d] + "\t" + fpr[d] + "\t" + fpr[d].length + "\t" + ftm[d].toFixed(1) + "\t" + fcg[d].toFixed(1) + "\t" + flc[d] + "\n";
        }

        let fppr = [];
        let fppn = [];
        let fptm = [];
        let fplc = [];
        let fpcg = [];
        let fppx = [];
        let fprn = -1;
        if (probedesign) {
            fprn = PrimerDesignForward(sq, plist, msk, minlc, probeminlen, probemaxlen, mintm + probetm, maxtm + probetm, x1, x2, primeroverlap, fppr, fppn, fptm, fplc, fppx, fpcg, (n + 1) + "F_" + probename, "h", "n", "", mplist);
            for (let d = 0; d < fprn + 1; d++) {
                if (multiplex) {
                    mplist.push(fppr[d]);
                }
                resultarea1 += fppn[d] + "\t" + fppr[d] + "\t" + fppr[d].length + "\t" + fptm[d].toFixed(1) + "\t" + fpcg[d].toFixed(1) + "\t" + fplc[d] + "\n";
            }
        }


        // reverse direction    
        if (maskrepeats) { msk = msk.reverse(); }
        let rpr = [];
        let rpn = [];
        let rtm = [];
        let rlc = [];
        let rpx = [];
        let rcg = [];
        x1 = lseqs - b[3];
        x2 = lseqs - b[2] + 1 - minlen;
        if (invertedPCR) {
            x1 = lseqs - b[1];
            x2 = lseqs - b[0] + 1 - minlen;
        }
        let rn = PrimerDesignReverse(cs, plist, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, (n + 1) + "R_", ends5, ends3, rtail, mplist);
        for (let d = rn; d > -1; d--) {
            if (multiplex) {
                mplist.push(rpr[d]);
            }
            rpx[d] = lseqs - (rpx[d] + rpr[d].length) + 1;
            resultarea1 += rpn[d] + "\t" + rpr[d] + "\t" + rpr[d].length + "\t" + rtm[d].toFixed(1) + "\t" + rcg[d].toFixed(1) + "\t" + rlc[d] + "\n";
        }

        let rppr = [];
        let rppn = [];
        let rptm = [];
        let rplc = [];
        let rppx = [];
        let rpcg = [];
        let rprn = -1;
        if (probedesign) {
            rprn = PrimerDesignReverse(cs, plist, msk, minlc, probeminlen, probemaxlen, mintm + probetm, maxtm + probetm, x1, x2, primeroverlap, rppr, rppn, rptm, rplc, rppx, rpcg, (n + 1) + "R_" + probename, "h", "n", "", mplist);
            for (let d = rprn; d > -1; d--) {
                if (multiplex) {
                    mplist.push(rppr[d]);
                }
                rppx[d] = lseqs - (rppx[d] + rppr[d].length) + 1;
                resultarea1 += rppn[d] + "\t" + rppr[d] + "\t" + rppr[d].length + "\t" + rptm[d].toFixed(1) + "\t" + rpcg[d].toFixed(1) + "\t" + rplc[d] + "\n";
            }
        }

 // combinations           
            if (fn > -1 && rn > -1) {
                if (multiplex) {
                    CombiningPairs(sq, invertedPCR, probedesign, minpcr, maxpcr, fn, rn, fpr, fpn, ftm, flc, fcg, fpx, fppr, fppn, fptm, fplc, fpcg, fppx, fprn, rpr, rpn, rtm, rlc, rpx, rcg, rppr, rppn, rptm, rplc, rppx, rpcg, rprn, multiplex, 'seq' + n);
                }
                else {
                    resultarea2 += name_seq[n] + ":\n";
                    resultarea2 += CombiningPairs(sq, invertedPCR, probedesign, minpcr, maxpcr, fn, rn, fpr, fpn, ftm, flc, fcg, fpx, fppr, fppn, fptm, fplc, fpcg, fppx, fprn, rpr, rpn, rtm, rlc, rpx, rcg, rppr, rppn, rptm, rplc, rppx, rpcg, rprn, multiplex, 'seq' + n);
                }
            }
    
        resultarea1 += "\n";
    }

    if (multiplex) {
        let z = 99;
        for (let n = 0; n < n_seq; n++) {
            const id = 'seq' + n;
            if (dataStore[id]) {
                let m = dataStore[id].length;
                if (z > m) { z = m; }
            }
        }
        for (let k = 0; k < z; k++) {
            for (let n = 0; n < n_seq; n++) {
                const id = 'seq' + n;
                if (dataStore[id]) {
                    let m = dataStore[id];
                    resultarea2 += m[k].fpn + "\t" + m[k].fpr + "\t" + m[k].fprl + "\t" + m[k].ftm.toFixed(1) + "\t" + m[k].fcg.toFixed(1) + "\t" + m[k].flc + "\n";
                    if ("fppn" in m[k]) {
                        resultarea2 += m[k].fppn + "\t" + m[k].fppr + "\t" + m[k].fppr.length + "\t" + m[k].fptm.toFixed(1) + "\t" + m[k].fpcg.toFixed(1) + "\t" + m[k].fplc + "\n";
                    }
                    if ("rppn" in m[k]) {
                        resultarea2 += m[k].rppn + "\t" + m[k].rppr + "\t" + m[k].rppr.length + "\t" + m[k].rptm.toFixed(1) + "\t" + m[k].rpcg.toFixed(1) + "\t" + m[k].rplc + "\n";
                    }
                    resultarea2 += m[k].rpn + "\t" + m[k].rpr + "\t" + m[k].rprl + "\t" + m[k].rtm.toFixed(1) + "\t" + m[k].rcg.toFixed(1) + "\t" + m[k].rlc + "\t" + m[k].prod + "/" + m[k].ta.toFixed(1) + "\n";
                }
            }
            resultarea2 += "\n";
        }
    }

    // End time
    let endTime = performance.now();
    let runtime = endTime - startTime;
    if (runtime > 999) {
        runtime = runtime / 1000;
        resultarea2 += "\nThe code took " + runtime.toFixed(0) + " seconds to run.\n\n";
    } else {
        resultarea2 += "\nThe code took " + runtime.toFixed(0) + " milliseconds to run.\n\n";
    }


    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

function CombiningPairs(sq, invertedPCR, probedesign, minpcr, maxpcr, fn, rn, fpr, fpn, ftm, flc, fcg, fpx, fppr, fppn, fptm, fplc, fpcg, fppx, fprn, rpr, rpn, rtm, rlc, rpx, rcg, rppr, rppn, rptm, rplc, rppx, rpcg, rprn, multiplex, id) {
    const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    const lseqs = sq.length;
    let resultarea = "";

    if (probedesign) {
        if (fprn > -1) {
            for (let f = fn; f > -1; f--) {
                let xf = fpx[f] + fpr[f].length;
                for (let m = 0; m <= fprn; m++) {
                    if (xf < fppx[m]) {
                        if (quickDimer(fpr[f], fppr[m]) === 0) {  // if (DimerLook2(fpr[f], fppr[m], min3) === -1) {
                            xf = fppx[m] + fppr[m].length;
                            for (let r = rn; r > -1; r--) {
                                let xr = rpx[r] + rpr[r].length - 1; //rpx[r] - rpr[r].length;
                                if (xf < rpx[r]) {
                                    let p = xr - fpx[f];
                                    if (minpcr === 0 || (minpcr > 0 && p >= minpcr && p <= maxpcr)) {
                                        if (quickDimer(fpr[f], rpr[r]) === -0) {      // if (DimerLook2(fpr[f], rpr[r], min3) === -1) {
                                            if (quickDimer(fppr[m], rpr[r]) === 0) { //  if (DimerLook2(fppr[m], rpr[r], min3) === -1) {
                                                let ta = Tm77(sq.substring(fpx[f], xr), KMg_M, 0);
                                                if (multiplex) {
                                                    addRecord(id, {
                                                        fpn: fpn[f], fpr: fpr[f], fprl: fpr[f].length, ftm: ftm[f], fcg: fcg[f], flc: flc[f],
                                                        fppn: fppn[m], fppr: fppr[m], fppl: fppr[m].length, fptm: fptm[m], fpcg: fpcg[m], fplc: fplc[m],
                                                        rpn: rpn[r], rpr: rpr[r], rprl: rpr[r].length, rtm: rtm[r], rcg: rcg[r], rlc: rlc[r], prod: p, ta: ta
                                                    });
                                                } else {
                                                    resultarea += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                                                    resultarea += fppn[m] + "\t" + fppr[m] + "\t" + fppr[m].length + "\t" + fptm[m].toFixed(1) + "\t" + fpcg[m].toFixed(1) + "\t" + fplc[m] + "\n";
                                                    resultarea += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n\n";
                                                }
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (rprn > -1) {
            if (!multiplex) resultarea += "\n";
            for (let f = fn; f > -1; f--) {
                let xf = fpx[f] + fpr[f].length;
                for (let m = 0; m <= rprn; m++) {
                    let xp = rppx[m] - rppr[m].length;
                    if (xf < xp) {
                        if (quickDimer(fpr[f], rppr[m]) === 0) {// if (DimerLook2(fpr[f], rppr[m], min3) === -1) {
                            for (let r = rn; r > -1; r--) {
                                let xr = rpx[r] + rpr[r].length - 1; //rpx[r] - rpr[r].length;
                                if (rppx[m] < rpx[r]) {
                                    let p = xr - fpx[f];
                                    if (minpcr === 0 || (minpcr > 0 && p >= minpcr && p <= maxpcr)) {
                                        if (quickDimer(fpr[f], rpr[r]) === 0) {     // if (DimerLook2(fpr[f], rpr[r], min3) === -1) {
                                            if (quickDimer(rppr[m], rpr[r]) === 0) {//  if (DimerLook2(rppr[m], rpr[r], min3) === -1) {
                                                let ta = Tm77(sq.substring(fpx[f], xr), KMg_M, 0);
                                                if (multiplex) {
                                                    addRecord(id, {
                                                        fpn: fpn[f], fpr: fpr[f], fprl: fpr[f].length, ftm: ftm[f], fcg: fcg[f], flc: flc[f],
                                                        rppn: rppn[m], rppr: rppr[m], rppl: rppr[m].length, rptm: rptm[m], rpcg: rpcg[m], rplc: rplc[m],
                                                        rpn: rpn[r], rpr: rpr[r], rprl: rpr[r].length, rtm: rtm[r], rcg: rcg[r], rlc: rlc[r], prod: p, ta: ta
                                                    });
                                                } else {
                                                    resultarea += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                                                    resultarea += rppn[m] + "\t" + rppr[m] + "\t" + rppr[m].length + "\t" + rptm[m].toFixed(1) + "\t" + rpcg[m].toFixed(1) + "\t" + rplc[m] + "\n";
                                                    resultarea += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n\n";
                                                }
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (invertedPCR) {
            for (let f = 0; f <= fn; f++) {
                for (let r = rn; r > -1; r--) {
                    if (!multiplex || rlc[r] > 0) {
                        let p = lseqs - (fpx[f] - rpx[r]) + rpr[r].length - 2;
                        if (minpcr === 0 || (minpcr > 0 && p >= minpcr && p <= maxpcr)) {
                            if (quickDimer(fpr[f], rpr[r]) === 0) {// if (DimerLook2(fpr[f], rpr[r], min3) === -1) {
                                let ta = Tm77(sq.substring(fpx[f], 0, rpx[r] + rpr[r].length - 1) + sq.substring(fpx[f]), KMg_M, 0);
                                if (multiplex) {
                                    addRecord(id, {
                                        fpn: fpn[f], fpr: fpr[f], fprl: fpr[f].length, ftm: ftm[f], fcg: fcg[f], flc: flc[f],
                                        rpn: rpn[r], rpr: rpr[r], rprl: rpr[r].length, rtm: rtm[r], rcg: rcg[r], rlc: rlc[r], prod: p, ta: ta
                                    });
                                    rlc[r] = 0;
                                    break;
                                } else {
                                    resultarea += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                                    resultarea += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n\n";
                                }
                            }
                        }
                    }
                }
            }
        } else {
            for (let f = 0; f <= fn; f++) {
                let xf = fpx[f] + fpr[f].length - 1;
                for (let r = rn; r > -1; r--) {
                    if (!multiplex || rlc[r] > 0) {
                        let xr = rpx[r] + rpr[r].length - 1; //  rpx[r] - rpr[r].length;
                        if (xf < rpx[r]) {
                            let p = xr - fpx[f];
                            if (minpcr === 0 || (minpcr > 0 && p >= minpcr && p <= maxpcr)) {
                                if (quickDimer(fpr[f], rpr[r]) === 0) { //  if (DimerLook2(fpr[f], rpr[r], min3) === -1) {
                                    let ta = Tm77(sq.substring(fpx[f], xr), KMg_M, 0);
                                    if (multiplex) {
                                        addRecord(id, {
                                            fpn: fpn[f], fpr: fpr[f], fprl: fpr[f].length, ftm: ftm[f], fcg: fcg[f], flc: flc[f],
                                            rpn: rpn[r], rpr: rpr[r], rprl: rpr[r].length, rtm: rtm[r], rcg: rcg[r], rlc: rlc[r], prod: p, ta: ta
                                        });
                                        rlc[r] = 0;
                                        break;
                                    } else {
                                        resultarea += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                                        resultarea += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return resultarea;
}

function PrimerDesignForward(sq, plist, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, fpr, fpn, ftm, flc, fpx, fcg, nm, e5, e3, tail, multiplist) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    // const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    const le5 = e5.length;
    const s5 = GeneratorVariants(e5).join(' ');

    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;

    let fn = -1;
    let q = 0;

    for (let y = x1; y < x1 + minlen - 1; y++) {
        if (msk[y] > 0) { q++; }
    }
    for (let x = x1; x < x2; x++) {
        if (msk[x + minlen - 1] > 0) { q++; }
        if (q === 0) {
            let x5 = x;
            let tm = maxtm + 1;
            let lc = 0;
            let q1 = 1;
            let s1 = sq.substring(x5 + minlen - le3, x5 + minlen);
            if (s3.indexOf(s1) > -1) {
                for (let d = 0; d < 1 + maxlen - minlen; d++) {
                    x5 = x - d;
                    if (x5 >= x1) {
                        s1 = sq.substring(x5, x5 + minlen + d);
                        if (s5.indexOf(s1.substring(0, le5)) > -1) {
                            tm = Tm(s1, salt, Mg_M, p_mkM);
                        }
                        if (tm >= mintm) {
                            break;
                        }
                    } else { break; }
                }

                if (tm >= mintm && tm <= maxtm) {
                    if (ssrrepeat(s1) === -1) {
                        lc = LingComplexity2(s1);
                        if (lc >= minlc) {
                            s1 = tail + s1;
                            if (quickDimer(s1) === 0) {//  if (DimerLook2(s1, s1, min3) === -1) {
                                q1 = 0;
                                if (plist.length > 0) {
                                    for (let j = 0; j < plist.length; j++) {
                                        if (quickDimer(s1, plist[j]) > 0) {// if (DimerLook2(s1, plist[j], min3) > -1) {
                                            q1 = 1;
                                            break;
                                        }
                                    }
                                }

                                if (multiplist.length > 0) {
                                    for (let j = 0; j < multiplist.length; j++) {
                                        if (quickDimer(s1, multiplist[j]) > 0) {// if (DimerLook2(s1, multiplist[j], min3) > -1) {
                                            q1 = 1;
                                            break;
                                        }
                                    }
                                }



                            }
                        }
                    }
                }
                if (q1 === 0) {
                    if (fn > -1) {
                        let p2 = fpx[fn] + fpr[fn].length - primeroverlap;
                        if (x5 < p2) {
                            if (flc[fn] < lc) {
                                fn--; q1 = 1;
                            }
                        } else { q1 = 1; }
                    }
                    else { q1 = 1; }

                    if (q1 === 1) {
                        fn++;
                        fpr[fn] = s1;
                        ftm[fn] = tm;
                        flc[fn] = lc;
                        fpx[fn] = x5;
                        fcg[fn] = CG(s1);
                        fpn[fn] = nm + (x5 + 1) + "-" + (x5 + s1.length);
                    }
                }
            }
        }
        if (msk[x] > 0) { q--; }
    }
    return fn;
}

function PrimerDesignReverse(cs, plist, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, nm, e5, e3, tail, multiplist) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;

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
            let tm = maxtm + 1;
            let lc = 0;
            let s1 = cs.substring(x5 + minlen - le3, x5 + minlen);
            if (s3.indexOf(s1) > -1) {
                for (let d = 0; d < 1 + maxlen - minlen; d++) {
                    x5 = x - d;
                    if (x5 >= x1) {
                        s1 = cs.substring(x5, x5 + minlen + d);
                        if (s5.indexOf(s1.substring(0, le5)) > -1) {
                            tm = Tm(s1, salt, Mg_M, p_mkM);
                        }
                        if (tm >= mintm) {
                            break;
                        }
                    } else { break; }
                }

                if (tm >= mintm && tm <= maxtm) {
                    if (ssrrepeat(s1) === -1) {
                        lc = LingComplexity2(s1);
                        if (lc >= minlc) {
                            s1 = tail + s1;
                            if (quickDimer(s1) === 0) {// if (DimerLook2(s1, s1, min3) === -1) {
                                q1 = 0;
                                if (plist.length > 0) {
                                    for (let j = 0; j < plist.length; j++) {
                                        if (quickDimer(s1, plist[j]) > 0) {// if (DimerLook2(s1, plist[j], min3) > -1) {
                                            q1 = 1;
                                            break;
                                        }
                                    }
                                }
                                if (multiplist.length > 0) {
                                    for (let j = 0; j < multiplist.length; j++) {
                                        if (quickDimer(s1, multiplist[j]) > 0) {// if (DimerLook2(s1, multiplist[j], min3) > -1) {
                                            q1 = 1;
                                            break;
                                        }
                                    }
                                }
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
                        } else { q1 = 1; }
                    }
                    else { q1 = 1; }

                    if (q1 === 1) {
                        rn++;
                        rpr[rn] = s1;
                        rtm[rn] = tm;
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
        if (b.length === 2) { // [SNP]
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