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
        b[1] = returnString.length - 1;
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
    }
    return { returnString, b, bx };
}

function analysis() {
    let mingap = parseInt(document.getElementById('mingap').value);
    let minlc = parseInt(document.getElementById('minlc').value);
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);
    const minpcr = parseInt(document.getElementById('minpcr').value);
    const maxpcr = parseInt(document.getElementById('maxpcr').value);
    const overlapping = document.getElementById('overlapping').checked;
    const maskrepeats = document.getElementById('repeats').checked;
    const ctconvert = document.getElementById('ctconvert').checked;

    const ftail = DNA(document.getElementById('ftail').value.toLowerCase().trim()); // "tttcacacaggaaacagctatgac";
    const rtail = DNA(document.getElementById('rtail').value.toLowerCase().trim()); // "tttcacacaggaaacagctatgac";

    const ends5 = "n";
    const flankadd = 50;
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

    let resultarea1 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n";
    let resultarea2 = "PrimerID\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\tFragment_Size(bp)/Tm(°C)\n\n";
    let startTime = performance.now();// Start time

    let mplist = [];

    for (let n = 0; n < n_seq; n++) {
        let panel1list = []; //multiplex list
        let panel2list = []; //multiplex list

        let resultpanel1 = "";
        let resultpanel2 = "";

        resultarea1 += name_seq[n] + ":\n";

        let RestSeq = Seq(seqs[n]);
        const b = RestSeq.b;
        const bx = RestSeq.bx;
        const seq = RestSeq.returnString;
        const lseqs = seq.length;
        let sq = "";
        let cs = "";

        let msk = new Array(lseqs).fill(0);
        let msk2 = new Array(lseqs).fill(0);
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
        if (maskrepeats) {
            msk2 = msk.reverse();
        }

        if (ctconvert) {
            sq = DNActBisulfiteForward(seq);
            cs = DNActBisulfiteReverse(seq);
        }
        else {
            sq = seq;
            cs = ComplementDNA(seq);
        }

        let t = 0;
        for (let i = 0; i < b.length; i += 2) {
            t++;
            if (b.length > 2) { resultarea1 += "Target" + (t) + "\n"; }

            let fpr = [];
            let fpn = [];
            let ftm = [];
            let flc = [];
            let fcg = [];
            let fpx = [];
            let x1 = b[i] - flankadd;
            let x2 = b[i + 1] + flankadd;

            if (x1 < 0) { x1 = 0; }
            if (x2 > lseqs) { x2 = lseqs - 1; }
            let fn = PrimerDesignForward(sq, plist, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, fpr, fpn, ftm, flc, fpx, fcg, (n + 1) + "F_", ends5, ends3, ftail, mplist);
            for (let d = 0; d < fn + 1; d++) {
                resultarea1 += fpn[d] + "\t" + fpr[d] + "\t" + fpr[d].length + "\t" + ftm[d].toFixed(1) + "\t" + fcg[d].toFixed(1) + "\t" + flc[d] + "\n";
            }

            // reverse direction    
            let rpr = [];
            let rpn = [];
            let rtm = [];
            let rlc = [];
            let rpx = [];
            let rcg = [];
            x1 = lseqs - b[i + 1] - flankadd;
            x2 = lseqs - b[i] + flankadd;

            if (x1 < 0) { x1 = 0; }
            if (x2 > lseqs) { x2 = lseqs - 1; }

            let rn = PrimerDesignReverse(cs, plist, msk2, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, (n + 1) + "R_", ends5, ends3, rtail, mplist);
            for (let d = rn; d > -1; d--) {
                rpx[d] = lseqs - (rpx[d] + rpr[d].length) + 1;
                resultarea1 += rpn[d] + "\t" + rpr[d] + "\t" + rpr[d].length + "\t" + rtm[d].toFixed(1) + "\t" + rcg[d].toFixed(1) + "\t" + rlc[d] + "\n";
            }
            const result = TilingPanels(sq, mingap, minlen, minpcr, maxpcr, fn, rn, fpr, fpn, ftm, flc, fcg, fpx, rpr, rpn, rtm, rlc, rpx, rcg, panel1list, panel2list);
            resultpanel1 += result.resultarea1;
            resultpanel2 += result.resultarea2;
        }

        resultarea1 += "\n";
        resultarea2 += "Panel A\n";
        resultarea2 += resultpanel1;
        resultarea2 += "Panel B\n";
        resultarea2 += resultpanel2;
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

function TilingPanels(sq, mingap, minlen, minpcr, maxpcr, fn, rn, fpr, fpn, ftm, flc, fcg, fpx, rpr, rpn, rtm, rlc, rpx, rcg, panel1list, panel2list) {
    const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    const mpcr = maxpcr - 100;
    if (minpcr < mpcr) { minpcr = mpcr; }
    const averpcr = mingap + ((1 + minpcr - minlen) / 2) | 0;

    // Panel 1
    let resultarea1 = "";
    let x1p1 = [];   // coordinates for aplicons
    let x2p1 = [];   // coordinates for aplicons
    let n = -1;      // PCR amplicon last x
    for (let f = 0; f <= fn; f++) {
        if (flc[f] < 1 || fpx[f] < n) continue;
        if (panel1list.some(d => quickDimer(fpx[f], d) > 0)) continue;
        const xf = fpx[f] + fpr[f].length - 1;

        for (let r = rn; r > -1; r--) {
            if (rlc[r] > 0 && rpx[r] > n) {
                const xr = rpx[r] + rpr[r].length - 1;
                if (xf < rpx[r]) {
                    const p = xr - fpx[f];

                    if (p < minpcr || p > maxpcr) continue;
                    if (quickDimer(fpr[f], rpr[r]) > 0) continue;
                    if (panel1list.some(d => quickDimer(rpr[r], d) > 0)) continue;
                    panel1list.push(fpr[f], rpr[r]);

                    const ta = Tm77(sq.substring(fpx[f], xr), KMg_M, 0);
                    x1p1.push(fpx[f]);
                    x2p1.push(xr);

                    resultarea1 += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                    resultarea1 += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n";

                    n = xr + mingap;
                    rlc[r] = 0;
                    flc[f] = 0;
                    break;
                }
            }
        }
    }

    // Panel 2
    let resultarea2 = "";
    let x1p2 = [];   // coordinates for aplicons
    let x2p2 = [];   // coordinates for aplicons 

    for (let t = 0; t < x1p1.length; t++) {
        const minStart = x1p1[t] + averpcr;

        for (let f = 0; f <= fn; f++) {
            if (flc[f] < 1 || fpx[f] < minStart) continue;
            if (panel2list.some(d => quickDimer(fpx[f], d) > 0)) continue;

            for (let r = rn; r > -1; r--) {
                if (rlc[r] > 0 && rpx[r] > x1p1[t] + averpcr) {

                    const xr = rpx[r] + rpr[r].length - 1;
                    const p = xr - fpx[f];

                    if (p < minpcr || p > maxpcr) continue;
                    if (quickDimer(fpr[f], rpr[r]) > 0) continue;
                    if (panel2list.some(d => quickDimer(rpr[r], d) > 0)) continue;
                    const ta = Tm77(sq.substring(fpx[f], xr), KMg_M, 0);
                    panel2list.push(fpr[f], rpr[r]);
                    x1p2.push(fpx[f]);
                    x2p2.push(xr);

                    resultarea2 += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                    resultarea2 += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n";

                    rlc[r] = 0;
                    flc[f] = 0;
                    f = fn;
                    break;
                }
            }
        }
    }
    return { resultarea1, resultarea2 };
}