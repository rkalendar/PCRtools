function analysis() {
    let minlc = parseInt(document.getElementById('minlc').value);
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);
    let snp = parseInt(document.getElementById('snp').value);

    let overlapping = document.getElementById('overlapping').checked;
    let maskrepeats = document.getElementById('repeats').checked;
    let ctconvert = document.getElementById('ctconvert').checked;

    let ends3 = document.getElementById('end3com').value.toLowerCase().trim();   // "swh ssw wsh sww www";
    if (DNA(ends3).length < 1) { ends3 = "n"; }


    if (minlc < 10) { minlc = 10; }
    if (minlc > 90) { minlc = 90; }
    if (minlen < 12) { minlen = 12; }
    if (maxlen < 12) { maxlen = 12; }
    if (minlen > 100) { minlen = 100; }
    if (maxlen > 500) { maxlen = 500; }
    if (minlen > maxlen) { maxlen = minlen; }
    if (mintm < 37) { mintm = 37; }
    if (mintm > 80) { mintm = 80; }
    if (snp < 1) { snp = 1; }

    let primeroverlap = 12;
    if (overlapping) { primeroverlap = minlen; }

    ReadResult = ReadingSeq(document.getElementById('inputTails').value);
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
    let resultarea1 = "";
    let resultarea2 = "";
    for (let n = 0; n < n_seq; n++) {
        resultarea1 += "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n" + name_seq[n] + ":\n";
        resultarea2 += "PrimerID\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\tFragment_Size(bp)/Tm(°C)\n" + name_seq[n] + ":\n";

        let resl1 = AsPCR(snp, seqs[n], tails, maskrepeats, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3);
        resultarea1 += resl1.resultarea1;
        resultarea2 += resl1.resultarea2;

        //REVERSE DIRECTION
        let resl2 = AsPCR(snp, ComplementDNA3(seqs[n]), tails, maskrepeats, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3);
        resultarea2 += "\n";
        resultarea1 += resl2.resultarea1;
        resultarea2 += resl2.resultarea2;
    }
    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

function ProbeDesignForward(snp, tails, s1, s2, pol, minlen, maxlen, mintm, nm) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 6;
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

    for (let i = 0; i < n; i++) {
        if (minpol > pol[i].length) { minpol = pol[i].length; }
        s[i] = (pol[i] + s2);
        t[i] = (s1 + pol[i] + s2);
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
        for (let i = minlen; i < maxlen + minlen; i++) {
            if (x1 - i < 0) { break; }
            let v = t[j].substring(x1 - i, x1);
            fpx[j] = x1 - i;
            fpr[j] = v;
            ftm[j] = Tm(fpr[j], salt, Mg_M, p_mkM);
            if (ftm[j] > mintm) {
                break;
            }
        }
        flc[j] = LingComplexity2(fpr[j]);
        fcg[j] = CG(fpr[j]);
        fpn[j] = nm + (fpx[j] + 1) + "-" + (fpx[j] + fpr[j].length + 1) + "_ASP" + (j + 1);
    }

    if (ntails > 0) {
        for (let i = 0; i < Math.min(ntails, n); i++) {
            fpr[i] = tails[i] + fpr[i];
        }
    }

    for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
            if (DimerLook2(fpr[i].toLowerCase(), fpr[j].toLowerCase(), min3) > 0) {
                flc[j] = 0;
                flc[i] = 0;
            }
        }
    }
    return { fpr, fpn, ftm, flc, fpx, fcg };
}

function PrimerDesignReverse(cs, msk, prlist, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, nm, e5, e3) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 6;


    const le5 = e5.length;
    const s5 = GeneratorVariants(e5).join(' ');

    s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;


    const lseqs = cs.length;
    const nprlist = prlist.length;

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
                    if (x5 >= x1) {
                        s1 = cs.substring(x5, x5 + minlen + d);
                        if (s5.indexOf(s1.substring(0, le5)) > -1) {
                            tm1 = Tm(s1, salt, Mg_M, p_mkM);
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

                    if (nprlist > 0 && q1 === 0) {
                        for (let j = 0; j < nprlist; j++) {
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

function Seqs(str) {
    let s1 = "";
    let s2 = "";
    let s3 = "";
    let s = str.split("[");
    if (s.length > 1) {
        s1 = DNA(s[0].trim());
        let p = s[1].split("]");
        if (p.length > 1) {
            s2 = p[0].trim();
            s3 = DNA(p[1].trim());
        }
    }
    else {
        return { s1, s2, s3 };
    }
    return { s1, s2, s3 };
}

function AsPCR(snp, inseq, tails, maskrepeats, minlc, minlen, maxlen, mintm, maxtm, primeroverlap, ctconvert, ends3) {
    const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    let resultarea1 = "";
    let resultarea2 = "";
    let rst = Seqs(inseq);
    let s1 = rst.s1;
    let s2 = rst.s3;
    let pol = rst.s2;
    let st = GeneratorVariants(pol);
    let seq = s1 + st.join('') + s2;
    const lseqs = seq.length;

    let sq = "";
    let cs = "";
    const end5 = "n";
    const end3 = ends3;


    let msk = new Array(lseqs).fill(0);
    if (maskrepeats) {
        msk = RepeatMask(seq);
    }

    if (ctconvert) {
        sq = DNActBisulfiteForward(seq);
        cs = DNActBisulfiteReverse(seq);
    }
    else {
        sq = seq;
        cs = ComplementDNA(seq);
    }

    //Forward direction - Probe design
    let result = ProbeDesignForward(snp, tails, s1, s2, st, minlen, maxlen, mintm, "F_");
    let fpr = result.fpr;
    let fpn = result.fpn;
    let ftm = result.ftm;
    let flc = result.flc;
    let fcg = result.fcg;
    let fpx = result.fpx;
    let fn = fpr.length;
    for (let d = 0; d < fn; d++) {
        resultarea1 += fpn[d] + "\t" + fpr[d] + "\t" + fpr[d].length + "\t" + ftm[d].toFixed(1) + "\t" + fcg[d].toFixed(1) + "\t" + flc[d] + "\n";
    }

    // reverse direction    
    if (maskrepeats) { msk = msk.reverse(); }
    let rpr = [];
    let rpn = [];
    let rtm = [];
    let rlc = [];
    let rpx = [];
    let rcg = [];
    x1 = 0;
    x2 = s2.length - minlen;
    let rn = PrimerDesignReverse(cs, msk, fpr, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, "R_", end5, end3);
    for (let d = rn; d > -1; d--) {
        rpx[d] = lseqs - (rpx[d] + rpr[d].length) + 1;
        resultarea1 += rpn[d] + "\t" + rpr[d] + "\t" + rpr[d].length + "\t" + rtm[d].toFixed(1) + "\t" + rcg[d].toFixed(1) + "\t" + rlc[d] + "\n";
    }
    resultarea1 += "\n\n";
    // combinations
    if (fn > -1 && rn > -1) {
        let x3asp = fpx[0] + fpr[0].length;
        for (let r = rn; r > -1; r--) {
            let x3r = rpx[r];
            if (x3asp < x3r) {
                let p = rpx[r] - fpx[0] + rpr[r].length - 1;
                let ta = Tm77(sq.substring(fpx[0], x3r + rpr[r].length - 1), KMg_M, 0);
                for (let f = 0; f < fn; f++) {
                    resultarea2 += fpn[f] + "\t" + fpr[f] + "\t" + fpr[f].length + "\t" + ftm[f].toFixed(1) + "\t" + fcg[f].toFixed(1) + "\t" + flc[f] + "\n";
                }
                resultarea2 += rpn[r] + "\t" + rpr[r] + "\t" + rpr[r].length + "\t" + rtm[r].toFixed(1) + "\t" + rcg[r].toFixed(1) + "\t" + rlc[r] + "\t" + p + "/" + ta.toFixed(1) + "\n\n";
            }
        }
    }
    return { resultarea1, resultarea2 };
}