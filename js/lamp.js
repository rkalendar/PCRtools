function analysis() {
    let minlc = parseInt(document.getElementById('minlc').value);
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);
    let maxf2b2 = parseInt(document.getElementById('maxpcr').value);
    let overlapping = document.getElementById('overlapping').checked;
    let maskrepeats = document.getElementById('repeats').checked;
    let ctconvert = document.getElementById('ctconvert').checked;
    let loopprimer = document.getElementById('loopprimer').checked;

    const ends5 = "n";
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

    let primeroverlap = 12;
    if (overlapping) { primeroverlap = 18; }

    const KMg_M = 0.055 + (3.795 * Math.sqrt(0.001));
    const min3 = 6;    // min3=5 makes LAMP no result
    const dd = 40;     // maximal distance between primers in group
    const f3_f2 = 20;  // maximal distance between F3-F2 primers 
    const lmg = 5000;   // maximal number groups for each side
    const fmg = 10;   // maximal number groups for each side
    /*
          The distance between 5' end of F2 and B2 is considered to be 120-200 bp and the distance between F2 and F3 as well as B2 and B3 is 0-20bp. 
          The distance for loop forming regions (5' of F2 to 3' of F1, 5' of B2 to 3' of B1) is 0-40bp.
   */
    ReadResult = ReadingSeq(document.getElementById('inputText').value);
    const name_seq = ReadResult.name_seq;
    const seqs = ReadResult.seqs;
    const n_seq = seqs.length;

    const result = ["", ""];
    let resultarea1 = "";
    let resultarea2 = "";
    for (let n = 0; n < n_seq; n++) {
        resultarea1 += "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n" + name_seq[n] + ":\n";

        ResSeq = Seq(seqs[n]);
        const seq = ResSeq.seq;
        const lseqs = seq.length;
        const b = ResSeq.b;
        const bx = ResSeq.bx;
        let sq = "";
        let cs = "";
        if (ctconvert) {
            sq = DNActBisulfiteForward(seq);
            cs = DNActBisulfiteReverse(seq);
        }
        else {
            sq = seq;
            cs = ComplementDNA(seq);
        }

        let msk = new Array(lseqs).fill(0);
        if (maskrepeats) {
            msk = RepeatMask(sq);
        }
        if (bx.length > 1) {
            for (let i = 0; i < bx.length; i += 2) {
                for (let j = bx[i]; j < bx[i + 1]; j++) {
                    msk[j]++;
                }
            }
        }

        let fpr = [];
        let fpn = [];
        let ftm = [];
        let flc = [];
        let fcg = [];
        let fpx = [];
        let x1 = b[0];
        let x2 = b[1] - minlen + 1;
        let fn = PrimerDesignForward(sq, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, fpr, fpn, ftm, flc, fpx, fcg, "F_", ends5, ends3);
        for (let d = 0; d < fn + 1; d++) {
            resultarea1 += fpn[d] + "\t" + fpr[d] + "\t" + fpr[d].length + "\t" + NumberToSeq(ftm[d], 1) + "\t" + NumberToSeq(fcg[d], 1) + "\t" + flc[d] + "\n";
        }

        let f1pr = [];
        let f1pn = [];
        let f1tm = [];
        let f1lc = [];
        let f1cg = [];
        let f1px = [];
        let f1n = PrimerDesignForward(sq, msk, minlc, minlen, maxlen, mintm + 5, maxtm + 5, x1, x2, primeroverlap, f1pr, f1pn, f1tm, f1lc, f1px, f1cg, "F1_", ends5, ends3);
        for (let d = 0; d < f1n + 1; d++) {
            resultarea1 += f1pn[d] + "\t" + f1pr[d] + "\t" + f1pr[d].length + "\t" + NumberToSeq(f1tm[d], 1) + "\t" + NumberToSeq(f1cg[d], 1) + "\t" + f1lc[d] + "\n";
        }

        let flooppr = [];
        let flooppn = [];
        let flooptm = [];
        let flooplc = [];
        let flooppx = [];
        let floopcg = [];
        let floopn = -1;
        if (loopprimer) {
            x1 = lseqs - x2;
            x2 = lseqs - x1 + 1 - minlen;
            floopn = PrimerDesignReverse(cs, msk, minlc, minlen, maxlen, mintm + 5, maxtm + 5, x1, x2, primeroverlap, flooppr, flooppn, flooptm, flooplc, flooppx, floopcg, "FL_", ends5, ends3);
            for (let d = floopn; d > -1; d--) {
                flooppx[d] = lseqs - flooppx[d] + 1;
                resultarea1 += flooppn[d] + "\t" + flooppr[d] + "\t" + flooppr[d].length + "\t" + NumberToSeq(flooptm[d], 1) + "\t" + NumberToSeq(floopcg[d], 1) + "\t" + flooplc[d] + "\n";
            }
        }

        //save Forward records
        let fset = [];
        let nfset = -1;
        if (loopprimer) {
            if (fn > 0 && f1n > 0 && floopn > 0) {

                //   F3 fpx[f3]---------> x3f3    F2 fpx[f2]---------> x3f2     FL  x3fl[fl] <---------------- flooppx[fl]       F1 f1px[f1]---------> x3f1
                // 5---------------------------------------------------------------------------------------------------------------------------------            
                outerLoop: for (let f3 = 0; f3 < fn; f3++) {
                    let x3f3 = fpx[f3] + fpr[f3].length - 1;
                    let f2count = 0;
                    Loopf2: for (let f2 = f3 + 1; f2 < fn; f2++) {
                        let x3f2 = fpx[f2] + fpr[f2].length - 1;
                        if (x3f3 + f3_f2 < fpx[f2]) { break; }
                        if (x3f3 < fpx[f2]) {
                            for (let d = floopn; d > -1; d--) {
                                let x3fl = flooppx[d] - flooppr[d].length + 1;
                                if (x3f2 + dd < x3fl) { break; }
                                if (x3f2 < x3fl) {
                                    for (let f1 = 0; f1 < f1n + 1; f1++) {
                                        if (flooppx[d] + dd < f1px[f1]) { break; }
                                        if (flooppx[d] < f1px[f1]) {
                                            let a1 = [];
                                            a1.push(fpr[f3]);
                                            a1.push(fpr[f2]);
                                            a1.push(f1pr[f1]);
                                            a1.push(flooppr[d]);
                                            if (DimerList(a1, min3) === -1) {
                                                fset.push(f3);
                                                fset.push(f2);
                                                fset.push(d);
                                                fset.push(f1);
                                                fset.push(fpx[f2]);
                                                fset.push(f1px[f1] + f1pr[f1].length);
                                                nfset++;
                                                f2count++;
                                                if (f2count > fmg) { break Loopf2; }
                                                //    resultarea2 += "F3_" + fpn[f3] + "\t" + fpr[f3] + "\t" + fpr[f3].length + "\t" + NumberToSeq(ftm[f3], 1) + "\t" + NumberToSeq(fcg[f3], 1) + "\t" + flc[f3] + "\n";
                                                //    resultarea2 += "F2_" + fpn[f2] + "\t" + fpr[f2] + "\t" + fpr[f2].length + "\t" + NumberToSeq(ftm[f2], 1) + "\t" + NumberToSeq(fcg[f2], 1) + "\t" + flc[f2] + "\n";
                                                //    resultarea2 += flooppn[d] + "\t" + flooppr[d] + "\t" + flooppr[d].length + "\t" + NumberToSeq(flooptm[d], 1) + "\t" + NumberToSeq(floopcg[d], 1) + "\t" + flooplc[d] + "\n";
                                                //    resultarea2 += f1pn[f1] + "\t" + f1pr[f1] + "\t" + f1pr[f1].length + "\t" + NumberToSeq(f1tm[f1], 1) + "\t" + NumberToSeq(f1cg[f1], 1) + "\t" + f1lc[f1] + "\n\n";
                                            }
                                            if (nfset > lmg) { break outerLoop; }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (fn > 0 && f1n > 0) {
                //   F3 fpx[f3]---------> x3f3    F2 fpx[f2]---------> x3f2          F1 fpx[f1]---------> x3f1
                // 5----------------------------------------------------------------------------------------------------                   
                outerLoop: for (let f3 = 0; f3 < fn; f3++) {
                    let x3f3 = fpx[f3] + fpr[f3].length - 1;
                    let f2count = 0;
                    Loopf2: for (let f2 = f3 + 1; f2 < fn; f2++) {
                        if (x3f3 + f3_f2 < fpx[f2]) { break; }
                        if (x3f3 < fpx[f2]) {
                            let x3f2 = fpx[f2] + fpr[f2].length - 1;
                            for (let f1 = 0; f1 < f1n + 1; f1++) {
                                if (f1px[f1] > x3f2 + dd) { break; }
                                if (f1px[f1] > x3f2) {
                                    let a1 = [];
                                    a1.push(fpr[f3]);
                                    a1.push(fpr[f2]);
                                    a1.push(f1pr[f1]);
                                    if (DimerList(a1, min3) === -1) {
                                        fset.push(f3);
                                        fset.push(f2);
                                        fset.push(f1);
                                        fset.push(fpx[f2]);
                                        fset.push(f1px[f1] + f1pr[f1].length);
                                        nfset++;
                                        f2count++;
                                        if (f2count > fmg) { break Loopf2; }
                                        //   resultarea2 += "F3_" + fpn[f3] + "\t" + fpr[f3] + "\t" + fpr[f3].length + "\t" + NumberToSeq(ftm[f3], 1) + "\t" + NumberToSeq(fcg[f3], 1) + "\t" + flc[f3] + "\n";
                                        //   resultarea2 += "F2_" + fpn[f2] + "\t" + fpr[f2] + "\t" + fpr[f2].length + "\t" + NumberToSeq(ftm[f2], 1) + "\t" + NumberToSeq(fcg[f2], 1) + "\t" + flc[f2] + "\n";
                                        //   resultarea2 += f1pn[f1] + "\t" + f1pr[f1] + "\t" + f1pr[f1].length + "\t" + NumberToSeq(f1tm[f1], 1) + "\t" + NumberToSeq(f1cg[f1], 1) + "\t" + f1lc[f1] + "\n\n";
                                    }
                                    if (nfset > lmg) { break outerLoop; }
                                }
                            }
                        }
                    }
                }
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
        let rn = PrimerDesignReverse(cs, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, "R_", ends5, ends3);
        for (let d = rn; d > -1; d--) {
            rpx[d] = lseqs - (rpx[d] + rpr[d].length) + 1;
            resultarea1 += rpn[d] + "\t" + rpr[d] + "\t" + rpr[d].length + "\t" + NumberToSeq(rtm[d], 1) + "\t" + NumberToSeq(rcg[d], 1) + "\t" + rlc[d] + "\n";
        }

        let r1pr = [];
        let r1pn = [];
        let r1tm = [];
        let r1lc = [];
        let r1px = [];
        let r1cg = [];
        x1 = lseqs - b[3];
        x2 = lseqs - b[2] + 1 - minlen;
        let r1n = PrimerDesignReverse(cs, msk, minlc, minlen, maxlen, mintm + 5, maxtm + 5, x1, x2, primeroverlap, r1pr, r1pn, r1tm, r1lc, r1px, r1cg, "B1_", ends5, ends3);
        for (let d = r1n; d > -1; d--) {
            r1px[d] = lseqs - (r1px[d] + r1pr[d].length) + 1;
            resultarea1 += r1pn[d] + "\t" + r1pr[d] + "\t" + r1pr[d].length + "\t" + NumberToSeq(r1tm[d], 1) + "\t" + NumberToSeq(r1cg[d], 1) + "\t" + r1lc[d] + "\n";
        }

        let rlooppr = [];
        let rlooppn = [];
        let rlooptm = [];
        let rlooplc = [];
        let rloopcg = [];
        let rlooppx = [];
        let rloopn = -1;
        if (loopprimer) {
            x1 = b[2];
            x2 = b[3] - minlen + 1;
            rloopn = PrimerDesignForward(sq, msk, minlc, minlen, maxlen, mintm + 5, maxtm + 5, x1, x2, primeroverlap, rlooppr, rlooppn, rlooptm, rlooplc, rlooppx, rloopcg, "BL_", ends5, ends3);
            for (let d = 0; d < rloopn + 1; d++) {
                resultarea1 += rlooppn[d] + "\t" + rlooppr[d] + "\t" + rlooppr[d].length + "\t" + NumberToSeq(rlooptm[d], 1) + "\t" + NumberToSeq(rloopcg[d], 1) + "\t" + rlooplc[d] + "\n";
            }
        }

        //save Reverse records
        let rset = [];
        let nrset = -1;
        if (loopprimer) {
            if (rn > 0 && r1n > 0 && rloopn > 0) {
                // 5---------------------------------------------------------------------------------------------------------------------------------  
                // R1 r1px[r1] <------------x5r1  RL rlooppx[rl] ---------->x3rl    rpx[r2]<------------x5r2 R2      rpx[r3]<------------x5r3 R3 
                outerLoop: for (let r3 = 0; r3 < rn; r3++) {
                    let f2count = 0;
                    Loopf2: for (let r2 = r3 + 1; r2 < rn; r2++) {
                        let x5r2 = rpx[r2] + rpr[r2].length - 1;
                        if (rpx[r3] > x5r2 + f3_f2) { break; }
                        if (rpx[r3] > x5r2) {
                            for (let d = rloopn; d > -1; d--) {
                                let x3rl = rlooppx[d] + rlooppr[d].length - 1;
                                if (rpx[r2] > x3rl + dd) { break; }
                                if (rpx[r2] > x3rl) {
                                    for (let r1 = 0; r1 <= r1n; r1++) {
                                        let x5r1 = r1px[r1] + r1pr[r1].length - 1;
                                        if (rlooppx[d] > x5r1 + dd) { break; }
                                        if (rlooppx[d] > x5r1) {
                                            let a1 = [];
                                            a1.push(rpr[r3]);
                                            a1.push(rpr[r2]);
                                            a1.push(r1pr[r1]);
                                            a1.push(rlooppr[d]);
                                            if (DimerList(a1, min3) === -1) {
                                                let x5r2 = rpx[r2] + rpr[r2].length - 1;
                                                rset.push(r3);
                                                rset.push(r2);
                                                rset.push(d);
                                                rset.push(r1);
                                                rset.push(r1px[r1]);
                                                rset.push(x5r2);
                                                nrset++;
                                                f2count++;
                                                if (f2count > fmg) { break Loopf2; }
                                                //       resultarea2 += "R3_" + rpn[r3] + "\t" + rpr[r3] + "\t" + rpr[r3].length + "\t" + NumberToSeq(rtm[r3], 1) + "\t" + NumberToSeq(rcg[r3], 1) + "\t" + rlc[r3] + "\n";
                                                //       resultarea2 += "R2_" + rpn[r2] + "\t" + rpr[r2] + "\t" + rpr[r2].length + "\t" + NumberToSeq(rtm[r2], 1) + "\t" + NumberToSeq(rcg[r2], 1) + "\t" + rlc[r2] + "\n";
                                                //       resultarea2 += rlooppn[d] + "\t" + rlooppr[d] + "\t" + rlooppr[d].length + "\t" + NumberToSeq(rlooptm[d], 1) + "\t" + NumberToSeq(rloopcg[d], 1) + "\t" + rlooplc[d] + "\n";
                                                //       resultarea2 += r1pn[r1] + "\t" + r1pr[r1] + "\t" + r1pr[r1].length + "\t" + NumberToSeq(r1tm[r1], 1) + "\t" + NumberToSeq(r1cg[r1], 1) + "\t" + r1lc[r1] + "\n\n";
                                            }
                                            if (nrset > lmg) { break outerLoop; }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (rn > 0 && r1n > 0) {
                outerLoop: for (let r3 = 0; r3 < rn; r3++) {
                    // 5----------------------------------------------------------------------------------------------  
                    // R1 r1px[r1] <------------x5r1    rpx[r2]<------------x5r2 R2      rpx[r3]<------------x5r3 R3                    
                    let f2count = 0;
                    Loopf2: for (let r2 = r3 + 1; r2 < rn; r2++) {
                        let x5r2 = rpx[r2] + rpr[r2].length - 1;
                        if (rpx[r3] > x5r2 + f3_f2) { break; }
                        if (rpx[r3] > x5r2) {
                            for (let r1 = 0; r1 <= r1n; r1++) {
                                let x5r1 = r1px[r1] + r1pr[r1].length - 1;
                                if (rpx[r2] > x5r1 + dd) { break; }
                                if (rpx[r2] > x5r1) {
                                    let a1 = [];
                                    a1.push(rpr[r3]);
                                    a1.push(rpr[r2]);
                                    a1.push(r1pr[r1]);
                                    if (DimerList(a1, min3) === -1) {
                                        let x5r2 = rpx[r2] + rpr[r2].length - 1;
                                        rset.push(r3);
                                        rset.push(r2);
                                        rset.push(r1);
                                        rset.push(r1px[r1]);
                                        rset.push(x5r2);
                                        nrset++;
                                        f2count++;
                                        if (f2count > fmg) { break Loopf2; }
                                        //    resultarea2 += "R3_" + rpn[r3] + "\t" + rpr[r3] + "\t" + rpr[r3].length + "\t" + NumberToSeq(rtm[r3], 1) + "\t" + NumberToSeq(rcg[r3], 1) + "\t" + rlc[r3] + "\n";
                                        //    resultarea2 += "R2_" + rpn[r2] + "\t" + rpr[r2] + "\t" + rpr[r2].length + "\t" + NumberToSeq(rtm[r2], 1) + "\t" + NumberToSeq(rcg[r2], 1) + "\t" + rlc[r2] + "\n";
                                        //    resultarea2 += r1pn[r1] + "\t" + r1pr[r1] + "\t" + r1pr[r1].length + "\t" + NumberToSeq(r1tm[r1], 1) + "\t" + NumberToSeq(r1cg[r1], 1) + "\t" + r1lc[r1] + "\n\n";
                                    }
                                    if (nrset > lmg) { break outerLoop; }
                                }
                            }
                        }
                    }
                }
            }
        }
        resultarea1 += "\n\n";

        if (nrset > -1 && nfset > -1) {
            resultarea2 += "PrimerID\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\tFragment_Size(bp)/Tm(°C)\n" + name_seq[n] + ":\n";

            let col = [];
            let ncol = -1;

            if (loopprimer) {
                for (let f = 0; f <= nfset; f++) {
                    let f1 = new Array(6).fill(0);
                    let v = 6 * f;
                    f1[0] = fset[v];
                    f1[1] = fset[v + 1];
                    f1[2] = fset[v + 2];
                    f1[3] = fset[v + 3];
                    f1[4] = fset[v + 4];// 5-F3
                    f1[5] = fset[v + 5];// F1->

                    for (let b = 0; b <= nrset; b++) {
                        let r1 = new Array(6).fill(0);
                        let v = 6 * b;
                        r1[0] = rset[v];
                        r1[1] = rset[v + 1];
                        r1[2] = rset[v + 2];
                        r1[3] = rset[v + 3];
                        r1[4] = rset[v + 4]; //<-B1
                        r1[5] = rset[v + 5]; //B3-5

                        let z = r1[5] - f1[4];

                        if (f1[5] <= r1[4] && z <= maxf2b2) {
                            let a1 = [];
                            a1.push(fpr[f1[0]]);
                            a1.push(fpr[f1[1]]);
                            a1.push(flooppr[f1[2]]);
                            a1.push(f1pr[f1[3]]);
                            let zlc = flc[f1[0]] + flc[f1[1]] + flooplc[f1[2]] + f1lc[f1[3]];
                            let a2 = [];
                            a2.push(r1pr[r1[3]]);
                            a2.push(rlooppr[r1[2]]);
                            a2.push(rpr[r1[1]]);
                            a2.push(rpr[r1[0]]);
                            zlc = zlc + rlc[r1[0]] + rlc[r1[1]] + rlooplc[f1[2]] + r1lc[r1[3]];

                            if (ncol > -1) {
                                if (col[col.length - 2] > z || col[col.length - 1] < zlc) {
                                    if (DimerList2(a1, a2, min3) === -1) {
                                        ncol++;
                                        col.push(f);
                                        col.push(b);
                                        col.push(z);
                                        col.push(zlc);
                                    }
                                }
                            } else {
                                if (DimerList2(a1, a2, min3) === -1) {
                                    ncol++;
                                    col.push(f);
                                    col.push(b);
                                    col.push(z);
                                    col.push(zlc);
                                }
                            }
                        }
                    }
                }

                // printing collections            
                for (let f = 0; f <= ncol; f++) {
                    let v = 4 * f;
                    let r = new Array(4).fill(0);
                    r[0] = col[v];     // f set
                    r[1] = col[v + 1]; // r set
                    r[2] = col[v + 2]; // size F3-B3
                    r[3] = col[v + 3]; // total LC

                    let f1 = new Array(6).fill(0);
                    let vf = 6 * r[0];
                    f1[0] = fset[vf];
                    f1[1] = fset[vf + 1];
                    f1[2] = fset[vf + 2];
                    f1[3] = fset[vf + 3];
                    f1[4] = fset[vf + 4];// 5-F3
                    f1[5] = fset[vf + 5];// F1->

                    let r1 = new Array(6).fill(0);
                    let vr = 6 * r[1];
                    r1[0] = rset[vr];
                    r1[1] = rset[vr + 1];
                    r1[2] = rset[vr + 2];
                    r1[3] = rset[vr + 3];
                    r1[4] = rset[vr + 4];//<-B1
                    r1[5] = rset[vr + 5];//B3-5

                    let a1 = [];
                    a1.push(fpr[f1[0]]);
                    a1.push(fpr[f1[1]]);
                    a1.push(flooppr[f1[2]]);
                    a1.push(f1pr[f1[3]]);
                    let a2 = [];
                    a2.push(r1pr[r1[3]]);
                    a2.push(rlooppr[r1[2]]);
                    a2.push(rpr[r1[1]]);
                    a2.push(rpr[r1[0]]);

                    resultarea2 += "F3_" + fpn[f1[0]] + "\t" + fpr[f1[0]] + "\t" + fpr[f1[0]].length + "\t" + NumberToSeq(ftm[f1[0]], 1) + "\t" + NumberToSeq(fcg[f1[0]], 1) + "\t" + flc[f1[0]] + "\n";
                    resultarea2 += "F2_" + fpn[f1[1]] + "\t" + fpr[f1[1]] + "\t" + fpr[f1[1]].length + "\t" + NumberToSeq(ftm[f1[1]], 1) + "\t" + NumberToSeq(fcg[f1[1]], 1) + "\t" + flc[f1[1]] + "\n";
                    resultarea2 += flooppn[f1[2]] + "\t" + flooppr[f1[2]] + "\t" + flooppr[f1[2]].length + "\t" + NumberToSeq(flooptm[f1[2]], 1) + "\t" + NumberToSeq(floopcg[f1[2]], 1) + "\t" + flooplc[f1[2]] + "\n";
                    resultarea2 += f1pn[f1[3]] + "\t" + f1pr[f1[3]] + "\t" + f1pr[f1[3]].length + "\t" + NumberToSeq(f1tm[f1[3]], 1) + "\t" + NumberToSeq(f1cg[f1[3]], 1) + "\t" + f1lc[f1[3]] + "\n";
                    resultarea2 += "FIP:" + f1pn[f1[3]] + "--F2_" + fpn[f1[1]] + "\t" + ComplementDNA(f1pr[f1[3]]) + "TTTTTT" + fpr[f1[1]] + "\n";
                    resultarea2 += r1pn[r1[3]] + "\t" + r1pr[r1[3]] + "\t" + r1pr[r1[3]].length + "\t" + NumberToSeq(r1tm[r1[3]], 1) + "\t" + NumberToSeq(r1cg[r1[3]], 1) + "\t" + r1lc[r1[3]] + "\n";
                    resultarea2 += rlooppn[r1[2]] + "\t" + rlooppr[r1[2]] + "\t" + rlooppr[r1[2]].length + "\t" + NumberToSeq(rlooptm[r1[2]], 1) + "\t" + NumberToSeq(rloopcg[r1[2]], 1) + "\t" + rlooplc[r1[2]] + "\n";
                    resultarea2 += "B2_" + rpn[r1[1]] + "\t" + rpr[r1[1]] + "\t" + rpr[r1[1]].length + "\t" + NumberToSeq(rtm[r1[1]], 1) + "\t" + NumberToSeq(rcg[r1[1]], 1) + "\t" + rlc[r1[1]] + "\n";
                    resultarea2 += "B3_" + rpn[r1[0]] + "\t" + rpr[r1[0]] + "\t" + rpr[r1[0]].length + "\t" + NumberToSeq(rtm[r1[0]], 1) + "\t" + NumberToSeq(rcg[r1[0]], 1) + "\t" + rlc[r1[0]] + "\n";
                    resultarea2 += "BIP:" + r1pn[r1[3]] + "--B2_" + rpn[r1[1]] + "\t" + ComplementDNA(r1pr[r1[3]]) + "TTTTTT" + rpr[r1[1]] + "\n";
                    let ta = Tm77(sq.substring(f1[4], r1[5]), KMg_M, 0);
                    resultarea2 += "F2-B2:" + fpn[f1[1]] + "-" + rpn[r1[1]] + "\t\t\t\t\t\t" + r[2] + "/" + NumberToSeq(ta, 1) + "\n\n";
                }
            }
            else {
                for (let f = 0; f <= nfset; f++) {
                    let f1 = new Array(5).fill(0);
                    let v = 5 * f;
                    f1[0] = fset[v];     //F1
                    f1[1] = fset[v + 1]; //F2 
                    f1[2] = fset[v + 2]; //F3
                    f1[3] = fset[v + 3]; // 5'-F2
                    f1[4] = fset[v + 4]; // F1->3'

                    for (let b = 0; b <= nrset; b++) {
                        let r1 = new Array(5).fill(0);
                        let v = 5 * b;
                        r1[0] = rset[v];
                        r1[1] = rset[v + 1];
                        r1[2] = rset[v + 2];
                        r1[3] = rset[v + 3]; //<-B1
                        r1[4] = rset[v + 4]; //B2-5'

                        let z = r1[4] - f1[3];

                        if (f1[4] <= r1[3] && z <= maxf2b2) {
                            let a1 = [];
                            a1.push(fpr[f1[0]]);
                            a1.push(fpr[f1[1]]);
                            a1.push(f1pr[f1[2]]);
                            let zlc = flc[f1[0]] + flc[f1[1]] + f1lc[f1[2]];
                            let a2 = [];
                            a2.push(r1pr[r1[2]]);
                            a2.push(rpr[r1[1]]);
                            a2.push(rpr[r1[0]]);
                            zlc = zlc + rlc[r1[0]] + rlc[r1[1]] + r1lc[r1[2]];

                            if (ncol > -1) {
                                if (col[col.length - 2] > z || col[col.length - 1] < zlc) {
                                    if (DimerList2(a1, a2, min3) === -1) {
                                        ncol++;
                                        col.push(f);
                                        col.push(b);
                                        col.push(z);
                                        col.push(zlc);
                                    }
                                }
                            } else {
                                if (DimerList2(a1, a2, min3) === -1) {
                                    ncol++;
                                    col.push(f);
                                    col.push(b);
                                    col.push(z);
                                    col.push(zlc);
                                }
                            }
                        }
                    }
                }

                // printing collections            
                for (let f = 0; f <= ncol; f++) {
                    let v = 4 * f;
                    let r = new Array(4).fill(0);
                    r[0] = col[v];     // F set
                    r[1] = col[v + 1]; // R set
                    r[2] = col[v + 2]; // size F3-B3
                    r[3] = col[v + 3]; // total LC

                    let f1 = new Array(5).fill(0);
                    let vf = 5 * r[0];
                    f1[0] = fset[vf];
                    f1[1] = fset[vf + 1];
                    f1[2] = fset[vf + 2];
                    f1[3] = fset[vf + 3]; // 5-F2
                    f1[4] = fset[vf + 4]; // F1->       

                    let r1 = new Array(5).fill(0);
                    let vr = 5 * r[1];
                    r1[0] = rset[vr];
                    r1[1] = rset[vr + 1];
                    r1[2] = rset[vr + 2];
                    r1[3] = rset[vr + 3]; //<-B1
                    r1[4] = rset[vr + 4]; //B2-5

                    let a1 = [];
                    a1.push(fpr[f1[0]]);
                    a1.push(fpr[f1[1]]);
                    a1.push(f1pr[f1[2]]);
                    let a2 = [];
                    a2.push(r1pr[r1[2]]);
                    a2.push(rpr[r1[1]]);
                    a2.push(rpr[r1[0]]);

                    resultarea2 += "F3_" + fpn[f1[0]] + "\t" + fpr[f1[0]] + "\t" + fpr[f1[0]].length + "\t" + NumberToSeq(ftm[f1[0]], 1) + "\t" + NumberToSeq(fcg[f1[0]], 1) + "\t" + flc[f1[0]] + "\n";
                    resultarea2 += "F2_" + fpn[f1[1]] + "\t" + fpr[f1[1]] + "\t" + fpr[f1[1]].length + "\t" + NumberToSeq(ftm[f1[1]], 1) + "\t" + NumberToSeq(fcg[f1[1]], 1) + "\t" + flc[f1[1]] + "\n";
                    resultarea2 += f1pn[f1[2]] + "\t" + f1pr[f1[2]] + "\t" + f1pr[f1[2]].length + "\t" + NumberToSeq(f1tm[f1[2]], 1) + "\t" + NumberToSeq(f1cg[f1[2]], 1) + "\t" + f1lc[f1[2]] + "\n";
                    resultarea2 += "FIP:" + f1pn[f1[2]] + "--F2_" + fpn[f1[1]] + "\t" + ComplementDNA(f1pr[f1[2]]) + "TTTTTT" + fpr[f1[1]] + "\n";
                    resultarea2 += r1pn[r1[2]] + "\t" + r1pr[r1[2]] + "\t" + r1pr[r1[2]].length + "\t" + NumberToSeq(r1tm[r1[2]], 1) + "\t" + NumberToSeq(r1cg[r1[2]], 1) + "\t" + r1lc[r1[2]] + "\n";
                    resultarea2 += "B2_" + rpn[r1[1]] + "\t" + rpr[r1[1]] + "\t" + rpr[r1[1]].length + "\t" + NumberToSeq(rtm[r1[1]], 1) + "\t" + NumberToSeq(rcg[r1[1]], 1) + "\t" + rlc[r1[1]] + "\n";
                    resultarea2 += "B3_" + rpn[r1[0]] + "\t" + rpr[r1[0]] + "\t" + rpr[r1[0]].length + "\t" + NumberToSeq(rtm[r1[0]], 1) + "\t" + NumberToSeq(rcg[r1[0]], 1) + "\t" + rlc[r1[0]] + "\n";
                    resultarea2 += "BIP:" + r1pn[r1[2]] + "--B2_" + rpn[r1[1]] + "\t" + ComplementDNA(r1pr[r1[2]]) + "TTTTTT" + rpr[r1[1]] + "\n";
                    let ta = Tm77(sq.substring(f1[3], r1[4]), KMg_M, 0);
                    resultarea2 += "F2-B2:" + fpn[f1[1]] + "-" + rpn[r1[1]] + "\t\t\t\t\t\t" + r[2] + "/" + NumberToSeq(ta, 1) + "\n\n";
                }
            }
        }
    }
    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

function DimerList(s, min3) {
    let n = s.length;
    if (n < 2) { return -1; }
    for (let i = 0; i < n - 1; i++) {
        for (let j = i + 1; j < n; j++) {
            if (DimerLook2(s[i], s[j], min3) > 0) { return 0; }
        }
    }
    return -1;
}

function DimerList2(s1, s2, min3) {
    let n1 = s1.length;
    let n2 = s2.length;
    if (n1 < 1 || n2 < 1) { return -1; }
    for (let i = 0; i < n1; i++) {
        for (let j = 0; j < n2; j++) {
            if (DimerLook2(s1[i], s2[j], min3) > 0) { return 0; }
        }
    }
    return -1;
}

function PrimerDesignReverse(cs, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, rpr, rpn, rtm, rlc, rpx, rcg, nm, e5, e3) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 6;
    const s5 = strgen(e5);

    s3data = GeneratorVariants2(e3);
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
            if (e3 === "n" || s3.indexOf(s1) > -1) {
                for (let d = 0; d < 1 + maxlen - minlen; d++) {
                    x5 = x - d;
                    if (x5 >= x1) {
                        s1 = cs.substring(x5, x5 + minlen + d);
                        if (e5 === "n" || s5.indexOf(s1.substring(0, 1)) > -1) {
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

function PrimerDesignForward(sq, msk, minlc, minlen, maxlen, mintm, maxtm, x1, x2, primeroverlap, fpr, fpn, ftm, flc, fpx, fcg, nm, e5, e3) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 6;
    const s5 = strgen(e5);

    s3data = GeneratorVariants2(e3);
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
            if (e3 === "n" || s3.indexOf(s1) > -1) {

                for (let d = 0; d < 1 + maxlen - minlen; d++) {
                    x5 = x - d;
                    if (x5 >= x1) {
                        s1 = sq.substring(x5, x5 + minlen + d);
                        if (e5 === "n" || s5.indexOf(s1.substring(0, 1)) > -1) {
                            tm = Tm(s1, salt, Mg_M, p_mkM);
                        }
                        if (tm >= mintm || !ssrrepeat(s1)) {
                            break;
                        }
                    } else { break; }
                }

                if (tm >= mintm && tm <= maxtm) {
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

function Seq(str) {
    let b = [];
    let b1 = [];
    let bx = [];
    let seq = "";
    for (let i = 0; i < str.length; i++) {
        let chr = str.charAt(i);
        if (chr === 'a' || chr === 't' || chr === 'c' || chr === 'g' || chr === 'r' || chr === 'y' || chr === 'm' || chr === 'k' || chr === 'w' || chr === 'b' || chr === 'd' || chr === 'v' || chr === 'h' || chr === 's' || chr === 'n') {
            seq += chr;
        }
        if (chr === 'u') {
            seq += 't';
        }
        if (chr === 'i') {
            seq += 'g';
        }
        if (chr === '[') {
            b1.push(seq.length);
        }
        if (chr === ']') {
            b1.push(-seq.length);
        }
        if (chr === '/') {
            bx.push(seq.length);
        }
    }
    if (b1.length === 0) {
        b.push(0);
        b.push(seq.length);
        b.push(0);
        b.push(seq.length);
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
        if (b.length === 2) {
            b.push(b[0]);
            b.push(b[1]);
        }
    }
    return { seq, b, bx };
}