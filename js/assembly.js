function analysis() {
    let minlc = parseInt(document.getElementById('minlc').value);
    let mintm = parseInt(document.getElementById('mintm').value);
    let maxtm = parseInt(document.getElementById('maxtm').value);
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);

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

    const min3 = 8;
    const result = ["", ""];
    let resultarea1 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n";
    let resultarea2 = "Location(ID)\tSequence(5'-3')\tLength(nt)\tTm(°C)\tCG(%)\tLinguistic_Complexity(%)\n";

    let fpr = new Array(n_seq).fill("");
    let fpn = new Array(n_seq).fill("");
    let ftm = new Array(n_seq).fill(0.0);
    let flc = new Array(n_seq).fill(0);
    let fcg = new Array(n_seq).fill(0.0);

    let rpr = new Array(n_seq).fill("");
    let rpn = new Array(n_seq).fill("");
    let rtm = new Array(n_seq).fill(0.0);
    let rlc = new Array(n_seq).fill(0);
    let rcg = new Array(n_seq).fill(0.0);


    let fn = 0;
    let rn = 0;

    for (let n = 0; n < n_seq; n++) {
        let RestSeq = Seq(seqs[n]);
        const seq = RestSeq.returnString;
        const cs = ComplementDNA(seq);

        if (PrimerDesignForward(n, seq, name_seq[n], plist, minlc, minlen, maxlen, mintm, fpr, fpn, ftm, flc, fcg, "F_", ends3, "")) {
            resultarea1 += fpn[n] + "\t" + fpr[n] + "\t" + fpr[n].length + "\t" + ftm[n].toFixed(1) + "\t" + fcg[n].toFixed(1) + "\t" + flc[n] + "\n";
        }
        if (PrimerDesignReverse(n, cs, name_seq[n], plist, minlc, minlen, maxlen, mintm, rpr, rpn, rtm, rlc, rcg, "R_", ends3, "")) {
            resultarea1 += rpn[n] + "\t" + rpr[n] + "\t" + rpr[n].length + "\t" + rtm[n].toFixed(1) + "\t" + rcg[n].toFixed(1) + "\t" + rlc[n] + "\n";
        }
    }


    // combinations
    if (DimerLook2(rpr[0], fpr[n_seq - 1], min3) > -1) {
        rlc[0] = 0;
        flc[n_seq - 1] = 0;
    }
    resultarea2 += rpn[0] + "\t" + rpr[0] + "\t" + rpr[0].length + "\t" + rtm[0].toFixed(1) + "\t" + rcg[0].toFixed(1) + "\t" + rlc[0] + "\n";
    resultarea2 += fpn[n_seq - 1] + "\t" + fpr[n_seq - 1] + "\t" + fpr[n_seq - 1].length + "\t" + ftm[n_seq - 1].toFixed(1) + "\t" + fcg[n_seq - 1].toFixed(1) + "\t" + flc[n_seq - 1] + "\n";
    

    for (let n = 1; n < n_seq - 1; n++) {
        if (fpr[n].length > 0 && rpr[n].length > 0 && rpr[n - 1].length > 0 && fpr[n + 1].length > 0) {
            let ftail = rpr[n - 1];
            let f1 = ComplementDNA(ftail).toUpperCase() + fpr[n];

            let rtail = fpr[n + 1];
            let r1 = ComplementDNA(rtail).toUpperCase() + rpr[n];

            if (DimerLook2(f1, f1, min3) > -1) {
                flc[n] = 0;
            }
            if (DimerLook2(r1, r1, min3) > -1) {
                rlc[n] = 0;
            }
            if (DimerLook2(f1, r1, min3) > -1) {
                rlc[n] = 0;
                flc[n] = 0;
            }

            resultarea2 += fpn[n] + "\t" + f1 + "\t" + fpr[n].length + "\t" + ftm[n].toFixed(1) + "\t" + fcg[n].toFixed(1) + "\t" + flc[n] + "\n";
            resultarea2 += rpn[n] + "\t" + r1 + "\t" + rpr[n].length + "\t" + rtm[n].toFixed(1) + "\t" + rcg[n].toFixed(1) + "\t" + rlc[n] + "\n";
        }
        else {
            resultarea2 += "\nProblem with " + name_seq[n] + "\n";
        }
    }

    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

function PrimerDesignForward(n, sq, seqname, plist, minlc, minlen, maxlen, mintm, fpr, fpn, ftm, flc, fcg, nm, e3, tail) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 8;

    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;
    const nplist = plist.length;

    for (let x = minlen; x < minlen + maxlen; x++) {
        let s1 = sq.substring(0, x);
        let s2 = s1.substring(x - le3);
        if (s3.indexOf(s2) > -1) {
            const tm1 = Tm(s1, salt, Mg_M, p_mkM);
            let q1 = 1;
            if (tm1 >= mintm) {
                if (ssrrepeat(s1) === -1) {
                    let lc = LingComplexity2(s1);
                    if (lc >= minlc) {
                        s1 = tail + s1;
                        if (DimerLook2(s1, s1, min3) === -1) {
                            q1 = 0;
                            if (nplist > 0) {
                                for (let j = 0; j < nplist; j++) {
                                    if (DimerLook2(s1, plist[j], min3) > -1) {
                                        q1 = 1;
                                        break;
                                    }
                                }
                            }

                            if (q1 === 0) {
                                fpr[n] = s1;
                                ftm[n] = tm1;
                                flc[n] = lc;
                                fcg[n] = CG(s1);
                                fpn[n] = seqname + "_" + nm + "1-" + s1.length;
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}


function PrimerDesignReverse(n, sq, seqname, plist, minlc, minlen, maxlen, mintm, fpr, fpn, ftm, flc, fcg, nm, e3, tail) {
    const salt = 0.055;
    const Mg_M = 0.001;
    const p_mkM = 0.2;
    const min3 = 8;

    const s3data = GeneratorVariants2(e3);
    const s3 = s3data.b.join(' ');
    const le3 = s3data.z;
    const nplist = plist.length;
    const lseqs = sq.length;

    for (let x = minlen; x < minlen + maxlen; x++) {
        let s1 = sq.substring(0, x);
        let s2 = s1.substring(x - le3);
        if (s3.indexOf(s2) > -1) {
            const tm1 = Tm(s1, salt, Mg_M, p_mkM);
            let q1 = 1;
            if (tm1 >= mintm) {
                if (ssrrepeat(s1) === -1) {
                    let lc = LingComplexity2(s1);
                    if (lc >= minlc) {
                        s1 = tail + s1;
                        if (DimerLook2(s1, s1, min3) === -1) {
                            q1 = 0;
                            if (nplist > 0) {
                                for (let j = 0; j < nplist; j++) {
                                    if (DimerLook2(s1, plist[j], min3) > -1) {
                                        q1 = 1;
                                        break;
                                    }
                                }
                            }

                            if (q1 === 0) {
                                fpr[n] = s1;
                                ftm[n] = tm1;
                                flc[n] = lc;
                                fcg[n] = CG(s1);
                                fpn[n] = seqname + "_" + nm + lseqs + "-" + (lseqs - s1.length + 1);
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
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