function analysis() {
    const salt_M = 0.055;
    const Mg_M = 0.001;
    const pr_mkM = 0.2;
    const kmer = 12;
    let minlen = parseInt(document.getElementById('minlen').value);
    let maxlen = parseInt(document.getElementById('maxlen').value);
    let minerror = parseInt(document.getElementById('minerror').value);
    let seqshow = document.getElementById('seqshow').checked;
    let allmatches = document.getElementById('allmatches').checked;
    let onlyamplicons = document.getElementById('onlyamplicons').checked;
    let ctconvert = document.getElementById('ctconvert').checked;

    if (minlen < 12) {
        minlen = 12;
    }
    if (maxlen < 12) {
        maxlen = 12;
    }
    if (minlen > maxlen) {
        maxlen = minlen;
    }

    ReadPrimers = ReadingSeq(document.getElementById('inputPrimers').value);
    let nm = ReadPrimers.name_seq;
    let pr = ReadPrimers.seqs;
    let primers = new Array(pr.length + 1).fill("");
    let name_primers = new Array(pr.length + 1).fill("");
    let n_primers = 1;

    for (let n = 0; n < pr.length; n++) {
        const s = DNASeq(pr[n]);
        if (s.length >= kmer) {
            primers[n_primers] = s;
            name_primers[n_primers] = nm[n];
            n_primers++;
        }
    }

    if (n_primers < 1) { return; }
    let tmspri = new Array(1 + n_primers).fill(0);

    const d = new Map();
    for (let n = 1; n < n_primers; n++) {
        const l = primers[n].length;
        tmspri[n] = Tm(primers[n], salt_M, Mg_M, pr_mkM);
        for (let y = 0; y < l - kmer + 1; y++) {
            const variants = GeneratorVariants(primers[n].substring(y, y + kmer));
            for (const variant of variants) {
                d.set(variant, (d.get(variant) || 0) + 1);
            }
        }
    }

    const p = new Map();
    for (let n = 1; n < n_primers; n++) {
        const l = primers[n].length;
        for (let y = 0; y < l - kmer + 1; y++) {
            const variants = GeneratorVariants(primers[n].substring(y, y + kmer));
            for (const variant of variants) {
                if (d.has(variant)) {
                    const k = d.get(variant);
                    const m = new Array(k + k + 1).fill(0);
                    m[0] = 1;
                    m[1] = n;
                    m[2] = y;
                    p.set(variant, m);
                    d.delete(variant);
                }
                else {
                    if (p.has(variant)) {
                        let m = p.get(variant);
                        m[0]++;
                        let z = m[0] * 2 - 1;
                        m[z] = n;
                        m[z + 1] = y;
                    }
                }
            }
        }
    }

    ReadResult = ReadingSeq(document.getElementById('inputText').value);
    let name_seq = ReadResult.name_seq;
    let seqs = ReadResult.seqs;
    let n_seq = seqs.length;
    let resultarea1 = "";
    //  ctconvert = true;
    for (let n = 0; n < n_seq; n++) {
        let startTime = performance.now();

        let seq = DNASeq(seqs[n]);
        let fseq = seq;
        let rseq = "";

        if (ctconvert) {
            fseq = DNActBisulfiteForward(seq);
            rseq = DNActBisulfiteReverse(seq);
        }
        else {
            rseq = ComplementDNA(seq);
        }

        resultarea1 += ">" + name_seq[n] + " : " + fseq.length + "nt\n\n";
        let fprimers = Searching(fseq, kmer, minerror, p, primers, name_primers, n_primers);

        // End time
        let endTime = performance.now();
        let runtime = endTime - startTime;
        if (runtime > 999) {
            runtime = runtime / 1000;
            resultarea1 += "The code took " + runtime.toFixed(0) + " seconds to run.\n\n";
        } else {
            resultarea1 += "The code took " + runtime.toFixed(0) + " milliseconds to run.\n\n";
        }
        if (fprimers.size > 0 && allmatches && !onlyamplicons) {
            resultarea1 += Print(0, fseq, fprimers, primers, name_primers, n_primers);
        }

        let rprimers = Searching(rseq, kmer, minerror, p, primers, name_primers, n_primers);
        if (rprimers.size > 0 && allmatches && !onlyamplicons) {
            resultarea1 += Print(2, rseq, rprimers, primers, name_primers, n_primers);
        }

        if (fprimers.size > 0 && rprimers.size > 0) {
            resultarea1 += "**** PCR primer pairs combinations:\n\n"
            resultarea1 += Combine(seq, fprimers, rprimers, primers, name_primers, tmspri, minlen, maxlen, seqshow, onlyamplicons);
        }

    }
    return resultarea1;
}

function Searching(seq, kmer, minerror, p, primers) {
    const fpr = new Map();
    const lseqs = seq.length;

    for (let y = 0; y < lseqs - kmer + 1; y++) {
        let s = seq.substring(y, y + kmer);

        if (p.has(s)) {
            let k = p.get(s);
            for (let x = 1; x < k.length; x += 2) {
                let n = k[x];
                let x1 = k[x + 1];
                let y1 = y - x1;    // location 5-end primer at target              
                let l = 0;
                let t = y + kmer;
                let e = 0;
                const lpr = primers[n].length;
                const aprimer = AntiSenseDNA(primers[n]);
                const v = lpr < 22 ? lpr - kmer - 2 - minerror : 20 - kmer - minerror; // test

                if (x1 + kmer <= lpr) {
                    for (let j = x1 + kmer; j < lpr; j++, t++) {
                        if (t > lseqs) { break; }
                        if (aDNA[aprimer.charCodeAt(j)][seq.charCodeAt(t)] > 24) {
                            l++;
                        } else {
                            e++;
                            if (e > minerror) { break; }
                        }
                    }
                }

                if (x1 > 0 && e <= minerror) {
                    t = y;
                    for (let j = x1 - 1; j > 0; j--) {
                        if (aDNA[aprimer.charCodeAt(j)][seq.charCodeAt(--t)] > 24) {
                            l++;
                        }
                    }
                }
                if (l > v && e <= minerror) {
                    if (fpr.has(y1)) {
                        let k = fpr.get(y1);
                        if (!k.includes(n)) {
                            k.push(n);
                            fpr.set(y1, k);
                        }
                    } else {
                        let k = [];
                        k.push(n);
                        fpr.set(y1, k);
                    }
                }
            }
        }
    }
    return fpr;
}

function Print(direction, seq, prim, primers, name_primers, n_primers) {
    let resultarea = "";
    const sign = " ";
    const salt_M = 0.055;
    const Mg_M = 0.001;
    const pr_mkM = 0.2;
    const l = seq.length;

    if (prim.size > 0) {
        for (let h = 1; h < n_primers; h++) {
            const s1 = "5-" + primers[h] + "->";
            const lpr = primers[h].length;
            let s2 = "";
            let s2a = "";
            let sb = "";
            let u = 0;

            for (const [key, value] of prim) {
                let y = key;   // primer location on target
                let k = value; // index primers 1.....n_primers

                for (let i = 0; i < k.length; i++) {
                    if (k[i] === h) {
                        if (u === 0) {
                            u = 1;
                            resultarea += ">" + name_primers[h] + "\n";
                            resultarea += primers[h] + "\n\n";
                        }
                        const y0 = y;
                        y = y - 2;
                        if (y > -1) {
                            s2 = seq.substring(y, y + lpr + 4);
                        }
                        else {
                            s2 = sign.repeat(- y) + seq.substring(0, lpr + 2);
                        }
                        s2a = AntiSenseDNA(s2);
                        const d = Tmelting2(s1, s2a, pr_mkM, salt_M, Mg_M);
                        const tm = d[0];
                        let f = [];
                        let sim = 0;
                        for (let v = 2; v < s1.length - 1; v++) {
                            const value = aDNA[s1.charCodeAt(v)][s2a.charCodeAt(v)];
                            if (value > 65) {
                                f.push("|");
                                sim++;
                            } else if (value > 24) {
                                f.push(":");
                                sim++;
                            } else {
                                f.push(" ");
                            }
                        }
                        sb = "  " + f.join("");
                        sim = (sim * 100) / lpr;

                        if (direction === 2) {
                            resultarea += (l - y0 - lpr + 1) + "<-" + (l - y0) + " Tm=" + tm.toFixed(1) + " " + sim.toFixed(1) + "%\n";
                            resultarea += "<-" + ReverseDNA(primers[h]) + "-5" + "\n";
                            resultarea += " " + ReverseDNA(sb) + "\n";
                            resultarea += ComplementDNA(s2.trim()) + "\n\n";
                        } else {
                            resultarea += (1 + y0) + "->" + (y0 + lpr) + " Tm=" + tm.toFixed(1) + " " + sim.toFixed(1) + "%\n";
                            resultarea += s1 + "\n";
                            resultarea += sb + "\n";
                            resultarea += s2 + "\n\n";
                        }
                    }
                }
            }
        }
    }
    return resultarea;
}

function Combine(seq, fpr, rpr, primers, name_primers, tmspri, minlen, maxlen, seqshow, onlyamplicons) {
    let resultarea = "";
    const l = seq.length;
    const rpr2 = new Map();
    for (const [key, value] of rpr) {
        let y = key;   // primer location on target
        let k = value; // index primers 1.....n_primers
        y = 2 + l - y; // new 5-end
        rpr2.set(y, k);
    }
    for (const [key, value] of fpr) {
        let x1 = key;   // forward primer location on target
        let f = value;  // index forward primers 1.....n_primers
        for (const [key, value] of rpr2) {
            let x2 = key;   // reverse primer location on target
            let r = value;  // index reverse primers 1.....n_primers
            let dx = x2 - x1;
            if (dx > minlen && dx < maxlen) {
                for (let i = 0; i < f.length; i++) {
                    if (!onlyamplicons) {
                        resultarea += ">" + name_primers[f[i]] + " " + (1 + x1) + "->" + (x1 + primers[f[i]].length) + "\n";
                        resultarea += primers[f[i]] + "\n";
                    }
                    for (let j = 0; j < r.length; j++) {
                        if (!onlyamplicons) {
                            const Ta = Tannealing(tmspri[f[i]], tmspri[r[j]], dx);
                            resultarea += ">" + name_primers[r[j]] + " " + (x2 - primers[r[j]].length - 1) + "<-" + (x2 - 2) + "\n";
                            resultarea += primers[r[j]] + "\n";
                            resultarea += "Amplicon size: " + (dx - 2) + "bp " + "Ta=" + Ta.toFixed(0) + "°C\n\n";
                            if (seqshow) {
                                resultarea += ">" + (1 + x1) + "-" + (x2 - 2) + "\n";
                                resultarea += seq.substring(x1, x2) + "\n\n";
                            }
                        } else {
                            resultarea += (dx - 2) + "\n";
                        }
                    }
                }
            }
        }
    }
    return resultarea;
}

function Tannealing(Tm1, Tm2, ln) {
    let t = (Tm1 < Tm2) ? Tm1 : Tm2;
    t += Math.log(ln);
    return t;
}