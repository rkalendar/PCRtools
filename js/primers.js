function analysis(sequence) {
    let salt_M = 0.055;
    let Mg_M = 0.001;
    let pr_mkM = 0.2;
    let Mg0 = NumToDouble(document.getElementById('mg_concentration').value);
    let salt0 = NumToDouble(document.getElementById('salt_concentration').value);
    let p0 = NumToDouble(document.getElementById('primer_concentration').value); //parseInt(document.getElementById('primer_concentration').value); // 
    let min3 = parseInt(document.getElementById('sensivity').value);

    const selectBox = document.getElementById("selectBox");
    // const selectedOption = selectBox.options[selectBox.selectedIndex];
    // const selectedValue = selectedOption.value;

    if (salt0 > 0) {
        salt_M = salt0 / 1000;
    }
    if (Mg0 > 0) {
        Mg_M = Mg0 / 1000;
    }
    if (p0 > 0) {
        pr_mkM = p0;
    }
    if (min3 < 3) {
        min3 = 3;
    }

    ReadResult = ReadingSeq(sequence.toLowerCase().trim());
    let name_seq = ReadResult.name_seq;
    let seqs = ReadResult.seqs;
    let n_seq = seqs.length;

    const result = ["", ""];
    let resultarea = "";
    let zname = 4;
    let zseq = 8;
    let Sequence_width_title = 0;
    for (let n = 0; n < n_seq; n++) {
        seqs[n] = PrimerSeq(seqs[n]);
        let lseqs = seqs[n].length;
        let lname = name_seq[n].length;
        if (Sequence_width_title < seqs[n].length)
            Sequence_width_title = seqs[n].length;
        if (zname < lname) {
            zname = lname;
        }
        if (zseq < lseqs) {
            zseq = lseqs;
        }
    }
    let a2 = "\t";
    if (zname > 4) {
        let value = new Array(zname - 4).fill(' ');
        a2 = value.join('') + "\t";
    }
    resultarea = "Name" + a2 + fitToWidth("Sequence", Sequence_width_title) + "\tTm°C\tCG%\tLC%\tnt\tA\tT\tC\tG\tExtinction coefficient(l/(mol·cm)\tMolecular weight(g/mol)\tnmol\tµg/OD260";
    console.time("Cross Primer Dimers: DimerLookX");

    for (let n = 0; n < n_seq; n++) {
        let an = 0, tn = 0, cn = 0, gn = 0, inosine = 0, un = 0;
        let lseq = seqs[n].length;
        for (let i = 0; i < lseq; i++) {
            let chr = seqs[n].charAt(i);
            if (chr === 'm') {
                an += 0.5;
                cn += 0.5;
            }
            if (chr === 'r') {
                an += 0.5;
                gn += 0.5;
            }
            if (chr === 'w') {
                an += 0.5;
                tn += 0.5;
            }
            if (chr === 's') {
                gn += 0.5;
                cn += 0.5;
            }
            if (chr === 'y') {
                tn += 0.5;
                cn += 0.5;
            }
            if (chr === 'k') {
                gn += 0.5;
                tn += 0.5;
            }
            if (chr === 'v') {
                an += 0.333;
                cn += 0.333;
                gn += 0.333;
            }
            if (chr === 'h') {
                an += 0.333;
                cn += 0.333;
                tn += 0.333;
            }
            if (chr === 'b') {
                tn += 0.333;
                cn += 0.333;
                gn += 0.333;
            }
            if (chr === 'd') {
                an += 0.333;
                tn += 0.333;
                gn += 0.333;
            }
            if (chr === 'n') {
                an += 0.25;
                cn += 0.25;
                gn += 0.25;
                tn += 0.25;
            }
            if (chr === 'x') {
                an += 0.25;
                cn += 0.25;
                gn += 0.25;
                tn += 0.25;
            }
            if (chr === 'a') {
                an++;
            }
            if (chr === 't') {
                tn++;
            }
            if (chr === 'c') {
                cn++;
            }
            if (chr === 'g') {
                gn++;
            }
            if (chr === 'i') {
                inosine++;
            }
            if (chr === 'u') {
                tn++;
            }
            if (chr == 'e') {
                an = an + 1;
                chr = 'a';
            }
            if (chr == 'f') {
                cn = cn + 1;
                chr = 'c';
            }
            if (chr == 'j') {
                gn = gn + 1;
                chr = 'g';
            }
            if (chr == 'l') {
                tn = tn + 1;
                chr = 't';
            }
        }
        if (lseq > 1) {
            let d1 = CG(seqs[n]);
            let w = MW(an, tn, cn, gn, inosine, un);
            let M = e260(seqs[n]);
            let nmol = 1000000 / M;
            let mg = (w * 1000) / M;
            if (name_seq[n] === null || name_seq[n] === "" || name_seq[n] === " ") {
                name_seq[n] = "seq_" + (n + 1);
            }
            let s2 = "\t";
            if (zname - name_seq[n].length > 0) {
                let value = new Array(zname - name_seq[n].length).fill(' ');
                s2 = value.join('') + "\t";
            }

            resultarea += "\n"
                + name_seq[n] + s2
                + fitToWidth(seqs[n], Sequence_width_title) + "\t"
                + NumberToSeq(Tm(seqs[n], salt_M, Mg_M, pr_mkM)) + "\t"
                + NumberToSeq(d1) + "\t"
                + LingComplexity(seqs[n]) + "\t"
                + lseq + "\t"
                + NumberToSeq(an) + "\t"
                + NumberToSeq(tn) + "\t"
                + NumberToSeq(cn) + "\t"
                + NumberToSeq(gn) + "\t"
                + fitToWidth(NumberToSeq(M), "Extinction coefficient(l/(molÂ·cm)") + "\t"
                + fitToWidth(NumberToSeq(w), "Molecular weight(g/mol)") + "\t"
                + fitToWidth(NumberToSeq(nmol), "nmol") + "\t"
                + NumberToSeq(mg);
        }
    }
    result[0] += resultarea;


    // Dimer showing
    if (n_seq > -1) {
        resultarea = "                Self-Dimers:\n";
        for (let n = 0; n < n_seq; n++) {
            let sl = seqs[n].length;
            sl = sl + sl + sl;
            let x0 = [];
            let x1 = [];
            let x2 = [];
            let x3 = [];
            k = DimerLookX(seqs[n], seqs[n], min3, min3, x0, x1, x2, x3); //   k = DimersLook(seqs[n], seqs[n], 3, 3, x0, x1, x2, x3);
            if (k > -1) {
                if (k == 0) { resultarea += "\nDimer for: " + name_seq[n]; }
                else { resultarea += "\n" + (k + 1) + " dimers for: " + name_seq[n]; }
                for (let w = 0; w < k + 1; w++) {
                    let d = new Array(4).fill(0);
                    d = Tmelting2(x1[w], x3[w], pr_mkM, salt_M, Mg_M);
                    if (x0[w] === 1) {
                        resultarea += "\n3'dimer: Tm=" + d[0].toFixed(1) + "°C; dG=" + d[1].toFixed(1) + " kcal/mol\n";
                    } else {
                        resultarea += "\nTm=" + d[0].toFixed(1) + "°C; dG=" + d[1].toFixed(1) + " kcal/mol\n";
                    }
                    resultarea += x1[w] + "\n";
                    resultarea += x2[w] + "\n";
                    resultarea += x3[w] + "\n";
                }
            }
        }
    }
    if (n_seq > 0) {
        resultarea += "\n               Cross Primer Dimers:\n";
        for (let n = 0; n < n_seq - 1; n++) {
            for (let f = n + 1; f < n_seq; f++) {
                let sl1 = seqs[n].length;
                let sl2 = seqs[f].length;
                let sl = sl1 + sl1 + sl2;
                if (sl1 > sl2) { sl = sl2 + sl2 + sl1; }
                let x0 = [];
                let x1 = [];
                let x2 = [];
                let x3 = [];
                k = DimerLookX(seqs[n], seqs[f], min3, min3, x0, x1, x2, x3);    //  k = DimersLook(seqs[n], seqs[f], 4, 5, x0, x1, x2, x3);                             
                if (k > -1) {
                    resultarea += "\n" + name_seq[n] + " with " + name_seq[f];
                    for (let w = 0; w < k + 1; w++) {
                        let d = new Array(4).fill(0);
                        d = Tmelting2(x1[w], x3[w], pr_mkM, salt_M, Mg_M);
                        if (x0[w] === 1) {
                            resultarea += "\n3'dimer: Tm=" + d[0].toFixed(1) + "°C dG=" + d[1].toFixed(1) + " kcal/mol\n";
                        } else {
                            resultarea += "\nTm=" + d[0].toFixed(1) + "°C dG=" + d[1].toFixed(1) + " kcal/mol\n";
                        }
                        resultarea += x1[w] + "\n";
                        resultarea += x2[w] + "\n";
                        resultarea += x3[w] + "\n";

                    }
                }
            }
        }
    }
    console.timeEnd("Cross Primer Dimers: DimerLookX");
    result[1] = resultarea;
    return result;
}