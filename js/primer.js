function analysis(sequence) {
        let salt = 0.05;
        let Mg_M = 0.0;
        let p = 0.0000002;
        var Mg0 = document.getElementById('mg_concentration');
        Mg0 = NumToDouble(Mg0.value);
        var salt0 = document.getElementById('salt_concentration');
        salt0 = NumToDouble(salt0.value);
        var p0 = document.getElementById('primer_concentration');
        p0 = NumToDouble(p0.value);
        var min3 = document.getElementById('sensivity');
        min3 = NumToDouble(min3.value);


        if (salt0 > 0) {
                salt = salt0 / 1000;
        }
        if (Mg0 > 0) {
                Mg_M = Mg0 / 1000;
        }
        if (p0 > 0) {
                p = p0 / 1000000;
        }


        const s = sequence.toLowerCase().trim();
        const l = s.length;
        var k = -1;
        var j = 0;
        var z = 0;
        var x = 0;
        var n_seq = -1;
        var name_seq = {};
        var seqs = {};

        if (l < 2)
                return;
        if (s.indexOf(">", 0) > -1) {
                do {
                        k = s.indexOf(">", k + 1);
                        if (k > -1) {
                                n_seq = n_seq + 1;
                        }
                } while (k > -1);
                name_seq = {};
                seqs = {};
                k = -1;
                j = 0;
                n_seq = -1;
                do {
                        k = s.indexOf(">", k + 1);
                        if (k > -1) {
                                j = s.indexOf("\r", k + 1);
                                z = s.indexOf("\n", k + 1);
                                x = s.indexOf(" ", k + 1);
                                if (j === -1 && z !== -1) {
                                        j = z;
                                }
                                if (j === -1 & z === -1 && x > -1) {
                                        j = x;
                                }
                                if (j > -1) {
                                        n_seq = n_seq + 1;
                                        name_seq[n_seq] = s.substring(k + 1, j).trim();
                                        k = s.indexOf(">", j + 1);
                                        if (k !== -1) {
                                                seqs[n_seq] = s.substring(j + 1, k - 1);
                                                k = k - 1;
                                        } else {
                                                seqs[n_seq] = s.substring(j + 1, l);
                                        }
                                }
                        }
                } while (k > -1);
        } //(s.indexOf(">",0)>-1)


        // Tab reading
        if (n_seq === -1 && s.indexOf("\t", 0) > -1) {
                const a1 = s.split("\n");
                for (let i = 0; i < a1.length; i++) {
                        if (a1[i].length > 3) {
                                const a2 = a1[i].split(" ");
                                if (a2.length > 0) {
                                        if (a2.length === 1) {
                                                n_seq = n_seq + 1;
                                                name_seq[n_seq] = "SEQ" + (1 + n_seq);
                                                seqs[n_seq] = a2[0].trim();
                                        } else {
                                                n_seq = n_seq + 1;
                                                name_seq[n_seq] = a2[0].trim();
                                                seqs[n_seq] = a2[1].trim();
                                                if (seqs[n_seq].length === 0) { seqs[n_seq] = name_seq[n_seq]; name_seq[n_seq] = "SEQ" + (1 + n_seq); }
                                                for (let j = 2; j < a2.length; j++) {
                                                        if (a2[j].trim.length !== 0) {
                                                                n_seq = n_seq + 1;
                                                                name_seq[n_seq] = a2[0].trim();
                                                                seqs[n_seq] = a2[j].trim();
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        // " " blank reading
        if (n_seq === -1 && s.indexOf(" ", 0) > -1) {
                const a1 = s.split("\n");
                for (let i = 0; i < a1.length; i++) {
                        if (a1[i].length > 3) {
                                const a2 = a1[i].split(" ");
                                if (a2.length > 0) {
                                        if (a2.length === 1) {
                                                n_seq = n_seq + 1;
                                                name_seq[n_seq] = "SEQ" + (1 + n_seq);
                                                seqs[n_seq] = a2[0].trim();
                                        } else {
                                                n_seq = n_seq + 1;
                                                name_seq[n_seq] = a2[0].trim();
                                                seqs[n_seq] = a2[1].trim();
                                                if (seqs[n_seq].length === 0) { seqs[n_seq] = name_seq[n_seq]; name_seq[n_seq] = "SEQ" + (1 + n_seq); }
                                                for (let j = 2; j < a2.length; j++) {
                                                        if (a2[j].trim.length !== 0) {
                                                                n_seq = n_seq + 1;
                                                                name_seq[n_seq] = a2[0].trim();
                                                                seqs[n_seq] = a2[j].trim();
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        if (n_seq === -1) {
                n_seq = 0;
                name_seq[n_seq] = "SEQ1";
                seqs[n_seq] = s;
        }

        // single sequence
        number_of_sequences.innerHTML = "Number of sequences: " + (1 + n_seq);


        const result = ["", ""];
        var resultarea = "";

        var zname = 4;
        var zseq = 8;
        var Sequence_width_title = 0;
        for (var n = 0; n < n_seq + 1; n++) {
                seqs[n] = PrimerSeq(seqs[n]);
                var lseqs = seqs[n].length;
                var lname = name_seq[n].length;
                if (Sequence_width_title < seqs[n].length)
                        Sequence_width_title = seqs[n].length;
                if (zname < lname) {
                        zname = lname;
                }
                if (zseq < lseqs) {
                        zseq = lseqs;
                }
        }
        var a2 = "\t";
        if (zname > 4) {
                var value = [];
                for (var q = 0; q < zname - 4; q++) {
                        value.push(' ');
                }
                a2 = value.join('') + "\t";
        }

        var a3 = "\t";
        if (zseq > 8) {
                var value = [];
                for (var q = 0; q < zseq - 8; q++) {
                        value.push(' ');
                }
                a3 = value.join('') + "\t";
        }

        resultarea = "Name" + a2 + fitToWidth("Sequence", Sequence_width_title) + "\tTm°C\tCG%\tLC%\tnt\tA\tT\tC\tG\tExtinction coefficient(l/(mol·cm)\tMolecular weight(g/mol)\tnmol\tµg/OD260";

        for (var n = 0; n < n_seq + 1; n++) {
                var an = 0;
                var tn = 0;
                var cn = 0;
                var gn = 0;
                var in_var = 0;
                var un = 0;
                const sq = seqs[n];
                var lseq = seqs[n].length;
                for (var i = 0; i < lseq; i++) {
                        var chr = seqs[n].charAt(i);
                        if (chr == 'm') {
                                an = an + 0.5;
                                cn = cn + 0.5;
                        }
                        if (chr == 'r') {
                                an = an + 0.5;
                                gn = gn + 0.5;
                        }
                        if (chr == 'w') {
                                an = an + 0.5;
                                tn = tn + 0.5;
                        }
                        if (chr == 's') {
                                gn = gn + 0.5;
                                cn = cn + 0.5;
                        }
                        if (chr == 'y') {
                                tn = tn + 0.5;
                                cn = cn + 0.5;
                        }
                        if (chr == 'k') {
                                gn = gn + 0.5;
                                tn = tn + 0.5;
                        }
                        if (chr == 'v') {
                                an = an + 0.333;
                                cn = cn + 0.333;
                                gn = gn + 0.333;
                        }
                        if (chr == 'h') {
                                an = an + 0.333;
                                cn = cn + 0.333;
                                tn = tn + 0.333;
                        }
                        if (chr == 'b') {
                                tn = tn + 0.333;
                                cn = cn + 0.333;
                                gn = gn + 0.333;
                        }
                        if (chr == 'd') {
                                an = an + 0.333;
                                tn = tn + 0.333;
                                gn = gn + 0.333;
                        }
                        if (chr == 'n') {
                                an = an + 0.25;
                                cn = cn + 0.25;
                                gn = gn + 0.25;
                                tn = tn + 0.25;
                        }
                        if (chr == 'x') {
                                an = an + 0.25;
                                cn = cn + 0.25;
                                gn = gn + 0.25;
                                tn = tn + 0.25;
                        }
                        if (chr == 'a') {
                                an = an + 1;
                        }
                        if (chr == 't') {
                                tn = tn + 1;
                        }
                        if (chr == 'c') {
                                cn = cn + 1;
                        }
                        if (chr == 'g') {
                                gn = gn + 1;
                        }
                        if (chr == 'i') {
                                in_var = in_var + 1;
                        }
                        if (chr == 'u') {
                                tn = tn + 1;
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
                        var d1 = CG(seqs[n]);
                        var w = MW(an, tn, cn, gn, in_var, un);
                        var M = e260(seqs[n]);
                        var nmol = 1000000 / M;
                        var mg = (w * 1000) / M;
                        if (name_seq[n] == null || name_seq[n] == "" || name_seq[n] == " ") {
                                name_seq[n] = "seq_" + (n + 1);
                        }
                        var s2 = "\t";
                        if (zname - name_seq[n].length > 0) {
                                var value = [];
                                for (var q = 0; q < zname - name_seq[n].length; q++) {
                                        value.push(' ');
                                }
                                s2 = value.join('') + "\t";
                        }
                        var s3 = "\t";
                        if (zseq - seqs[n].length > 0) {
                                var value = [];
                                for (var q = 0; q < zseq - seqs[n].length; q++) {
                                        value.push(' ');
                                }
                                s3 = value.join('') + "\t";
                        }
                        var s4 = "                    ";
                        if (w > 9999.999) {
                                s4 = "             ";
                        }
                        if (w > 99999.999) {
                                s4 = "      ";
                        }

                        resultarea += "\n"
                                + name_seq[n] + s2
                                + fitToWidth(seqs[n], Sequence_width_title) + "\t"
                                + NumberToSeq(Tm(seqs[n], salt, Mg_M, p)) + "\t"
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
                resultarea = "                Self-Dimers:\n\n";
                for (var n = 0; n < n_seq + 1; n++) {
                        var sl = seqs[n].length;
                        sl = sl + sl + sl;
                        var x0 = {};
                        var x1 = {};
                        var x2 = {};
                        var x3 = {};
                        k = DimerLook(seqs[n], seqs[n], min3, x0, x1, x2, x3);
                        if (k > -1) {
                                if (k == 0) { resultarea += "1 dimer for: " + name_seq[n] + "\n"; }
                                else { resultarea += (k + 1) + " dimers for: " + name_seq[n] + "\n"; }
                                for (var w = 0; w < k + 1; w++) {
                                        //Tmelting2(s5, s3, pr_conc, kcl_mM) 
                                        let d = new Array(4).fill(0);
                                        d = Tmelting2(x1[w], x3[w], p, salt);
                                        resultarea += "Tm=" + d[0].toFixed(1) + "°C; dG=" + d[1].toFixed(1) + " kcal/mol \n";
                                        resultarea += x1[w] + "\n";
                                        resultarea += x2[w] + "\n";
                                        resultarea += x3[w] + "\n\n";
                                }
                        }
                }
        }
        if (n_seq > 0) {
                resultarea += "\n               Cross Primer Dimers:\n\n";
                for (var n = 0; n < n_seq; n++) {
                        for (var f = n + 1; f < n_seq + 1; f++) {
                                var sl1 = seqs[n].length;
                                var sl2 = seqs[f].length;
                                var sl = sl1 + sl1 + sl2;
                                if (sl1 > sl2) { sl = sl2 + sl2 + sl1; }

                                var x0 = {};
                                var x1 = {};
                                var x2 = {};
                                var x3 = {};
                                k = DimerLook(seqs[n], seqs[f], min3, x0, x1, x2, x3);
                                if (k > -1) {
                                        resultarea += name_seq[n] + " with " + name_seq[f] + "\n";
                                        resultarea += name_seq[n] + "\n";
                                        for (var w = 0; w < k + 1; w++) {
                                                //Tmelting2(s5, s3, pr_conc, kcl_mM) 
                                                let d = new Array(4).fill(0);
                                                d = Tmelting2(x1[w], x3[w], p, salt);
                                                resultarea += "Tm=" + d[0].toFixed(1) + "°C; dG=" + d[1].toFixed(1) + " kcal/mol \n";
                                                resultarea += x1[w] + "\n";
                                                resultarea += x2[w] + "\n";
                                                resultarea += x3[w] + "\n\n";

                                        }
                                }
                        }
                }
        }

        result[1] = resultarea;
        return result;
}