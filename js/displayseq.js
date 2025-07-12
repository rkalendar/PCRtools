function Display(s) {
    let k = -1;
    let j = 0;
    let z = 0;
    let x = 0;
    let d = "";
    let n = -1;
    const l = s.length;
    s = s.toLowerCase();
    if (l < 2)
        return;
    let name_seq = [];
    let seqs = [];
    if (s.indexOf(">", 0) > -1) {
        do {
            k = s.indexOf(">", k + 1);
            if (k > -1) {
                n++;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
        k = -1;
        j = 0;
        n = -1;
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
                    n++;
                    name_seq[n] = s.substring(k + 1, j).trim();
                    k = s.indexOf(">", j + 1);
                    if (k !== -1) {
                        seqs[n] = s.substring(j + 1, k - 1);
                        k = k - 1;
                    } else {
                        seqs[n] = s.substring(j + 1, l);
                    }
                }
            }
        } while (k > -1);
    } //(s.indexOf(">",0)>-1)
    // Tab reading
    if (n === -1 && s.indexOf("\t", 0) > -1) {
        k = -1;
        do {
            k = s.indexOf("\t", k + 1);
            if (k > -1) {
                n++;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
        k = -1;
        j = 0;
        x = 0;
        n = -1;
        do {
            k = s.indexOf("\t", k + 1);
            if (k > -1) {
                j = s.indexOf("\r", k + 3);
                z = s.indexOf("\n", k + 3);
                if (j === -1 && z !== -1) {
                    j = z;
                }
                if (j > -1) {
                    n++;
                    name_seq[n] = s.substring(x, k).trim();
                    seqs[n] = s.substring(k, j);
                    x = j + 1;
                    k = j;
                } else {
                    n++;
                    name_seq[n] = s.substring(x, k).trim();
                    seqs[n] = s.substring(k, l);
                    break;
                }
            }
        } while (k > -1);
    }
    // " " blank reading
    if (n === -1 && s.indexOf(" ", 0) > -1) {
        k = -1;
        do {
            k = s.indexOf(" ", k + 1);
            if (k > -1) {
                n++;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
        k = -1;
        j = 0;
        x = 0;
        n = -1;
        do {
            k = s.indexOf(" ", k + 1);
            if (k > -1) {
                j = s.indexOf("\r", k + 3);
                z = s.indexOf("\n", k + 3);
                if (j === -1 && z !== -1)
                    j = z;
                if (j > -1) {
                    n++;
                    name_seq[n] = s.substring(x, k).trim();
                    seqs[n] = s.substring(k, j);
                    x = j + 1;
                    k = j;
                } else {
                    n++;
                    name_seq[n] = s.substring(x, k).trim();
                    seqs[n] = s.substring(k, l);
                    break;
                }
            }
        } while (k > -1);
    }

    if (n === -1) {
        name_seq.push("Seq1");
        seqs.push(s);
    }

    n = 0;
    if (name_seq.length > 0) {
        let s = "";
        for (let j = 0; j < name_seq.length; j++) {
            let z = DNA(seqs[j]);
            if (z.length > 0) {
                s += z;
                n++;
            }
        }
        d = " " + n + " : " + s.length + " nt (CG=" + CG(s).toFixed(1) + "%)";
    }
    return d;
}
