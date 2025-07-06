function ReadingSeq(s) {
    let k = -1;
    let j = 0;
    let z = 0;
    let x = 0;
    let n_seq = -1;
    const l = s.length;
    s = s.toLowerCase();

    let name_seq = [];
    let seqs = [];

    if (l < 2) { return { name_seq, seqs }; }

    if (s.indexOf(">", 0) > -1) {
        do {
            k = s.indexOf(">", k + 1);
            if (k > -1) {
                n_seq = n_seq + 1;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
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
        k = -1;
        do {
            k = s.indexOf("\t", k + 1);
            if (k > -1) {
                n_seq = n_seq + 1;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
        k = -1;
        j = 0;
        x = 0;
        n_seq = -1;
        do {
            k = s.indexOf("\t", k + 1);
            if (k > -1) {
                j = s.indexOf("\r", k + 3);
                z = s.indexOf("\n", k + 3);
                if (j === -1 && z !== -1) {
                    j = z;
                }
                if (j > -1) {
                    n_seq = n_seq + 1;
                    name_seq[n_seq] = s.substring(x, k).trim();
                    seqs[n_seq] = s.substring(k, j);
                    x = j + 1;
                    k = j;
                } else {
                    n_seq = n_seq + 1;
                    name_seq[n_seq] = s.substring(x, k).trim();
                    seqs[n_seq] = s.substring(k, l);
                    break;
                }
            }
        } while (k > -1);
    }
    // " " blank reading
    if (n_seq === -1 && s.indexOf(" ", 0) > -1) {
        k = -1;
        do {
            k = s.indexOf(" ", k + 1);
            if (k > -1) {
                n_seq = n_seq + 1;
            }
        } while (k > -1);
        name_seq = [];
        seqs = [];
        k = -1;
        j = 0;
        x = 0;
        n_seq = -1;
        do {
            k = s.indexOf(" ", k + 1);
            if (k > -1) {
                j = s.indexOf("\r", k + 3);
                z = s.indexOf("\n", k + 3);
                if (j === -1 && z !== -1)
                    j = z;
                if (j > -1) {
                    n_seq = n_seq + 1;
                    name_seq[n_seq] = s.substring(x, k).trim();
                    seqs[n_seq] = s.substring(k, j);
                    x = j + 1;
                    k = j;
                } else {
                    n_seq = n_seq + 1;
                    name_seq[n_seq] = s.substring(x, k).trim();
                    seqs[n_seq] = s.substring(k, l);
                    break;
                }
            }
        } while (k > -1);
    }
    if (n_seq === -1) {
        n_seq = 0;
        name_seq[n_seq] = "Seq1";
        seqs[n_seq] = s;
    }
    return { name_seq, seqs };
}