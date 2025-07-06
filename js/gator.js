function analysis() {
    let probe = NumToDouble(document.getElementById('probe').value);
    let distance = NumToDouble(document.getElementById('distance').value);
    if (probe < 35) {
        probe = 35;
    }
    if (probe > 40) {
        probe = 40;
    }
    if (distance < 5) {
        distance = 5;
    }
    let sequence = document.getElementById('inputText').value;
    const s = sequence.toLowerCase().trim();
    const l = s.length;
    const min3 = 5;
    let k = -1;
    let n_seq = -1;
    let j = 0, z = 0, x = 0;

    if (l < 2)
        return;
    let name_seq = [];
    let seqs = [];
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
        name_seq[n_seq] = "SEQ1";
        seqs[n_seq] = s;

    }

    const result = ["", ""];
    let resultarea = "";
    for (let n = 0; n < n_seq + 1; n++) {
        resultarea = "Location(ID)\tLinguistic_Complexity(%)\tSequence(5'-3')\n" + name_seq[n].toUpperCase() + ":\n";
        var pr1 = [];
        var pn1 = [];
        var lc1 = [];
        var lc2 = [];
        var p1x1 = [];
        var p1x2 = [];
        var pr2 = [];
        var pn2 = [];
        var p2x1 = [];
        var p2x2 = [];


        let b = [];
        let bx = [];
        const sq = Seq(seqs[n], b, bx);
        const lseqs = sq.length;

        let x1 = b[0];
        let x2 = b[1] - probe + 1;
        let x3 = b[2];
        let x4 = b[3] - probe + 1;

        let msk = RepeatMask(sq);
        if (bx.length > 1) {
            for (let y = 0; y < bx.length; y += 2) {
                msk = msk.fill(1, bx[y], bx[y + 1]);
            }
        }

        let q = 0;
        for (let y = x1; y < x1 + probe - 1; y++) {
            if (msk[y] > 0) { q++; }
        }
        for (let x = x1; x < x2; x++) {
            if (msk[x + probe - 1] > 0) { q++; }
            if (q === 0) {
                const s1 = sq.substring(x, x + probe);
                if (DimerLook(s1, s1, min3) === -1) {
                    if (repeat(s1) === -1) {
                        lc1.push(LingComplexity(s1));
                        pr1.push(s1);
                        p1x1.push(x + 1);
                        p1x2.push(x + probe + 1);
                        pn1.push("F_" + (x + 1) + "-" + (x + probe + 1));
                    }
                }
            }
            if (msk[x] > 0) { q--; }
        }

        q = 0;
        for (let y = x3; y < x3 + probe - 1; y++) {
            if (msk[y] > 0) { q++; }
        }
        for (let x = x3; x < x4; x++) {
            if (msk[x + probe - 1] > 0) { q++; }
            if (q === 0) {
                const s1 = sq.substring(x, x + probe);
                if (DimerLook(s1, s1, min3) === -1) {
                    if (repeat(s1) === -1) {
                        lc2.push(LingComplexity(s1));
                        pr2.push(s1);
                        p2x1.push(x + 1);
                        p2x2.push(x + probe + 1);
                        pn2.push("F_" + (x + 1) + "-" + (x + probe + 1));
                    }
                }
            }
            if (msk[x] > 0) { q--; }
        }


        for (let i = 0; i < pr1.length; i++) {
            for (let j = 0; j < pr2.length; j++) {
                if (x1 === x3 & x2 === x4) {
                    if (p2x1[j] + 1 > p1x2[i] + distance) {
                        if (DimerLook(pr1[i], pr2[j], min3) === -1) {
                            resultarea += pn1[i] + "\t" + lc1[i] + "\t" + pr1[i] + "\n";
                            resultarea += pn2[j] + "\t" + lc2[j] + "\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + pr2[j] + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" + "\n\n";
                            break;

                        }
                    }
                } else {
                    if (DimerLook(pr1[i], pr2[j], min3) === -1) {
                        resultarea += pn1[i] + "\t" + lc2[i] + "\t" + pr1[i] + "\n";
                        resultarea += pn2[j] + "\t" + lc2[j] + "\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + pr2[j] + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" + "\n\n";
                        break;
                    }
                }
            }
        }
        result[0] += resultarea;

        // reverse direction    
        resultarea = "Location(ID)\tLinguistic_Complexity(%)\tSequence(5'-3')\n" + name_seq[n].toUpperCase() + ":\n";
        var pr1 = [];
        var pn1 = [];
        var lc1 = [];
        var lc2 = [];
        var p1x1 = [];
        var p1x2 = [];
        var pr2 = [];
        var pn2 = [];
        var p2x1 = [];
        var p2x2 = [];

        const cs = ComplementDNA(sq);
        msk = msk.reverse();

        x1 = lseqs - b[3];
        x2 = lseqs - b[2] + 1 - probe;

        q = 0;
        for (let y = x1; y < x1 + probe - 1; y++) {
            if (msk[y] > 0) { q++; }
        }
        for (let x = x1; x < x2; x++) {
            if (msk[x + probe - 1] > 0) { q++; }
            if (q === 0) {
                const s1 = cs.substring(x, x + probe);
                if (DimerLook(s1, s1, min3) === -1) {
                    if (repeat(s1) === -1) {
                        lc1.push(LingComplexity(s1));
                        pr1.push(s1);
                        p1x1.push(x + 1);
                        p1x2.push(x + probe + 1);
                        pn1.push("R_" + (lseqs - x) + "-" + (lseqs - (x + probe) + 1));
                    }
                }
            }
            if (msk[x] > 0) { q--; }
        }

        x3 = lseqs - b[1];
        x4 = lseqs - b[0] + 1 - probe;
        q = 0;
        for (let y = x3; y < x3 + probe - 1; y++) {
            if (msk[y] > 0) { q++; }
        }
        for (let x = x3; x < x4; x++) {
            if (msk[x + probe - 1] > 0) { q++; }
            if (q === 0) {
                const s1 = cs.substring(x, x + probe);
                if (DimerLook(s1, s1, min3) === -1) {
                    if (repeat(s1) === -1) {
                        lc2.push(LingComplexity(s1));
                        pr2.push(s1);
                        p2x1.push(x + 1);
                        p2x2.push(x + probe + 1);
                        pn2.push("R_" + (lseqs - x) + "-" + (lseqs - (x + probe) + 1));
                    }
                }
            }
            if (msk[x] > 0) { q--; }
        }

        for (let i = 0; i < pr1.length; i++) {
            for (let j = 0; j < pr2.length; j++) {
                if (x1 === x3 & x2 === x4) {
                    if (p2x1[j] + 1 > p1x2[i] + distance) {
                        if (DimerLook(pr1[i], pr2[j], min3) === -1) {
                            resultarea += pn2[j] + "\t" + lc2[j] + "\t" + pr2[j] + "\n";
                            resultarea += pn1[i] + "\t" + lc1[i] + "\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + pr1[i] + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" + "\n\n";
                            break;
                        }
                    }
                } else {
                    if (DimerLook(pr1[i], pr2[j], min3) === -1) {
                        resultarea += pn2[j] + "\t" + lc2[j] + "\t" + pr2[j] + "\n";
                        resultarea += pn1[i] + "\t" + lc1[i] + "\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + pr1[i] + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t" + "\n\n";
                        break;
                    }

                }
            }
        }
        result[1] += resultarea;
    }
    return result;
}

function NumToDouble(val) {
    val = val.replace(',', '.');
    val = val.replace(/[^0-9.-]/g, "");
    let out = parseFloat(val);
    if (isNaN(out))
        return 0;
    else
        return out;
}

function Seq(str, b, bx) {
    let returnString = "";
    let b1 = [];
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
        b.push(0);
        b.push(returnString.length);
        b.push(0);
        b.push(returnString.length);
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
    return returnString;
}

function ReverseDNA(source) {
    return source.split("").reverse().join("");
}

function repeat(str) {
    const t = ['gggg', 'ccccc', 'aaaaa', 'ttttt', 'gcgcgc', 'cgcgcg', 'atatat', 'tatata'];
    for (let i = 0; i < t.length - 1; i++) {
        if (str.indexOf(t[i]) > -1) {
            return 0;
        }
    }
    return -1;
}

function AntiSenseDNA(source) {
    let result = [];
    const d = new Array(128).fill(0);
    d[97] = 116;  //'t' a
    d[98] = 118;  //'v' b
    d[99] = 103;  //'g' c
    d[100] = 104; //'h' d
    d[103] = 99;  //'c' g
    d[104] = 100; //'d' h
    d[105] = 99;  //'i' i->c
    d[107] = 109; //'m' k
    d[109] = 107; //'k' m
    d[110] = 110; //'n' n
    d[114] = 121; //'y' r
    d[115] = 115; //'s' s
    d[116] = 97;  //'a' t
    d[117] = 97;  //'a' u
    d[118] = 98;  //'b' v
    d[119] = 119; //'w' w
    d[121] = 114; //'r' y
    for (let i = 0; i < source.length; i++) {
        if (d[source.charCodeAt(i)]) {
            result.push(String.fromCharCode(d[source.charCodeAt(i)]));
        }
    }
    return result.join('');
}

function ComplementDNA(source) {
    return AntiSenseDNA(ReverseDNA(source));
}

function DimerLook(r1, r2, n_base) {
    let z1 = "";
    let z2 = "";
    let z3 = "";
    let j = 0;
    let x = 0;
    let y = 0;
    let r = 0;
    let w = 0;
    let l = 0;
    let l1 = r1.length;
    let l2 = r2.length;

    if (n_base < 3) {
        n_base = 3;
    }
    if (n_base > 10) {
        n_base = 10;
    }
    let initK = n_base - 1;
    if (l1 < l2) {
        z1 = r2;
        z3 = ReverseDNA(r1);
        l = l1;
        l1 = l2;
        l2 = l;
    } else {
        z1 = r1;
        z3 = ReverseDNA(r2);
    }
    z2 = AntiSenseDNA(z3);

    let r1x = z1;
    let r2y = z2;
    for (y = 0; y < l2 - initK; y++) {
        x = -1;
        while (true) {
            x = z1.indexOf(z2.substring(y, y + initK), x + 1);
            if (x === -1) {
                break;
            }
            w = l2 + x - y;
            z = y + l1 - x - 1;
            j = 0;
            l = initK;
            r = initK;

            while (x + l <= l1 - 1) {
                if (y + l > l2 - 1) {
                    break;
                }
                if (r1x[(x + l)] === r2y[(y + l)]) {
                    l += 1;
                    r += 1;
                } else {
                    if (x + l + 1 > l1 - 1) {
                        break;
                    }
                    if (y + l + 1 > l2 - 1) {
                        break;
                    }
                    if (r1x[(x + l + 1)] !== r2y[(y + l + 1)]) {
                        break;
                    }
                    l += 2;
                    r += 1;
                }
            }

            while (x - j - 1 >= 0) {
                if (y - j - 1 < 0) {
                    break;
                }
                if (r1x[(x - j - 1)] === r2y[(y - j - 1)]) {
                    l += 1;
                    r += 1;
                    j += 1;
                } else {
                    if (x - j - 2 < 0) {
                        break;
                    }
                    if (y - j - 2 < 0) {
                        break;
                    }
                    if (r1x[(x - j - 2)] !== r2y[(y - j - 2)]) {
                        break;
                    }
                    l += 2;
                    r += 1;
                    j += 2;
                }
            }
            if (r > n_base) {
                return 0;
            }

        }
    }
    return -1;
}

function RepeatMask(r1) {
    const h = 12;
    const l = r1.length;
    const r2 = ComplementDNA(r1);
    const m = new Array(l).fill(0);
    const d = new Map();
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r1.substring(y, y + h);
        if (!s.includes('n')) {
            if (d.has(s)) {
                let k = d.get(s);
                d.set(s, ++k);
            } else {
                d.set(s, 1);
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r2.substring(y, y + h);
        if (!s.includes('n')) {
            if (d.has(s)) {
                let k = d.get(s);
                d.set(s, ++k);
            } else {
                d.set(s, 1);
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r1.substring(y, y + h);
        if (!s.includes('n')) {
            const k = d.get(s);
            if (k > 1) {
                m.fill(1, y, y + h);
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r2.substring(y, y + h);
        if (!s.includes('n')) {
            const k = d.get(s);
            if (k > 1) {
                m.fill(1, l - y - h, l - y); // ( value, start, end)
            }
        }
    }
    return m;
}

function LingComplexity(source) {
    let r = 0;
    let k = 2;
    let n1 = 0;
    let n2 = 0;
    let n3 = 0;
    let n4 = 0;
    const l = source.length;
    if (l < 4) {
        return 100;
    }

    const dx = new Array(128).fill(0);
    dx[97] = 0; // a = 0
    dx[119] = 0; // at                   
    dx[101] = 0; // a  
    dx[108] = 1; // t = 1            
    dx[116] = 1; // t 
    dx[117] = 1; // u  
    dx[98] = 2; // tgc         
    dx[99] = 2; // c = 2
    dx[102] = 2; // c        
    dx[104] = 2; // atc   
    dx[121] = 2; // tc        
    dx[109] = 2; // ac
    dx[115] = 2; // gc        
    dx[103] = 3; // g = 3
    dx[105] = 3; // g 
    dx[106] = 3; // g      
    dx[107] = 3; // gt   
    dx[114] = 3; // ag         
    dx[118] = 3; // agc   
    dx[100] = 3; // atg 
    dx[110] = 4; // atgc

    const y1 = new Array(5).fill(0);
    const y2 = new Array(5).fill(0).map(() => new Array(5).fill(0));
    const y3 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0)));
    let t = 4 + (l > 16 ? 16 : l - 1);
    if (l > 12) {
        t += (l > 65 ? 64 : l - 2);
        k = 3;
    }
    if (l > 48) {
        t += (l > 258 ? 256 : l - 3);
        k = 4;
    }

    if (k === 4) {
        const y4 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0))));
        for (let i = 0; i < l - 3; i++) {
            n1 = dx[source.charCodeAt(i)];
            n2 = dx[source.charCodeAt(i + 1)];
            n3 = dx[source.charCodeAt(i + 2)];
            n4 = dx[source.charCodeAt(i + 3)];
            if (y1[n1] === 0) {
                r++;
                y1[n1] = 1;
            }
            if (y2[n1][n2] === 0) {
                r++;
                y2[n1][n2] = 1;
            }
            if (y3[n1][n2][n3] === 0) {
                r++;
                y3[n1][n2][n3] = 1;
            }
            if (y4[n1][n2][n3][n4] === 0) {
                r++;
                y4[n1][n2][n3][n4] = 1;
            }
        }
        if (y3[n2][n3][n4] === 0) {
            r++;
            y3[n2][n3][n4] = 1;
        }
        if (y2[n3][n4] === 0) {
            r++;
            y2[n3][n4] = 1;
        }
        if (y2[n2][n3] === 0) {
            r++;
            y2[n2][n3] = 1;
        }
        if (y1[n2] === 0) {
            r++;
            y1[n2] = 1;
        }
        if (y1[n3] === 0) {
            r++;
            y1[n3] = 1;
        }
        if (y1[n4] === 0) {
            r++;
            y1[n4] = 1;
        }
    }

    if (k === 3) {
        for (let i = 0; i < l - 2; i++) {
            n1 = dx[source.charCodeAt(i)];
            n2 = dx[source.charCodeAt(i + 1)];
            n3 = dx[source.charCodeAt(i + 2)];
            if (y1[n1] === 0) {
                r++;
                y1[n1] = 1;
            }
            if (y2[n1][n2] === 0) {
                r++;
                y2[n1][n2] = 1;
            }
            if (y3[n1][n2][n3] === 0) {
                r++;
                y3[n1][n2][n3] = 1;
            }
        }
        if (y2[n2][n3] === 0) {
            r++;
            y2[n2][n3] = 1;
        }
        if (y1[n2] === 0) {
            r++;
            y1[n2] = 1;
        }
        if (y1[n3] === 0) {
            r++;
            y1[n3] = 1;
        }
    }

    if (k === 2) {
        for (let i = 0; i < l - 1; i++) {
            n1 = dx[source.charCodeAt(i)];
            n2 = dx[source.charCodeAt(i + 1)];
            if (y1[n1] === 0) {
                r++;
                y1[n1] = 1;
            }
            if (y2[n1][n2] === 0) {
                r++;
                y2[n1][n2] = 1;
            }
        }
        if (y1[n2] === 0) {
            r++;
            y1[n2] = 1;
        }
    }
    return Math.floor((100 * r) / t);
}