function DNASeq(str) {
        const validChars = new Set(['a', 't', 'c', 'g', 'r', 'y', 'm', 'k', 'w', 'b', 'd', 'v', 'h', 's', 'n']);
        let returnString = "";
        for (let i = 0; i < str.length; i++) {
                let chr = str.charAt(i);
                if (validChars.has(chr)) {
                        returnString += chr;
                } else {
                        switch (chr) {
                                case 'i':
                                        returnString += 'g';
                                        break;
                                case 'u':
                                        returnString += 't';
                                        break;
                                case 'e':
                                        returnString += 'a';
                                        break;
                                case 'f':
                                        returnString += 'c';
                                        break;
                                case 'j':
                                        returnString += 'g';
                                        break;
                                case 'l':
                                        returnString += 't';
                                        break;
                                default:
                                        // Optionally handle other cases
                                        break;
                        }
                }
        }
        return returnString;
}
function DNA(str) {
        var returnString = "";
        for (var i = 0; i < str.length; i++) {
                var chr = str.charAt(i);
                if (chr === 'a' || chr === 't' || chr === 'c' || chr === 'g' || chr === 'r' || chr === 'y' || chr === 'm' || chr === 'k' || chr === 'w' || chr === 'b' || chr === 'd' || chr === 'v' || chr === 'h' || chr === 's' || chr === 'n') {
                        returnString += chr;
                }
        }
        return returnString;
}

function PrimerSeq(str) {
        var returnString = "";
        for (var i = 0; i < str.length; i++) {
                var chr = str.charAt(i);
                if (chr === 'a' || chr === 't' || chr === 'c' || chr === 'g' || chr === 'r' || chr === 'y' || chr === 'm' || chr === 'k' || chr === 'w' || chr === 'b' || chr === 'd' || chr === 'v' || chr === 'h' || chr === 's' || chr === 'n' || chr === 'u' || chr === 'i' || chr === 'j' || chr === 'l' || chr === 'f' || chr === 'e') {
                        returnString += chr;
                }
        }
        return returnString;
}

function ReverseDNA(source) {
        return source.split("").reverse().join("");
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
        d[105] = 99;  //'i'->c
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

        d[40] = 41;  //'(' 
        d[41] = 40;  //')' 

        //  d[47] = 92;  //'/' 
        //  d[92] = 47;  //'\' 

        d[91] = 93;  //'[' 
        d[93] = 91;  //']' 

        d[123] = 125;  //'{' 
        d[125] = 123;  //'}' 

        for (var i = 0; i < source.length; i++) {
                if (d[source.charCodeAt(i)] > 0) {
                        result.push(String.fromCharCode(d[source.charCodeAt(i)]));
                }
                else {
                        result.push(String.fromCharCode(source.charCodeAt(i)));
                }
        }
        return result.join('');
}
function ComplementDNA(source) {
        const d = new Array(128).fill(0);
        d[97] = 116;  //'t' a
        d[98] = 118;  //'v' b
        d[99] = 103;  //'g' c
        d[100] = 104; //'h' d
        d[103] = 99;  //'c' g
        d[104] = 100; //'d' h
        d[105] = 99;  //'i'->c
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
        let b = source.split('');  // Convert string to array
        let l = source.length;
        let n = Math.floor(l / 2);

        for (let i = 0; i < n; i++) {
                let v = d[source.charCodeAt(i)];
                if (v > 0) {
                        let t = d[source.charCodeAt(l - i - 1)];
                        b[l - i - 1] = String.fromCharCode(v);
                        b[i] = String.fromCharCode(t);
                }
        }
        if ((l % 2) === 1) {
                let v = d[source.charCodeAt(n)];
                if (v > 0) {
                        b[n] = String.fromCharCode(v);
                }
        }
        return b.join('');          // Convert array back to string
}

function ComplementDNA3(source) {
        return AntiSenseDNA(source).split("").reverse().join("");
}
function ComplementDNA2(source) {
        const cdna = { 'a': 't', 'b': 'v', 'c': 'g', 'd': 'h', 'g': 'c', 'h': 'd', 'k': 'm', 'm': 'k', 'n': 'n', 'y': 'r', 's': 's', 't': 'a', 'u': 'a', 'v': 'b', 'w': 'w', 'r': 'y' };
        let b = source.split('');  // Convert string to array
        let l = source.length;
        let n = Math.floor(l / 2);

        for (let i = 0; i < n; i++) {
                if (cdna[b[i]]) {
                        let t = cdna[b[l - i - 1]];
                        b[l - i - 1] = cdna[b[i]];
                        b[i] = t;
                }
        }
        if ((l % 2) === 1) {
                if (cdna[b[n]]) {
                        b[n] = cdna[b[n]];
                }
        }
        return b.join('');          // Convert array back to string
}

function MW(an, tn, cn, gn, ni, un) {
        return (an * 313.21) + ((tn - un) * 304.2) + (cn * 289.18) + ((gn - ni) * 329.21) + (ni * 314) + (un * 290.169) - 61.96;
}

function CG(str) {
        const sdn = new Array(128).fill(0);
        //M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T
        sdn[98] = 0.66666;   //B=(C/G/T)
        sdn[99] = 1;         //c
        sdn[100] = 0.33333;  //D=(A/G/T)
        sdn[103] = 1;        //g
        sdn[104] = 0.33333;  //H=(A/C/T)
        sdn[105] = 1;        //i=G
        sdn[107] = 0.5;      //K=(G/T)
        sdn[109] = 0.5;      //M=(A/C)
        sdn[110] = 0.5;      //N=(A/G/C/T)
        sdn[114] = 0.5;      //R=(A/G)
        sdn[115] = 1;        //S=(G/C)
        sdn[118] = 0.66666;  //V=(A/G/C)
        sdn[121] = 0.5;      //Y=(C/T)
        let dr = parseFloat(0);
        const l = str.length;
        for (var i = 0; i < l; i++) {
                dr += sdn[str.charCodeAt(i)];
        }
        return ((100 * dr) / l);
}

var LingComplexity = function (source) {
        let r = 0;
        let k = 2;
        let n1 = 0;
        let n2 = 0;
        let n3 = 0;
        let n4 = 0;
        let n5 = 0;
        let n6 = 0;
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
        if (l > 192) {
                t += (l > 1027 ? 1024 : l - 4);
                k = 5;
        }
        if (l > 768) {
                t += (l > 5000 ? 4096 : l - 5);
                k = 6;
        }

        if (k === 6) {
                const y4 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0))));
                const y5 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0)))));
                const y6 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0))))));
                for (let i = 0; i < l - 5; i++) {
                        n1 = dx[source.charCodeAt(i)];
                        n2 = dx[source.charCodeAt(i + 1)];
                        n3 = dx[source.charCodeAt(i + 2)];
                        n4 = dx[source.charCodeAt(i + 3)];
                        n5 = dx[source.charCodeAt(i + 4)];
                        n6 = dx[source.charCodeAt(i + 5)];
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
                        if (y5[n1][n2][n3][n4][n5] === 0) {
                                r++;
                                y5[n1][n2][n3][n4][n5] = 1;
                        }
                        if (y6[n1][n2][n3][n4][n5][n6] === 0) {
                                r++;
                                y6[n1][n2][n3][n4][n5][n6] = 1;
                        }
                }
                if (y5[n2][n3][n4][n5][n6] === 0) {
                        r++;
                        y5[n2][n3][n4][n5][n6] = 1;
                }
                if (y4[n3][n4][n5][n6] === 0) {
                        r++;
                        y4[n3][n4][n5][n6] = 1;
                }
                if (y4[n2][n3][n4][n5] === 0) {
                        r++;
                        y4[n2][n3][n3][n5] = 1;
                }
                if (y3[n4][n5][n6] === 0) {
                        r++;
                        y3[n4][n5][n6] = 1;
                }
                if (y3[n3][n4][n5] === 0) {
                        r++;
                        y3[n3][n4][n5] = 1;
                }
                if (y3[n2][n3][n4] === 0) {
                        r++;
                        y3[n2][n3][n4] = 1;
                }
                if (y2[n5][n6] === 0) {
                        r++;
                        y2[n5][n6] = 1;
                }
                if (y2[n4][n5] === 0) {
                        r++;
                        y2[n4][n5] = 1;
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
                if (y1[n5] === 0) {
                        r++;
                        y1[n5] = 1;
                }
                if (y1[n6] === 0) {
                        r++;
                        y1[n6] = 1;
                }
        }

        if (k === 5) {
                const y4 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0))));
                const y5 = new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0).map(() => new Array(5).fill(0)))));
                for (let i = 0; i < l - 4; i++) {
                        n1 = dx[source.charCodeAt(i)];
                        n2 = dx[source.charCodeAt(i + 1)];
                        n3 = dx[source.charCodeAt(i + 2)];
                        n4 = dx[source.charCodeAt(i + 3)];
                        n5 = dx[source.charCodeAt(i + 4)];
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
                        if (y5[n1][n2][n3][n4][n5] === 0) {
                                r++;
                                y5[n1][n2][n3][n4][n5] = 1;
                        }
                }
                if (y4[n2][n3][n4][n5] === 0) {
                        r++;
                        y4[n2][n3][n3][n5] = 1;
                }
                if (y3[n3][n4][n5] === 0) {
                        r++;
                        y3[n3][n4][n5] = 1;
                }
                if (y3[n2][n3][n4] === 0) {
                        r++;
                        y3[n2][n3][n4] = 1;
                }
                if (y2[n4][n5] === 0) {
                        r++;
                        y2[n4][n5] = 1;
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
                if (y1[n5] === 0) {
                        r++;
                        y1[n5] = 1;
                }
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


function LingComplexity2(source) {
        let r = 0;
        let k = 2;
        let n1 = 0;
        let n2 = 0;
        let n3 = 0;
        let n4 = 0;
        let n5 = 0;
        let n6 = 0;
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

function RepeatMask(r1) {
        const h = 12;
        const l = r1.length;
        const r2 = ComplementDNA(r1);
        const m = new Array(l).fill(0);
        const d = new Map();
        for (var y = 0; y < r1.length - h + 1; y++) {
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
        for (var y = 0; y < r1.length - h + 1; y++) {
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
        for (var y = 0; y < r1.length - h + 1; y++) {
                let s = r1.substring(y, y + h);
                if (!s.includes('n')) {
                        const k = d.get(s);
                        if (k > 1) {
                                m.fill(1, y, y + h);
                        }
                }
        }
        for (var y = 0; y < r1.length - h + 1; y++) {
                let s = r2.substring(y, y + h);
                if (!s.includes('n')) {
                        const k = d.get(s);
                        if (k > 1) {
                                m.fill(2, l - y - h, l - y); // ( value, start, end)
                        }
                }
        }
        return m;
}

function Tm65(str, dmso) {
        let cx = 0;
        let rx = 0;
        let hx = 0;
        const lseq = str.length
        for (var i = 0; i < lseq; i++) {
                var chr = str.charAt(i);
                if (chr === 'c' || chr === 'g' || chr === 's' || chr === 'i') {
                        cx++;
                }
                if (chr === 'r' || chr === 'y' || chr === 'k' || chr === 'n') {
                        rx++;
                }
                if (chr === 'd' || chr === 'h') {
                        hx++;
                }
                if (chr === 'b' || chr === 'v') {
                        hx += 2;
                }

        }
        return (lseq < 1) ? 0 : (64.9 - dmso + (41 * (cx + (rx / 2) + (hx / 3) - 16.4)) / lseq);
}
function Tm77(str, Salt_M, dmso) {
        let cx = 0;
        let rx = 0;
        let hx = 0;
        const lseq = str.length
        for (var i = 0; i < lseq; i++) {
                var chr = str.charAt(i);
                if (chr === 'c' || chr === 'g' || chr === 's' || chr === 'i') {
                        cx++;
                }
                if (chr === 'r' || chr === 'y' || chr === 'k' || chr === 'n') {
                        rx++;
                }
                if (chr === 'd' || chr === 'h') {
                        hx++;
                }
                if (chr === 'b' || chr === 'v') {
                        hx += 2;
                }
        }
        return (lseq < 1) ? 0 : (77.1 - dmso + 11.7 * Math.log10(Salt_M) + (41 * (cx + (rx / 2) + (hx / 3)) - 528) / lseq);       //"°C (Tm = 77.1 + 11.7Log[K+] + (41*(G + C) - 528)/L)\n\n"
}
function Tm85(str, Salt_M, dmso) {
        let cx = 0;
        let rx = 0;
        let hx = 0;
        const lseq = str.length
        for (var i = 0; i < lseq; i++) {
                var chr = str.charAt(i);
                if (chr === 'c' || chr === 'g' || chr === 's' || chr === 'i') {
                        cx++;
                }
                if (chr === 'r' || chr === 'y' || chr === 'k' || chr === 'n') {
                        rx++;
                }
                if (chr === 'd' || chr === 'h') {
                        hx++;
                }
                if (chr === 'b' || chr === 'v') {
                        hx += 2;
                }
        }
        return (lseq < 1) ? 0 : (81.5 - dmso + 16.6 * Math.log10(Salt_M) + (41 * (cx + (rx / 2) + (hx / 3)) - 675) / lseq);
}
function DNActBisulfiteForward(source) {
        let s = AntiSenseDNA(source);
        s = s.replace(/gc/g, "11");
        s = s.replace(/c/g, "t");
        s = s.replace(/11/g, "gc");
        s = AntiSenseDNA(s);
        return s;
}
function DNActBisulfiteReverse(source) {
        let s = source.replace(/cg/g, "11");
        s = s.replace(/c/g, "t");
        s = s.replace(/11/g, "cg");
        s = ComplementDNA(s);
        return s;
}
function ssrrepeat(str) {
        const substrings = ['gggg', 'ccccc', 'aaaaa', 'ttttt', 'gcgcgc', 'cgcgcg', 'atatat', 'tatata'];
        for (const sub of substrings) {
                if (str.includes(sub)) {
                        return 0;
                }
        }
        return -1;
}
function strgen(str) {
        const mappings = { 'b': "t c g", 'd': "a t g", 'h': "a t c", 'k': "t g", 'm': "a c", 'n': "a t c g", 'r': "a g", 's': "c g", 'v': "a c g", 'w': "a t", 'y': "t c" };
        return mappings[str] || str;
}

function GeneratorVariants(pol) {
        if (pol === undefined) { return "a t c g" }
        //pol = "swh ssw wsh sww www";
        let k = 1;
        let n = pol.indexOf('/');
        let b = [];

        if (n > -1) {
                b = pol.split("/");
                for (let i = 0; i < b.length; i++) {
                        b[i] = DNA(b[i]);
                }
        }
        else {
                pol = DNA(pol);
                for (let i = 0; i < pol.length; i++) {
                        let s = strgen(pol.charAt(i));
                        let c = s.split(" ");
                        let t = c.length;
                        if (t > 1) {
                                k = k * t;
                        }
                }
                let d = new Array(k).fill("");

                let f = 1;
                for (let i = 0; i < pol.length; i++) {
                        let s = strgen(pol.charAt(i));
                        let c = s.split(" ");
                        let t = c.length;
                        if (t > 1) {
                                let z = [];
                                for (let j = 0; j < f; j++) {
                                        for (let i = 0; i < t; i++) {
                                                z.push(d[j] + c[i]);
                                        }
                                }
                                f = z.length;
                                for (let j = 0; j < f; j++) {
                                        d[j] = z[j];
                                }
                                z.length = 0;
                        }
                        else {
                                for (let j = 0; j <= f; j++) {
                                        d[j] += pol.charAt(i);
                                }
                        }
                }
                for (let i = 0; i < k; i++) {
                        b.push(d[i]);
                }
        }
        return b;
}
function findCommonString(arr) {
        const frequencyMap = new Map();
        for (const str of arr) {
                if (frequencyMap.has(str)) {
                        frequencyMap.set(str, frequencyMap.get(str) + 1);
                } else {
                        frequencyMap.set(str, 1);
                }
        }
        for (const [key, value] of frequencyMap) {
                if (value > 1) {
                        return key; // Return the first common string found
                }
        }
        return -1;
}
function countUniqueSNPs(arr) {
        const n = arr.length;
        const l = arr[0].length;
        let uniqueSNPCount = 0;

        // Loop through each character position
        for (let i = 0; i < l; i++) {
                const seenChars = new Set();
                let allUnique = true;

                // Check each string in the array at position i
                for (const str of arr) {
                        const chr = str.charAt(i);
                        if (seenChars.has(chr)) {
                                allUnique = false;
                                break; // Stop checking this position if duplicate found
                        }
                        seenChars.add(chr);
                }

                // If all characters at this position are unique, increment the counter
                if (allUnique) {
                        uniqueSNPCount++;
                }
        }
        return uniqueSNPCount;
}

function GeneratorVariants2(pol) {
        if (pol === undefined) { return "a t c g", 1; }
        //  pol = "swh ssw wsh sww www";
        let b = pol.split(" ");
        const n = b.length;
        let z = 0;

        for (let j = 0; j < n; j++) {
                b[j] = DNA(b[j].trim());
                let k = 1;
                z = b[j].length;

                for (let i = 0; i < z; i++) {
                        let s = strgen(b[j].charAt(i));
                        let c = s.split(" ");
                        let t = c.length;
                        if (t > 1) {
                                k = k * t;
                        }
                }
                let d = new Array(k).fill("");

                let f = 1;
                for (let i = 0; i < z; i++) {
                        let s = strgen(b[j].charAt(i));
                        let c = s.split(" ");
                        let t = c.length;
                        if (t > 1) {
                                let z = [];
                                for (let j = 0; j < f; j++) {
                                        for (let i = 0; i < t; i++) {
                                                z.push(d[j] + c[i]);
                                        }
                                }
                                f = z.length;
                                for (let r = 0; r < f; r++) {
                                        d[r] = z[r];
                                }
                                z.length = 0;
                        }
                        else {
                                for (let r = 0; r <= f; r++) {
                                        d[r] += b[j].charAt(i);
                                }
                        }
                }
                for (let i = 0; i < k; i++) {
                        b.push(d[i]);
                }
        }
        return { b, z };
}
