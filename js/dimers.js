// "M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T and I"
//        a  CGT  c  AGT  a   c   g   ACT g   g   GT  t   AC  ATGC AG  GC  t  t   AGC AT   CT
//        a   b   c   d   e   f   g   h   i   j   k   l   m    n   r   s   t  u   v   w    y
//        97  98  99  100 101 102 103 104 105 106 107 108 109 110 114 115 116 117 118 119 121
//         G-C > A-T > G·G > G·T > G·A         
const aDNA = new Array(128).fill(0).map(() => new Array(128).fill(0));
aDNA[97][98] = 33;  //a-b
aDNA[97][100] = 33; //a-d
aDNA[97][104] = 33; //a-h
aDNA[97][105] = 50; //a-i
aDNA[97][107] = 50; //a-k
aDNA[97][108] = 100; //a-l
aDNA[97][110] = 25;  //a-n
aDNA[97][116] = 100; //a-t
aDNA[97][117] = 100; //a-u
aDNA[97][119] = 50; //a-w
aDNA[97][121] = 50; //a-y
aDNA[98][97] = 33;  //b-a
aDNA[98][98] = 67;  //b-b
aDNA[98][99] = 33;  //b-c
aDNA[98][100] = 67; //b-d
aDNA[98][101] = 33; //b-e
aDNA[98][102] = 33; //b-f
aDNA[98][103] = 33; //b-g
aDNA[98][104] = 67; //b-h
aDNA[98][105] = 33; //b-i
aDNA[98][106] = 33; //b-j
aDNA[98][107] = 33; //b-k
aDNA[98][109] = 67; //b-m
aDNA[98][110] = 75; //b-n
aDNA[98][114] = 67; //b-r
aDNA[98][115] = 67; //b-s
aDNA[98][118] = 100; //b-v
aDNA[98][119] = 33;  //b-w
aDNA[98][121] = 33;  //b-y
aDNA[99][98] = 33;   //c-b
aDNA[99][100] = 33;  //c-d
aDNA[99][103] = 100; //c-g
aDNA[99][105] = 100; //c-i
aDNA[99][106] = 100; //c-j
aDNA[99][107] = 50;  //c-k
aDNA[99][110] = 25;  //c-n
aDNA[99][114] = 50;  //c-r
aDNA[99][115] = 50;  //c-s
aDNA[99][118] = 33;  //c-v
aDNA[100][97] = 33;  //d-a
aDNA[100][98] = 67;  //d-b
aDNA[100][99] = 33;  //d-c
aDNA[100][100] = 67; //d-d
aDNA[100][101] = 33; //d-e
aDNA[100][102] = 33; //d-f
aDNA[100][104] = 100; //d-h
aDNA[100][107] = 33;  //d-k
aDNA[100][108] = 33;  //d-l
aDNA[100][109] = 67;  //d-m
aDNA[100][110] = 75;  //d-n
aDNA[100][114] = 33; //d-r
aDNA[100][115] = 33; //d-s
aDNA[100][116] = 33; //d-t
aDNA[100][117] = 33; //d-u
aDNA[100][118] = 67; //d-v
aDNA[100][119] = 67; //d-w
aDNA[100][121] = 67; //d-y
aDNA[101][98] = 33; //e-b
aDNA[101][100] = 33; //e-d
aDNA[101][104] = 33; //e-h
aDNA[101][105] = 50; //e-i
aDNA[101][107] = 50; //e-k
aDNA[101][108] = 100; //e-l
aDNA[101][110] = 25; //e-n
aDNA[101][116] = 100; //e-t
aDNA[101][117] = 100; //e-u
aDNA[101][119] = 50; //e-w
aDNA[101][121] = 50; //e-y
aDNA[102][98] = 33; //f-b
aDNA[102][100] = 33; //f-d
aDNA[102][103] = 100; //f-g
aDNA[102][105] = 100; //f-i
aDNA[102][106] = 100; //f-j
aDNA[102][107] = 50; //f-k
aDNA[102][110] = 25; //f-n
aDNA[102][114] = 50; //f-r
aDNA[102][115] = 50; //f-s
aDNA[102][118] = 33; //f-v
aDNA[103][98] = 33; //g-b
aDNA[103][99] = 100; //g-c
aDNA[103][102] = 100; //g-f
aDNA[103][104] = 33; //g-h
aDNA[103][105] = 50; //g-i
aDNA[103][109] = 50; //g-m
aDNA[103][110] = 25; //g-n
aDNA[103][115] = 50; //g-s
aDNA[103][118] = 33; //g-v
aDNA[103][121] = 50; //g-y
aDNA[104][97] = 33; //h-a
aDNA[104][98] = 67; //h-b
aDNA[104][100] = 100; //h-d
aDNA[104][101] = 33; //h-e
aDNA[104][103] = 33; //h-g
aDNA[104][104] = 67; //h-h
aDNA[104][105] = 33; //h-i
aDNA[104][106] = 33; //h-j
aDNA[104][107] = 67; //h-k
aDNA[104][108] = 33; //h-l
aDNA[104][109] = 33; //h-m
aDNA[104][110] = 75; //h-n
aDNA[104][114] = 67; //h-r
aDNA[104][115] = 33; //h-s
aDNA[104][116] = 33; //h-t
aDNA[104][117] = 33; //h-u
aDNA[104][118] = 67; //h-v
aDNA[104][119] = 67; //h-w
aDNA[104][121] = 33; //h-y
aDNA[105][97] = 50; //i-a
aDNA[105][98] = 33; //i-b
aDNA[105][99] = 100; //i-c
aDNA[105][101] = 50; //i-e
aDNA[105][102] = 100; //i-f
aDNA[105][103] = 50; //i-g
aDNA[105][104] = 33; //i-h
aDNA[105][106] = 50; //i-j
aDNA[105][108] = 50; //i-l
aDNA[105][109] = 50; //i-m
aDNA[105][110] = 25; //i-n
aDNA[105][115] = 50; //i-s
aDNA[105][116] = 50; //i-t
aDNA[105][117] = 50; //i-u
aDNA[105][118] = 33; //i-v
aDNA[105][121] = 50; //i-y
aDNA[106][98] = 33; //j-b
aDNA[106][99] = 100; //j-c
aDNA[106][102] = 100; //j-f
aDNA[106][104] = 33; //j-h
aDNA[106][105] = 50; //j-i
aDNA[106][109] = 50; //j-m
aDNA[106][110] = 25; //j-n
aDNA[106][115] = 50; //j-s
aDNA[106][118] = 33; //j-v
aDNA[106][121] = 50; //j-y
aDNA[107][97] = 50; //k-a
aDNA[107][98] = 33; //k-b
aDNA[107][99] = 50; //k-c
aDNA[107][100] = 33; //k-d
aDNA[107][101] = 50; //k-e
aDNA[107][102] = 50; //k-f
aDNA[107][104] = 67; //k-h
aDNA[107][109] = 100; //k-m
aDNA[107][110] = 50; //k-n
aDNA[107][114] = 50; //k-r
aDNA[107][115] = 50; //k-s
aDNA[107][118] = 67; //k-v
aDNA[107][119] = 50; //k-w
aDNA[107][121] = 50; //k-y
aDNA[108][97] = 100; //l-a
aDNA[108][100] = 33; //l-d
aDNA[108][101] = 100; //l-e
aDNA[108][104] = 33; //l-h
aDNA[108][105] = 50; //l-i
aDNA[108][109] = 50; //l-m
aDNA[108][110] = 25; //l-n
aDNA[108][114] = 50; //l-r
aDNA[108][118] = 33; //l-v
aDNA[108][119] = 50; //l-w
aDNA[109][98] = 67; //m-b
aDNA[109][100] = 67; //m-d
aDNA[109][103] = 50; //m-g
aDNA[109][104] = 33; //m-h
aDNA[109][105] = 50; //m-i
aDNA[109][106] = 50; //m-j
aDNA[109][107] = 100; //m-k
aDNA[109][108] = 50; //m-l
aDNA[109][110] = 50; //m-n
aDNA[109][114] = 50; //m-r
aDNA[109][115] = 50; //m-s
aDNA[109][116] = 50; //m-t
aDNA[109][117] = 50; //m-u
aDNA[109][118] = 33; //m-v
aDNA[109][119] = 50; //m-w
aDNA[109][121] = 50; //m-y
aDNA[110][97] = 25; //n-a
aDNA[110][98] = 75; //n-b
aDNA[110][99] = 25; //n-c
aDNA[110][100] = 75; //n-d
aDNA[110][101] = 25; //n-e
aDNA[110][102] = 25; //n-f
aDNA[110][103] = 25; //n-g
aDNA[110][104] = 75; //n-h
aDNA[110][105] = 25; //n-i
aDNA[110][106] = 25; //n-j
aDNA[110][107] = 50; //n-k
aDNA[110][108] = 25; //n-l
aDNA[110][109] = 50; //n-m
aDNA[110][110] = 100; //n-n
aDNA[110][114] = 50; //n-r
aDNA[110][115] = 50; //n-s
aDNA[110][116] = 25; //n-t
aDNA[110][117] = 25; //n-u
aDNA[110][118] = 75; //n-v
aDNA[110][119] = 50; //n-w
aDNA[110][121] = 50; //n-y
aDNA[114][98] = 67; //r-b
aDNA[114][99] = 50; //r-c
aDNA[114][100] = 33; //r-d
aDNA[114][102] = 50; //r-f
aDNA[114][104] = 67; //r-h
aDNA[114][107] = 50; //r-k
aDNA[114][108] = 50; //r-l
aDNA[114][109] = 50; //r-m
aDNA[114][110] = 50; //r-n
aDNA[114][115] = 50; //r-s
aDNA[114][116] = 50; //r-t
aDNA[114][117] = 50; //r-u
aDNA[114][118] = 33; //r-v
aDNA[114][119] = 50; //r-w
aDNA[114][121] = 100; //r-y
aDNA[115][98] = 67; //s-b
aDNA[115][99] = 50; //s-c
aDNA[115][100] = 33; //s-d
aDNA[115][102] = 50; //s-f
aDNA[115][103] = 50; //s-g
aDNA[115][104] = 33; //s-h
aDNA[115][105] = 50; //s-i
aDNA[115][106] = 50; //s-j
aDNA[115][107] = 50; //s-k
aDNA[115][109] = 50; //s-m
aDNA[115][110] = 50; //s-n
aDNA[115][114] = 50; //s-r
aDNA[115][115] = 100; //s-s
aDNA[115][118] = 67; //s-v
aDNA[115][121] = 50; //s-y
aDNA[116][97] = 100; //t-a
aDNA[116][100] = 33; //t-d
aDNA[116][101] = 100; //t-e
aDNA[116][104] = 33; //t-h
aDNA[116][105] = 50; //t-i
aDNA[116][109] = 50; //t-m
aDNA[116][110] = 25; //t-n
aDNA[116][114] = 50; //t-r
aDNA[116][118] = 33; //t-v
aDNA[116][119] = 50; //t-w
aDNA[117][97] = 100; //u-a
aDNA[117][100] = 33; //u-d
aDNA[117][101] = 100; //u-e
aDNA[117][104] = 33; //u-h
aDNA[117][105] = 50; //u-i
aDNA[117][109] = 50; //u-m
aDNA[117][110] = 25; //u-n
aDNA[117][114] = 50; //u-r
aDNA[117][118] = 33; //u-v
aDNA[117][119] = 50; //u-w
aDNA[118][98] = 100; //v-b
aDNA[118][99] = 33; //v-c
aDNA[118][100] = 67; //v-d
aDNA[118][102] = 33; //v-f
aDNA[118][103] = 33; //v-g
aDNA[118][104] = 67; //v-h
aDNA[118][105] = 33; //v-i
aDNA[118][106] = 33; //v-j
aDNA[118][107] = 67; //v-k
aDNA[118][108] = 33; //v-l
aDNA[118][109] = 33; //v-m
aDNA[118][110] = 75; //v-n
aDNA[118][114] = 33; //v-r
aDNA[118][115] = 67; //v-s
aDNA[118][116] = 33; //v-t
aDNA[118][117] = 33; //v-u
aDNA[118][118] = 67; //v-v
aDNA[118][119] = 33; //v-w
aDNA[118][121] = 67; //v-y
aDNA[119][97] = 50; //w-a
aDNA[119][98] = 33; //w-b
aDNA[119][100] = 67; //w-d
aDNA[119][101] = 50; //w-e
aDNA[119][104] = 67; //w-h
aDNA[119][107] = 50; //w-k
aDNA[119][108] = 50; //w-l
aDNA[119][109] = 50; //w-m
aDNA[119][110] = 50; //w-n
aDNA[119][114] = 50; //w-r
aDNA[119][116] = 50; //w-t
aDNA[119][117] = 50; //w-u
aDNA[119][118] = 33; //w-v
aDNA[119][119] = 100; //w-w
aDNA[119][121] = 50; //w-y
aDNA[121][97] = 50; //y-a
aDNA[121][98] = 33; //y-b
aDNA[121][100] = 67; //y-d
aDNA[121][101] = 50; //y-e
aDNA[121][103] = 50; //y-g
aDNA[121][104] = 33; //y-h
aDNA[121][105] = 50; //y-i
aDNA[121][106] = 50; //y-j
aDNA[121][107] = 50; //y-k
aDNA[121][109] = 50; //y-m
aDNA[121][110] = 50; //y-n
aDNA[121][114] = 100; //y-r
aDNA[121][115] = 50; //y-s
aDNA[121][118] = 67; //y-v
aDNA[121][119] = 50; //y-w

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
const cdna = new Array(128).fill(0);
cdna[65] = 116;//A
cdna[66] = 118;//B
cdna[67] = 103;//C
cdna[68] = 104;//D
cdna[71] = 99; //G
cdna[72] = 100;//H
cdna[73] = 99; //I
cdna[75] = 109;//K
cdna[77] = 107;//M
cdna[78] = 110;//N
cdna[82] = 121;//R
cdna[83] = 115;//S
cdna[84] = 97; //T
cdna[85] = 97; //U
cdna[86] = 98; //V
cdna[87] = 119;//W
cdna[89] = 114;//Y        
cdna[97] = 116;//  t <- a
cdna[98] = 118;//  v <- b
cdna[99] = 103;//  g <- c
cdna[100] = 104;// h <- d
cdna[103] = 99; // c <- g
cdna[104] = 100;// d <- h
cdna[105] = 99; // i <- g
cdna[107] = 109;// m <- k
cdna[109] = 107;// k <- m
cdna[110] = 110;// n <- n
cdna[114] = 121;// y <- r
cdna[115] = 115;// s <- s
cdna[116] = 97; // a <- t
cdna[117] = 97; // a <- u
cdna[118] = 98; // b <- v
cdna[119] = 119;// w <- w
cdna[121] = 114;// r <- y    

function DimerLook(r1, r2, n_base, m0, m1, m2, m3) {
    let z1 = "";
    let z2 = "";
    let z3 = "";
    let j = 0;
    let n = 0;
    let x = 0;
    let y = 0;
    let z = 0;
    let r = 0;
    let f = - 1;
    let w = 0;
    let l = 0;
    let l1 = r1.length;
    let l2 = r2.length;
    r1 = r1.toLowerCase();
    r2 = r2.toLowerCase();

    if (n_base < 1) n_base = 1;
    if (n_base > 10) n_base = 10;
    let initK = n_base + 1;
    let StrLen = n_base + 2;
    let FilterSize = n_base + 3;
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
    let px = new Array(l1 + l2 + l2).fill(0);
    let r1x = z1;
    let r2y = z2;
    for (y = 0; y < l2 - initK; y++) {
        x = - 1;
        while (true) {
            x = z1.indexOf(z2.substring(y, y + initK), x + 1);
            if (x == - 1) {
                break;
            }
            w = l2 + x - y;
            z = y + l1 - x - 1;
            j = 0;
            l = initK;
            r = initK;
            if (px[w] == 0) {
                while (x + l <= l1 - 1) {
                    if (y + l > l2 - 1) {
                        break;
                    }
                    if (r1x[(x + l)] == r2y[(y + l)]) {
                        l += 1;
                        r += 1;
                    } else {
                        if (x + l + 1 > l1 - 1) {
                            break;
                        }
                        if (y + l + 1 > l2 - 1) {
                            break;
                        }
                        if (r1x[(x + l + 1)] != r2y[(y + l + 1)]) { break; }
                        l += 2;
                        r += 1;
                    }

                }

                while (x - j - 1 >= 0) {
                    if (y - j - 1 < 0) {
                        break;
                    }
                    if (r1x[(x - j - 1)] == r2y[(y - j - 1)]) {
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
                        if (r1x[(x - j - 2)] != r2y[(y - j - 2)]) { break; }
                        l += 2;
                        r += 1;
                        j += 2;
                    }

                }

                if (r > StrLen) {
                    f++;
                    m0[f] = 0;
                    if (r > FilterSize) {
                        m0[f] = 3;
                    }

                    if ((y - j < 2) && (x - j > 1)) {
                        m0[f] = 1;
                    }

                    if ((z < l2) && (x + 1 - j > l1 - 1)) {
                        m0[f] = 2;
                    }

                    if (m0[f] > 0) {
                        px[w] = 1;
                        if (w > l2 - 1) {
                            w -= l2;
                            m1[f] = ("5-" + z1 + "->");
                            let v = [];
                            for (j = 0; j < w; j++) {
                                v.push(' ');
                            }
                            m3[f] = (v.join('') + "<-" + z3 + "-5");
                            v = [];
                            for (j = 0; j < w + 2; j++) {
                                v.push(' ');
                            }
                            m2[f] = v.join('');
                            n = l1 - w;
                            if (n > l2 - 1) {
                                n = l2;
                            }
                            for (j = 0; j < n; j++) {
                                if (r1x[(w + j)] == r2y[j])
                                    m2[f] = (m2[f] + "|");
                                else {
                                    m2[f] = (m2[f] + " ");
                                }
                            }
                        }
                        else {
                            w = l2 - w;
                            m3[f] = ("<-" + z3 + "-5");
                            let v = [];
                            for (j = 0; j < w + 2; j++) {
                                v.push(' ');
                            }
                            m1[f] = (v.join('') + "5-" + z1 + "->");
                            v = [];
                            for (j = 0; j < w + 2; j++) {
                                v.push(' ');
                            }
                            m2[f] = v.join('');
                            for (j = 0; j < l2 - w; j++) {
                                if (r1x[j] == r2y[(j + w)])
                                    m2[f] = (m2[f] + "|");
                                else
                                    m2[f] = (m2[f] + " ");
                            }
                        }
                    }
                    else {
                        f--;
                    }
                }
            }
        }
    }
    return f;
}
function DimerLook3(r1, r2, n_base) {
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

    r1 = r1.toLowerCase();
    r2 = r2.toLowerCase();

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



function DimerLookX(r1, r2, n_base, k, m0, m1, sb, m2) {
    let z = 0, f = 0, w = 0, x = 0, y = 0, l = 0, h = 0, j = 0, e = 0, x1 = 0, y1 = 0, fn = 0, x0 = 0, x2 = 0, y0 = 0, b1 = 0, b2 = 0;
    let n = -1;
    let z1 = "", z2 = "";
    let l1 = r1.length;
    let l2 = r2.length;
    let lmin = l2;
    let dmin = n_base + 1;
    r1 = r1.toLowerCase();
    r2 = r2.toLowerCase();

    if (dmin < 4) { dmin = 4; }
    if (dmin > 9) { dmin = 9; }
    if (lmin > l1) {
        lmin = l1;
        let temp = l1;
        l1 = l2;
        l2 = temp;
        z1 = r2;
        z2 = ReverseDNA(r1);
    } else {
        z1 = r1;
        z2 = ReverseDNA(r2);
    }
    if (lmin < 6) {
        return 0;
    }

    const lx = new Array(l1 + l2 + l2).fill(0);
    const u = new Array(l1 - 4).fill(0);
    const px = Array.from({ length: 5 }, () => Array.from({ length: 5 }, () => Array.from({ length: 5 }, () => Array.from({ length: 5 }, () => Array.from({ length: 5 }, () => Array(2).fill(0))))));
    const ax = new Array(5).fill(0);

    ax[0] = dx[z1.charCodeAt(0)];
    ax[1] = dx[z1.charCodeAt(1)];
    ax[2] = dx[z1.charCodeAt(2)];
    ax[3] = dx[z1.charCodeAt(3)];

    for (j = 4; j < l1; j++) {
        ax[4] = dx[z1.charCodeAt(j)];
        px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]++;
        ax[0] = ax[1];
        ax[1] = ax[2];
        ax[2] = ax[3];
        ax[3] = ax[4];
    }

    ax[0] = dx[z1.charCodeAt(0)];
    ax[1] = dx[z1.charCodeAt(1)];
    ax[2] = dx[z1.charCodeAt(2)];
    ax[3] = dx[z1.charCodeAt(3)];

    for (j = 4; j < l1; j++) {
        ax[4] = dx[z1.charCodeAt(j)];
        y = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];
        if (y > 0) {
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] = fn;
            u[fn] = j - 4;
            fn += y;
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0] = -1;
        } else {
            u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] - y] = j - 4;
            px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0]--;
        }
        ax[0] = ax[1];
        ax[1] = ax[2];
        ax[2] = ax[3];
        ax[3] = ax[4];
    }

    ax[0] = dx[cdna[z2.charCodeAt(0)]];
    ax[1] = dx[cdna[z2.charCodeAt(1)]];
    ax[2] = dx[cdna[z2.charCodeAt(2)]];
    ax[3] = dx[cdna[z2.charCodeAt(3)]];
    for (y = 0; y < l2 - 4; y++) {
        ax[4] = dx[cdna[z2.charCodeAt(y + 4)]];
        f = px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][0];

        if (f < 0) {
            for (h = 1; h <= -f; h++) {
                x = u[px[ax[0]][ax[1]][ax[2]][ax[3]][ax[4]][1] + h - 1];
                w = x - y;
                x2 = l2 + w - 1;

                if (w < 1) {
                    x0 = 0;
                    y0 = -w;
                } else {
                    x0 = w;
                    y0 = 0;
                }
                z = l2 - y0;
                if (l1 - x0 < z) {
                    z = l1 - x0;
                }

                if (lx[x2] === 0) {
                    l = 5;
                    e = 0;
                    for (j = 0; j <= k + k; j++) {
                        if (x + 6 + j > l1 || y + 6 + j > l2) {
                            break;
                        }
                        const value = aDNA[z1.charCodeAt(x + 5 + j)][z2.charCodeAt(y + 5 + j)];
                        if (value > 24) {
                            l++;
                        } else {
                            if (e === 1) {
                                break;
                            }
                            e = 1;
                            if (value > 24) {
                                l++;
                            }
                        }
                    }
                    x1 = x;
                    y1 = y;
                    for (j = 0; j <= k + k; j++) {
                        if (x - j - 1 < 0 || y - j - 1 < 0) {
                            break;
                        }
                        const value = aDNA[z1.charCodeAt(x - j - 1)][z2.charCodeAt(y - j - 1)];
                        if (value > 24) {
                            l++;
                        } else {
                            if (e === 1) {
                                break;
                            }
                            e = 1;
                        }
                        x1--;
                        y1--;
                    }

                    if (l > dmin) {
                        lx[x2] = l;
                        n++;
                        m0[n] = 0;
                        m1[n] = "<-" + z2 + "-5";
                        sb[n] = "  ";
                        m2[n] = "5-" + z1 + "->";
                        let value = [];
                        if (w > 0) {
                            value = new Array(w).fill(' ');
                            m0[n] = 1;
                            m1[n] = value.join('') + "<-" + z2 + "-5";
                            sb[n] = value.join('') + '  ';
                            m2[n] = "5-" + z1 + "->";
                        }
                        if (w < 0) {
                            m0[n] = 2;
                            m1[n] = "<-" + z2 + "-5";
                            value = new Array(2 - w).fill(' ');
                            sb[n] = value.join('');
                            value = new Array(-w).fill(' ');
                            m2[n] = value.join('') + "5-" + z1 + "->";
                        }
                        value = new Array(z).fill(' ');
                        for (j = 0; j < z; j++) {
                            const h1 = aDNA[z2.charCodeAt(y0 + j)][z1.charCodeAt(x0 + j)];
                            if (h1 > 24) {
                                value[j] = ':';
                            }
                            if (h1 > 70) {
                                value[j] = '|';
                            }
                        }
                        sb[n] = sb[n] + value.join('');
                    }
                }
            }
        }
        ax[0] = ax[1];
        ax[1] = ax[2];
        ax[2] = ax[3];
        ax[3] = ax[4];
    }
    return n;
}


function DimersLook(s1, s2, n_base, f_base, m0, m1, sb, m2) {
    // s1 = "atatatatatatatatatatatatat";
    // s2 = "atatatatatatatatatatatatat";
    /**
           <-cgcgcgctaaggggcgc-5
             ||||||
5-cgcggggaatcgcgcgc->     
     */
    //  s1 = "cgcggggaatcgcgcgc";
    //  s2 = "cgcggggaatcgcgcgc";

    s1 = s1.toLowerCase();
    s2 = s2.toLowerCase();

    let n = -1;
    let l1 = s1.length;
    let l2 = s2.length;
    let r1 = s1;
    let r2 = s2;
    if (l1 > l2) {
        r1 = s2;
        r2 = s1
    }
    l1 = r1.length;
    l2 = r2.length;
    let h = n_base;
    let f = f_base - 1;
    if (h < 3) { h = 3; }
    if (h > f) { f = h; }

    const z2 = ComplementDNA(r2);
    const z1 = ReverseDNA(r2);
    const d1 = new Array(1 + l1 - h).fill("");
    let lx = new Array(l1 + l2).fill(0);

    for (let y = 0; y <= l1 - h; y++) {
        d1[y] = r1.substring(y, y + h);
    }
    for (let y = 0; y < r2.length - h + 1; y++) {
        let g = z2.substring(y, y + h);
        let x = -1;
        do {
            x = d1.indexOf(g, x + 1);
            if (x > -1) {
                let w = x - y;
                let x2 = l2 + w - 1;
                if (lx[x2] === 0) {
                    let x0 = 0;
                    let y0 = 0;
                    let z = 0;
                    if (w < 1) {
                        x0 = 0;
                        y0 = -w;
                    } else {
                        x0 = w;
                        y0 = 0;
                    }
                    z = l2 - y0;
                    if (l1 - x0 < z) {
                        z = l1 - x0;
                    }

                    let l = h;
                    let k = Math.min(l1 - x - h, l2 - y - h);
                    for (j = h; j < h + k; j++) {
                        if (r1.charAt(x + j) === z2.charAt(y + j)) {
                            l++;
                        } else { break; }
                    }
                    if (l > f) {
                        lx[x2] = l;
                        n++;
                        m0[n] = 0;
                        m1[n] = "<-" + z1 + "-5";
                        sb[n] = "  ";
                        m2[n] = "5-" + r1 + "->";
                        let value = [];
                        if (w > 0) {
                            value = new Array(w).fill(' ');
                            m0[n] = 1;
                            m1[n] = value.join('') + "<-" + z1 + "-5";
                            sb[n] = value.join('') + '  ';
                            m2[n] = "5-" + r1 + "->";
                        }
                        if (w < 0) {
                            m0[n] = 2;
                            m1[n] = "<-" + z1 + "-5";
                            value = new Array(2 - w).fill(' ');
                            sb[n] = value.join('');
                            value = new Array(-w).fill(' ');
                            m2[n] = value.join('') + "5-" + r1 + "->";
                        }
                        value = new Array(z).fill(' ');
                        for (j = 0; j < z; j++) {
                            const h1 = aDNA[z1.charCodeAt(y0 + j)][r1.charCodeAt(x0 + j)];
                            if (h1 > 49) {
                                value[j] = ':';
                            }
                            if (h1 > 70) {
                                value[j] = '|';
                            }
                        }
                        sb[n] = sb[n] + value.join('');
                    }
                }
            }
        } while (x > -1);
    }
    return n;
}


function DimerLook2(r1, r2, n_base) {
    //    r1="tgacactagtgcagttgct";
    //    r2="cgctgacactagtgcaatagac";
    let h = n_base; //  dimer3 = >5
    if (h < 6) {
        h = 6;
    }
    if (h > 8) {
        h = 8;
    }
    return qDimer(r1, r2, h);
}

function qDimer(p1, p2, h) {
    const dimer5 = 8;  // >7
    const dimer3 = h;  // =7
    const n = 5;       // initial kmer=5

    let r1 = p1;
    let r2 = p2;
    if (p1.length > p2.length) {
        r1 = p2;
        r2 = p1;
    }

    /*
                     x[i] x0      
    5-gaaccgaagcgtacagtcgcc-> r1
                |    ||||||         
             <-cctctacagcggaaacctca-5        
                     y    y0       
                             */

    const r1Map = new Map();
    for (let y = 0; y <= r1.length - n; y++) {
        const sub = r1.substring(y, y + n);
        if (!r1Map.has(sub)) {
            r1Map.set(sub, new Set());
        }
        r1Map.get(sub).add(y);
    }

    const z2 = ComplementDNA(r2);
    for (let y = 0; y <= r2.length - n; y++) {
        const sub = z2.substring(y, y + n);
        if (r1Map.has(sub)) {
            const x2 = Array.from(r1Map.get(sub));

            for (let i = 0; i < x2.length; i++) {
                const x0 = r1.length - x2[i] - n;
                const y0 = z2.length - y - n;
                let l = n;
                let e = 0;

                if (x0 === 1 && y0 > 1) {
                    for (let k = 0; k < Math.min(y, x2[i]); k++) {
                        const s1 = r1.substring(x2[i] - k, x2[i] - k + 1);
                        const s2 = z2.substring(y - k, y - k + 1);
                        if (s1 === s2) {
                            l++;
                        } else {
                            if (++e > 1) {
                                break;
                            }
                        }
                    }
                } else {
                    for (let k = 0; k < Math.min(x0, y0); k++) {
                        const s1 = r1.substring(x2[i] + k + n, x2[i] + k + n + 1);
                        const s2 = z2.substring(y + k + n, y + k + n + 1);
                        if (s1 === s2) {
                            l++;
                        } else {
                            if (++e > 1) {
                                break;
                            }
                        }
                    }
                }

                if (l > dimer3) {
                    if (l > dimer5) {
                        return 1; // 5'-end dimer
                    }
                    if (y === 0 && x2[i] > 0) {
                        return 1; // 3'-end dimer
                    }
                    if (x0 === 1 && y0 > 1) {
                        return 1; // 3'-end dimer
                    }
                }
            }
        }
    }
    return 0;
}

function quickDimer(r1, r2 = null) {
    let h = 7;
    if (r2 === null) {
        return qDimer(r1, r1, h);
    }
    return qDimer(r1, r2, h);
}