class LowComplexitySequence {
  constructor(seq, telomer, minlenseq) {
    let b = Array.from(seq, char => char.charCodeAt(0));
    this.mapb = new Array(seq.length).fill(0);
    b = this.simpleRepeatsMasking(b, telomer, minlenseq);
    const b1 = [];
    const h = b.length - 1;
    for (let i = 0; i <= h; i++) {
      if (b[i] === 5) {
        let e;
        for (e = i + 1; e <= h; e++) {
          if (b[e] < 5) {
            break;
          }
        }
        if (minlenseq < e - i) {
          b1.push([i, e - i]);
        }
        i = e;
      }
    }

    this.ibloks = new Array(2 * b1.length).fill(0);
    let e = -1;
    this.totalssr = 0;
    for (let i = 0; i < b1.length; i++) {
      const z7 = b1[i];
      this.ibloks[++e] = z7[0];
      this.ibloks[++e] = z7[1];
      this.totalssr += z7[1];
      for (let j = z7[0]; j < z7[0] + z7[1]; j++) {
        this.mapb[j] = 4;
      }
    }
  }

  simpleRepeatsMasking(b, telomer, ssrlen) {
    const lmer = 17;
    const l = b.length;
    if (l < lmer + 3) {
      return b;
    }
    const Kmax = telomer;
    const dx2 = new Array(128).fill(4);
    dx2[97] = 0;  // a
    dx2[101] = 0; // e
    dx2[108] = 1; // l
    dx2[116] = 1; // t
    dx2[117] = 1; // u
    dx2[99] = 2;  // c
    dx2[102] = 2; // f
    dx2[103] = 3; // g
    dx2[105] = 3; // i
    dx2[106] = 3; // j
    dx2[119] = 4; // w
    dx2[98] = 4;  // b
    dx2[104] = 4; // h
    dx2[121] = 4; // y
    dx2[109] = 4; // m
    dx2[115] = 4; // s
    dx2[107] = 4; // k
    dx2[114] = 4; // r
    dx2[118] = 4; // v
    dx2[100] = 4; // d
    dx2[110] = 4; // n

    const v2 = Array.from({ length: 5 }, () => Array(5).fill(0));
    const v3 = Array.from({ length: 5 }, () => Array.from({ length: 5 }, () => Array(5).fill(0)));
    const k2 = new Array(l - lmer).fill(0);

    for (let i = 0; i < l; i++) {
      b[i] = dx2[b[i]];
    }

    for (let i = 0; i < lmer; i++) {
      const n1 = b[i];
      const n2 = b[i + 1];
      v2[n1][n2]++;
      if (v2[n1][n2] === 1) {
        k2[0]++;
      }
    }

    for (let i = 0; i < lmer - 1; i++) {
      const n1 = b[i];
      const n2 = b[i + 1];
      const n3 = b[i + 2];
      v3[n1][n2][n3]++;
      if (v3[n1][n2][n3] === 1) {
        k2[0]++;
      }
    }

    for (let i = 1; i < l - lmer; i++) {
      k2[i] = k2[i - 1];
      const n1 = b[i - 1];
      const n2 = b[i];
      const n3 = b[i + 1];
      v2[n1][n2]--;
      if (v2[n1][n2] === 0) {
        k2[i]--;
      }

      v3[n1][n2][n3]--;
      if (v3[n1][n2][n3] === 0) {
        k2[i]--;
      }

      const n1Next = b[i + lmer - 2];
      const n2Next = b[i + lmer - 1];
      const n3Next = b[i + lmer];
      v2[n2Next][n3Next]++;
      if (v2[n2Next][n3Next] === 1) {
        k2[i]++;
      }

      v3[n1Next][n2Next][n3Next]++;
      if (v3[n1Next][n2Next][n3Next] === 1) {
        k2[i]++;
      }
    }

    for (let i = 0; i < l - lmer; i++) {
      if (k2[i] < Kmax && b[i] < 4) {
        let x = i;
        let u = 0;
        let e = i + 1;
        while (e < l - lmer) {
          if (k2[e] > Kmax || b[e] === 4) {
            u++;
            if (u > lmer) {
              break;
            }
          }
          e++;
        }
        i = e - 1;
        if (e - x > ssrlen) {
          if (x > 0) {
            x--;
          }
          for (let h = x; h < e + lmer; h++) {
            b[h] = 5;
          }
        }
      }
    }
    return b;
  }

  getBlocks() {
    return [this.ibloks];
  }

  getMapBytes() {
    return this.mapb;
  }

  getIntBlocks() {
    return this.ibloks;
  }

  getTotalRepeats() {
    return this.totalssr;
  }
}