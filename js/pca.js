'use strict';
/**
 * Polymerase cycling assembly (PCA) — building a gene from overlapping oligos.
 *
 * The construct is tiled by oligonucleotides that alternate between the top and
 * the bottom strand and overlap their neighbours. In the first ~30 cycles the
 * overlaps anneal and the polymerase extends each oligo along its partner, so
 * fragments grow towards full length; end primers are then added for ~23 more
 * cycles to amplify the finished construct away from the incomplete pieces
 * (Stemmer et al., Gene 1995; Wikipedia, Polymerase cycling assembly).
 *
 * The geometry, written out once because everything below depends on it. With
 * overlaps O(1)..O(k-1), where O(i) covers [s(i), e(i)):
 *
 *     oligo 1     = [0,      e(1))
 *     oligo i     = [s(i-1), e(i))        2 <= i <= k-1
 *     oligo k     = [s(k-1), L)
 *
 * and the strands alternate, oligo 1 on top. Note what that does to the 3'
 * ends: a top oligo ends at the RIGHT edge of its overlap, a bottom oligo at
 * the LEFT edge — so for every overlap BOTH partners have their 3' end inside
 * it, which is exactly what has to happen for each to prime on the other. It
 * also means both ends of an overlap are 3' termini, and both want a G or C.
 *
 * What the design actually optimises: junctions are placed evenly and each
 * overlap is then chosen, within a window, for a melting temperature close to
 * the target — because one cycling program has to serve every junction at once.
 * Oligo length follows from the junctions rather than being imposed. This is a
 * greedy per-junction search, not the global optimisation DNAWorks runs over the
 * whole set; it is fast and predictable, and the resulting Tm spread is reported
 * so a bad junction is visible rather than hidden.
 *
 * There is deliberately no GC window on the overlap. At a fixed Tm target GC is
 * not free to vary independently — a low-GC overlap is selected longer and a
 * high-GC one shorter — so a GC range would charge for the same property twice.
 * What GC would be a proxy for is scored directly instead: hairpin stem, longest
 * homopolymer run, linguistic complexity, and a G or C at each 3' terminus.
 *
 * Reused from js/dna.js:      ComplementDNA, CG, LingComplexity6
 * Reused from js/DNAtm.js:    DNA_Tm, at the fixed conditions used site-wide
 *                             (55 mM Na+/K+, 1 mM Mg2+, 0.2 uM oligo) so every
 *                             Tm here is comparable with the other tools
 * Reused from js/dimers4.js:  not needed — the hairpin scan below is local
 */

/* eslint-disable no-unused-vars */

var PCA_SALT_M = 0.055, PCA_MG_M = 0.001, PCA_OLIGO_UM = 0.2;

function pcaTm(s) { return DNA_Tm(s, PCA_SALT_M, PCA_MG_M, PCA_OLIGO_UM); }
function pcaRC(s) { return ComplementDNA(String(s || '').toLowerCase()); }
function pcaUp(s) { return String(s == null ? '' : s).toUpperCase(); }
function pcaFmt(x, d) { return isFinite(x) ? Number(x).toFixed(d == null ? 1 : d) : 'n/a'; }
function pcaPad(s, w) { s = String(s == null ? '' : s); return s + ' '.repeat(Math.max(0, w - s.length)); }

/** Longest homopolymer run. */
function pcaMaxRun(s) {
  var best = 0, run = 0, prev = '';
  for (var i = 0; i < s.length; i++) {
    if (s.charAt(i) === prev) run++; else { run = 1; prev = s.charAt(i); }
    if (run > best) best = run;
  }
  return best;
}

var PCA_COMP1 = { a: 't', t: 'a', g: 'c', c: 'g' };

/** Longest hairpin stem inside `s`, anchored on an exact k-mer. */
function pcaHairpin(s, minStem, minLoop) {
  minStem = minStem || 4; minLoop = minLoop == null ? 3 : minLoop;
  var n = s.length, best = 0;
  for (var i = 0; i + minStem <= n; i++) {
    var rc = pcaRC(s.substr(i, minStem));
    var at = s.indexOf(rc, i + minStem + minLoop);
    if (at === -1) continue;
    var stem = minStem, a = i + minStem, b = at - 1;
    while (b - a >= minLoop && PCA_COMP1[s.charAt(a)] === s.charAt(b)) { stem++; a++; b--; }
    if (stem > best) best = stem;
  }
  return best;
}

// ── Input ────────────────────────────────────────────────────────────────────

/** First record only; keep A/C/G/T, fold U to T, drop everything else. */
function pcaParse(raw) {
  var lines = String(raw || '').split(/\r?\n/);
  var name = 'construct', body = [], seenBody = false;
  for (var i = 0; i < lines.length; i++) {
    if (/^\s*>/.test(lines[i])) {
      if (seenBody) break;
      name = lines[i].replace(/^\s*>\s*/, '').split(/[\s|]/)[0] || 'construct';
      continue;
    }
    body.push(lines[i]);
    if (lines[i].trim()) seenBody = true;
  }
  var joined = body.join('').toLowerCase().replace(/\s/g, '');
  var seq = joined.replace(/u/g, 't').replace(/[^acgt]/g, '');
  return { name: name, seq: seq, dropped: joined.length - seq.length };
}

// ── Overlap scoring ──────────────────────────────────────────────────────────

/**
 * How good is this stretch as an annealing overlap? Lower is better; the Tm
 * term dominates, because one cycling program has to serve every junction.
 */
function pcaScoreOverlap(ov, opt, ctx) {
  var tm = pcaTm(ov), gc = CG(ov), why = [];
  var cost = Math.abs(tm - opt.tmTarget) * 2.0;

  // A junction on repeated or low-complexity sequence is the one failure that
  // does not announce itself: the two oligos meeting there anneal just as
  // happily to the other copy, and the construct assembles in the wrong order.
  if (ctx && ctx.mask && ctx.at != null) {
    var mc = pcaMaskedCount(ctx.mask, ctx.at, ctx.at + ov.length);
    if (mc) { cost += 6 + mc * 1.5; why.push(mc + ' masked nt'); }
  }
  // The precise version of the same worry, inside this block: an overlap whose
  // sequence occurs twice cannot tell its partners apart at all.
  if (ctx && ctx.seq) {
    var firstAt = ctx.seq.indexOf(ov);
    if (firstAt !== -1 && ctx.seq.indexOf(ov, firstAt + 1) !== -1) {
      cost += 1000; why.push('sequence occurs twice in this block');
    } else if (ctx.rc && ctx.rc.indexOf(ov) !== -1) {
      cost += 1000; why.push('sequence also present inverted in this block');
    }
  }

  // Both ends of an overlap are 3' termini (see the header), so both want a
  // G or C to prime cleanly.
  var first = ov.charAt(0), last = ov.charAt(ov.length - 1);
  if (first !== 'g' && first !== 'c') { cost += 1.2; why.push('no G/C at the left edge'); }
  if (last !== 'g' && last !== 'c') { cost += 1.2; why.push('no G/C at the right edge'); }

  var run = pcaMaxRun(ov);
  if (run >= 5) { cost += 2.5 * (run - 4); why.push(run + '-base run'); }

  var hp = pcaHairpin(ov, 4, 3);
  if (hp >= 4) { cost += 2.0 * (hp - 3); why.push(hp + '-bp hairpin'); }

  if (typeof LingComplexity6 === 'function') {
    var lc = LingComplexity6(ov);
    if (lc < 60) { cost += (60 - lc) * 0.10; why.push('complexity ' + lc + '%'); }
  }
  return { cost: cost, tm: tm, gc: gc, hairpin: hp, run: run, why: why };
}

// ── Junction placement ───────────────────────────────────────────────────────

/**
 * Tile [0, L) with oligos. Junctions start evenly spaced — which keeps the
 * oligos even and, unlike walking left to right, cannot strand a runt at the
 * far end — and each is then searched within `window` for its best overlap.
 */
function pcaPlan(seq, opt, ctx) {
  var L = seq.length;
  ctx = ctx || {};
  if (L < opt.oligoMin) {
    return { error: 'Only ' + L + ' bases — shorter than one oligo. Order it as a single oligo instead.' };
  }
  if (L <= opt.oligoMax) {
    // Two oligos is the smallest thing that is still an assembly.
    var len2 = Math.min(opt.ovlMax, Math.max(opt.ovlMin, Math.round(L / 2)));
    var s2 = Math.floor((L - len2) / 2);
    var ov2 = { s: s2, e: s2 + len2, len: len2,
                sc: pcaScoreOverlap(seq.substr(s2, len2), opt, { mask: ctx.mask, at: s2, seq: seq, rc: ctx.rc }) };
    return {
      k: 2, overlaps: [ov2],
      oligos: [pcaMakeOligo(seq, 1, 0, ov2.e, true), pcaMakeOligo(seq, 2, ov2.s, L, false)]
    };
  }

  var ovlAvg = (opt.ovlMin + opt.ovlMax) / 2;
  var aimLen = (opt.oligoMin + opt.oligoMax) / 2;            // the range is the aim
  var step = Math.max(1, aimLen - ovlAvg);                   // advance per oligo
  var k = Math.max(2, Math.round((L - ovlAvg) / step));      // number of oligos
  // Keep every oligo inside the length range if the arithmetic allows it: the
  // mean oligo length is (L + (k-1)*overlap) / k.
  var guard = 0;
  while (k > 2 && (L + (k - 1) * opt.ovlMax) / k > opt.oligoMax && guard++ < 500) k++;
  guard = 0;
  while (k > 2 && (L + (k - 1) * opt.ovlMin) / k < opt.oligoMin && guard++ < 500) k--;

  var overlaps = [], prevEnd = 0, i;
  var spacing = (L - ovlAvg) / k;
  for (i = 1; i <= k - 1; i++) {
    // Equal oligo lengths put junction i STARTING at i*(L-v)/k. Centring it on
    // i*L/k instead leaves both terminal oligos about half an overlap short,
    // because each of them has a junction on one side only.
    var ideal = Math.round(i * (L - ovlAvg) / k);
    // The window is where we would LIKE to look; where we may actually look
    // starts after the previous overlap ends. Taking the lower bound from
    // prevEnd rather than clamping the window globally is what keeps a wide
    // window usable: two neighbours can still be pushed together, but the later
    // one is never squeezed out of having any legal position at all.
    // Everything after this junction still has to fit, so that bound is applied
    // to the search rather than left to reject candidates one junction too late.
    var maxS = L - ((k - i) * opt.oligoMin - (k - 1 - i) * opt.ovlMax);
    var maxEnd = (i < k - 1)
      ? L - ((k - i - 1) * opt.oligoMin - (k - 2 - i) * opt.ovlMax)
      : L;
    var lo = Math.max(0, Math.min(ideal - opt.window, L), prevEnd);
    var hi = Math.min(Math.max(lo, ideal + opt.window), maxS);
    if (hi < lo) {
      return { error: 'Junction ' + i + ' of ' + (k - 1) + ' has nowhere legal to go: the overlaps ahead of ' +
               'it need ' + (L - maxS) + ' nt and only ' + (L - lo) + ' nt are left. Shorten the overlap range, ' +
               'or raise the oligo length so the oligos advance further each step.' };
    }
    var best = null;
    for (var len = opt.ovlMin; len <= opt.ovlMax; len++) {
      for (var s = lo; s <= hi; s++) {
        if (s + len > L) break;
        if (s + len > maxEnd) break;      // would leave the next junction homeless
        var sc = pcaScoreOverlap(seq.substr(s, len), opt, { mask: ctx.mask, at: s, seq: seq, rc: ctx.rc });
        // A junction that wanders costs oligo-length evenness, so it has to earn
        // its place: 0.15/nt is worth about 0.075 degC of Tm improvement per nt.
        var cost = sc.cost + Math.abs(s - ideal) * 0.15 + Math.abs(len - ovlAvg) * 0.05;
        if (!best || cost < best.cost) best = { s: s, e: s + len, len: len, cost: cost, sc: sc };
      }
    }
    if (!best) {
      return { error: 'Could not place junction ' + i + ' of ' + (k - 1) + '. The oligos advance about ' +
               Math.round(spacing) + ' nt each, which leaves no room for an overlap of ' + opt.ovlMin +
               '-' + opt.ovlMax + ' nt. Shorten the overlap, or raise the oligo length.' };
    }
    overlaps.push(best);
    prevEnd = best.e;
  }

  // Oligos follow from the junctions.
  var oligos = [];
  for (i = 0; i < k; i++) {
    var a = i === 0 ? 0 : overlaps[i - 1].s;
    var b = i === k - 1 ? L : overlaps[i].e;
    oligos.push(pcaMakeOligo(seq, i + 1, a, b, (i % 2) === 0));
  }
  return { oligos: oligos, overlaps: overlaps, k: k, spacing: spacing };
}

function pcaMakeOligo(seq, idx, a, b, top) {
  var foot = seq.substring(a, b);
  return {
    idx: idx, start: a, end: b, top: top,
    seq: top ? foot : pcaRC(foot),
    len: b - a, gc: CG(foot), tm: pcaTm(foot),
    // a top oligo's 3' end is at its right edge, a bottom oligo's at its left
    threeAt: top ? b : a
  };
}

/** The oligo set must rebuild the target exactly — orientation included. */
function pcaVerify(seq, oligos) {
  var buf = new Array(seq.length).fill(null), i, j;
  for (i = 0; i < oligos.length; i++) {
    var o = oligos[i];
    var top = o.top ? o.seq : pcaRC(o.seq);            // back to top-strand sense
    if (top.length !== o.end - o.start) {
      return { ok: false, why: 'oligo ' + o.idx + ' length does not match its footprint' };
    }
    for (j = 0; j < top.length; j++) {
      var p = o.start + j;
      if (buf[p] !== null && buf[p] !== top.charAt(j)) {
        return { ok: false, why: 'oligos disagree at position ' + (p + 1) };
      }
      buf[p] = top.charAt(j);
    }
  }
  for (i = 0; i < seq.length; i++) {
    if (buf[i] === null) return { ok: false, why: 'position ' + (i + 1) + ' is not covered by any oligo' };
    if (buf[i] !== seq.charAt(i)) {
      return { ok: false, why: 'position ' + (i + 1) + ' rebuilds as ' + buf[i] + ', not ' + seq.charAt(i) };
    }
  }
  return { ok: true };
}

// ── Misannealing ─────────────────────────────────────────────────────────────

/**
 * The failure that ruins a PCA reaction is an oligo priming on the wrong
 * partner, so every 3'-terminal seed is searched against the whole construct,
 * both strands, and any match other than its own junction is reported.
 */
function pcaCrossCheck(seq, oligos, opt) {
  var hits = [], rc = pcaRC(seq), L = seq.length, seed = opt.seed;
  for (var i = 0; i < oligos.length; i++) {
    var o = oligos[i];
    if (o.seq.length < seed) continue;
    var tail = o.seq.slice(-seed);                     // 3'-terminal seed, 5'->3'
    // Where this oligo is meant to prime — its own 3' end, in each frame.
    var ownTop = o.top ? o.threeAt - seed : -1;
    var ownBot = o.top ? -1 : L - o.threeAt - seed;
    var at;
    for (at = seq.indexOf(tail); at !== -1; at = seq.indexOf(tail, at + 1)) {
      if (o.top && at === ownTop) continue;
      hits.push({ oligo: o.idx, strand: '+', at: at, seed: tail });
    }
    for (at = rc.indexOf(tail); at !== -1; at = rc.indexOf(tail, at + 1)) {
      if (!o.top && at === ownBot) continue;
      hits.push({ oligo: o.idx, strand: '-', at: L - at - seed, seed: tail });
    }
  }
  return hits;
}

// ── End primers ──────────────────────────────────────────────────────────────

/** Grow inwards from a construct end until the Tm target is met. */
function pcaEndPrimer(seq, fromStart, opt) {
  var best = null;
  for (var len = opt.primerMin; len <= opt.primerMax; len++) {
    if (len > seq.length) break;
    var s = fromStart ? seq.substr(0, len) : pcaRC(seq.slice(-len));
    var tm = pcaTm(s), gc = CG(s);
    var cost = Math.abs(tm - opt.primerTm) * 2 + (gc < 35 || gc > 65 ? 4 : 0) +
               (pcaMaxRun(s) > 4 ? 3 : 0);
    var last = s.charAt(s.length - 1);
    if (last !== 'g' && last !== 'c') cost += 1;
    if (!best || cost < best.cost) best = { seq: s, len: len, tm: tm, gc: gc, cost: cost };
  }
  return best;
}

function pcaTa(tm1, tm2, ampLen) { return Math.min(tm1, tm2) + Math.log(ampLen); }

// ── Repeats and low complexity ───────────────────────────────────────────────
//
// The other design modules mask the template before they place anything on it,
// using RepeatMask (repeated 15-mers, either strand) and LowComplexitySequence
// (SSR and telomere-like blocks); PCA masks the same way, for the same reason
// and with a sharper consequence. A primer that sits on a repeat merely risks a
// spurious product. A PCA *junction* that sits on a repeat is worse: the two
// oligos meeting there can just as well anneal to the other copy, and the
// assembly scrambles rather than fails, which is far harder to notice.
//
// So the mask is used twice. Junction overlaps are kept off masked positions,
// and — the part that actually solves the problem rather than reporting it —
// the construct is CUT between the copies of a long repeat, so that no single
// assembly reaction ever contains two copies of the same sequence. Inside one
// block the tiling is then unambiguous, and the blocks are joined afterwards
// through their shared unique ends, the way Gibson assembly joins fragments.

// The same settings pcr.js uses, so a stretch masked there is masked here:
// Kmax = 14 for telomere-like repeats (11 would target SSR), 21 nt minimum block.
var PCA_LC_KMAX = 14, PCA_LC_MINLEN = 21;

/** Per-position mask, in the convention pcr.js and qpcr.js use: >0 = suspect. */
function pcaMask(seq, opt) {
  var n = seq.length, msk = new Array(n).fill(0), i, j;
  var sources = [], missing = [];

  if (opt.maskRepeats) {
    if (typeof RepeatMask === 'function') {
      var rm = RepeatMask(seq);
      var covered = 0;
      for (i = 0; i < n; i++) if (rm[i] > 0) { msk[i]++; covered++; }
      sources.push({ what: 'repeated 15-mers', covered: covered });
    } else {
      // A filter that is asked for and cannot run must not read as a filter that
      // found nothing: that is indistinguishable from clean sequence, and it is
      // how a missing <script> tag hides.
      missing.push({ what: 'repeated 15-mers', needs: 'js/dna.js' });
    }
  }
  if (opt.maskLowComplexity) {
    if (typeof LowComplexitySequence === 'function') {
      var lc = new LowComplexitySequence(seq, PCA_LC_KMAX, PCA_LC_MINLEN);
      // getBlocks()[0] is a flat [start, length, start, length, ...] list.
      var blocks = lc.getBlocks()[0] || [], lcCovered = 0;
      for (i = 0; i < blocks.length; i += 2) {
        for (j = blocks[i]; j < blocks[i] + blocks[i + 1] && j < n; j++) { msk[j]++; lcCovered++; }
      }
      sources.push({ what: 'low complexity (SSR, telomere-like)', covered: lcCovered });
    } else {
      missing.push({ what: 'low complexity (SSR, telomere-like)', needs: 'js/LowComplexitySequence.js' });
    }
  }
  var total = 0;
  for (i = 0; i < n; i++) if (msk[i] > 0) total++;
  return { msk: msk, masked: total, sources: sources, missing: missing };
}

/** Is any position of [a, b) masked? */
function pcaMaskedAny(mask, a, b) {
  if (!mask) return false;
  for (var i = a; i < b; i++) if (mask.msk[i] > 0) return true;
  return false;
}
/** How many positions of [a, b) are masked. */
function pcaMaskedCount(mask, a, b) {
  if (!mask) return 0;
  var c = 0;
  for (var i = a; i < b; i++) if (mask.msk[i] > 0) c++;
  return c;
}

/**
 * Long repeats as families of copies, which is what block splitting needs and
 * what a per-position mask cannot give: the mask says "this is repeated
 * somewhere", the families say "these two stretches are copies of each other".
 *
 * Seeded on exact k-mers of length `minLen` on both strands, extended to the
 * maximal match, then merged so that overlapping copies collapse into one.
 */
function pcaFindRepeats(seq, minLen) {
  var n = seq.length, k = Math.max(8, minLen | 0), i;
  if (n < 2 * k) return [];
  var rc = pcaRC(seq);
  var index = new Map(), pairs = [];

  // forward k-mers
  for (i = 0; i + k <= n; i++) {
    var key = seq.substr(i, k);
    var at = index.get(key);
    if (at === undefined) index.set(key, [i]); else at.push(i);
  }
  // direct repeats: any k-mer seen more than once
  var CAP = 20000, tooMany = false;
  index.forEach(function (positions) {
    if (positions.length < 2 || pairs.length > CAP) { if (positions.length >= 2) tooMany = true; return; }
    // Chain consecutive copies only: n copies give n-1 seeds, not n squared.
    for (var a = 0; a < positions.length - 1; a++) {
      pairs.push({ i: positions[a], j: positions[a + 1], len: k, inverted: false });
    }
  });
  // inverted repeats: a forward k-mer that also occurs on the other strand
  for (i = 0; i + k <= n && pairs.length <= CAP; i++) {
    var revKey = rc.substr(i, k);
    var hit = index.get(revKey);
    if (!hit || hit.length > 50) continue;      // ignore a k-mer that is everywhere
    var fwdOfRev = n - i - k;                 // where this piece sits on the top strand
    for (var h = 0; h < hit.length; h++) {
      if (hit[h] >= fwdOfRev) continue;       // count each inverted pair once
      pairs.push({ i: hit[h], j: fwdOfRev, len: k, inverted: true });
    }
  }
  if (!pairs.length) return [];

  // Extend each seed to its maximal match, then drop the ones swallowed by a
  // longer neighbour so a 60 nt repeat is one family, not forty overlapping ones.
  var ext = [];
  pairs.forEach(function (p) {
    var a = p.i, b = p.j, len = p.len;
    if (!p.inverted) {
      // Direct repeat: seq[a+x] === seq[b+x]. Both copies grow the same way.
      while (a + len < n && b + len < n && seq.charAt(a + len) === seq.charAt(b + len)) len++;
      while (a > 0 && b > 0 && seq.charAt(a - 1) === seq.charAt(b - 1)) { a--; b--; len++; }
    } else {
      // Inverted repeat: seq[a+x] pairs with the COMPLEMENT of seq[b+len-1-x],
      // so as one copy grows to the right the other grows to the left.
      while (a + len < n && b > 0 && PCA_COMP1[seq.charAt(a + len)] === seq.charAt(b - 1)) { b--; len++; }
      while (a > 0 && b + len < n && PCA_COMP1[seq.charAt(a - 1)] === seq.charAt(b + len)) { a--; len++; }
    }
    var tandem = false;
    if (b < a + len) { len = b - a; tandem = true; }   // copies run into each other
    if (len < k) return;
    ext.push({ a: a, b: b, len: len, inverted: p.inverted, tandem: tandem });
  });
  ext.sort(function (x, y) { return y.len - x.len || x.a - y.a; });

  var kept = [];
  ext.forEach(function (e) {
    for (var i2 = 0; i2 < kept.length; i2++) {
      var q = kept[i2];
      // already inside a longer family with the same offset between copies
      if (e.a >= q.a && e.a + e.len <= q.a + q.len && e.b >= q.b && e.b + e.len <= q.b + q.len) return;
    }
    kept.push(e);
  });

  // Group copies into families: two entries belong together when they share a copy.
  var families = [];
  kept.forEach(function (e) {
    var fam = null;
    for (var f = 0; f < families.length; f++) {
      var copies = families[f].copies;
      for (var c = 0; c < copies.length; c++) {
        if (pcaOverlaps(copies[c], e.a, e.len) || pcaOverlaps(copies[c], e.b, e.len)) { fam = families[f]; break; }
      }
      if (fam) break;
    }
    if (!fam) { fam = { copies: [], len: e.len, inverted: e.inverted, tandem: e.tandem }; families.push(fam); }
    pcaAddCopy(fam.copies, e.a, e.len);
    pcaAddCopy(fam.copies, e.b, e.len);
    if (e.len > fam.len) fam.len = e.len;
    if (e.inverted) fam.inverted = true;
    if (e.tandem) fam.tandem = true;
  });

  families.forEach(function (f) { f.copies.sort(function (x, y) { return x.start - y.start; }); });
  families.sort(function (x, y) { return y.len - x.len || x.copies[0].start - y.copies[0].start; });
  var out = families.filter(function (f) { return f.copies.length > 1; });
  if (tooMany || pairs.length > CAP) out.truncated = true;
  return out;
}

function pcaOverlaps(copy, start, len) {
  return copy.start < start + len && start < copy.start + copy.len;
}
function pcaAddCopy(copies, start, len) {
  for (var i = 0; i < copies.length; i++) {
    if (pcaOverlaps(copies[i], start, len)) {
      var end = Math.max(copies[i].start + copies[i].len, start + len);
      copies[i].start = Math.min(copies[i].start, start);
      copies[i].len = end - copies[i].start;
      return;
    }
  }
  copies.push({ start: start, len: len });
}

/**
 * Where the construct may be cut into blocks. A cut is only legal where the
 * whole shared stretch — `blockOvl` bases that both neighbouring blocks will
 * carry — is unique, because that stretch is what joins the two blocks
 * afterwards. Joining on a repeat would put the ambiguity back in.
 */
function pcaCutAllowed(mask, seq, at, opt) {
  var half = Math.ceil(opt.blockOvl / 2);
  var a = Math.max(0, at - half), b = Math.min(seq.length, at + half);
  return pcaMaskedCount(mask, a, b) === 0;
}

/**
 * Split into blocks so that (a) no block is longer than blockMax and (b) no
 * block contains two copies of the same repeat family — the copies are pushed
 * into separate assembly reactions, which is the only way the tiling inside a
 * reaction can be unambiguous. Blocks share `blockOvl` bases so they can be
 * fused afterwards.
 */
function pcaBlocksByRepeats(seq, opt, mask, families) {
  var L = seq.length, notes = [];
  if (!opt.blockMax && !(opt.splitOnRepeats && families.length)) {
    return { blocks: [{ start: 0, end: L, seq: seq }], notes: notes };
  }

  // For each position, which families have a copy starting at or covering it.
  var forced = [];              // positions where a cut is REQUIRED before continuing
  if (opt.splitOnRepeats) {
    families.forEach(function (f, fi) {
      for (var c = 1; c < f.copies.length; c++) {
        // a cut has to fall between the end of copy c-1 and the start of copy c
        var lo = f.copies[c - 1].start + f.copies[c - 1].len;
        var hi = f.copies[c].start;
        if (hi > lo) forced.push({ lo: lo, hi: hi, fam: fi, len: f.len });
      }
    });
    forced.sort(function (a, b) { return a.hi - b.hi; });
  }

  var blocks = [], start = 0, guard = 0;
  while (start < L && guard++ < 1000) {
    var hardMax = opt.blockMax ? Math.min(L, start + opt.blockMax) : L;
    // the first repeat gap that must be cut inside this block
    var need = null;
    for (var i = 0; i < forced.length; i++) {
      var f = forced[i];
      if (f.lo <= start) continue;                       // its earlier copy is behind us
      if (f.hi <= start) continue;
      if (f.hi > hardMax && f.lo > hardMax) continue;    // beyond this block anyway
      // the second copy would land in this block: cut before it
      if (f.hi <= hardMax) { need = f; break; }
    }

    var target = need ? Math.min(need.hi, hardMax) : hardMax;
    if (target >= L) { blocks.push({ start: start, end: L, seq: seq.substring(start, L) }); break; }

    // Pick the cut: as late as possible, but in unique sequence, and not so
    // early that the block becomes a stub.
    var minEnd = Math.min(L, start + Math.max(opt.oligoMin * 2, opt.blockOvl * 2));
    var cut = -1, from = Math.max(minEnd, need ? need.lo : minEnd);
    for (var p = target; p >= from; p--) {
      if (pcaCutAllowed(mask, seq, p, opt)) { cut = p; break; }
    }
    if (cut < 0) {
      // Nowhere unique to cut. Take the least-masked position rather than give up.
      var bestP = target, bestScore = Infinity;
      for (var q = target; q >= from; q--) {
        var sc = pcaMaskedCount(mask, Math.max(0, q - opt.blockOvl / 2), Math.min(L, q + opt.blockOvl / 2));
        if (sc < bestScore) { bestScore = sc; bestP = q; }
      }
      cut = bestP;
      notes.push('The join at position ' + cut + ' had to be placed on masked sequence — there was no unique ' +
                 'window of ' + opt.blockOvl + ' bp to cut in. That block junction is the weak point of the assembly.');
    }
    if (need && cut < need.lo) {
      notes.push('A repeat of ' + need.len + ' bp could not be split between blocks: there is no room to cut ' +
                 'between its copies. Both copies stay in one reaction, so the tiling there is ambiguous.');
    }
    blocks.push({ start: start, end: cut, seq: seq.substring(start, cut) });
    start = Math.max(cut - opt.blockOvl, blocks[blocks.length - 1].start + 1);
  }
  if (start < L && blocks.length) {
    var last = blocks[blocks.length - 1];
    if (last.end < L) blocks.push({ start: start, end: L, seq: seq.substring(start, L) });
  }
  return { blocks: blocks, notes: notes };
}

/** Which repeat families still have two copies inside one block. */
function pcaFamiliesInBlock(families, blk) {
  var out = [];
  families.forEach(function (f, fi) {
    var inside = f.copies.filter(function (c) {
      return c.start >= blk.start && c.start + c.len <= blk.end;
    });
    if (inside.length > 1) out.push({ fam: fi, len: f.len, copies: inside });
  });
  return out;
}

// ── Blocks ───────────────────────────────────────────────────────────────────

/** Split a long construct into separately assembled blocks that share ends. */
function pcaBlocks(seq, opt) {
  var L = seq.length;
  if (!opt.blockMax || L <= opt.blockMax) return [{ start: 0, end: L, seq: seq }];
  var n = Math.ceil((L - opt.blockOvl) / (opt.blockMax - opt.blockOvl));
  var span = Math.ceil((L + (n - 1) * opt.blockOvl) / n);
  var out = [], at = 0;
  for (var i = 0; i < n && at < L; i++) {
    var end = i === n - 1 ? L : Math.min(L, at + span);
    out.push({ start: at, end: end, seq: seq.substring(at, end) });
    if (end >= L) break;
    at = end - opt.blockOvl;
  }
  return out;
}

// ── Rendering ────────────────────────────────────────────────────────────────

function pcaRenderOligos(target, plans, opt, warn) {
  var t = [];
  t.push('POLYMERASE CYCLING ASSEMBLY — oligonucleotide set');
  t.push('');
  t.push('Construct : ' + target.name + '   ' + target.seq.length + ' bp   GC ' + pcaFmt(CG(target.seq), 1) + '%');
  t.push('Design    : oligos ' + opt.oligoMin + '–' + opt.oligoMax + ' nt, overlaps ' +
         opt.ovlMin + '–' + opt.ovlMax + ' nt at Tm ' + pcaFmt(opt.tmTarget, 0) + ' °C');
  t.push('Blocks    : ' + plans.length +
         (plans.length > 1 ? '  (assembled separately, sharing ' + opt.blockOvl + ' bp ends)' : ''));
  var totalOligos = plans.reduce(function (s, p) { return s + p.plan.oligos.length; }, 0);
  t.push('Oligos    : ' + totalOligos + ' in total');
  t.push('');

  plans.forEach(function (blk, bi) {
    var pl = blk.plan;
    if (plans.length > 1) {
      t.push('═'.repeat(104));
      t.push('BLOCK ' + (bi + 1) + ' of ' + plans.length + '   construct ' + (blk.start + 1) + '..' + blk.end +
             '   ' + blk.seq.length + ' bp   ' + pl.oligos.length + ' oligos');
      t.push('═'.repeat(104));
    }
    var head = pcaPad('name', 12) + pcaPad('strand', 11) + pcaPad('position', 14) +
               pcaPad('len', 5) + pcaPad('Tm full', 9) + pcaPad('GC%', 6) + 'sequence (5\'->3\')';
    t.push(head);
    t.push('-'.repeat(100));
    pl.oligos.forEach(function (o) {
      var name = (plans.length > 1 ? 'B' + (bi + 1) + '-' : '') + (o.top ? 'F' : 'R') + o.idx;
      t.push(pcaPad(name, 12) + pcaPad(o.top ? 'top (+)' : 'bottom (−)', 11) +
             pcaPad((blk.start + o.start + 1) + '-' + (blk.start + o.end), 14) +
             pcaPad(o.len, 5) + pcaPad(pcaFmt(o.tm, 1), 9) + pcaPad(pcaFmt(o.gc, 0), 6) +
             pcaUp(o.seq));
    });
    t.push('');
    t.push('Junction overlaps — these all have to anneal under one program:');
    t.push('  ' + pcaPad('#', 5) + pcaPad('position', 14) + pcaPad('len', 5) + pcaPad('Tm', 7) +
           pcaPad('GC%', 6) + pcaPad('sequence', 28) + 'notes');
    pl.overlaps.forEach(function (ov, j) {
      t.push('  ' + pcaPad(j + 1, 5) + pcaPad((blk.start + ov.s + 1) + '-' + (blk.start + ov.e), 14) +
             pcaPad(ov.len, 5) + pcaPad(pcaFmt(ov.sc.tm, 1), 7) + pcaPad(pcaFmt(ov.sc.gc, 0), 6) +
             pcaPad(pcaUp(blk.seq.substring(ov.s, ov.e)), 28) + (ov.sc.why.join(', ') || '–'));
    });
    var tms = pl.overlaps.map(function (o) { return o.sc.tm; });
    if (tms.length) {
      var lo = Math.min.apply(null, tms), hi = Math.max.apply(null, tms);
      t.push('');
      t.push('  Overlap Tm ' + pcaFmt(lo, 1) + '–' + pcaFmt(hi, 1) + ' °C, spread ' + pcaFmt(hi - lo, 1) + ' °C' +
             (hi - lo > opt.tmSpread ? '   ← wider than the ' + opt.tmSpread + ' °C you asked for' : ''));
      if (hi - lo > opt.tmSpread) {
        warn.push('Block ' + (bi + 1) + ': the overlap Tm spread is ' + pcaFmt(hi - lo, 1) +
                  ' °C. Widen the overlap length range so the search has more room, or split the construct into more blocks.');
      }
    }
    var bad = pl.oligos.filter(function (o) { return o.len < opt.oligoMin || o.len > opt.oligoMax; });
    if (bad.length) {
      t.push('  ! ' + bad.length + ' oligo(s) outside ' + opt.oligoMin + '–' + opt.oligoMax + ' nt: ' +
             bad.map(function (o) { return '#' + o.idx + ' (' + o.len + ')'; }).join(', '));
      warn.push('Block ' + (bi + 1) + ': ' + bad.length +
                ' oligo(s) are outside the requested length range — the junction arithmetic could not do better at this construct length.');
    }
    t.push('');
    t.push('  Reassembly check: ' + (pl.verified.ok
      ? 'the set rebuilds this block exactly, orientation included.'
      : 'FAILED — ' + pl.verified.why));
    if (!pl.verified.ok) warn.push('Block ' + (bi + 1) + ' does not rebuild the target — do not order these oligos.');
    t.push('');
  });

  t.push('Naming: F = top strand, R = bottom strand, numbered along the construct.');
  t.push('"Tm full" is the whole oligo against its complement. What the cycling program follows is');
  t.push('the junction overlap Tm in the second table.');
  t.push('Order them unpurified at the smallest scale — PCA needs very little of each, and it is the');
  t.push('end primers, not oligo purity, that select full-length product at the end.');
  return t.join('\n');
}

function pcaRenderMap(target, plans, opt) {
  var t = [], WIDTH = 92;
  t.push('ASSEMBLY MAP');
  t.push('');
  t.push('Each oligo drawn against the construct: > runs 5\'->3\' on the top strand, < on the bottom.');
  t.push('Where two neighbours meet they share an overlap — that is where they prime on each other.');
  t.push('');
  plans.forEach(function (blk, bi) {
    var pl = blk.plan, L = blk.seq.length;
    var col = function (p) { return Math.round(p * WIDTH / L); };
    if (plans.length > 1) {
      t.push('── Block ' + (bi + 1) + '   construct ' + (blk.start + 1) + '..' + blk.end + ', ' + L + ' bp ──');
    }
    t.push('  ' + pcaPad('', 7) + '1' + ' '.repeat(Math.max(0, WIDTH - 2)) + String(L));
    pl.oligos.forEach(function (o) {
      var a = col(o.start), b = Math.max(a + 1, col(o.end));
      var name = (o.top ? 'F' : 'R') + o.idx;
      t.push('  ' + pcaPad(name, 7) + ' '.repeat(a) + (o.top ? '>' : '<').repeat(b - a));
    });
    t.push('');
  });
  t.push('The map is ' + WIDTH + ' columns wide whatever the construct length, so a column is a relative');
  t.push('measure — read exact coordinates from the oligo table.');
  return t.join('\n');
}

function pcaRenderProtocol(target, plans, opt, warn) {
  var t = [];
  t.push('END PRIMERS AND CYCLING');
  t.push('');
  plans.forEach(function (blk, bi) {
    var f = pcaEndPrimer(blk.seq, true, opt), r = pcaEndPrimer(blk.seq, false, opt);
    if (!f || !r) { t.push('Block ' + (bi + 1) + ': too short for an end primer.'); return; }
    var pre = plans.length > 1 ? 'B' + (bi + 1) + '-' : '';
    t.push((plans.length > 1 ? 'Block ' + (bi + 1) + ' — ' : '') +
           'end primers, added after the assembly cycles to amplify full-length product:');
    t.push('  ' + pcaPad('name', 12) + pcaPad('len', 5) + pcaPad('Tm', 7) + pcaPad('GC%', 6) + 'sequence (5\'->3\')');
    t.push('  ' + pcaPad(pre + 'PCA-F', 12) + pcaPad(f.len, 5) + pcaPad(pcaFmt(f.tm, 1), 7) +
           pcaPad(pcaFmt(f.gc, 0), 6) + pcaUp(f.seq));
    t.push('  ' + pcaPad(pre + 'PCA-R', 12) + pcaPad(r.len, 5) + pcaPad(pcaFmt(r.tm, 1), 7) +
           pcaPad(pcaFmt(r.gc, 0), 6) + pcaUp(r.seq));
    t.push('  product ' + blk.seq.length + ' bp   Ta ≈ ' + pcaFmt(pcaTa(f.tm, r.tm, blk.seq.length), 0) +
           ' °C   ΔTm ' + pcaFmt(Math.abs(f.tm - r.tm), 1) + ' °C');
    if (Math.abs(f.tm - r.tm) > 4) {
      warn.push('Block ' + (bi + 1) + ': the end primers differ by ' + pcaFmt(Math.abs(f.tm - r.tm), 1) +
                ' °C. The construct ends fix these sequences, so add a few bases of vector or tag sequence if the pair has to be matched.');
    }
    t.push('');
  });

  var longest = Math.max.apply(null, plans.map(function (b) { return b.seq.length; }));
  var ext = Math.max(30, Math.ceil(longest / 1000 * 60));
  var allOv = [].concat.apply([], plans.map(function (b) {
    return b.plan.overlaps.map(function (o) { return o.sc.tm; });
  }));
  var ovMin = allOv.length ? Math.min.apply(null, allOv) : opt.tmTarget;

  t.push('─'.repeat(96));
  t.push('SUGGESTED CYCLING — the reference protocol, to be optimised for your construct');
  t.push('─'.repeat(96));
  t.push('');
  t.push('  1. Pool the oligos equimolar, to roughly 20–50 nM each. Nothing else goes in yet:');
  t.push('     at this stage the oligos are both template and primer.');
  t.push('  2. Assembly, ~30 cycles:  95 °C 30 s → ' + pcaFmt(Math.max(45, ovMin - 5), 0) +
         ' °C 30 s → 72 °C ' + ext + ' s');
  t.push('     The annealing step follows the WEAKEST overlap (' + pcaFmt(ovMin, 1) + ' °C) less a few degrees,');
  t.push('     because every junction has to anneal, not just the average one.');
  t.push('  3. Add the end primers to ~0.4 µM each and run ~23 further cycles at the Ta above.');
  t.push('     This is what pulls full-length product away from the incomplete fragments.');
  t.push('  4. Gel-purify the band at the expected size before cloning, and sequence it.');
  t.push('');
  t.push('  Use a proof-reading polymerase — but note what it can and cannot do here: every base of');
  t.push('  the construct comes from a synthetic oligo, and a proof-reader will not correct an oligo');
  t.push('  that was synthesised wrong. Synthesis error, not polymerase error, is the thing to expect,');
  t.push('  which is why the product is sequenced rather than trusted.');
  if (plans.length > 1) {
    t.push('');
    t.push('  Blocks: assemble and amplify each separately, gel-purify, then join them. Neighbouring');
    t.push('  blocks share ' + opt.blockOvl + ' bp, so they can be fused by overlap-extension PCR using the outermost');
    t.push('  primers, or by Gibson assembly — the Gibson Assembly tool designs that route.');
  }
  return t.join('\n');
}

function pcaRenderRepeats(target, plans, opt, mask, families, warn) {
  var t = [], L = target.seq.length;
  t.push('REPEATS, MASKING AND HOW THE BLOCKS WERE CUT');
  t.push('');

  // ── masking ───────────────────────────────────────────────────────────────
  t.push('Masking — the same two filters the primer-design modules use:');
  (mask.missing || []).forEach(function (m) {
    t.push('  ' + pcaPad(m.what, 40) + 'NOT RUN — ' + m.needs + ' is not loaded on this page');
    warn.push('The ' + m.what + ' filter was switched on but could not run: ' + m.needs +
              ' is not loaded on this page, so nothing was masked by it. This is a page fault, not a ' +
              'property of your sequence — the design below was made without that filter.');
  });
  if (!mask.sources.length && !(mask.missing || []).length) {
    t.push('  Both masks are switched off. Junctions are being placed without knowing which parts of the');
    t.push('  construct are repeated, which is exactly how an assembly ends up scrambled rather than failed.');
    warn.push('Repeat and low-complexity masking are both off, so junctions may sit on repeated sequence.');
  } else {
    mask.sources.forEach(function (s) {
      t.push('  ' + pcaPad(s.what, 40) + pcaPad(s.covered + ' nt', 12) + pcaFmt(100 * s.covered / L, 1) + '% of the construct');
    });
    t.push('  ' + pcaPad('masked in total', 40) + pcaPad(mask.masked + ' nt', 12) + pcaFmt(100 * mask.masked / L, 1) + '%');
    t.push('');
    t.push('  A masked position is one a junction should not sit on. Where a junction had to be placed on');
    t.push('  masked sequence anyway, its row in the Oligos tab says so.');
  }
  t.push('');

  // ── repeat families ───────────────────────────────────────────────────────
  t.push('─'.repeat(96));
  t.push('Long repeats — copies of the same sequence, at least ' + opt.repeatMin + ' bp');
  t.push('─'.repeat(96));
  if (!families.length) {
    t.push('  None. Nothing in this construct repeats itself over ' + opt.repeatMin + ' bp or more, so the tiling is');
    t.push('  unambiguous wherever the junctions fall.');
  } else {
    t.push('  ' + pcaPad('#', 5) + pcaPad('length', 8) + pcaPad('kind', 10) + pcaPad('copies at', 40) + 'separated?');
    families.forEach(function (f, fi) {
      var where = f.copies.map(function (c) { return (c.start + 1) + '-' + (c.start + c.len); }).join(', ');
      // which block each copy landed in
      var inBlocks = f.copies.map(function (c) {
        for (var b = 0; b < plans.length; b++) {
          if (c.start >= plans[b].start && c.start + c.len <= plans[b].end) return b + 1;
        }
        return '?';
      });
      var uniqueBlocks = inBlocks.filter(function (v, i2) { return inBlocks.indexOf(v) === i2; });
      var separated = uniqueBlocks.length === inBlocks.length && inBlocks.indexOf('?') === -1;
      t.push('  ' + pcaPad(fi + 1, 5) + pcaPad(f.len + ' bp', 8) +
             pcaPad(f.tandem ? 'tandem' : (f.inverted ? 'inverted' : 'direct'), 10) +
             pcaPad(where.length > 38 ? where.slice(0, 37) + '…' : where, 40) +
             (separated ? 'yes — blocks ' + inBlocks.join(', ')
                        : (f.tandem ? 'no — tandem, cannot be cut apart'
                                    : 'NO — blocks ' + inBlocks.join(', ') + ', same reaction')));
    });
    t.push('');
    t.push('');
    if (families.some(function (f) { return f.tandem; })) {
      t.push('  A tandem family is copies lying back to back, with nothing between them to cut at. Those');
      t.push('  cannot be separated by any block boundary — the junctions are simply kept off them, and the');
      t.push('  oligos spanning them are flagged in Warnings.');
      t.push('');
    }
  t.push('  "Separated" is the whole point. Ambiguity only matters inside one reaction: if the two copies');
    t.push('  of a repeat are assembled in different tubes, no oligo in either tube can anneal to the wrong');
    t.push('  one, and the tiling is unambiguous even though the finished construct still contains both.');
  }
  t.push('');

  // ── blocks ────────────────────────────────────────────────────────────────
  t.push('─'.repeat(96));
  t.push('Blocks');
  t.push('─'.repeat(96));
  t.push('  ' + pcaPad('#', 5) + pcaPad('span', 16) + pcaPad('length', 9) + pcaPad('oligos', 8) +
         pcaPad('masked', 9) + 'repeats still inside');
  plans.forEach(function (blk, bi) {
    var inside = pcaFamiliesInBlock(families, blk);
    var mc = pcaMaskedCount(mask, blk.start, blk.end);
    t.push('  ' + pcaPad(bi + 1, 5) + pcaPad((blk.start + 1) + '-' + blk.end, 16) +
           pcaPad(blk.seq.length + ' bp', 9) + pcaPad(blk.plan.oligos.length, 8) +
           pcaPad(pcaFmt(100 * mc / blk.seq.length, 0) + '%', 9) +
           (inside.length ? inside.map(function (x) { return '#' + (x.fam + 1) + ' (' + x.len + ' bp ×' + x.copies.length + ')'; }).join(', ')
                          : 'none'));
    if (inside.length) {
      var f0 = families[inside[0].fam];
      warn.push('Block ' + (bi + 1) + ' still contains ' + inside[0].copies.length + ' copies of repeat #' +
                (inside[0].fam + 1) + ' (' + inside[0].len + ' bp). The oligos there cannot tell the copies apart. ' +
                (f0 && f0.tandem
                  ? 'They are tandem, so no block boundary can separate them — this stretch has to be built another way, or edited if synonymous codons allow it.'
                  : 'Shorten the block length so the cut falls between them, or edit the construct if synonymous codons allow it.'));
    }
  });

  if (plans.length > 1) {
    t.push('');
    t.push('  Joins: neighbouring blocks share ' + opt.blockOvl + ' bp, and the cut was placed so that shared stretch is');
    t.push('  unique sequence. Assemble and amplify each block separately, gel-purify, then fuse them through');
    t.push('  those shared ends — by overlap-extension PCR with the outermost primers, or by Gibson assembly.');
    t.push('  The Gibson Assembly tool designs that second step from the same block boundaries.');
    t.push('');
    t.push('  ' + pcaPad('join', 8) + pcaPad('shared stretch', 18) + pcaPad('masked?', 10) + 'sequence');
    for (var b = 1; b < plans.length; b++) {
      var a0 = plans[b].start, a1 = Math.min(plans[b - 1].end, plans[b].start + opt.blockOvl);
      var shared = target.seq.substring(a0, a1);
      var dirty = pcaMaskedCount(mask, a0, a1);
      t.push('  ' + pcaPad(b + '→' + (b + 1), 8) + pcaPad((a0 + 1) + '-' + a1, 18) +
             pcaPad(dirty ? dirty + ' nt' : 'clean', 10) +
             pcaUp(shared.length > 46 ? shared.slice(0, 45) + '…' : shared));
    }
  } else {
    t.push('');
    t.push('  One block: the whole construct is assembled in a single reaction, so there is nothing to join.');
  }
  return t.join('\n');
}

function pcaRenderWarnings(target, plans, opt, extra) {
  var t = [], w = [];
  t.push('WARNINGS AND LIMITS');
  t.push('');

  plans.forEach(function (blk, bi) {
    var hits = pcaCrossCheck(blk.seq, blk.plan.oligos, opt);
    if (hits.length) {
      var which = [];
      hits.forEach(function (h) { if (which.indexOf(h.oligo) === -1) which.push(h.oligo); });
      w.push('Block ' + (bi + 1) + ': ' + hits.length + ' case(s) where an oligo\'s 3\'-terminal ' + opt.seed +
             ' nt also match somewhere other than its own junction, so it can prime in the wrong place. Oligos ' +
             which.map(function (n) { return '#' + n; }).join(', ') +
             '. This is the classic PCA failure — a repeat in the construct that the tiling cannot avoid.');
    }
  });

  plans.forEach(function (blk, bi) {
    blk.plan.oligos.forEach(function (o) {
      var hp = pcaHairpin(o.seq, 5, 3);
      if (hp >= 6) {
        w.push('Block ' + (bi + 1) + ' oligo #' + o.idx + ' folds back on itself over ' + hp +
               ' bp — it may not be available to anneal.');
      }
    });
  });

  if (target.dropped) {
    w.push(target.dropped + ' non-ACGT character(s) were dropped from the input. Ambiguity codes cannot be built as one oligo set — pick a base first.');
  }
  if (target.seq.length < 100) {
    w.push('The construct is only ' + target.seq.length + ' bp. Below about 100 bp a single pair of long overlapping oligos, annealed and extended once, is simpler than cycling assembly.');
  }
  if (plans.length === 1 && target.seq.length > 800) {
    w.push('Assembling ' + target.seq.length + ' bp in one reaction is ambitious: PCA is most reliable up to roughly 500–600 bp per reaction. Set a block length so the construct is built in pieces and joined.');
  }
  (extra || []).forEach(function (m) { w.push(m); });

  if (!w.length) t.push('• Nothing flagged.');
  else w.forEach(function (m) { t.push('• ' + m); });

  t.push('');
  t.push('Scope and limits:');
  t.push('  1) Junctions are placed evenly and each overlap is then chosen, within a window, for a Tm');
  t.push('     near the target. That is a greedy per-junction search, not a global optimisation of the');
  t.push('     whole set the way DNAWorks anneals over it, so a construct with a hostile stretch may');
  t.push('     have a better solution this does not find. The Tm spread is reported so you can judge.');
  t.push('  2) Tm is the site-wide nearest-neighbour value (55 mM Na⁺/K⁺, 1 mM Mg²⁺, 0.2 µM oligo), so');
  t.push('     it is comparable with every other tool here. A real assembly mix is not at those');
  t.push('     conditions — treat the numbers as a design-time reference, not a prediction.');
  t.push('  3) The misannealing check is an exact ' + opt.seed + '-nt match of each 3\' end against the construct.');
  t.push('     It cannot see a mismatched interaction that still primes, and it knows nothing about the');
  t.push('     vector the construct is going into.');
  t.push('  4) Nothing here optimises codons. If the construct is a coding sequence and you are free to');
  t.push('     choose synonymous codons, a codon-shuffling design can flatten the Tm spread and break up');
  t.push('     repeats that this tool can only report.');
  t.push('  5) Repeats are handled by separating them, not by solving them. Two copies of the same');
  t.push('     sequence cannot be told apart by any choice of junctions, so the construct is cut between');
  t.push('     them and each copy assembled in its own reaction; the Repeats tab says which families were');
  t.push('     separated and which, for want of room to cut, were not.');
  return t.join('\n');
}

// ── Orchestration ────────────────────────────────────────────────────────────

function pcaNum(id, dflt) {
  var el = document.getElementById(id);
  if (!el) return dflt;
  var v = parseFloat(el.value);
  return isFinite(v) ? v : dflt;
}
function pcaVal(id, dflt) { var el = document.getElementById(id); return el && el.value != null ? el.value : dflt; }
function pcaChk(id, dflt) { var el = document.getElementById(id); return el ? !!el.checked : !!dflt; }

function pcaSetPanes(map) {
  for (var id in map) { var el = document.getElementById(id); if (el) el.value = map[id]; }
}
function pcaPending(note) {
  pcaSetPanes({ analysisResult1: note, analysisResult2: '', analysisResult3: '',
                analysisResult4: '', analysisResult5: '' });
}

function runPCA() {
  var raw = pcaVal('inputText', '');
  if (!String(raw).trim()) {
    pcaPending('Paste the sequence you want to build, 5\'->3\', with or without a FASTA header.\n' +
               'The tool returns the overlapping oligos that assemble it, the end primers that\n' +
               'amplify the finished construct, and a cycling program.\n\n' +
               'Click "Design oligo set" to run, or "Load Example" for a worked case.');
    return;
  }
  var target = pcaParse(raw);
  if (!target.seq) { pcaPending('[input] No A/C/G/T bases found. Paste DNA, with or without a FASTA header.'); return; }

  var opt = {
    oligoMin: pcaNum('oligo_min', 40), oligoMax: pcaNum('oligo_max', 60),
    ovlMin: pcaNum('ovl_min', 16), ovlMax: pcaNum('ovl_max', 26),
    tmTarget: pcaNum('tm_target', 60), tmSpread: pcaNum('tm_spread', 5),
    window: Math.max(0, Math.min(30, pcaNum('junction_window', 8))),
    blockMax: Math.max(0, pcaNum('block_max', 500)),
    blockOvl: Math.max(20, pcaNum('block_ovl', 30)),
    primerTm: pcaNum('primer_tm', 62), primerMin: pcaNum('primer_min', 18), primerMax: pcaNum('primer_max', 30),
    seed: Math.max(6, Math.min(20, pcaNum('cross_seed', 10))),
    maskRepeats: pcaChk('mask_repeats', true),
    maskLowComplexity: pcaChk('mask_lowcomplexity', true),
    splitOnRepeats: pcaChk('split_repeats', true),
    // The page offers 50-600 with a default of 120, and the report quotes this
    // value back ("copies of the same sequence, at least N bp"), so the ceiling
    // has to be the page's: a lower one silently made every setting above it
    // mean the same thing. The floor stays permissive rather than matching the
    // spinner, because a smaller value typed by hand is legitimate - it only
    // finds shorter families - and pcaFindRepeats floors the seed at 8 anyway.
    // A longer seed is cheaper, not dearer: it yields fewer pairs to extend.
    repeatMin: Math.max(12, Math.min(600, pcaNum('repeat_min', 120)))
  };
  if (opt.oligoMin > opt.oligoMax) { pcaPending('[input] The oligo length range is inverted.'); return; }
  if (opt.ovlMin > opt.ovlMax) { pcaPending('[input] The overlap length range is inverted.'); return; }
  if (opt.ovlMax >= opt.oligoMin) {
    pcaPending('[input] The longest overlap (' + opt.ovlMax + ' nt) must be shorter than the shortest oligo (' +
               opt.oligoMin + ' nt), or an oligo would be nothing but overlap.');
    return;
  }
  if (2 * opt.ovlMax > opt.oligoMax) {
    pcaPending('[input] An oligo in the middle of the construct carries an overlap at BOTH ends, so the ' +
               'longest oligo (' + opt.oligoMax + ' nt) has to hold two overlaps of up to ' + opt.ovlMax +
               ' nt. Shorten the longest overlap to ' + Math.floor(opt.oligoMax / 2) + ' nt or less, or raise the oligo length.');
    return;
  }
  if (opt.blockMax && opt.blockMax <= opt.blockOvl) {
    pcaPending('[input] The block length must be longer than the block overlap.');
    return;
  }

  try {
    var warn = [], failed = null;
    // Mask the whole construct first, find its long repeats, and let both decide
    // where the blocks are cut — a repeat is only harmless once its copies are
    // in different reactions.
    var mask = pcaMask(target.seq, opt);
    var families = (opt.splitOnRepeats || opt.maskRepeats) ? pcaFindRepeats(target.seq, opt.repeatMin) : [];
    var split = pcaBlocksByRepeats(target.seq, opt, mask, families);
    var blocks = split.blocks;
    split.notes.forEach(function (n) { warn.push(n); });

    blocks.forEach(function (b) {
      if (failed) return;
      // Junctions are placed against the block's OWN mask: a repeat whose other
      // copy sits in a different block cannot confuse this reaction.
      var bmask = pcaMask(b.seq, opt);
      var pl = pcaPlan(b.seq, opt, { mask: bmask, rc: pcaRC(b.seq) });
      if (pl.error) { failed = pl.error; return; }
      pl.verified = pcaVerify(b.seq, pl.oligos);
      b.plan = pl;
      b.mask = bmask;
    });
    if (failed) { pcaPending('[design] ' + failed); return; }

    var tab1 = pcaRenderOligos(target, blocks, opt, warn);
    var tab2 = pcaRenderMap(target, blocks, opt);
    var tab3 = pcaRenderProtocol(target, blocks, opt, warn);
    var tab4 = pcaRenderRepeats(target, blocks, opt, mask, families, warn);
    var tab5 = pcaRenderWarnings(target, blocks, opt, warn);
    pcaSetPanes({ analysisResult1: tab1, analysisResult2: tab2, analysisResult3: tab3,
                  analysisResult4: tab4, analysisResult5: tab5 });
  } catch (err) {
    pcaPending('[error] ' + (err && err.message ? err.message : err));
    if (window.console) console.error('PCA error:', err);
  }
}

// ── Example, pending state, reset ────────────────────────────────────────────

var PCA_EXAMPLES = {
  cds: {
    note: 'Example loaded — a 456 bp coding fragment to build from scratch.\nClick "Design oligo set" to run it.',
    seq:
      '>demo_cds — 456 bp coding fragment for de novo assembly\n' +
      'atgagcaaaggcgaagaactgtttaccggcgttgtgccgattctggttgaactggatggcgatgtgaacggccataaattcagcgtgagcggcgaaggcgaa\n' +
      'ggcgatgcgacctatggcaaactgaccctgaaatttatttgcaccaccggcaaactgccggtgccgtggccgaccctggtgaccaccctgacctatggcgtg\n' +
      'cagtgctttagccgctatccggatcatatgaaacagcatgattttttcaaaagcgcgatgccggaaggctatgtgcaggaacgcaccatttttttcaaagat\n' +
      'gatggcaactataaaacccgcgcggaagtgaaatttgaaggcgataccctggtgaaccgcattgaactgaaaggcattgattttaaagaagatggcaacatt\n' +
      'ctgggccataaactggaatataactataacagccataacgtgtatatt'
  }
};

function pcaClearFeedback() {
  if (window.SeqValidate) window.SeqValidate.clear('inputFeedback');
}

/** A design belongs to the sequence it was run on, so an edit retires it. */
function pcaAwaitRun() {
  var el = document.getElementById('inputText');
  if (!el || !el.value.trim()) { runPCA(); return; }
  pcaPending('Click "Design oligo set" to run the design for the sequence in the box.');
}

function loadExamplePCA(name) {
  var ex = PCA_EXAMPLES[name] || PCA_EXAMPLES.cds;
  var el = document.getElementById('inputText');
  if (!el) return;
  el.value = ex.seq;
  pcaClearFeedback();
  pcaPending(ex.note);
}

function clearPCA() {
  var el = document.getElementById('inputText');
  if (el) el.value = '';
  pcaClearFeedback();
  runPCA();
}

if (typeof window !== 'undefined') {
  window.runPCA = runPCA;
  window.pcaAwaitRun = pcaAwaitRun;
  window.pcaPending = pcaPending;
  window.loadExamplePCA = loadExamplePCA;
  window.clearPCA = clearPCA;
}
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    pcaParse: pcaParse, pcaPlan: pcaPlan, pcaVerify: pcaVerify, pcaBlocks: pcaBlocks,
    pcaCrossCheck: pcaCrossCheck, pcaEndPrimer: pcaEndPrimer, pcaScoreOverlap: pcaScoreOverlap,
    pcaMask: pcaMask, pcaFindRepeats: pcaFindRepeats, pcaBlocksByRepeats: pcaBlocksByRepeats
  };
}
