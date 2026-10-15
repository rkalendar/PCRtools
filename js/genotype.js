'use strict';
/**
 * Insertion and knock-in genotyping — assays that answer whether an engineered
 * cassette is present at a locus, and in how many copies.
 *
 * Two assays are designed from one context sequence, and they are complementary:
 *
 *   KNOCK-IN assay — one primer sits INSIDE the cassette, its partner in the
 *   flanking genome. The product therefore exists only on the engineered allele;
 *   wild-type DNA has nothing for the cassette primer to anneal to.
 *
 *       Fwd ---------->              Probe ---->
 *     5'------[==== cassette ====]----------------------------------
 *                                            <---------- Rev
 *
 *   EMPTY-SITE assay — both primers sit in the flanks, spanning the insertion
 *   point, so the product is short on the wild-type allele and longer by the
 *   whole cassette on the engineered one. With a large cassette that second
 *   product does not amplify under qPCR conditions, which is what makes the
 *   assay wild-type specific.
 *
 *       Fwd ---------->   |   <---------- Rev        | = insertion point
 *     5'-------------------------------------------------
 *
 * Run together, the pair reads zygosity: wild type gives the empty-site signal
 * only, a homozygote the knock-in signal only, a heterozygote both. That is the
 * whole point of designing them as a pair rather than one at a time.
 *
 * What this tool does that a plain primer picker cannot: it knows which bases
 * belong to the cassette and which to the genome. Every rule below that mentions
 * a position — "the primer must not be able to bind wild type", "the probe must
 * straddle the insertion point", "the 3' end lying in the flank must be weak" —
 * needs that distinction, and it is derived once, either from the markup or by
 * comparing the two alleles.
 *
 * Tm conventions. Full-length Tm is the site-wide nearest-neighbour value
 * (55 mM Na+/K+, 1 mM Mg2+, 0.2 uM oligo) so it is comparable with every other
 * tool here. Tm_3'end is the Tm of the 3'-terminal GT_TM3_WINDOW nucleotides;
 * it is reported because assay-design rules for this class of work are written
 * in it, and a long primer with a weak 3' end behaves very differently from a
 * short one with the same overall Tm.
 *
 * Reused from js/dna.js                     ComplementDNA, CG, LingComplexity6,
 *                                           RepeatMask
 * Reused from js/DNAtm.js                   DNA_Tm
 * Reused from js/LowComplexitySequence.js   LowComplexitySequence
 */

/* eslint-disable no-unused-vars */

var GT_SALT_M = 0.055, GT_MG_M = 0.001, GT_OLIGO_UM = 0.2;

/* The 3'-terminal window whose Tm is reported as Tm_3'end. Ten nucleotides is
   what reproduces the Tm_3'end column of the FastPCR assay-design reports this
   tool's rules come from (mean absolute error 2.2 C over their tabulated
   oligos), so a threshold quoted against that column can be used here as
   written. It is a short-oligo duplex Tm, nothing more. */
var GT_TM3_WINDOW = 10;

/* Same settings pcr.js uses, so a stretch masked there is masked here. */
var GT_LC_KMAX = 14, GT_LC_MINLEN = 21;

var GT_COMP = { a: 't', t: 'a', g: 'c', c: 'g' };

function gtTm(s)   { return DNA_Tm(String(s).toLowerCase(), GT_SALT_M, GT_MG_M, GT_OLIGO_UM); }
function gtTm3(s)  { s = String(s); return gtTm(s.slice(-GT_TM3_WINDOW)); }
function gtRC(s)   { return ComplementDNA(String(s || '').toLowerCase()); }
function gtUp(s)   { return String(s == null ? '' : s).toUpperCase(); }
function gtFmt(x, d) { return isFinite(x) ? Number(x).toFixed(d == null ? 1 : d) : 'n/a'; }
function gtPad(s, w) { s = String(s == null ? '' : s); return s + ' '.repeat(Math.max(0, w - s.length)); }
function gtPadL(s, w) { s = String(s == null ? '' : s); return ' '.repeat(Math.max(0, w - s.length)) + s; }

/** Longest homopolymer run. */
function gtMaxRun(s) {
  var best = 0, run = 0, prev = '';
  for (var i = 0; i < s.length; i++) {
    if (s.charAt(i) === prev) run++; else { run = 1; prev = s.charAt(i); }
    if (run > best) best = run;
  }
  return best;
}

/** Longest hairpin stem inside `s`, anchored on an exact k-mer. */
function gtHairpin(s, minStem, minLoop) {
  minStem = minStem || 4; minLoop = minLoop == null ? 3 : minLoop;
  var n = s.length, best = 0;
  for (var i = 0; i + minStem <= n; i++) {
    var rc = gtRC(s.substr(i, minStem));
    var at = s.indexOf(rc, i + minStem + minLoop);
    if (at === -1) continue;
    var stem = minStem, a = i + minStem, b = at - 1;
    while (b - a >= minLoop && GT_COMP[s.charAt(a)] === s.charAt(b)) { stem++; a++; b--; }
    if (stem > best) best = stem;
  }
  return best;
}

// ── IUPAC 3'-end patterns ────────────────────────────────────────────────────

var GT_IUPAC = {
  a: 'a', c: 'c', g: 'g', t: 't',
  r: 'ag', y: 'ct', s: 'gc', w: 'at', k: 'gt', m: 'ac',
  b: 'cgt', d: 'agt', h: 'act', v: 'acg', n: 'acgt'
};

/**
 * Does the oligo's 3' end match one of the slash- or space-separated patterns?
 * A pattern is read 5'->3' and compared against the same number of terminal
 * bases, so "ca" means ...CA-3'. Patterns are IUPAC, as on the other pages.
 */
function gtEnd3Match(oligo, spec) {
  if (!spec) return true;
  var pats = String(spec).toLowerCase().split(/[\s/,]+/).filter(Boolean);
  if (!pats.length) return true;
  var s = String(oligo).toLowerCase();
  for (var p = 0; p < pats.length; p++) {
    var pat = pats[p];
    if (pat.length > s.length) continue;
    var tail = s.slice(-pat.length), ok = true;
    for (var i = 0; i < pat.length; i++) {
      var set = GT_IUPAC[pat.charAt(i)];
      if (!set || set.indexOf(tail.charAt(i)) === -1) { ok = false; break; }
    }
    if (ok) return true;
  }
  return false;
}

// ── Input ────────────────────────────────────────────────────────────────────

/**
 * Keep the bases and nothing else, case preserved. U is read as T. An IUPAC
 * ambiguity code is kept as N rather than dropped: dropping it would shift
 * every position after it, so the coordinates reported would no longer be
 * those of the sequence pasted, and an oligo could be laid over a base nobody
 * knows. An N stays in the template and no oligo may cover it.
 */
function gtClean(s) {
  return String(s || '').replace(/\s/g, '')
    .replace(/u/g, 't').replace(/U/g, 'T')
    .replace(/[rywsmkbdhv]/g, 'n').replace(/[RYWSMKBDHV]/g, 'N')
    .replace(/[^acgtnACGTN]/g, '');
}

function gtRecords(raw) {
  var lines = String(raw || '').split(/\r?\n/), recs = [], cur = null;
  for (var i = 0; i < lines.length; i++) {
    if (/^\s*>/.test(lines[i])) {
      cur = { name: lines[i].replace(/^\s*>\s*/, '').split(/[\s|]/)[0] || ('seq' + (recs.length + 1)), body: '' };
      recs.push(cur);
    } else if (lines[i].trim()) {
      if (!cur) { cur = { name: 'seq1', body: '' }; recs.push(cur); }
      cur.body += lines[i];
    }
  }
  return recs;
}

/**
 * Parse the markup form. The conventions are the ones used in this class of
 * work, and only the ones that mean something here:
 *
 *   [CASSETTE/]   the bracketed sequence is present on one allele and absent on
 *   [/CASSETTE]   the other — an insertion, which is what this tool designs for
 *   /.../         a region no oligo may overlap
 *
 * Everything outside the brackets is flanking genome.
 */
function gtParseMarkup(body) {
  var warn = [];
  var excluded = [], out = '', i, ch;
  var ins = null, insAt = -1, insOnFirst = true;

  // one pass, tracking which marker we are inside
  var buf = '', state = 'flank', slashRun = null;
  for (i = 0; i < body.length; i++) {
    ch = body.charAt(i);
    if (ch === '[') {
      if (state !== 'flank') { warn.push('Nested "[" at position ' + (i + 1) + ' was ignored.'); continue; }
      state = 'bracket'; buf = ''; continue;
    }
    if (ch === ']') {
      if (state !== 'bracket') { warn.push('Unmatched "]" at position ' + (i + 1) + ' was ignored.'); continue; }
      state = 'flank';
      var parts = buf.split('/');
      if (parts.length === 1) {
        // a plain [region] — a placement hint, not an allele; keep the bases
        out += gtClean(parts[0]);
      } else if (parts.length === 2) {
        var left = gtClean(parts[0]), right = gtClean(parts[1]);
        if (ins !== null && (left || right)) {
          // one run designs for one cassette; a second group would silently
          // replace the first, so it is kept as flank and said so
          warn.push('A second bracket group, [' + gtUp(left || right).slice(0, 12) + (left ? '…/]' : '…]') +
                    ', was found after the cassette; only the first group marks the cassette and the ' +
                    'bases of this one were kept as flank. If they are one cassette, give it as one group.');
          out += left + right;
        } else if (left && !right) {
          ins = left; insAt = out.length; insOnFirst = true;
          out += ins;                                   // the long allele is the reference
        } else if (!left && right) {
          ins = right; insAt = out.length; insOnFirst = false;
        } else if (left && right) {
          warn.push('[' + gtUp(left).slice(0, 12) + '/' + gtUp(right).slice(0, 12) +
                    '] is a substitution, not an insertion. This tool designs insertion and ' +
                    'empty-site assays; for a substitution or SNP use the qPCR probe tool.');
          ins = left; insAt = out.length; insOnFirst = true;
          out += ins;
        } else {
          warn.push('An empty [/] carries no cassette, so there is nothing to design against.');
        }
      } else {
        warn.push('More than one "/" inside brackets is not supported here.');
        out += gtClean(buf.replace(/\//g, ''));
      }
      buf = ''; continue;
    }
    if (ch === '/' && state === 'flank') {
      if (slashRun === null) slashRun = out.length;
      else { if (out.length > slashRun) excluded.push([slashRun, out.length]); slashRun = null; }
      continue;
    }
    if (state === 'bracket') { buf += ch; continue; }
    out += gtClean(ch);
  }
  if (state === 'bracket') warn.push('An unclosed "[" was ignored.');
  if (slashRun !== null) warn.push('An unclosed "/" region was ignored.');

  // Excluded regions were recorded against `out`, which lacks the cassette when
  // it was given as [/CASSETTE]. They are returned on the engineered allele,
  // where the design runs, so a region names the same stretch of flank
  // whichever way the cassette was marked.
  if (ins !== null && !insOnFirst && ins.length) {
    var shifted = [];
    for (i = 0; i < excluded.length; i++) {
      var r0 = excluded[i][0], r1 = excluded[i][1];
      if (r1 <= insAt) shifted.push([r0, r1]);
      else if (r0 >= insAt) shifted.push([r0 + ins.length, r1 + ins.length]);
      else {                      // straddles the point: the flank either side, not the cassette
        shifted.push([r0, insAt]);
        shifted.push([insAt + ins.length, r1 + ins.length]);
      }
    }
    excluded = shifted;
  }

  return { seq: out, insert: ins, insAt: insAt, insOnFirst: insOnFirst, excluded: excluded, warn: warn };
}

/**
 * Derive the cassette from two alleles instead of from markup: the insert is
 * what is left when the shared prefix and the shared suffix are taken away.
 *
 * This is exact for a clean insertion, which is the case this tool is for, and
 * it degrades honestly: if the two sequences also differ elsewhere, the reported
 * insert widens to cover everything between the first and last difference, and
 * the caller is told.
 */
function gtDeriveInsert(longSeq, shortSeq) {
  var L = longSeq.toLowerCase(), S = shortSeq.toLowerCase();
  if (L.length <= S.length) return { error: 'The two sequences are the same length (or the engineered allele was given second), so there is no insertion to find. Alleles that differ by a substitution are the qPCR probe tool\'s case.' };
  var p = 0;
  while (p < S.length && L.charAt(p) === S.charAt(p)) p++;
  var q = 0;
  while (q < S.length - p && L.charAt(L.length - 1 - q) === S.charAt(S.length - 1 - q)) q++;
  var insert = longSeq.substring(p, L.length - q);
  var expected = L.length - S.length;
  var note = null;
  if (insert.length !== expected) {
    note = 'The alleles differ by ' + expected + ' nt in length but by ' + insert.length +
           ' nt of sequence, so there are changes beside the insertion. The cassette is reported as ' +
           'everything between the first and the last difference; check the boundaries before ordering.';
  }
  return { up: longSeq.substring(0, p), insert: insert, down: longSeq.substring(L.length - q), note: note };
}

/** One parsed locus, whichever way it was given. */
function gtParse(raw) {
  var recs = gtRecords(raw);
  if (!recs.length) return { error: 'Paste the locus: either one sequence with the cassette in [brackets/], or two sequences — engineered allele first, wild type second.' };

  if (recs.length >= 2) {
    var a = gtClean(recs[0].body), b = gtClean(recs[1].body);
    if (!a || !b) return { error: 'One of the two records has no A/C/G/T bases.' };
    var swapped = false;
    if (b.length > a.length) { var t = a; a = b; b = t; swapped = true; }
    var d = gtDeriveInsert(a, b);
    if (d.error) return { error: d.error };
    if (!d.insert) return { error: 'The two sequences are identical once aligned, so there is no cassette to design against.' };
    var w = [];
    if (swapped) w.push('The second record was longer, so it was taken as the engineered allele.');
    if (d.note) w.push(d.note);
    if (recs.length > 2) w.push('Only the first two records were used; ' + (recs.length - 2) + ' more were ignored.');
    return {
      name: recs[0].name, source: 'two alleles (' + recs[0].name + ' vs ' + recs[1].name + ')',
      up: d.up, insert: d.insert, down: d.down, excluded: [], warn: w
    };
  }

  var m = gtParseMarkup(recs[0].body);
  if (!m.seq) return { error: 'No A/C/G/T bases found.' };
  if (m.insert === null) {
    return { error: 'No cassette was marked. Put it in brackets with a slash — [CASSETTE/] if the bracketed sequence is the insertion — or paste the two alleles as two FASTA records instead.' };
  }
  return {
    name: recs[0].name, source: 'markup',
    up: m.seq.substring(0, m.insAt),
    insert: m.insert,
    down: m.seq.substring(m.insAt + (m.insOnFirst ? m.insert.length : 0)),
    excluded: m.excluded, warn: m.warn
  };
}

/** The two allele sequences the assays will be judged against. */
function gtAlleles(loc) {
  return {
    ki: (loc.up + loc.insert + loc.down).toLowerCase(),
    wt: (loc.up + loc.down).toLowerCase(),
    insStart: loc.up.length,                    // in KI coordinates
    insEnd: loc.up.length + loc.insert.length,
    point: loc.up.length                        // the insertion point, in WT coordinates
  };
}

/**
 * Excluded regions are kept in engineered-allele coordinates. On wild type the
 * cassette is gone, so a region downstream of it moves up by the cassette's
 * length, a region inside it vanishes, and one that straddles a junction keeps
 * only its flank part.
 */
function gtExcludedWT(excluded, alleles) {
  var out = [], a = alleles.insStart, b = alleles.insEnd, L = b - a;
  for (var i = 0; i < (excluded || []).length; i++) {
    var r0 = excluded[i][0], r1 = excluded[i][1];
    if (r0 < a) out.push([r0, Math.min(r1, a)]);
    if (r1 > b) out.push([Math.max(r0, b) - L, r1 - L]);
  }
  return out;
}

// ── Masking ──────────────────────────────────────────────────────────────────

/** Per-position mask, in the convention pcr.js and qpcr.js use: >0 = suspect. */
function gtMask(seq, opt) {
  var n = seq.length, msk = new Array(n).fill(0), i, j;
  var sources = [], missing = [];
  if (opt.maskRepeats) {
    if (typeof RepeatMask === 'function') {
      var rm = RepeatMask(seq), covered = 0;
      for (i = 0; i < n; i++) if (rm[i] > 0) { msk[i]++; covered++; }
      sources.push({ what: 'repeated 15-mers', covered: covered });
    } else {
      missing.push({ what: 'repeated 15-mers', needs: 'js/dna.js' });
    }
  }
  if (opt.maskLowComplexity) {
    if (typeof LowComplexitySequence === 'function') {
      var lc = new LowComplexitySequence(seq, GT_LC_KMAX, GT_LC_MINLEN);
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

function gtMaskedCount(mask, a, b) {
  if (!mask) return 0;
  var c = 0;
  for (var i = Math.max(0, a); i < Math.min(mask.msk.length, b); i++) if (mask.msk[i] > 0) c++;
  return c;
}

function gtOverlapsAny(ranges, a, b) {
  for (var i = 0; i < (ranges || []).length; i++) {
    if (ranges[i][0] < b && a < ranges[i][1]) return true;
  }
  return false;
}

// ── Oligo candidates ─────────────────────────────────────────────────────────

/**
 * Every oligo that fits entirely inside [lo, hi) and passes the filters.
 *
 * `strand` is the oligo's own direction: '+' reads the template left to right
 * and has its 3' end at the right edge of its footprint, '-' is the reverse
 * complement and has its 3' end at the LEFT edge. Which end is the 3' end
 * decides everything downstream, so it is carried explicitly rather than
 * inferred.
 */
function gtOligos(seq, lo, hi, strand, opt, mask, excluded) {
  var out = [], mid = (opt.tmMin + opt.tmMax) / 2;
  lo = Math.max(0, lo); hi = Math.min(seq.length, hi);
  for (var len = opt.lnMin; len <= opt.lnMax; len++) {
    for (var s = lo; s + len <= hi; s++) {
      if (gtOverlapsAny(excluded, s, s + len)) continue;
      var foot = seq.substring(s, s + len);
      if (!/^[acgt]+$/.test(foot)) continue;
      var oligo = strand === '+' ? foot : gtRC(foot);
      // the cheap tests first; the nearest-neighbour Tm is the costly one
      if (!gtEnd3Match(oligo, opt.end3)) continue;
      var run = gtMaxRun(oligo);
      if (run > opt.maxRun) continue;
      var lc = typeof LingComplexity6 === 'function' ? LingComplexity6(oligo) : 100;
      if (lc < opt.lcMin) continue;
      var tm = gtTm(oligo);
      if (tm < opt.tmMin || tm > opt.tmMax) continue;
      var tm3 = gtTm3(oligo);
      if (tm3 < opt.tm3Min || tm3 > opt.tm3Max) continue;
      var hp = gtHairpin(oligo, 4, 3);
      var mc = gtMaskedCount(mask, s, s + len);
      var cost = Math.abs(tm - mid) * 2
               + (hp >= 4 ? 2.0 * (hp - 3) : 0)
               + (run >= 4 ? 1.0 * (run - 3) : 0)
               + (100 - lc) * 0.05
               + mc * 1.5;
      out.push({
        seq: oligo, start: s, end: s + len, len: len, strand: strand,
        threeAt: strand === '+' ? s + len : s,          // 3'-terminal base index
        tm: tm, tm3: tm3, gc: CG(oligo), lc: lc, run: run, hairpin: hp,
        masked: mc, cost: cost
      });
    }
  }
  out.sort(function (a, b) { return a.cost - b.cost; });
  return out;
}

// ── Probes ───────────────────────────────────────────────────────────────────

/**
 * Split the probe at a template position and melt each half separately.
 *
 * This is what "balanced over the site" means in practice: a probe centred on a
 * junction is not simply a probe that covers it, it is one whose two halves have
 * a deliberate temperature difference, so that the half carrying the junction
 * dominates binding and a mismatch there costs more than it would in the middle
 * of a long, uniformly stable probe.
 */
function gtSplitTm(probe, offset) {
  // offset = how many bases of the probe lie 5' of the split, in probe
  // coordinates and in the probe's own direction
  var left = probe.slice(0, offset), right = probe.slice(offset);
  return {
    left: left, right: right,
    tmLeft: left.length ? gtTm(left) : 0,
    tmRight: right.length ? gtTm(right) : 0
  };
}

/**
 * The balance rule as it is written in this field: the half carrying the site
 * melts near 40 C, the other near 30, and the two differ by a margin. All three
 * numbers are parameters, because they are a convention rather than a
 * measurement; the defaults are the conventional ones.
 */
function gtBalanceOK(sp, opt) {
  var d = sp.tmLeft - sp.tmRight;
  return d >= opt.balDelta && sp.tmLeft <= opt.balLeftMax && sp.tmRight >= opt.balRightMin;
}

/**
 * Probe candidates covering template position `point`, on either strand.
 * `mode` is 'balance' for the split-Tm rule above, or 'cover' for a probe that
 * merely has to contain the point with a margin.
 */
function gtProbesAt(seq, point, opt, mask, excluded, mode) {
  var out = [], pmid = (opt.ptmMin + opt.ptmMax) / 2;
  for (var len = opt.plnMin; len <= opt.plnMax; len++) {
    for (var s = Math.max(0, point - len + 1); s + len <= seq.length && s <= point; s++) {
      if (gtOverlapsAny(excluded, s, s + len)) continue;
      var foot = seq.substring(s, s + len);
      if (!/^[acgt]+$/.test(foot)) continue;
      // a footprint's homopolymer run and duplex Tm are the same on either
      // strand, so they are tested once, and the cheap one first
      if (gtMaxRun(foot) > opt.maxRun) continue;
      var tm = gtTm(foot);
      if (tm < opt.ptmMin || tm > opt.ptmMax) continue;
      var mc = gtMaskedCount(mask, s, s + len);
      for (var si = 0; si < 2; si++) {
        var strand = si === 0 ? '+' : '-';
        var probe = strand === '+' ? foot : gtRC(foot);
        // no G at the 5' end of a reporter, the usual constraint
        if (opt.probeNo5G && probe.charAt(0) === 'g') continue;
        // where the point falls inside the probe, in the probe's own direction
        var offset = strand === '+' ? (point - s) : (s + len - point);
        var margin = Math.min(offset, len - offset);
        if (margin < opt.probeMargin) continue;
        var lc = typeof LingComplexity6 === 'function' ? LingComplexity6(probe) : 100;
        if (lc < opt.lcMin) continue;
        var sp = gtSplitTm(probe, offset);
        var balanced = gtBalanceOK(sp, opt);
        if (mode === 'balance' && !balanced) continue;
        var cost = Math.abs(tm - pmid) * 2
                 + (mode === 'balance' ? Math.abs((sp.tmLeft - sp.tmRight) - opt.balDelta) * 0.5
                                       : Math.abs(offset - len / 2) * 0.6)
                 + (gtHairpin(probe, 4, 3) >= 4 ? 2 : 0)
                 + (100 - lc) * 0.05 + mc * 1.5;
        out.push({
          seq: probe, start: s, end: s + len, len: len, strand: strand,
          tm: tm, tm3: gtTm3(probe), gc: CG(probe), lc: lc, masked: mc,
          offset: offset, margin: margin, split: sp, balanced: balanced, cost: cost
        });
      }
    }
  }
  out.sort(function (a, b) { return a.cost - b.cost; });
  return out;
}

/** Probe candidates anywhere strictly inside [lo, hi), for the knock-in assay. */
function gtProbesIn(seq, lo, hi, opt, mask, excluded, preferAt) {
  var out = [], pmid = (opt.ptmMin + opt.ptmMax) / 2;
  lo = Math.max(0, lo); hi = Math.min(seq.length, hi);
  for (var len = opt.plnMin; len <= opt.plnMax; len++) {
    for (var s = lo; s + len <= hi; s++) {
      if (gtOverlapsAny(excluded, s, s + len)) continue;
      var foot = seq.substring(s, s + len);
      if (!/^[acgt]+$/.test(foot)) continue;
      // strand-independent tests once per footprint, the cheap one first
      if (gtMaxRun(foot) > opt.maxRun) continue;
      var tm = gtTm(foot);
      if (tm < opt.ptmMin || tm > opt.ptmMax) continue;
      var mc = gtMaskedCount(mask, s, s + len);
      // a probe that straddles the junction is more specific than one beside
      // it, so it is preferred — but not required, because the junction is
      // often the hardest sequence in the locus
      var spans = preferAt != null && s < preferAt && preferAt < s + len;
      for (var si = 0; si < 2; si++) {
        var strand = si === 0 ? '+' : '-';
        var probe = strand === '+' ? foot : gtRC(foot);
        if (opt.probeNo5G && probe.charAt(0) === 'g') continue;
        var lc = typeof LingComplexity6 === 'function' ? LingComplexity6(probe) : 100;
        if (lc < opt.lcMin) continue;
        var cost = Math.abs(tm - pmid) * 2 + (spans ? -3 : 0)
                 + (gtHairpin(probe, 4, 3) >= 4 ? 2 : 0)
                 + (100 - lc) * 0.05 + mc * 1.5;
        out.push({
          seq: probe, start: s, end: s + len, len: len, strand: strand,
          tm: tm, tm3: gtTm3(probe), gc: CG(probe), lc: lc, masked: mc,
          spansJunction: spans, cost: cost
        });
      }
    }
  }
  out.sort(function (a, b) { return a.cost - b.cost; });
  return out;
}

// ── Pairing ──────────────────────────────────────────────────────────────────

function gtAmpliconLen(fwd, rev) { return rev.end - fwd.start; }

/**
 * Longest stretch over which two oligos can anneal to each other.
 *
 * Oligo A pairs with oligo B wherever a piece of A equals the reverse
 * complement of a piece of B, so the longest complementary run is the longest
 * common substring of A and revcomp(B).
 */
function gtMaxComplement(a, b) {
  var A = String(a).toLowerCase(), B = gtRC(b);
  var m = A.length, n = B.length, best = 0, i, j;
  var prev = new Array(n + 1).fill(0), cur;
  for (i = 1; i <= m; i++) {
    cur = new Array(n + 1).fill(0);
    for (j = 1; j <= n; j++) {
      if (A.charAt(i - 1) === B.charAt(j - 1)) {
        cur[j] = prev[j - 1] + 1;
        if (cur[j] > best) best = cur[j];
      }
    }
    prev = cur;
  }
  return best;
}

/**
 * The stretch of template an oligo can usefully occupy, given that the product
 * has to cross `junction` and be no longer than ampMax.
 *
 * Deriving this before generating candidates rather than filtering afterwards is
 * what keeps the search both correct and quick. Filtering afterwards means
 * ranking thousands of oligos the geometry had already excluded, and then -- the
 * way this first went wrong -- truncating that ranked list to the cheapest few
 * hundred, which are not the ones near the junction at all.
 */
function gtWindow(len, junction, side, opt) {
  return side === 'before'
    ? [Math.max(0, junction - opt.ampMax), junction]
    : [junction, Math.min(len, junction + opt.ampMax)];
}

/**
 * Best forward/reverse/probe combination, by total cost, with the amplicon in
 * range and the probe strictly between the primers.
 *
 * Candidates arrive already restricted to the windows above and sorted by cost,
 * so the heads of the lists are both near the junction and good. The join takes
 * the best GT_JOIN of each side; for a given pair the first probe that fits is
 * by construction the cheapest one that fits, so the probe scan stops there.
 */
var GT_JOIN = 80, GT_JOIN_P = 500;

function gtAssemble(fwds, revs, probes, opt, want) {
  var F = fwds.slice(0, GT_JOIN), R = revs.slice(0, GT_JOIN), P = probes.slice(0, GT_JOIN_P);
  var best = null, tried = 0, dimers = 0;
  var ampMid = (opt.ampMin + opt.ampMax) / 2;
  // P is cost-ordered, so its head is the most any probe can take off a pair's
  // cost. That bound is not zero: the junction bonus makes some probes negative.
  var pFloor = P.length ? P[0].cost : 0;
  // A probe's pairing with a primer does not depend on which pair it is tried
  // with, so each (probe, primer) check is made once and remembered.
  var dimF = [], dimR = [];
  function dimer(cache, k, i, primer) {
    var row = cache[k] || (cache[k] = []);
    if (row[i] === undefined) row[i] = gtMaxComplement(P[k].seq, primer.seq) >= opt.dimerMax;
    return row[i];
  }
  for (var i = 0; i < F.length; i++) {
    for (var j = 0; j < R.length; j++) {
      if (R[j].start < F[i].end) continue;                 // primers must not overlap
      var len = gtAmpliconLen(F[i], R[j]);
      if (len < opt.ampMin || len > opt.ampMax) continue;
      tried++;
      var pairCost = F[i].cost + R[j].cost
                   + Math.abs(F[i].tm - R[j].tm) * 1.5     // a matched pair anneals together
                   + Math.abs(len - ampMid) * 0.02;
      if (best && pairCost + pFloor >= best.cost) continue; // cannot win even with the best probe
      // Two oligos that anneal to each other consume each other instead of the
      // template. In a repetitive locus the cheapest forward and reverse can
      // even come out as exact reverse complements, which is a perfect dimer
      // and no assay at all, so this is a rejection rather than a penalty.
      if (gtMaxComplement(F[i].seq, R[j].seq) >= opt.dimerMax) { dimers++; continue; }
      for (var k = 0; k < P.length; k++) {
        var p = P[k];
        if (p.start < F[i].end || p.end > R[j].start) continue;
        if (dimer(dimF, k, i, F[i]) || dimer(dimR, k, j, R[j])) continue;
        var cost = pairCost + p.cost;
        if (!best || cost < best.cost) {
          best = { fwd: F[i], rev: R[j], probe: p, amplicon: len, cost: cost, want: want };
        }
        break;                                              // cost-ordered: this is the best fit
      }
    }
  }
  return {
    best: best, pairsTried: tried, dimersRejected: dimers,
    capped: fwds.length > GT_JOIN || revs.length > GT_JOIN || probes.length > GT_JOIN_P
  };
}

// ── In-silico check ──────────────────────────────────────────────────────────

/**
 * Where an oligo can prime on a sequence. An exact match is reported as such; a
 * near match is reported only when the 3'-terminal GT_ANCHOR bases are exact,
 * because a polymerase extends from the 3' end and a mismatch there is what
 * actually prevents priming.
 *
 * This is a local check against the two allele sequences given, not a genome
 * search. It can only say that an assay distinguishes THESE two alleles.
 */
var GT_ANCHOR = 5;

function gtFindSites(oligo, seq, maxMism) {
  var hits = [], n = seq.length, L = oligo.length;
  var fwd = oligo.toLowerCase(), rev = gtRC(oligo);
  function scan(pat, orient) {
    for (var s = 0; s + L <= n; s++) {
      var mism = 0, ok = true;
      for (var i = 0; i < L; i++) {
        if (seq.charAt(s + i) !== pat.charAt(i)) {
          mism++;
          if (mism > maxMism) { ok = false; break; }
          // the 3'-terminal anchor must be exact; for the forward orientation
          // that is the right-hand end of the match, for the reverse the left
          var fromThree = orient === '+' ? (L - 1 - i) : i;
          if (fromThree < GT_ANCHOR) { ok = false; break; }
        }
      }
      if (ok) hits.push({ at: s, end: s + L, orient: orient, mism: mism });
    }
  }
  scan(fwd, '+');
  scan(rev, '-');
  return hits;
}

/** Amplicons an oligo pair would give on one sequence. */
function gtAmplicons(fwdSeq, revSeq, seq, maxMism, maxLen) {
  var f = gtFindSites(fwdSeq, seq, maxMism).filter(function (h) { return h.orient === '+'; });
  var r = gtFindSites(revSeq, seq, maxMism).filter(function (h) { return h.orient === '-'; });
  var out = [];
  for (var i = 0; i < f.length; i++) {
    for (var j = 0; j < r.length; j++) {
      // the two footprints have to be separate stretches of template; when they
      // coincide the "product" is the primer pair annealed to each other
      if (r[j].at < f[i].end) continue;
      var len = r[j].end - f[i].at;
      if (len > maxLen) continue;
      out.push({ at: f[i].at, end: r[j].end, len: len, mism: f[i].mism + r[j].mism });
    }
  }
  out.sort(function (a, b) { return a.len - b.len; });
  return out;
}

/** What each assay should do on each allele. */
function gtInSilico(assays, alleles, opt) {
  var rows = [];
  var ceiling = Math.max(opt.ampMax * 6, alleles.ki.length);
  ['ki', 'wt'].forEach(function (which) {
    var seq = alleles[which];
    assays.forEach(function (a) {
      if (!a.design) return;
      var d = a.design;
      var amps = gtAmplicons(d.fwd.seq, d.rev.seq, seq, opt.verifyMism, ceiling);
      var probeHits = gtFindSites(d.probe.seq, seq, opt.verifyMism);
      // a probe only reports if it lies inside a product
      var inProduct = amps.some(function (m) {
        return probeHits.some(function (h) { return h.at >= m.at && h.end <= m.end; });
      });
      // amplicons are sorted by length, so the shortest decides whether anything
      // amplifies at all; only then does it matter where the reporter sits
      rows.push({
        assay: a.key, allele: which,
        amplicons: amps, probeHits: probeHits.length, inProduct: inProduct,
        verdict: !amps.length ? 'no product'
               : (amps[0].len > opt.ampMax ? 'product too long to amplify'
                 : (!inProduct ? 'product, but the reporter does not sit in it' : 'signal'))
      });
    });
  });
  return rows;
}

// ── Design ───────────────────────────────────────────────────────────────────

/**
 * Knock-in assay. One primer inside the cassette, its partner in the flank.
 * `orient` 'fwd-in-cassette' puts the forward primer on the cassette and reads
 * out into the downstream flank; 'rev-in-cassette' is the mirror image.
 */
function gtDesignKI(loc, alleles, opt, mask, excluded, orient) {
  var seq = alleles.ki, a = alleles.insStart, b = alleles.insEnd;
  var fwds, revs, probes, region, wF, wR;
  if (orient === 'rev-in-cassette') {
    // ---- Fwd in upstream flank ---->   <---- Rev in cassette
    // the product crosses the cassette's 5' junction, at a
    wF = gtWindow(seq.length, a, 'before', opt);
    wR = gtWindow(seq.length, a, 'after', opt); wR[1] = Math.min(wR[1], b);
    fwds = gtOligos(seq, wF[0], wF[1], '+', opt, mask, excluded);
    revs = gtOligos(seq, wR[0], wR[1], '-', opt, mask, excluded);
    probes = gtProbesIn(seq, wF[0], wR[1], opt, mask, excluded, a);
    region = 'forward primer in the upstream flank, reverse primer on the cassette';
  } else {
    //        Fwd in cassette ---->   probe ---->      <---- Rev in flank
    // the product crosses the cassette's 3' junction, at b
    wF = gtWindow(seq.length, b, 'before', opt); wF[0] = Math.max(wF[0], a);
    wR = gtWindow(seq.length, b, 'after', opt);
    fwds = gtOligos(seq, wF[0], wF[1], '+', opt, mask, excluded);
    revs = gtOligos(seq, wR[0], wR[1], '-', opt, mask, excluded);
    probes = gtProbesIn(seq, wF[0], wR[1], opt, mask, excluded, b);
    region = 'forward primer on the cassette, reverse primer in the downstream flank';
  }
  var r = gtAssemble(fwds, revs, probes, opt, 'ki');
  return {
    key: 'knock-in', orient: orient, region: region, design: r.best,
    counts: { fwd: fwds.length, rev: revs.length, probe: probes.length },
    pairsTried: r.pairsTried, dimersRejected: r.dimersRejected, capped: r.capped,
    coords: 'engineered allele'
  };
}

/**
 * Empty-site assay. Both primers in the flanks, the probe straddling the
 * insertion point, designed on the wild-type allele where that point is a
 * single position rather than a cassette.
 */
function gtDesignWT(loc, alleles, opt, mask, excluded, probeMode) {
  var seq = alleles.wt, pt = alleles.point;
  // the mask was built on the engineered allele; rebuild it for wild type so
  // positions line up
  var wmask = gtMask(seq, opt);
  // the excluded regions, likewise, are moved onto wild-type coordinates
  var wex = gtExcludedWT(excluded, alleles);
  var wF = gtWindow(seq.length, pt, 'before', opt), wR = gtWindow(seq.length, pt, 'after', opt);
  var fwds = gtOligos(seq, wF[0], wF[1], '+', opt, wmask, wex);
  var revs = gtOligos(seq, wR[0], wR[1], '-', opt, wmask, wex);
  var probes = gtProbesAt(seq, pt, opt, wmask, wex, probeMode);
  var r = gtAssemble(fwds, revs, probes, opt, 'wt');
  return {
    key: 'empty site', orient: 'both primers in the flanks', probeMode: probeMode,
    region: 'forward and reverse primers either side of the insertion point, reporter across it',
    design: r.best, mask: wmask,
    counts: { fwd: fwds.length, rev: revs.length, probe: probes.length },
    pairsTried: r.pairsTried, dimersRejected: r.dimersRejected, capped: r.capped,
    coords: 'wild-type allele'
  };
}

// ── Rendering ────────────────────────────────────────────────────────────────

function gtOligoRow(name, o, offset) {
  offset = offset || 0;
  return gtPad(name, 16) +
         gtPad(o.strand === '+' ? 'top (+)' : 'bottom (-)', 12) +
         gtPad((offset + o.start + 1) + '-' + (offset + o.end), 14) +
         gtPadL(o.len, 4) + '  ' +
         gtPadL(gtFmt(o.tm, 1), 7) + '  ' +
         gtPadL(gtFmt(o.tm3, 1), 9) + '  ' +
         gtPadL(gtFmt(o.gc, 0), 5) + '  ' +
         gtUp(o.seq);
}

function gtRenderAssays(loc, alleles, assays, opt, warn) {
  var t = [];
  t.push('INSERTION / KNOCK-IN GENOTYPING ASSAYS');
  t.push('');
  t.push('Locus      : ' + loc.name + '   read from ' + loc.source);
  t.push('Cassette   : ' + loc.insert.length + ' bp, at position ' + (loc.up.length + 1) +
         ' of the engineered allele');
  t.push('Flanks     : ' + loc.up.length + ' bp upstream, ' + loc.down.length + ' bp downstream');
  t.push('Alleles    : engineered ' + alleles.ki.length + ' bp, wild type ' + alleles.wt.length + ' bp');
  t.push('');

  assays.forEach(function (a) {
    t.push('═'.repeat(100));
    t.push(a.key === 'knock-in' ? 'KNOCK-IN ASSAY — signal only when the cassette is present'
                                : 'EMPTY-SITE ASSAY — signal from the unmodified locus');
    t.push('═'.repeat(100));
    t.push('  ' + a.region);
    t.push('  coordinates below are on the ' + a.coords);
    t.push('');
    if (!a.design) {
      t.push('  No set found. Candidates: ' + a.counts.fwd + ' forward, ' + a.counts.rev +
             ' reverse, ' + a.counts.probe + ' reporter' +
             (a.pairsTried ? '; ' + a.pairsTried + ' primer pairs were in the amplicon range' : '') +
             (a.dimersRejected ? '; ' + a.dimersRejected + ' pairs were rejected as dimers' : '') + '.');
      var why = [];
      if (!a.counts.fwd) why.push('no forward primer passed the filters');
      if (!a.counts.rev) why.push('no reverse primer passed the filters');
      if (!a.counts.probe) why.push('no reporter passed the filters' +
        (a.probeMode === 'balance' ? ' — the balance rule is the usual reason, try "cover the point" instead' : ''));
      if (a.counts.fwd && a.counts.rev && !a.pairsTried) why.push('no primer pair fell inside the amplicon range');
      if (a.pairsTried && a.counts.probe) why.push('no reporter fitted between a usable primer pair');
      if (a.dimersRejected && a.dimersRejected >= a.pairsTried) why.push('every pair in range annealed to itself');
      t.push('  Why: ' + (why.join('; ') || 'the combination could not be satisfied') + '.');
      warn.push(a.key + ' assay: no oligo set was found. ' + (why[0] || '') );
      t.push('');
      return;
    }
    var d = a.design, off = 0;
    t.push('  ' + gtPad('name', 16) + gtPad('strand', 12) + gtPad('position', 14) +
           gtPadL('len', 4) + '  ' + gtPadL('Tm', 7) + '  ' + gtPadL("Tm_3'end", 9) + '  ' +
           gtPadL('GC%', 5) + '  sequence (5\'->3\')');
    t.push('  ' + '-'.repeat(96));
    var tag = a.key === 'knock-in' ? 'KI' : 'WT';
    t.push('  ' + gtOligoRow(tag + '-F', d.fwd, off));
    t.push('  ' + gtOligoRow(tag + '-R', d.rev, off));
    t.push('  ' + gtOligoRow(tag + '-P', d.probe, off));
    t.push('');
    t.push('  amplicon ' + d.amplicon + ' bp   primer dTm ' + gtFmt(Math.abs(d.fwd.tm - d.rev.tm), 1) +
           ' C   Ta ~ ' + gtFmt(Math.min(d.fwd.tm, d.rev.tm) + Math.log(d.amplicon), 0) + ' C' +
           '   longest primer-primer pairing ' + gtMaxComplement(d.fwd.seq, d.rev.seq) + ' bp');
    if (a.key === 'empty site') {
      var sp = d.probe.split;
      t.push('');
      t.push('  Reporter balance across the insertion point:');
      t.push('    ' + gtPad('half carrying the point (5\' side)', 38) + gtPad(gtUp(sp.left), 34) +
             gtFmt(sp.tmLeft, 1) + ' C');
      t.push('    ' + gtPad('other half (3\' side)', 38) + gtPad(gtUp(sp.right), 34) +
             gtFmt(sp.tmRight, 1) + ' C');
      t.push('    ' + gtPad('difference', 38) + gtPad('', 34) + gtFmt(sp.tmLeft - sp.tmRight, 1) + ' C');
      t.push('    rule: difference >= ' + opt.balDelta + ', 5\' half <= ' + opt.balLeftMax +
             ', 3\' half >= ' + opt.balRightMin + '  ->  ' + (d.probe.balanced ? 'satisfied' : 'NOT satisfied'));
      if (!d.probe.balanced) {
        warn.push('The empty-site reporter does not satisfy the balance rule; it was accepted because the ' +
                  'reporter mode was set to cover the point rather than balance across it.');
      }
    } else if (d.probe.spansJunction) {
      t.push('  the reporter straddles the cassette/genome junction, which is as specific as this assay gets');
    } else {
      t.push('  the reporter lies beside the junction, not across it — specificity rests on the cassette primer');
    }
    t.push('');
  });

  t.push('Naming: F forward primer, R reverse primer, P reporter. Tm is the site-wide');
  t.push('nearest-neighbour value (55 mM Na+/K+, 1 mM Mg2+, 0.2 uM); Tm_3\'end is the Tm of');
  t.push('the 3\'-terminal ' + GT_TM3_WINDOW + ' nucleotides, the quantity assay-design rules for this work');
  t.push('are usually written in.');
  return t.join('\n');
}

function gtRenderZygosity(loc, alleles, assays, rows, opt, warn) {
  var t = [];
  t.push('ZYGOSITY — WHAT THE PAIR READS OUT');
  t.push('');
  var have = {};
  assays.forEach(function (a) { if (a.design) have[a.key] = true; });
  if (!have['knock-in'] || !have['empty site']) {
    t.push('  Both assays are needed to call zygosity, and ' +
           (!have['knock-in'] && !have['empty site'] ? 'neither was designed.'
            : (!have['knock-in'] ? 'the knock-in assay was not designed.' : 'the empty-site assay was not designed.')));
    t.push('');
  } else {
    t.push('  ' + gtPad('genotype', 18) + gtPad('knock-in assay', 20) + gtPad('empty-site assay', 20) + 'reading');
    t.push('  ' + '-'.repeat(84));
    t.push('  ' + gtPad('wild type +/+', 18) + gtPad('no signal', 20) + gtPad('signal', 20) + 'no cassette at this locus');
    t.push('  ' + gtPad('heterozygote +/-', 18) + gtPad('signal', 20) + gtPad('signal', 20) + 'one allele carries it');
    t.push('  ' + gtPad('homozygote -/-', 18) + gtPad('signal', 20) + gtPad('no signal', 20) + 'both alleles carry it');
    t.push('');
    t.push('  Run them on the same plate. The empty-site assay is the internal control that');
    t.push('  separates a true homozygote from a failed reaction: a well with no signal in');
    t.push('  either assay is a failed well, not a genotype.');
    t.push('');
  }

  t.push('─'.repeat(100));
  t.push('In-silico check against the two allele sequences');
  t.push('─'.repeat(100));
  t.push('  ' + gtPad('assay', 14) + gtPad('allele', 20) + gtPad('products', 26) + 'verdict');
  t.push('  ' + '-'.repeat(96));
  var expect = {
    'knock-in|ki': 'signal', 'knock-in|wt': 'no product',
    'empty site|wt': 'signal', 'empty site|ki': 'any'
  };
  var trouble = [];
  rows.forEach(function (r) {
    var amps = r.amplicons.length
      ? r.amplicons.slice(0, 2).map(function (m) { return m.len + ' bp' + (m.mism ? ' (' + m.mism + ' mism)' : ''); }).join(', ') +
        (r.amplicons.length > 2 ? ' +' + (r.amplicons.length - 2) + ' more' : '')
      : 'none';
    var key = r.assay + '|' + r.allele;
    var want = expect[key] || 'any';
    var ok = want === 'any' || want === r.verdict;
    t.push('  ' + gtPad(r.assay, 14) +
           gtPad(r.allele === 'ki' ? 'engineered' : 'wild type', 20) +
           gtPad(amps, 26) + r.verdict + (ok ? '' : '   <- expected ' + want));
    if (!ok) trouble.push(r.assay + ' on the ' + (r.allele === 'ki' ? 'engineered' : 'wild-type') +
                          ' allele: ' + r.verdict + ', expected ' + want);
  });
  t.push('');
  if (trouble.length) {
    t.push('  ' + trouble.length + ' result(s) are not what the design intends — see Warnings.');
    trouble.forEach(function (m) { if (warn) warn.push('In-silico check: ' + m + '.'); });
  } else {
    t.push('  Every assay behaves as intended on both alleles.');
  }
  t.push('');
  var es = null;
  assays.forEach(function (a) { if (a.key === 'empty site' && a.design) es = a; });
  if (es) {
    t.push('  The empty-site assay is expected to give a product on the engineered allele too —');
    t.push('  longer by the ' + loc.insert.length + ' bp cassette. Whether that product appears in a real');
    t.push('  reaction depends on the extension time; with a cassette of this size the product');
    t.push('  would be ' + (es.design.amplicon + loc.insert.length) + ' bp, which is why the assay reads');
    t.push('  as wild-type specific.');
    t.push('');
  }
  t.push('  Scope: this is a check against the two sequences you supplied, both strands, with');
  t.push('  up to ' + opt.verifyMism + ' mismatch(es) outside the 3\'-terminal ' + GT_ANCHOR + ' bases. It is not a genome');
  t.push('  search — it cannot tell you the assay is unique in the genome, only that it');
  t.push('  distinguishes these two alleles. Use In Silico PCR against your own sequence set,');
  t.push('  or a genome BLAT, for the specificity question.');
  return t.join('\n');
}

function gtRenderContext(loc, alleles, assays, opt, mask) {
  var t = [], W = 100;
  t.push('LOCUS MAP AND MASKING');
  t.push('');
  var ki = alleles.ki, a = alleles.insStart, b = alleles.insEnd;
  var col = function (p) { return Math.round(p * W / ki.length); };
  t.push('  Engineered allele, ' + ki.length + ' bp, drawn to ' + W + ' columns:');
  var bar = '';
  for (var i = 0; i < W; i++) {
    var p = i * ki.length / W;
    bar += (p >= a && p < b) ? '#' : '-';
  }
  t.push('  ' + bar);
  t.push('  ' + ' '.repeat(col(a)) + '^' + (b - a > 2 ? ' '.repeat(Math.max(0, col(b) - col(a) - 2)) + '^' : ''));
  t.push('  ' + ' '.repeat(col(a)) + 'cassette ' + (b - a) + ' bp');
  t.push('');
  assays.forEach(function (as) {
    if (!as.design) return;
    var d = as.design, onKI = as.coords === 'engineered allele';
    var L = onKI ? ki.length : alleles.wt.length;
    var c = function (p) { return Math.round(p * W / L); };
    var line = new Array(W).fill(' ');
    function put(o, chv) {
      var s = c(o.start), e = Math.max(s + 1, c(o.end));
      for (var x = s; x < e && x < W; x++) line[x] = chv;
    }
    put(d.fwd, '>'); put(d.rev, '<'); put(d.probe, '=');
    t.push('  ' + gtPad(as.key, 12) + line.join(''));
  });
  t.push('  ' + gtPad('', 12) + '> forward primer   < reverse primer   = reporter');
  t.push('');

  t.push('─'.repeat(100));
  t.push('Masking — the same two filters the primer-design modules use');
  t.push('─'.repeat(100));
  (mask.missing || []).forEach(function (m) {
    t.push('  ' + gtPad(m.what, 40) + 'NOT RUN — ' + m.needs + ' is not loaded on this page');
  });
  if (!mask.sources.length && !(mask.missing || []).length) {
    t.push('  Both masks are switched off, so oligos may sit on repeated or low-complexity sequence.');
  } else {
    mask.sources.forEach(function (s) {
      t.push('  ' + gtPad(s.what, 40) + gtPad(s.covered + ' nt', 12) +
             gtFmt(100 * s.covered / ki.length, 1) + '% of the engineered allele');
    });
    t.push('  ' + gtPad('masked in total', 40) + gtPad(mask.masked + ' nt', 12) +
           gtFmt(100 * mask.masked / ki.length, 1) + '%');
    t.push('');
    t.push('  A masked position is one an oligo should not sit on. Where an oligo had to be placed');
    t.push('  on masked sequence anyway, its row shows a non-zero masked count in the candidate');
    t.push('  tables and the cost of that choice is already in its ranking.');
  }
  t.push('');
  if (loc.excluded && loc.excluded.length) {
    t.push('  Regions you excluded with /.../ : ' +
           loc.excluded.map(function (r) { return (r[0] + 1) + '-' + r[1]; }).join(', ') +
           '. No oligo was allowed to overlap them.');
    t.push('');
  }
  t.push('─'.repeat(100));
  t.push('Candidates that passed the filters');
  t.push('─'.repeat(100));
  t.push('  ' + gtPad('assay', 14) + gtPadL('forward', 9) + gtPadL('reverse', 9) +
         gtPadL('reporter', 10) + gtPadL('pairs in range', 16) + gtPadL('dimers rejected', 18));
  assays.forEach(function (as) {
    t.push('  ' + gtPad(as.key, 14) + gtPadL(as.counts.fwd, 9) + gtPadL(as.counts.rev, 9) +
           gtPadL(as.counts.probe, 10) + gtPadL(as.pairsTried, 16) +
           gtPadL(as.dimersRejected || 0, 18) +
           (as.capped ? '   (list capped)' : ''));
  });
  return t.join('\n');
}

function gtRenderWarnings(loc, assays, opt, extra, mask) {
  var t = [], w = [];
  t.push('WARNINGS AND LIMITS');
  t.push('');
  (mask.missing || []).forEach(function (m) {
    w.push('The ' + m.what + ' filter was switched on but could not run: ' + m.needs +
           ' is not loaded on this page, so nothing was masked by it. This is a page fault, not a ' +
           'property of your sequence.');
  });
  (loc.warn || []).forEach(function (m) { w.push(m); });
  (extra || []).forEach(function (m) { w.push(m); });

  assays.forEach(function (a) {
    if (!a.design) return;
    var d = a.design;
    [['forward primer', d.fwd], ['reverse primer', d.rev], ['reporter', d.probe]].forEach(function (pair) {
      var o = pair[1];
      if (o.masked) {
        w.push(a.key + ' ' + pair[0] + ' covers ' + o.masked + ' masked nt — it sits partly on repeated or ' +
               'low-complexity sequence and may prime elsewhere.');
      }
      var hp = gtHairpin(o.seq, 5, 3);
      if (hp >= 6) w.push(a.key + ' ' + pair[0] + ' folds back on itself over ' + hp + ' bp.');
    });
    if (Math.abs(d.fwd.tm - d.rev.tm) > 4) {
      w.push(a.key + ' assay: the primers differ by ' + gtFmt(Math.abs(d.fwd.tm - d.rev.tm), 1) +
             ' C. One annealing temperature has to serve both, so the weaker primer sets the Ta.');
    }
  });

  if (loc.insert.length < 30) {
    w.push('The cassette is only ' + loc.insert.length + ' bp. A primer sitting inside it has little room, and ' +
           'the empty-site assay cannot separate the two alleles by product length — this is a small ' +
           'insertion, better handled as an InDel with allele-specific probes in the qPCR tool.');
  }
  if (loc.up.length < 60 || loc.down.length < 60) {
    w.push('Flanking sequence is short (' + loc.up.length + ' bp upstream, ' + loc.down.length +
           ' bp downstream). Give at least a few hundred bases either side or the search has nowhere to put ' +
           'the outer primer.');
  }

  if (!w.length) t.push('• Nothing flagged.');
  else w.forEach(function (m) { t.push('• ' + m); });

  t.push('');
  t.push('Scope and limits:');
  t.push('  1) The in-silico check compares the two alleles you supplied and nothing else. It says');
  t.push('     whether the assay separates them, not whether it is unique in a genome.');
  t.push('  2) Tm is the site-wide nearest-neighbour value at fixed conditions, so it is comparable');
  t.push('     with every other tool here but is not the Tm of your reaction. Where a rule in this');
  t.push('     field is quoted at a different convention — vendor software commonly runs several');
  t.push('     degrees higher — compare the Tm_3\'end column rather than the full-length Tm, which is');
  t.push('     where the two conventions agree most closely for short oligos.');
  t.push('  3) The reporter is placed and scored, but no dye or quencher chemistry is modelled: no');
  t.push('     MGB stabilisation, no dye-dependent Tm shift, no quencher. A reporter that looks');
  t.push('     right here still has to be ordered in the chemistry you actually use.');
  t.push('  4) Nothing here checks copy number beyond presence or absence. A multi-copy integrant');
  t.push('     reads as a homozygote on this pair of assays; separating those needs a quantitative');
  t.push('     copy-number assay against a reference gene, which this tool does not design.');
  t.push('  5) The primer pair search is greedy over a capped candidate list, not exhaustive. The');
  t.push('     candidate counts are reported so you can see when the cap bites.');
  return t.join('\n');
}

// ── Orchestration ────────────────────────────────────────────────────────────

function gtNum(id, dflt) {
  var el = document.getElementById(id);
  if (!el) return dflt;
  var v = parseFloat(el.value);
  return isFinite(v) ? v : dflt;
}
function gtVal(id, dflt) { var el = document.getElementById(id); return el && el.value != null ? el.value : dflt; }
function gtChk(id, dflt) { var el = document.getElementById(id); return el ? !!el.checked : !!dflt; }

function gtSetPanes(map) {
  for (var id in map) { var el = document.getElementById(id); if (el) el.value = map[id]; }
}
function gtPending(note) {
  gtSetPanes({ analysisResult1: note, analysisResult2: '', analysisResult3: '', analysisResult4: '' });
}

function runGenotype() {
  var raw = gtVal('inputText', '');
  if (!String(raw).trim()) {
    gtPending('Paste the locus, 5\'->3\'. Two ways to give it:\n\n' +
              '  one sequence   ...flank[CASSETTE/]flank...   the bracketed sequence is the insertion\n' +
              '  two sequences  engineered allele first, wild type second — the cassette is found by\n' +
              '                 comparing them\n\n' +
              'The tool returns a knock-in assay, an empty-site assay, and what the pair reads out\n' +
              'for each genotype.\n\n' +
              'Click "Design assays" to run, or "Load Example" for a worked locus.');
    return;
  }

  var loc = gtParse(raw);
  if (loc.error) { gtPending('[input] ' + loc.error); return; }

  var opt = {
    lnMin: gtNum('ln_min', 18), lnMax: gtNum('ln_max', 28),
    tmMin: gtNum('tm_min', 58), tmMax: gtNum('tm_max', 64),
    tm3Min: gtNum('tm3_min', 20), tm3Max: gtNum('tm3_max', 40),
    end3: gtVal('end3', ''),
    plnMin: gtNum('pln_min', 16), plnMax: gtNum('pln_max', 28),
    ptmMin: gtNum('ptm_min', 58), ptmMax: gtNum('ptm_max', 68),
    probeNo5G: gtChk('probe_no5g', true),
    probeMargin: Math.max(1, gtNum('probe_margin', 4)),
    probeMode: gtVal('probe_mode', 'balance'),
    balDelta: gtNum('bal_delta', 9),
    balLeftMax: gtNum('bal_left_max', 41),
    balRightMin: gtNum('bal_right_min', 20),
    ampMin: gtNum('amp_min', 80), ampMax: gtNum('amp_max', 220),
    maxRun: Math.max(3, gtNum('max_run', 4)),
    lcMin: gtNum('lc_min', 60),
    verifyMism: Math.max(0, Math.min(3, gtNum('verify_mism', 1))),
    dimerMax: Math.max(4, gtNum('dimer_max', 9)),
    maskRepeats: gtChk('mask_repeats', true),
    maskLowComplexity: gtChk('mask_lowcomplexity', true),
    wantKI: gtChk('assay_ki', true),
    wantWT: gtChk('assay_wt', true),
    orient: gtVal('orient', 'fwd-in-cassette')
  };
  if (opt.lnMin > opt.lnMax) { gtPending('[input] The primer length range is inverted.'); return; }
  if (opt.tmMin > opt.tmMax) { gtPending('[input] The primer Tm range is inverted.'); return; }
  if (opt.tm3Min > opt.tm3Max) { gtPending('[input] The 3\'-end Tm range is inverted.'); return; }
  if (opt.plnMin > opt.plnMax) { gtPending('[input] The reporter length range is inverted.'); return; }
  if (opt.ptmMin > opt.ptmMax) { gtPending('[input] The reporter Tm range is inverted.'); return; }
  if (opt.ampMin > opt.ampMax) { gtPending('[input] The amplicon size range is inverted.'); return; }
  if (!opt.wantKI && !opt.wantWT) { gtPending('[input] Switch on at least one of the two assays.'); return; }
  if (!loc.insert || !loc.insert.length) { gtPending('[input] The cassette is empty, so there is nothing to design against.'); return; }

  try {
    var warn = [];
    var alleles = gtAlleles(loc);
    var mask = gtMask(alleles.ki, opt);

    var assays = [];
    if (opt.wantKI) assays.push(gtDesignKI(loc, alleles, opt, mask, loc.excluded, opt.orient));
    if (opt.wantWT) assays.push(gtDesignWT(loc, alleles, opt, mask, loc.excluded, opt.probeMode));

    var rows = gtInSilico(assays, alleles, opt);

    // the first two panes add to `warn` as they render, so the warnings pane is built last
    var pane1 = gtRenderAssays(loc, alleles, assays, opt, warn);
    var pane2 = gtRenderZygosity(loc, alleles, assays, rows, opt, warn);
    var pane3 = gtRenderContext(loc, alleles, assays, opt, mask);
    var pane4 = gtRenderWarnings(loc, assays, opt, warn, mask);
    gtSetPanes({ analysisResult1: pane1, analysisResult2: pane2, analysisResult3: pane3, analysisResult4: pane4 });
  } catch (err) {
    gtPending('[error] ' + (err && err.message ? err.message : err));
    if (window.console) console.error('genotype error:', err);
  }
}

// ── Examples, pending state, reset ───────────────────────────────────────────

/* Synthetic loci. The sequences are generated, not taken from any real project:
   a worked example on a public page should not carry someone's locus. */
var GT_EXAMPLES = {
  floxed: {
    note: 'Example loaded — a floxed cassette in a synthetic genomic context, given as markup.\n' +
          'Click "Design assays" to run it.',
    seq:
      '>floxed_locus synthetic - 1248 bp context, a 608 bp loxP-flanked cassette\n' +
      'ccgggcctcgccagttgaggacttcggctgaccctggtagtcccattttttcctgtgataagtgggaatgagtccaaaatcgatgacgtaagtcggagaa\n' +
      'cgttactatttctaacaattttttatcaagttatcctacactttcatggttttcgtctaattcgtcgagtgttactgtggttcgtttcgagtgcctctta\n' +
      'cgccgatgtcgctttttcacgagcgtttagaaaatctgtttcagcatctcaaggccgaaaaaaagaccagatttttgaactcatggaatcagtacacaat\n' +
      'ggctcaggtccaaagctact[ataacttcgtatagcatacattatacgaagttatgctttcagacccagaacgggacacaacggaagtcaacgcgtgggt\n' +
      'ttaggactctatgttgtttgcatctatcgtatcgctgtgtactaccggggattagtaggcggagacacctagatccgtccaggacctagcaccttagttg\n' +
      'ctaatgcactgggccgaactgttaggcaacgggcgaaaaccgaccactaacggaagcagcccggatgtgcactaccggggcctgaacatcaaccacgttg\n' +
      'ctagggtatcgagctgacatctagccctcagccgccttaatgaggctgtagtgcccgggtgtccggaaaatattttccatgctttatgcaccgcggaagg\n' +
      'tatgacggcaatatatttgctggctatataatcatatacacaactatcatggcaagggtattaacgtgggatggcaagtgccggcaagtgaatgcctggt\n' +
      'aagaattccggccgcgtgtgtcacatttttttgacaacaggagaaggcgcgacgtcggtatggctggacacttcgtataacaaaggaatgggtaaataac\n' +
      'ttcgtatagcatacattatacgaagttat/]gtacactggtaaagttccgggacgtggctaatataatacggtttttgggtagtcatactcgataatcat\n' +
      'aatgtccggtgccacctatgatatgtacttgacggcattaaggattccacgaaatcggcaatcagcgtattagtcaaaccgagatgaatctagaatatta\n' +
      'ttatgtctaattgcgaatacctagcgaaattgcctcattggagcttgacttggtaagtcatgaacctaatgtctcgaagtaggattagaaagccagccct\n' +
      'cctatctattgttatttgacgctagaatgcaaggtcaattacactcgtatc'
  },
  twoalleles: {
    note: 'Example loaded — the same locus given as two alleles instead of markup, so the cassette\n' +
          'is found by comparing them. Click "Design assays" to run it.',
    seq: null         // filled in below from the markup example
  }
};

(function () {
  // Build the two-allele example from the markup one, so the two examples are
  // demonstrably the same locus read two different ways.
  var m = GT_EXAMPLES.floxed.seq;
  var body = m.split('\n').slice(1).join('');
  var p = gtParseMarkup(body);
  if (!p.insert) return;
  var up = p.seq.substring(0, p.insAt), ins = p.insert;
  var down = p.seq.substring(p.insAt + ins.length);
  function wrap(s) {
    var out = [];
    for (var i = 0; i < s.length; i += 100) out.push(s.substr(i, 100));
    return out.join('\n');
  }
  GT_EXAMPLES.twoalleles.seq =
    '>engineered_allele synthetic\n' + wrap(up + ins + down) + '\n' +
    '>wild_type_allele synthetic\n' + wrap(up + down);
})();

function gtClearFeedback() {
  if (window.SeqValidate) window.SeqValidate.clear('inputFeedback');
}

/** A design belongs to the sequence it was run on, so an edit retires it. */
function gtAwaitRun() {
  var el = document.getElementById('inputText');
  if (!el || !el.value.trim()) { runGenotype(); return; }
  gtPending('Click "Design assays" to run the design for the sequence in the box.');
}

function loadExampleGenotype(name) {
  var ex = GT_EXAMPLES[name] || GT_EXAMPLES.floxed;
  var el = document.getElementById('inputText');
  if (!el) return;
  el.value = ex.seq;
  gtClearFeedback();
  gtPending(ex.note);
}

function clearGenotype() {
  var el = document.getElementById('inputText');
  if (el) el.value = '';
  gtClearFeedback();
  runGenotype();
}

if (typeof window !== 'undefined') {
  window.runGenotype = runGenotype;
  window.gtAwaitRun = gtAwaitRun;
  window.gtPending = gtPending;
  window.loadExampleGenotype = loadExampleGenotype;
  window.clearGenotype = clearGenotype;
  window.GT_EXAMPLES = GT_EXAMPLES;
}
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    gtClean: gtClean, gtParse: gtParse, gtParseMarkup: gtParseMarkup, gtDeriveInsert: gtDeriveInsert,
    gtAlleles: gtAlleles, gtExcludedWT: gtExcludedWT,
    gtOligos: gtOligos, gtProbesAt: gtProbesAt, gtProbesIn: gtProbesIn,
    gtSplitTm: gtSplitTm, gtBalanceOK: gtBalanceOK, gtFindSites: gtFindSites,
    gtAmplicons: gtAmplicons, gtInSilico: gtInSilico, gtAssemble: gtAssemble,
    gtMaxComplement: gtMaxComplement, gtDesignKI: gtDesignKI, gtDesignWT: gtDesignWT,
    gtTm: gtTm, gtTm3: gtTm3, gtEnd3Match: gtEnd3Match, gtMask: gtMask
  };
}
