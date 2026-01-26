/* =====================================================================================
 * - FASTA parse/format
 * - DNA sanitize (ACGTN only)
 * - k-mer Jaccard distances (fast)
 * - UPGMA guide tree
 * - Progressive profile alignment:
 * - exact DP for small profiles
 * - banded DP for larger profiles
 * ===================================================================================== */

(function () {
  "use strict";

  // -------------------------- UI helpers --------------------------
  const elIn = document.getElementById("inputText");
  const elOut = document.getElementById("outFasta");
  const elClustal = document.getElementById("outClustal");
  const elLog = document.getElementById("outLog") || document.getElementById("log");
  const btnRun = document.getElementById("btnRun");
  const btnDownload = document.getElementById("btnDownload");

  const elClustalW = document.getElementById("clustalW");
  const elK = document.getElementById("kmerK");
  const elMatch = document.getElementById("scMatch");
  const elMismatch = document.getElementById("scMismatch");
  const elGap = document.getElementById("scGap");
  const elBandW = document.getElementById("bandW");
  const elExactMax = document.getElementById("exactMax");

  function log(msg, cls) {
    // Support both <textarea> and <div> log targets
    if (!elLog) return;
    if (elLog.tagName === "TEXTAREA") {
      elLog.value += (msg + "\n");
      elLog.scrollTop = elLog.scrollHeight;
      return;
    }
    const div = document.createElement("div");
    div.textContent = msg;
    if (cls) div.className = cls;
    elLog.appendChild(div);
  }
  function clearLog() {
    if (!elLog) return;
    if (elLog.tagName === "TEXTAREA") { elLog.value = ""; return; }
    elLog.innerHTML = "";
  }

  function downloadText(filename, text) {
    const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }


  // -------------------------- UI helpers (extended) --------------------------

  const $ = (id) => document.getElementById(id);

  function isTextArea(el) {
    return el && el.tagName === "TEXTAREA";
  }

  function readText(el, fallback = "") {
    if (!el) return fallback;
    if (isTextArea(el) || ("value" in el)) return (el.value ?? fallback);
    return (el.textContent ?? fallback);
  }

  function writeText(el, text) {
    if (!el) return;
    if (isTextArea(el) || ("value" in el)) el.value = (text ?? "");
    else el.textContent = (text ?? "");
  }

  function appendText(el, text, addNewline = true) {
    if (!el) return;
    const t = (text ?? "") + (addNewline ? "\n" : "");
    if (isTextArea(el)) {
      el.value += t;
      el.scrollTop = el.scrollHeight;
      return;
    }
    const div = document.createElement("div");
    div.textContent = (text ?? "");
    el.appendChild(div);
  }

  function writeHTML(el, html) {
    if (!el) return;
    el.innerHTML = (html ?? "");
  }

  function readNumber(el, def, min = -Infinity, max = Infinity) {
    const raw = readText(el, "");
    let v = Number.parseFloat(raw);
    if (!Number.isFinite(v)) v = def;
    if (Number.isFinite(min)) v = Math.max(min, v);
    if (Number.isFinite(max)) v = Math.min(max, v);
    return v;
  }

  function setDisabled(el, disabled) {
    if (!el) return;
    el.disabled = !!disabled;
  }

  function setButtonBusy(btn, busy, labelBusy = "Running", labelIdle = "Generate") {
    if (!btn) return;
    if (busy) {
      btn.dataset._origText = btn.textContent;
      btn.dataset._origBg = btn.style.backgroundColor;
      btn.textContent = labelBusy;
      btn.style.backgroundColor = "#f31221";
      btn.disabled = true;
    } else {
      btn.textContent = btn.dataset._origText || labelIdle;
      btn.style.backgroundColor = btn.dataset._origBg || "";
      btn.disabled = false;
    }
  }

  // Unified input/output helpers for this page
  const ui = {
    getInputFasta: () => readText(elIn, ""),
    setOutputFasta: (s) => writeText(elOut, s),
    setOutputClustal: (s) => writeText(elClustal, s),
    clearOutputClustal: () => { if (elClustal) writeText(elClustal, ""); },
    logLine: (s) => appendText(elLog, s, true),
    clearLog: () => {
      if (!elLog) return;
      if (isTextArea(elLog)) writeText(elLog, "");
      else elLog.innerHTML = "";
    }
  };

  // -------------------------- FASTA parse/format --------------------------

  function parseFasta(text) {
    const lines = text.split(/\r?\n/);
    const recs = [];
    let id = null;
    let seq = [];
    for (const raw of lines) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith(">")) {
        if (id !== null) recs.push({ id, seq: seq.join("") });
        id = line.slice(1).trim() || `seq${recs.length + 1}`;
        seq = [];
      } else {
        seq.push(line.replace(/\s+/g, ""));
      }
    }
    if (id !== null) recs.push({ id, seq: seq.join("") });
    if (recs.length === 0) throw new Error("No FASTA records found.");
    return recs;
  }

  function formatFasta(records, width = 80) {
    const out = [];
    for (const r of records) {
      out.push(">" + r.id);
      const s = r.seq;
      for (let i = 0; i < s.length; i += width) out.push(s.slice(i, i + width));
    }
    return out.join("\n");
  }

  // -------------------------- DNA sanitize --------------------------

  function sanitizeDNA(seq) {
    // Keep A C G T I N only; convert U->T;  I->G;
    const s = seq.toUpperCase().replace(/U/g, "T").replace(/I/g, "G");
    let out = "";
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (c === "A" || c === "C" || c === "G" || c === "T" || c === "N") out += c;
    }
    return out;
  }

  // -------------------------- Scoring --------------------------

  function scorePair(a, b, P) {
    // a,b: char (A,C,G,T,N or '-') ; simple linear gap penalty
    if (a === "-" && b === "-") return 0;
    if (a === "-" || b === "-") return P.gap;
    return (a === b) ? P.match : P.mismatch;
  }

  // -------------------------- k-mer distance (Jaccard) --------------------------

  function kmerSet(seq, k) {
    // Using hashed strings; for k<=12 and moderate data OK in browser
    const set = new Set();
    if (seq.length < k) return set;
    for (let i = 0; i <= seq.length - k; i++) {
      set.add(seq.slice(i, i + k));
    }
    return set;
  }

  function jaccardDistance(setA, setB) {
    if (setA.size === 0 && setB.size === 0) return 0;
    let inter = 0;
    // iterate smaller
    const [small, big] = setA.size < setB.size ? [setA, setB] : [setB, setA];
    for (const x of small) if (big.has(x)) inter++;
    const union = setA.size + setB.size - inter;
    return union === 0 ? 0 : 1.0 - (inter / union);
  }

  function computeDistancesByKmer(records, k) {
    const N = records.length;
    const kmers = new Array(N);
    for (let i = 0; i < N; i++) kmers[i] = kmerSet(records[i].seq, k);

    const dist = Array.from({ length: N }, () => new Float64Array(N));
    for (let i = 0; i < N; i++) {
      dist[i][i] = 0;
      for (let j = i + 1; j < N; j++) {
        const d = jaccardDistance(kmers[i], kmers[j]);
        dist[i][j] = d;
        dist[j][i] = d;
      }
    }
    return dist;
  }

  // -------------------------- UPGMA --------------------------

  function upgma(distMatrix) {
    const N = distMatrix.length;
    let nextId = N;

    let clusters = [];
    for (let i = 0; i < N; i++) {
      clusters.push({
        id: i,
        members: [i],
        node: { leaf: true, index: i, members: [i], height: 0 }
      });
    }

    function clusterDistance(c1, c2) {
      let sum = 0, cnt = 0;
      for (const i of c1.members) for (const j of c2.members) { sum += distMatrix[i][j]; cnt++; }
      return sum / cnt;
    }

    while (clusters.length > 1) {
      let bestI = 0, bestJ = 1;
      let bestD = clusterDistance(clusters[0], clusters[1]);

      for (let i = 0; i < clusters.length; i++) {
        for (let j = i + 1; j < clusters.length; j++) {
          const d = clusterDistance(clusters[i], clusters[j]);
          if (d < bestD) { bestD = d; bestI = i; bestJ = j; }
        }
      }

      const A = clusters[bestI], B = clusters[bestJ];
      const members = A.members.concat(B.members);
      const node = { left: A.node, right: B.node, members: members.slice(), height: bestD / 2 };
      const merged = { id: nextId++, members, node };
      const jIdx = Math.max(bestI, bestJ);
      const iIdx = Math.min(bestI, bestJ);
      clusters.splice(jIdx, 1);
      clusters.splice(iIdx, 1);
      clusters.push(merged);
    }
    return clusters[0].node;
  }

  // -------------------------- Profiles --------------------------

  function profileFromSingle(record) {
    return { ids: [record.id], aligned: [record.seq] };
  }

  function profileLen(P) { return P.aligned[0].length; }

  // Column scores for profile-profile
  function colScore(PA, ia, PB, ib, sc) {
    let s = 0;
    const A = PA.aligned, B = PB.aligned;
    for (let i = 0; i < A.length; i++) {
      const a = A[i][ia];
      for (let j = 0; j < B.length; j++) {
        const b = B[j][ib];
        s += scorePair(a, b, sc);
      }
    }
    return s;
  }

  function gapScoreLeft(PA, ia, PB, sc) {
    // column ia in A vs gap-column in B
    let s = 0;
    const A = PA.aligned, B = PB.aligned;
    for (let i = 0; i < A.length; i++) {
      const a = A[i][ia];
      for (let j = 0; j < B.length; j++) s += scorePair(a, "-", sc);
    }
    return s;
  }

  function gapScoreUp(PA, PB, ib, sc) {
    // gap-column in A vs column ib in B
    let s = 0;
    const A = PA.aligned, B = PB.aligned;
    for (let i = 0; i < A.length; i++) {
      for (let j = 0; j < B.length; j++) s += scorePair("-", B[j][ib], sc);
    }
    return s;
  }

  // Exact DP profile-profile (O(LA*LB))
  function alignProfilesExact(PA, PB, sc) {
    const LA = profileLen(PA), LB = profileLen(PB);
    const dp = Array.from({ length: LA + 1 }, () => new Int32Array(LB + 1));
    const tb = Array.from({ length: LA + 1 }, () => new Uint8Array(LB + 1)); // 1 diag, 2 up, 3 left

    dp[0][0] = 0;
    for (let i = 1; i <= LA; i++) { dp[i][0] = dp[i - 1][0] + gapScoreLeft(PA, i - 1, PB, sc); tb[i][0] = 2; }
    for (let j = 1; j <= LB; j++) { dp[0][j] = dp[0][j - 1] + gapScoreUp(PA, PB, j - 1, sc); tb[0][j] = 3; }

    for (let i = 1; i <= LA; i++) {
      for (let j = 1; j <= LB; j++) {
        const d = dp[i - 1][j - 1] + colScore(PA, i - 1, PB, j - 1, sc);
        const u = dp[i - 1][j] + gapScoreLeft(PA, i - 1, PB, sc);
        const l = dp[i][j - 1] + gapScoreUp(PA, PB, j - 1, sc);
        let best = d, dir = 1;
        if (u > best) { best = u; dir = 2; }
        if (l > best) { best = l; dir = 3; }
        dp[i][j] = best; tb[i][j] = dir;
      }
    }

    // Traceback ops
    let i = LA, j = LB;
    const ops = [];
    while (i > 0 || j > 0) {
      const dir = tb[i][j];
      if (dir === 1) { ops.push("D"); i--; j--; }
      else if (dir === 2) { ops.push("U"); i--; }
      else { ops.push("L"); j--; }
    }
    ops.reverse();
    return applyOpsToProfiles(PA, PB, ops);
  }

  // Banded DP profile-profile (O((LA+LB)*bandW))
  function alignProfilesBanded(PA, PB, sc, bandW) {
    const LA = profileLen(PA), LB = profileLen(PB);

    // band centered around diagonal: j ~ i * (LB/LA) ; but we keep simple diagonal with offset
    // We'll use j in [i - w, i + w] with w = bandW/2, but scaled if LA != LB.
    const w = Math.max(32, Math.floor(bandW / 2));

    // dp rows as Maps over j within band
    const NEG = -2147483640; // "minus infinity" within int32
    const dp = Array.from({ length: LA + 1 }, () => new Int32Array(0)); // placeholders
    const tb = Array.from({ length: LA + 1 }, () => new Uint8Array(0));
    const jMinArr = new Int32Array(LA + 1);
    const jMaxArr = new Int32Array(LA + 1);

    function bandForRow(i) {
      // diagonal estimate: i * LB/LA
      const diag = (LA === 0) ? 0 : Math.round(i * (LB / LA));
      const jMin = Math.max(0, diag - w);
      const jMax = Math.min(LB, diag + w);
      return [jMin, jMax];
    }

    for (let i = 0; i <= LA; i++) {
      const [jMin, jMax] = bandForRow(i);
      jMinArr[i] = jMin; jMaxArr[i] = jMax;
      const len = jMax - jMin + 1;
      dp[i] = new Int32Array(len);
      tb[i] = new Uint8Array(len);
      for (let t = 0; t < len; t++) dp[i][t] = NEG;
    }

    function get(dpRow, jMin, j) {
      const idx = j - jMin;
      if (idx < 0 || idx >= dpRow.length) return NEG;
      return dpRow[idx];
    }
    function set(dpRow, jMin, j, val) {
      const idx = j - jMin;
      dpRow[idx] = val;
    }
    function setTB(tbRow, jMin, j, val) {
      const idx = j - jMin;
      tbRow[idx] = val;
    }
    function getTB(tbRow, jMin, j) {
      const idx = j - jMin;
      return tbRow[idx];
    }

    // init dp[0][0] if in band
    {
      const jMin0 = jMinArr[0];
      if (0 >= jMin0 && 0 <= jMaxArr[0]) set(dp[0], jMin0, 0, 0);
    }

    // Initialize first row within band: dp[0][j]
    for (let j = 1; j <= LB; j++) {
      const jMin0 = jMinArr[0];
      if (j < jMin0 || j > jMaxArr[0]) continue;
      const prev = get(dp[0], jMin0, j - 1);
      if (prev === NEG) continue;
      const v = prev + gapScoreUp(PA, PB, j - 1, sc);
      set(dp[0], jMin0, j, v);
      setTB(tb[0], jMin0, j, 3);
    }

    // Initialize first column: dp[i][0]
    for (let i = 1; i <= LA; i++) {
      const jMin = jMinArr[i];
      if (0 < jMin || 0 > jMaxArr[i]) continue;
      const prev = get(dp[i - 1], jMinArr[i - 1], 0);
      if (prev === NEG) continue;
      const v = prev + gapScoreLeft(PA, i - 1, PB, sc);
      set(dp[i], jMin, 0, v);
      setTB(tb[i], jMin, 0, 2);
    }

    // fill
    for (let i = 1; i <= LA; i++) {
      const jMin = jMinArr[i], jMax = jMaxArr[i];
      for (let j = jMin; j <= jMax; j++) {
        if (i === 0 && j === 0) continue;

        const dPrev = get(dp[i - 1], jMinArr[i - 1], j - 1);
        const uPrev = get(dp[i - 1], jMinArr[i - 1], j);
        const lPrev = get(dp[i], jMin, j - 1);

        let best = NEG, dir = 0;

        if (dPrev !== NEG && j - 1 >= 0) {
          const v = dPrev + colScore(PA, i - 1, PB, j - 1, sc);
          if (v > best) { best = v; dir = 1; }
        }
        if (uPrev !== NEG) {
          const v = uPrev + gapScoreLeft(PA, i - 1, PB, sc);
          if (v > best) { best = v; dir = 2; }
        }
        if (lPrev !== NEG && j - 1 >= 0) {
          const v = lPrev + gapScoreUp(PA, PB, j - 1, sc);
          if (v > best) { best = v; dir = 3; }
        }

        if (best !== NEG) {
          set(dp[i], jMin, j, best);
          setTB(tb[i], jMin, j, dir);
        }
      }
    }

    // Traceback from (LA, LB)
    const endRowMin = jMinArr[LA];
    if (LB < endRowMin || LB > jMaxArr[LA]) {
      // if target is outside band, fallback to exact
      return alignProfilesExact(PA, PB, sc);
    }
    if (get(dp[LA], endRowMin, LB) === NEG) {
      // unreachable within band, fallback to exact
      return alignProfilesExact(PA, PB, sc);
    }

    let i = LA, j = LB;
    const ops = [];
    while (i > 0 || j > 0) {
      const rowMin = jMinArr[i];
      const dir = getTB(tb[i], rowMin, j);
      if (dir === 1) { ops.push("D"); i--; j--; }
      else if (dir === 2) { ops.push("U"); i--; }
      else if (dir === 3) { ops.push("L"); j--; }
      else {
        // safety fallback
        return alignProfilesExact(PA, PB, sc);
      }
    }
    ops.reverse();
    return applyOpsToProfiles(PA, PB, ops);
  }

  function applyOpsToProfiles(PA, PB, ops) {
    const newIds = PA.ids.concat(PB.ids);
    const buildersA = PA.aligned.map(() => []);
    const buildersB = PB.aligned.map(() => []);

    let colA = 0, colB = 0;
    for (const op of ops) {
      if (op === "D") {
        for (let s = 0; s < PA.aligned.length; s++) buildersA[s].push(PA.aligned[s][colA]);
        for (let s = 0; s < PB.aligned.length; s++) buildersB[s].push(PB.aligned[s][colB]);
        colA++; colB++;
      } else if (op === "U") {
        for (let s = 0; s < PA.aligned.length; s++) buildersA[s].push(PA.aligned[s][colA]);
        for (let s = 0; s < PB.aligned.length; s++) buildersB[s].push("-");
        colA++;
      } else { // "L"
        for (let s = 0; s < PA.aligned.length; s++) buildersA[s].push("-");
        for (let s = 0; s < PB.aligned.length; s++) buildersB[s].push(PB.aligned[s][colB]);
        colB++;
      }
    }

    const aligned = [];
    for (const b of buildersA) aligned.push(b.join(""));
    for (const b of buildersB) aligned.push(b.join(""));

    return { ids: newIds, aligned };
  }

  function alignProfilesAuto(PA, PB, sc, exactMaxCols, bandW) {
    const LA = profileLen(PA), LB = profileLen(PB);
    const maxCols = Math.max(LA, LB);

    if (maxCols <= exactMaxCols) return alignProfilesExact(PA, PB, sc);

    // banded for larger cases
    return alignProfilesBanded(PA, PB, sc, bandW);
  }

  // -------------------------- Build MSA from guide tree --------------------------

  function buildMSA(node, records, sc, exactMaxCols, bandW) {
    if (node.leaf) return profileFromSingle(records[node.index]);
    const L = buildMSA(node.left, records, sc, exactMaxCols, bandW);
    const R = buildMSA(node.right, records, sc, exactMaxCols, bandW);
    return alignProfilesAuto(L, R, sc, exactMaxCols, bandW);
  }

  // -------------------------- Stats --------------------------

  function consensus(alignedRecords) {
    const L = alignedRecords[0].seq.length;
    let out = "";
    for (let k = 0; k < L; k++) {
      const counts = new Map();
      for (const r of alignedRecords) {
        const c = r.seq[k];
        if (c === "-") continue;
        counts.set(c, (counts.get(c) || 0) + 1);
      }
      if (counts.size === 0) { out += "-"; continue; }
      let bestC = "N", bestN = -1;
      for (const [c, n] of counts.entries()) {
        if (n > bestN) { bestN = n; bestC = c; }
      }
      out += bestC;
    }
    return out;
  }


  // -------------------------- CLUSTAL formatting --------------------------

  function safeClustalName(name, used) {
    // CLUSTAL traditionally shows up to 60 chars; ensure uniqueness
    const mcut = 60;
    let base = (name || "seq").trim().replace(/\s+/g, "_");
    if (base.length > mcut) base = base.slice(0, mcut);
    let out = base;
    let n = 2;
    while (used.has(out)) {
      const suffix = "_" + n;
      const cut = Math.max(1, mcut - suffix.length);
      out = base.slice(0, cut) + suffix;
      n++;
    }
    used.add(out);
    return out;
  }

  function clustalConsensusLine(blockSeqs) {
    // blockSeqs: array of strings (same length) for the block
    const L = blockSeqs[0].length;
    let line = "";
    for (let k = 0; k < L; k++) {
      let c0 = null;
      let nonGap = 0;
      let ok = true;
      for (const s of blockSeqs) {
        const c = s[k];
        if (c === "-") continue;
        nonGap++;
        if (c0 === null) c0 = c;
        else if (c !== c0) { ok = false; break; }
      }
      // conservative: mark '*' if all non-gap are identical and at least 2 sequences contribute
      line += (ok && nonGap >= 2) ? "*" : " ";
    }
    return line;
  }

  function toClustal(alignedRecords, blockSize = 220) {
    if (!alignedRecords || alignedRecords.length === 0) return "CLUSTAL\n\n";
    const L = alignedRecords[0].seq.length;
    for (const r of alignedRecords) {
      if (r.seq.length !== L) throw new Error("Aligned sequences have different lengths; cannot format CLUSTAL.");
    }

    const used = new Set();
    const names = alignedRecords.map(r => safeClustalName(r.id, used));
    const nameW = Math.max(...names.map(n => n.length), 10) + 2;

    const lines = [];
    lines.push("CLUSTAL W multiple sequence alignment");
    lines.push("");

    for (let pos = 0; pos < L; pos += blockSize) {
      const block = alignedRecords.map(r => r.seq.slice(pos, pos + blockSize));
      for (let i = 0; i < alignedRecords.length; i++) {
        const nm = names[i].padEnd(nameW, " ");
        lines.push(nm + block[i]);
      }
      const cons = clustalConsensusLine(block);
      lines.push("".padEnd(nameW, " ") + cons);
      lines.push("");
    }

    return lines.join("\n");
  }
  // -------------------------- Main run --------------------------

  async function runMSA() {
    ui.clearLog();
    ui.setOutputFasta("");
    ui.clearOutputClustal();
    setDisabled(btnDownload, true);

    const sc = {
      match: readNumber(elMatch, 2),
      mismatch: readNumber(elMismatch, -1),
      gap: readNumber(elGap, -2)
    };
    const k = Math.max(3, Math.min(12, readNumber(elK, 6)));
    const bandW = Math.max(64, readNumber(elBandW, 512));
    const exactMaxCols = Math.max(500, readNumber(elExactMax, 4000));
    const clustalW = Math.max(60, Math.min(300, readNumber(elClustalW, 110)));
    const input = ui.getInputFasta();
    const inputBytes = new Blob([input]).size;

    log(`Input size: ${inputBytes.toLocaleString()} bytes`);

    if (inputBytes > 1_200_000) {
      log("Warning: input exceeds ~1.2 MB. Browser MSA may be slow or fail due to memory/time constraints.", "warn");
    }

    let records = parseFasta(input);
    log(`Sequences: ${records.length}`);

    // sanitize DNA
    let totalLen = 0;
    records = records.map(r => {
      const s = sanitizeDNA(r.seq);
      totalLen += s.length;
      return { id: r.id, seq: s };
    });

    log(`Total DNA length: ${totalLen.toLocaleString()} nt`);
    if (records.length < 2) throw new Error("Need at least 2 sequences.");

    // Basic guardrails
    // Pairwise distance uses k-mer sets: memory grows with unique k-mers.
    // For very long sequences, choose higher k to reduce collision? not necessarily.
    if (totalLen > 1_000_000) {
      log("Warning: total length near 1Mb. Consider fewer/shorter sequences for reliable runtime in browser.", "warn");
    }

    // Compute k-mer distances
    log(`Computing k-mer distances (k=${k})...`);
    const t0 = performance.now();
    const dist = computeDistancesByKmer(records, k);
    const t1 = performance.now();
    log(`Distance matrix built in ${(t1 - t0).toFixed(1)} ms`);

    // UPGMA
    log("Building UPGMA guide tree...");
    const t2 = performance.now();
    const tree = upgma(dist);
    const t3 = performance.now();
    log(`Guide tree built in ${(t3 - t2).toFixed(1)} ms`);

    // Progressive alignment
    log(`Progressive alignment (exactMaxCols=${exactMaxCols}, bandW=${bandW})...`);
    const t4 = performance.now();
    const profile = buildMSA(tree, records, sc, exactMaxCols, bandW);
    let t5 = performance.now() - t4;

    if (t5 > 999) {
      t5 = t5 / 1000;
      log(`Alignment completed in ${(t5).toFixed(1)} s`);
    } else {
      log(`Alignment completed in ${(t5).toFixed(1)} ms`);
    }

    // Output
    const aligned = profile.ids.map((id, i) => ({ id, seq: profile.aligned[i] }));
    const alnLen = aligned[0].seq.length;
    log(`Aligned columns: ${alnLen.toLocaleString()}`);

    const cons = consensus(aligned);
    const outText = formatFasta(aligned) + "\n\n>CONSENSUS\n" + cons + "\n";
    ui.setOutputFasta(outText);

    // CLUSTAL output (aligned sequences only)
    if (elClustal) {
      try {
        ui.setOutputClustal(toClustal(aligned, clustalW));
      } catch (e) {
        ui.setOutputClustal("CLUSTAL formatting error: " + (e && e.message ? e.message : String(e)));
      }
    }

    renderClustalColored(aligned, clustalW);

    const colorTabBtn = document.querySelector('.tablink[onclick*="res_color"]');
    if (colorTabBtn) colorTabBtn.click();
 openTab(null, 'res_color');
 
    setDisabled(btnDownload, false);
    btnDownload.onclick = () => downloadText("msa_alignment.fasta", outText);
  }

  // -------------------------- Bindings --------------------------

  btnRun.addEventListener("click", () => {
    setButtonBusy(btnRun, true, "Running", "Generate");

    setTimeout(async () => {
      try {
        await runMSA(); // запуск вашего выравнивания
      } catch (err) {
        log("Error: " + (err && err.message ? err.message : String(err)), "err");
      } finally {
        setButtonBusy(btnRun, false, "Running", "Generate");
      }
    }, 0);
  });

  function columnConservation(chars) {
    // chars: массив символов в колонке (включая '-'?)
    const counts = new Map();
    let nonGap = 0;

    for (const c of chars) {
      if (c === "-") continue;
      nonGap++;
      counts.set(c, (counts.get(c) || 0) + 1);
    }
    if (nonGap === 0) return 0; // все gaps

    let best = 0;
    for (const v of counts.values()) best = Math.max(best, v);

    return best / nonGap; // 0..1
  }

  function consToClass(consRatio) {
    // маппинг в 0..4
    if (consRatio >= 1.0) return "c4";
    if (consRatio >= 0.85) return "c3";
    if (consRatio >= 0.65) return "c2";
    if (consRatio >= 0.45) return "c1";
    return "c0";
  }

  function escapeHtml(s) {
    return s.replace(/[&<>"']/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]));
  }

  function renderClustalColored(alignedRecords, blockSize = 110) {
    const host = document.getElementById("clustalColor");
    if (!host) return;

    if (!alignedRecords || alignedRecords.length === 0) {
      host.innerHTML = "";
      return;
    }

    const L = alignedRecords[0].seq.length;
    for (const r of alignedRecords) {
      if (r.seq.length !== L) {
        host.innerHTML = `<div class="warn">Aligned sequences have different lengths; cannot render.</div>`;
        return;
      }
    }

    const used = new Set();
    const names = alignedRecords.map(r => safeClustalName(r.id, used));
    const nameW = Math.max(...names.map(n => n.length), 10) + 2;

    let pre = "";
    pre += `${escapeHtml("CLUSTAL W multiple sequence alignment")}\n\n`;

    for (let pos = 0; pos < L; pos += blockSize) {
      const block = alignedRecords.map(r => r.seq.slice(pos, pos + blockSize));

      // conservation class for each column in this block
      const classes = [];
      for (let k = 0; k < block[0].length; k++) {
        const colChars = block.map(s => s[k]);
        const ratio = columnConservation(colChars);
        classes.push(consToClass(ratio));
      }

      // sequence lines
      for (let i = 0; i < alignedRecords.length; i++) {
        const nm = escapeHtml(names[i].padEnd(nameW, " "));
        const s = block[i];

        let lineHtml = "";
        for (let k = 0; k < s.length; k++) {
          const ch = s[k];
          const cls = (ch === "-") ? "c0" : classes[k];
          lineHtml += `<span class="${cls}">${escapeHtml(ch)}</span>`;
        }

        pre += `${nm}${lineHtml}\n`;
      }

      // consensus line (same logic as in toClustal)
      const cons = clustalConsensusLine(block);
      pre += `${"".padEnd(nameW, " ")}${escapeHtml(cons)}\n\n`;
    }

    host.innerHTML = `
      <div class="clustal-panel">
        <div class="clustal-head">
          <div class="clustal-title">Colored CLUSTAL</div>
          <div class="clustal-hint">blockSize=${blockSize}, columns=${L}</div>
        </div>
        <pre class="clustal-pre">${pre}</pre>
      </div>
    `;
  }



})();
