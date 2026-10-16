/* ────────────────────────────────────────────────────────────────────────────
   example-link.js — deep-link into a tool with its worked example preloaded.

   Opening any module with "?example=1" (or "?example=<name>" where the module's
   loadExample() accepts a name, as MSA does) fills the input with the example
   data the module already ships, exactly as clicking "Load Example" would.

   This is what the DigitalGens landing page and the documentation link to, so a
   first-time user reaches a working input in one click. Modules that design on
   demand (TSDR) still need their run button pressed once the example is in.
   ──────────────────────────────────────────────────────────────────────────── */
(function () {
  "use strict";

  function run() {
    var value;
    try {
      value = new URLSearchParams(window.location.search).get("example");
    } catch (e) {
      return;                       // very old browser: silently skip
    }
    if (value === null) return;     // parameter absent — nothing to do

    if (typeof window.loadExample !== "function") {
      console.warn("[example-link] This page has no loadExample() to call.");
      return;
    }

    var isDefault = (value === "" || value === "1" || value === "true");
    var name = isDefault ? null : value;

    // Some pages (MSA) offer several named examples and their loadExample()
    // requires the name — calling it bare would assign examples[undefined],
    // i.e. the literal string "undefined". When no name was given, take the
    // one from the page's own first example button.
    if (isDefault && window.loadExample.length > 0) {
      name = firstDeclaredExampleName();
      if (name === null) {
        console.warn("[example-link] This page needs a named example; " +
                     'use "?example=<name>".');
        return;
      }
    }

    var input = document.getElementById("inputText");
    var before = input ? input.value : null;

    try {
      if (name === null) window.loadExample();
      else window.loadExample(name);
    } catch (err) {
      console.error("[example-link] Could not load the example:", err);
      return;
    }

    // Never leave a bad assignment sitting in the sequence box: it would be
    // treated as sequence on the next run.
    if (input && (input.value === "undefined" || input.value === "null")) {
      console.error('[example-link] loadExample("' + name + '") produced no data.');
      input.value = before === null ? "" : before;
    }
  }

  /** Read the argument out of the page's first onclick="loadExample('x')". */
  function firstDeclaredExampleName() {
    var nodes = document.querySelectorAll("[onclick*='loadExample']");
    for (var i = 0; i < nodes.length; i++) {
      var m = /loadExample\(\s*['"]([^'"]+)['"]\s*\)/.exec(nodes[i].getAttribute("onclick"));
      if (m) return m[1];
    }
    return null;
  }

  if (document.readyState === "complete") {
    run();
  } else {
    window.addEventListener("load", run);
  }
})();
