/* ════════════════════════════════════════════════════════════════════════════
   Version stamp — single source of truth.

   Bump these two lines when you cut a release, and tag the repository to match.
   Nothing else needs editing: every page carries <span id="version"></span> in
   its footer, and the value is also exposed as
   document.documentElement.dataset.version for anyone inspecting the page.

   Kept separate from the theme logic below, which returns early on pages with
   no navbar.
   ════════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  var DIGITALGENS_VERSION = "1.0";
  var DIGITALGENS_UPDATED = "2026-10-15";      // ISO 8601

  var MONTHS = ["January", "February", "March", "April", "May", "June",
                "July", "August", "September", "October", "November", "December"];

  function pretty(iso) {
    var m = /^(\d{4})-(\d{2})-(\d{2})$/.exec(iso);
    if (!m) return iso;
    return Number(m[3]) + " " + MONTHS[Number(m[2]) - 1] + " " + m[1];
  }

  document.documentElement.setAttribute("data-version", DIGITALGENS_VERSION);

  document.addEventListener("DOMContentLoaded", function () {
    // Most pages set the year in their own inline script; references.html does not.
    // Idempotent: only fills the span if nothing else already has.
    var y = document.getElementById("year");
    if (y && !y.textContent) y.textContent = new Date().getFullYear();

    var el = document.getElementById("version");
    if (!el) return;
    el.textContent = "DigitalGens v" + DIGITALGENS_VERSION +
                     ", updated " + pretty(DIGITALGENS_UPDATED);
  });
})();

/* ════════════════════════════════════════════════════════════════════════════
   theme.js — light / dark switch for the navbar (Bootstrap 5.3 native theming)
   Place in: /assets/js/theme.js
   Load in <head> (without `defer`) AFTER the stylesheet links, e.g.

       <link href="assets/css/site-theme.css" rel="stylesheet">
       <script src="assets/js/theme.js"></script>

   What it does
   ────────────
   1. Sets <html data-bs-theme="…"> immediately (before paint, so no flash).
      • If the visitor chose a theme before  → that choice (localStorage).
      • Otherwise                            → follows the OS setting (auto).
   2. After the DOM loads, drops a sun/moon button into the navbar that flips
      the theme and remembers the choice.
   3. While the visitor has NOT made an explicit choice, it keeps following the
      OS as it changes. Once they click the button, their choice wins.

   No libraries, no markup edits required — just this file and the <script> tag.
   ════════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  var KEY = "theme";          // localStorage key
  var root = document.documentElement;

  function stored() {
    try {
      var v = localStorage.getItem(KEY);
      return v === "light" || v === "dark" ? v : null;
    } catch (e) { return null; }
  }

  function osPrefersDark() {
    return !!(window.matchMedia &&
              window.matchMedia("(prefers-color-scheme: dark)").matches);
  }

  function current() {
    return root.getAttribute("data-bs-theme") === "dark" ? "dark" : "light";
  }

  function apply(theme) {
    root.setAttribute("data-bs-theme", theme === "dark" ? "dark" : "light");
  }

  /* 1 ─ apply the right theme as early as possible (runs while <head> parses) */
  apply(stored() || (osPrefersDark() ? "dark" : "light"));

  /* sun / moon glyphs (Bootstrap Icons paths, MIT) — inherit button colour */
  var SUN  = '<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" fill="currentColor" viewBox="0 0 16 16" aria-hidden="true"><path d="M8 11a3 3 0 1 1 0-6 3 3 0 0 1 0 6m0 1a4 4 0 1 0 0-8 4 4 0 0 0 0 8M8 0a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 0m0 13a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 13m8-5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1 0-1h2a.5.5 0 0 1 .5.5M3 8a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1 0-1h2A.5.5 0 0 1 3 8m10.657-5.657a.5.5 0 0 1 0 .707l-1.414 1.415a.5.5 0 1 1-.707-.708l1.414-1.414a.5.5 0 0 1 .707 0m-9.193 9.193a.5.5 0 0 1 0 .707L3.05 13.657a.5.5 0 0 1-.707-.707l1.414-1.414a.5.5 0 0 1 .707 0m9.193 2.121a.5.5 0 0 1-.707 0l-1.414-1.414a.5.5 0 0 1 .707-.708l1.414 1.415a.5.5 0 0 1 0 .707M4.464 4.465a.5.5 0 0 1-.707 0L2.343 3.05a.5.5 0 1 1 .707-.707l1.414 1.414a.5.5 0 0 1 0 .708"/></svg>';
  var MOON = '<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18" fill="currentColor" viewBox="0 0 16 16" aria-hidden="true"><path d="M6 .278a.77.77 0 0 1 .08.858 7.2 7.2 0 0 0-.878 3.46c0 4.021 3.278 7.277 7.318 7.277q.792-.001 1.533-.16a.79.79 0 0 1 .81.316.73.73 0 0 1-.031.893A8.35 8.35 0 0 1 8.344 16C3.734 16 0 12.286 0 7.71 0 4.266 2.114 1.312 5.124.06A.75.75 0 0 1 6 .278"/></svg>';

  document.addEventListener("DOMContentLoaded", function () {
    var btn = document.getElementById("theme-toggle");

    /* 2 ─ if the page doesn't already provide a button, inject one into the nav */
    if (!btn) {
      var nav = document.querySelector(".navbar .navbar-nav");
      if (!nav) return;                       // no navbar on this page → nothing to do

      var li = document.createElement("li");
      li.className = "nav-item d-flex align-items-center ms-lg-2";

      btn = document.createElement("button");
      btn.id = "theme-toggle";
      btn.type = "button";
      btn.className = "btn btn-sm btn-outline-light border-0 px-2 d-inline-flex align-items-center";

      li.appendChild(btn);
      nav.appendChild(li);
    }

    function render() {
      var dark = current() === "dark";
      // show the icon of the mode you'll switch TO
      btn.innerHTML = dark ? SUN : MOON;
      var label = dark ? "Switch to light theme" : "Switch to dark theme";
      btn.setAttribute("aria-label", label);
      btn.setAttribute("title", label);
    }

    btn.addEventListener("click", function () {
      var next = current() === "dark" ? "light" : "dark";
      apply(next);
      try { localStorage.setItem(KEY, next); } catch (e) {}
      render();
    });

    render();

    /* 3 ─ follow OS changes until the visitor makes an explicit choice */
    if (window.matchMedia) {
      var mq = window.matchMedia("(prefers-color-scheme: dark)");
      var onChange = function (e) {
        if (stored() === null) { apply(e.matches ? "dark" : "light"); render(); }
      };
      if (mq.addEventListener)      mq.addEventListener("change", onChange);
      else if (mq.addListener)      mq.addListener(onChange);   // Safari < 14
    }
  });
})();
