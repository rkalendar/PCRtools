/* ════════════════════════════════════════════════════════════════════════════
   Embed mode — lets a tool page be shown inside the combined FastPCR workbench
   (fastpcr.html) without any change to the tool's own logic.
   ----------------------------------------------------------------------------
   Every tool page defines the SAME global functions (analysis(), Seq(), …) and
   the SAME element ids (inputText, mintm, maxtm, repeats, mainTabs, …), so the
   tools can never share one document — they would overwrite each other. The
   workbench therefore loads each tool in its own <iframe>, which keeps every
   tool's globals and ids isolated exactly as they are today.

   This script only adapts the PRESENTATION when a page is loaded with ?embed=1:
     • hides the site chrome (navbar, page hero, footer) so the tool's panel
       blends into the parent page,
     • reports its content height to the parent so the iframe can grow to fit
       (no inner scrollbar).
   Loaded WITHOUT defer and doing nothing at all on a normal (non-embedded)
   page load, so standalone tool pages are completely unaffected.
   ════════════════════════════════════════════════════════════════════════════ */
(function () {
    'use strict';
    if (!/[?&]embed=1(?:&|$)/.test(window.location.search)) return;   // normal page load → no-op

    // Hide the chrome as early as possible (before first paint) to avoid a flash
    // of the duplicated navbar/footer inside the frame.
    const style = document.createElement('style');
    style.textContent =
        'nav.navbar, header.page-hero, footer.site-footer { display: none !important; }' +
        'main#main { padding-top: 0 !important; padding-bottom: 0 !important; }' +
        '.content-surface { border: 0 !important; box-shadow: none !important; }' +
        'body { background: transparent !important; overflow-x: hidden; }';
    (document.head || document.documentElement).appendChild(style);

    // Tell the parent how tall we are, so it can size the iframe and there is a
    // single page-level scrollbar instead of one per tool.
    //
    // Measure the BODY, never documentElement: documentElement.scrollHeight is at
    // least the viewport height, and the viewport IS the iframe we are resizing —
    // so using it feeds back on itself and the frame grows without bound. The body
    // height reflects the real content and stays stable as the frame is resized.
    let last = 0;
    const report = () => {
        const b = document.body;
        if (!b) return;
        const h = Math.ceil(Math.max(b.scrollHeight, b.offsetHeight, b.getBoundingClientRect().height));
        if (h > 0 && Math.abs(h - last) > 1) {          // ignore sub-pixel jitter
            last = h;
            window.parent.postMessage({ type: 'fastpcr:height', height: h }, '*');
        }
    };

    document.addEventListener('DOMContentLoaded', report);
    window.addEventListener('load', report);
    // The tools rewrite big <textarea>/result panes on every run, and Bootstrap tab
    // switches change the visible height — observe the body instead of polling.
    if (window.ResizeObserver) {
        document.addEventListener('DOMContentLoaded', () => {
            new ResizeObserver(report).observe(document.body);
        });
    }
    // Bootstrap tab shows inside the tool (Input / Report / Results) resize it too.
    document.addEventListener('shown.bs.tab', report);
})();
