/* ============================================================================
   paginate.js — Paged.js configuration + Table-of-Contents helper.

   Inlined by build-inline-css.sh into the <head>, BEFORE the Paged.js polyfill
   (Paged.js reads window.PagedConfig at startup). It does two small things:

     1. after():  once Paged.js has laid the document into A4 pages, fill each
                  TOC entry's page number from the page it actually landed on
                  (accurate — no estimation), and bind smooth-scroll on links.
     2. fallback: if the polyfill is absent (file opened without it), still make
                  TOC links clickable.

   What this file NO LONGER does (vs. the old .sheet-era version):
     - No page-count estimation. Page numbers + the per-page footer are produced
       by CSS @page margin-boxes (see Section 9 of report-linear.css), resolved
       by Paged.js via counter(page) / counter(pages).
     - No overflow detection. CSS break-* rules handle pagination.

   NO DEPENDENCIES beyond Paged.js. Runs once.
   ============================================================================ */

(function () {
  "use strict";

  /* Fill TOC page numbers from the pages Paged.js produced. */
  function fillTocNumbers() {
    document.querySelectorAll('a.toc-entry[href^="#"]').forEach(function (link) {
      var id = link.getAttribute('href').slice(1);
      var target = document.getElementById(id);
      var numEl = link.querySelector('.toc-page-num');
      if (!target || !numEl) return;
      var pageEl = target.closest('.pagedjs_page');
      var n = pageEl && pageEl.getAttribute('data-page-number');
      if (n) numEl.textContent = n;
    });
  }

  /* On-screen page chrome: grey backdrop, centered white A4 sheets with a shadow
     and a gap between them, so the browser view matches the printed pagination.
     Injected AFTER Paged.js finishes because Paged.js strips @media screen rules
     from the author CSS. Scoped to @media screen, so print is never affected. */
  function injectPreviewChrome() {
    if (document.getElementById("pagedjs-preview-chrome")) return;
    var s = document.createElement("style");
    s.id = "pagedjs-preview-chrome";
    s.textContent =
      "@media screen{" +
        "html{background:#e9ecef}" +
        ".pagedjs_pages{padding:8mm 0}" +
        ".pagedjs_page{margin:0 auto 8mm;" +
          "box-shadow:0 0 0 1px #dcdfe3,0 6px 18px rgba(0,0,0,.14)}" +
      "}";
    document.head.appendChild(s);
  }

  /* Make TOC entries smooth-scroll to their section. */
  function bindTocLinks() {
    document.querySelectorAll('a.toc-entry[href^="#"]').forEach(function (link) {
      link.addEventListener('click', function (e) {
        var target = document.getElementById(this.getAttribute('href').slice(1));
        if (target) {
          e.preventDefault();
          target.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }
      });
    });
  }

  /* Paged.js reads this object before it paginates. */
  window.PagedConfig = {
    auto: true,
    after: function () {
      injectPreviewChrome();
      fillTocNumbers();
      bindTocLinks();
    }
  };

  /* Graceful fallback: if the polyfill never loads, still bind TOC links. */
  if (typeof window.PagedPolyfill === "undefined") {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", bindTocLinks);
    } else {
      bindTocLinks();
    }
  }
})();
