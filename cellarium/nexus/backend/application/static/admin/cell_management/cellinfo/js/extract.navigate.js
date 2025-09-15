(function () {
  'use strict';

  function onClickExtract(e, btn) {
    try {
      var ns = window.CellInfoFilters || {};
      var dataset = (typeof ns.getDataset === 'function') ? ns.getDataset() : '';
      var filters = (typeof ns.collectFiltersPayload === 'function') ? ns.collectFiltersPayload() : {};
      var href = btn.getAttribute('href') || window.location.href;
      var url = new URL(href, window.location.origin);
      if (dataset) url.searchParams.set('dataset', dataset);
      if (filters && typeof filters === 'object' && Object.keys(filters).length > 0) {
        url.searchParams.set('filters', JSON.stringify(filters));
      }
      window.location.href = url.toString();
    } catch (_) {
      // Fallback: let the default link work
      return;
    }
    e.preventDefault();
  }

  function init() {
    var btn = document.getElementById('extract-btn');
    if (!btn) return;
    btn.addEventListener('click', function (e) { onClickExtract(e, btn); });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
