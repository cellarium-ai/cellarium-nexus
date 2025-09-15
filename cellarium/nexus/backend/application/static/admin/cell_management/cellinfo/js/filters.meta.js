(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  let _fieldsMeta = null;
  let _fieldsPromise = null;

  ns.readyFieldsMeta = function () {
    if (_fieldsPromise) return _fieldsPromise;
    // Prefer embedded JSON from the template for zero-latency init
    try {
      const script = document.getElementById('filters-fields-meta');
      if (script && script.textContent) {
        const data = JSON.parse(script.textContent);
        _fieldsMeta = Array.isArray(data) ? data : (data.fields || []);
        _fieldsPromise = Promise.resolve(_fieldsMeta);
        console.debug('[Cell Info Filters][meta] using embedded fields meta');
        return _fieldsPromise;
      }
    } catch (_) { /* noop */ }

    // Fallback to API if embedding is not present
    _fieldsPromise = fetch(ns.URLS.fields, { credentials: 'same-origin' })
      .then(r => r.json())
      .then(data => {
        _fieldsMeta = Array.isArray(data) ? data : (data.fields || []);
        console.debug('[Cell Info Filters][meta] fetched fields meta from API');
        return _fieldsMeta;
      })
      .catch((e) => {
        console.warn('[Cell Info Filters][meta] failed to load fields meta from API and no embedded JSON present', e);
        _fieldsMeta = [];
        return _fieldsMeta;
      });
    return _fieldsPromise;
  };

  ns.getFieldsMeta = function () {
    return Array.isArray(_fieldsMeta) ? _fieldsMeta : [];
  };

  ns.getDataset = function () {
    const selected = document.querySelector('#dataset-selector .dataset-card[aria-selected="true"]');
    return selected ? selected.getAttribute('data-dataset') || '' : '';
  };

  console.debug('[Cell Info Filters][meta] init');
})();
