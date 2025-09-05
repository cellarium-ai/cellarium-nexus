(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  let _fieldsMeta = null;
  let _fieldsPromise = null;

  ns.readyFieldsMeta = function () {
    if (_fieldsPromise) return _fieldsPromise;
    _fieldsPromise = fetch(ns.URLS.fields, { credentials: 'same-origin' })
      .then(r => r.json())
      .then(data => {
        _fieldsMeta = Array.isArray(data) ? data : (data.fields || []);
        console.debug('[Cell Info Filters][meta] fetched fields meta from API');
        return _fieldsMeta;
      })
      .catch((e) => {
        console.warn('[Cell Info Filters][meta] failed to fetch fields meta from API, falling back to embedded JSON if present', e);
        const script = document.getElementById('filters-fields-meta');
        if (!script) { _fieldsMeta = []; return _fieldsMeta; }
        try {
          const data = JSON.parse(script.textContent);
          _fieldsMeta = Array.isArray(data) ? data : (data.fields || []);
        } catch (err) {
          _fieldsMeta = [];
        }
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
