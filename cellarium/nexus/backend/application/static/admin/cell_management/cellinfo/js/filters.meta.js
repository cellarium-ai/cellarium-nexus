(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  let _fieldsMeta = null;
  let _fieldsMetaAll = null;
  let _fieldsPromise = null;

  // Load all fields metadata from embedded JSON
  function _loadFieldsMetaAll() {
    if (_fieldsMetaAll !== null) return;
    try {
      const script = document.getElementById('filters-fields-meta-all');
      if (script && script.textContent) {
        _fieldsMetaAll = JSON.parse(script.textContent) || {};
        console.debug('[Cell Info Filters][meta] loaded fields meta for all datasets');
      }
    } catch (_) {
      _fieldsMetaAll = {};
    }
  }

  ns.readyFieldsMeta = function () {
    if (_fieldsPromise) return _fieldsPromise;
    // Load all datasets metadata first
    _loadFieldsMetaAll();
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

  /**
   * Get fields metadata for a specific dataset.
   * @param {string} datasetName - Name of the dataset
   * @returns {Array} - Array of field metadata objects
   */
  ns.getFieldsMetaForDataset = function (datasetName) {
    _loadFieldsMetaAll();
    if (_fieldsMetaAll && datasetName && _fieldsMetaAll[datasetName]) {
      return _fieldsMetaAll[datasetName];
    }
    return ns.getFieldsMeta();
  };

  /**
   * Update the active fields metadata when switching datasets.
   * @param {string} datasetName - Name of the newly selected dataset
   * @returns {Array} - The updated fields metadata
   */
  ns.updateFieldsMeta = function (datasetName) {
    _loadFieldsMetaAll();
    if (_fieldsMetaAll && datasetName && _fieldsMetaAll[datasetName]) {
      _fieldsMeta = _fieldsMetaAll[datasetName];
      console.debug('[Cell Info Filters][meta] updated fields meta for dataset:', datasetName);
    } else {
      console.warn('[Cell Info Filters][meta] no fields meta found for dataset:', datasetName);
    }
    return _fieldsMeta || [];
  };

  ns.getDataset = function () {
    const selected = document.querySelector('#dataset-selector .dataset-card[aria-selected="true"]');
    return selected ? selected.getAttribute('data-dataset') || '' : '';
  };

  console.debug('[Cell Info Filters][meta] init');
})();
