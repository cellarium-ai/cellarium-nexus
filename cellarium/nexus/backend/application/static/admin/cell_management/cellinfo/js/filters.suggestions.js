(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Cache of suggestions keyed by dataset then field
  // Structure: { [dataset: string]: { [fieldKey: string]: string[] } }
  const _cache = {};
  let _pending = null;

  ns.readySuggestions = function (dataset) {
    const ds = dataset || ns.getDataset() || '';
    if (ds && _cache[ds]) return Promise.resolve(_cache[ds]);
    if (_pending) return _pending;
    const url = ns.URLS.suggestions_all + '?dataset=' + encodeURIComponent(ds);
    _pending = fetch(url, { credentials: 'same-origin' })
      .then(r => r.json())
      .then(data => {
        const map = (data && data.suggestions) || {};
        _cache[ds] = map;
        _pending = null;
        console.debug('[Cell Info Filters][suggestions] loaded for dataset', ds, map);
        return map;
      })
      .catch((e) => {
        _pending = null;
        console.warn('[Cell Info Filters][suggestions] failed to load', e);
        _cache[ds] = {};
        return _cache[ds];
      });
    return _pending;
  };

  ns.getSuggestions = function (fieldKey, dataset) {
    const ds = dataset || ns.getDataset() || '';
    const map = _cache[ds] || {};
    const arr = map[fieldKey] || [];
    return Array.isArray(arr) ? arr : [];
  };

  console.debug('[Cell Info Filters][suggestions] init');
})();
