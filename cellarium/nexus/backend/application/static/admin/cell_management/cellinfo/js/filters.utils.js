(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Build <option> list from allowed operator keys
  ns.buildOperatorOptions = function (allowed) {
    const frag = document.createDocumentFragment();
    for (const op of allowed) {
      const opt = document.createElement('option');
      opt.value = op;
      opt.textContent = ns.OP_LABELS[op] || op;
      frag.appendChild(opt);
    }
    return frag;
  };

  // Find field meta by key
  ns.findFieldMeta = function (fieldsMeta, key) {
    return fieldsMeta.find(f => f.key === key) || null;
  };

  // Return first non-empty option value from a select
  ns.firstNonEmptyOption = function (selectEl) {
    if (!(selectEl instanceof HTMLSelectElement)) return '';
    for (let i = 0; i < selectEl.options.length; i += 1) {
      const opt = selectEl.options[i];
      if (opt && opt.value !== '') return opt.value;
    }
    return '';
  };

  // Ensure the first option is an empty placeholder for value widgets only
  ns.ensureEmptyPlaceholder = function (selectEl, text) {
    if (!(selectEl instanceof HTMLSelectElement)) return;
    const first = selectEl.options[0];
    if (!first || first.value !== '') {
      const opt = document.createElement('option');
      opt.value = '';
      opt.textContent = text || 'Select value';
      selectEl.insertBefore(opt, selectEl.firstChild);
    } else if (!first.textContent) {
      first.textContent = text || 'Select value';
    }
  };

  console.debug('[Cell Info Filters][utils] init');
})();
