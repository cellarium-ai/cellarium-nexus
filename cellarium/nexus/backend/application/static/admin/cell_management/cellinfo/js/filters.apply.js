(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  function getCookie(name) {
    const cookies = document.cookie ? document.cookie.split('; ') : [];
    for (const c of cookies) {
      const [k, ...rest] = c.split('=');
      if (decodeURIComponent(k) === name) {
        return decodeURIComponent(rest.join('='));
      }
    }
    return null;
  }

  ns.collectFiltersPayload = function () {
    const rows = Array.from(document.querySelectorAll('#filters-rows .filter-row'));
    const fieldsMeta = ns.getFieldsMeta();
    const payload = [];
    for (const row of rows) {
      const index = row.getAttribute('data-index');
      const fieldSel = row.querySelector(`select[name="filters-${index}-field"]`);
      const operSel = row.querySelector(`select[name="filters-${index}-operator"]`);
      const valueName = `filters-${index}-value`;
      const valueEl = row.querySelector(`[name="${valueName}"]`);
      if (!fieldSel || !operSel || !valueEl) continue;

      const fieldKey = fieldSel.value;
      const operator = operSel.value;
      const meta = ns.findFieldMeta(fieldsMeta, fieldKey);
      const type = meta ? meta.type : 'string';

      let value = null;
      if (type === 'number') {
        const v = parseFloat(valueEl.value);
        if (!Number.isFinite(v)) continue; // skip invalid
        value = v;
      } else if (type === 'boolean') {
        const v = String(valueEl.value).toLowerCase();
        if (v !== 'true' && v !== 'false') continue;
        value = v === 'true';
      } else {
        // categorical multi-select
        if (valueEl instanceof HTMLSelectElement) {
          value = Array.from(valueEl.selectedOptions).map(o => o.value);
        } else {
          value = valueEl.value ? [valueEl.value] : [];
        }
        if (!Array.isArray(value)) value = [];
      }

      payload.push({ field: fieldKey, operator, value });
    }
    return payload;
  };

  ns.postCount = function (payload) {
    const csrftoken = getCookie('csrftoken');
    return fetch(ns.URLS.count, {
      method: 'POST',
      credentials: 'same-origin',
      headers: {
        'Content-Type': 'application/json',
        ...(csrftoken ? { 'X-CSRFToken': csrftoken } : {}),
      },
      body: JSON.stringify(payload),
    }).then((r) => {
      if (!r.ok) throw new Error('Count request failed with status ' + r.status);
      return r.json();
    });
  };

  console.debug('[Cell Info Filters][apply] init');
})();
