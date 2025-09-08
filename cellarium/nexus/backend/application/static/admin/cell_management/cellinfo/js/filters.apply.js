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
    const filters = {};
    for (const row of rows) {
      const index = row.getAttribute('data-index');
      const fieldSel = row.querySelector(`select[name="filters-${index}-field"]`);
      const operSel = row.querySelector(`select[name="filters-${index}-operator"]`);
      const valueName = `filters-${index}-value`;
      const valueEl = row.querySelector(`[name="${valueName}"]`);
      if (!fieldSel || !operSel || !valueEl) continue;

      let fieldKey = fieldSel.value;
      let operator = operSel.value;
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
        // string/categorical
        if (valueEl instanceof HTMLSelectElement) {
          const arr = Array.from(valueEl.selectedOptions).map((o) => o.value).filter((v) => v !== '');
          value = arr;
        } else {
          const v = valueEl.value || '';
          value = v;
        }
      }

      // Normalize to backend format "column__operator": value
      // Adjust operator/value for categorical lists under eq/not_eq
      if (operator === 'eq' || operator === 'not_eq') {
        if (Array.isArray(value)) {
          if (value.length === 0) {
            continue; // nothing to apply
          } else if (value.length === 1) {
            value = value[0];
          } else {
            operator = operator === 'eq' ? 'in' : 'not_in';
          }
        }
      }

      // Ensure membership operators have array values
      if ((operator === 'in' || operator === 'not_in') && !Array.isArray(value)) {
        value = value === '' || value === null || value === undefined ? [] : [String(value)];
      }

      const key = `${fieldKey}__${operator}`;
      filters[key] = value;
    }
    return filters;
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
