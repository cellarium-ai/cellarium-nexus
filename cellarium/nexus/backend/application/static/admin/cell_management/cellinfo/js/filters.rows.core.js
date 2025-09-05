(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Handle field change: updates operator options and value widget
  ns.onFieldChange = function (rowEl, fieldsMeta) {
    const index = parseInt(rowEl.getAttribute('data-index') || '0', 10);
    const fieldSel = rowEl.querySelector('select[name^="filters-"][name$="-field"]');
    const operSel = rowEl.querySelector('select[name^="filters-"][name$="-operator"]');
    if (!fieldSel || !operSel) return;
    const fieldKey = fieldSel.value;
    const meta = ns.findFieldMeta(fieldsMeta, fieldKey);
    let fType = meta ? meta.type : 'string';

    // If backend marks string with suggest=true, treat as categorical (multi-select with AJAX)
    if (meta && (fType === 'string' || fType === 'text') && meta.suggest) {
      fType = 'categorical';
    }
    // Normalize types: default 'string' -> plain text; only special case 'categorical'
    if (fType === 'string' || fType === 'text' || !fType) fType = 'text';

    const allowedOps = meta && Array.isArray(meta.operators) && meta.operators.length ? meta.operators : Object.keys(ns.OP_LABELS);
    console.debug('[Cell Info Filters][rows.core] onFieldChange', { index, fieldKey, resolvedType: fType, allowedOps, meta });

    // Rebuild operator options and refresh Select2; default to 'eq' or first allowed
    const $ = window.jQuery || (window.django && window.django.jQuery);
    if ($ && $(operSel).data('select2')) {
      console.debug('[Cell Info Filters][rows.core] rebuilding operator options (select2 present)');
      $(operSel).empty();
      const frag = ns.buildOperatorOptions(allowedOps);
      const tmp = document.createElement('div');
      tmp.appendChild(frag);
      const opts = Array.from(tmp.childNodes);
      for (const opt of opts) {
        operSel.appendChild(opt);
      }
      const defaultOp = allowedOps.includes('eq') ? 'eq' : (allowedOps[0] || '');
      $(operSel).val(defaultOp).trigger('change.select2');
    } else {
      console.debug('[Cell Info Filters][rows.core] rebuilding operator options');
      operSel.innerHTML = '';
      operSel.appendChild(ns.buildOperatorOptions(allowedOps));
      const defaultOp = allowedOps.includes('eq') ? 'eq' : (allowedOps[0] || '');
      operSel.value = defaultOp;
      ns.initSelect2IfAvailable(operSel.parentElement);
    }

    ns.replaceValueWidget(rowEl, fType, index, meta, fieldKey);
  };

  // Reindex inputs after add/remove
  ns.reindexRows = function (rowsContainer) {
    const rows = Array.from(rowsContainer.querySelectorAll('.filter-row'));
    rows.forEach((row, i) => {
      const oldIndex = row.getAttribute('data-index');
      if (String(oldIndex) === String(i)) return;
      row.setAttribute('data-index', String(i));
      row.querySelectorAll('input, select, textarea').forEach((el) => {
        if (el.name) el.name = el.name.replace(/filters-(\d+)-/g, `filters-${i}-`);
        if (el.id) el.id = el.id.replace(/id_filters-(\d+)-/g, `id_filters-${i}-`);
      });
    });
  };

  // Bind events for a single row
  ns.bindRowEvents = function (rowEl, fieldsMeta) {
    const fieldSel = rowEl.querySelector('select[name^="filters-"][name$="-field"]');
    if (fieldSel) {
      const handler = () => ns.onFieldChange(rowEl, fieldsMeta);
      fieldSel.addEventListener('change', handler);
      const $ = window.jQuery || (window.django && window.django.jQuery);
      if ($) {
        try {
          $(fieldSel).on('select2:select select2:clear', handler);
        } catch (e) { /* ignore */ }
      }
      const idx = rowEl.getAttribute('data-index');
      console.debug('[Cell Info Filters][rows.core] bound field change events', { index: idx, name: fieldSel.name });
    }

    const removeBtn = rowEl.querySelector('.filter-remove');
    if (removeBtn) {
      removeBtn.addEventListener('click', () => {
        const rowsContainer = rowEl.parentElement;
        const formEl = document.getElementById('filters-formset');
        rowEl.remove();
        ns.reindexRows(rowsContainer);
        ns.setTotalForms(formEl, rowsContainer.querySelectorAll('.filter-row').length);
        ns.updateApplyVisibility();
      });
    }
  };

  console.debug('[Cell Info Filters][rows.core] init');
})();
