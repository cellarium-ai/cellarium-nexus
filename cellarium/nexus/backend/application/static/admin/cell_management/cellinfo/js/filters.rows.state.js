(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Initialize already-rendered rows, bind events, and init Select2
  ns.initExistingRows = function (fieldsMeta) {
    const rowsContainer = document.getElementById('filters-rows');
    const rows = rowsContainer ? Array.from(rowsContainer.querySelectorAll('.filter-row')) : [];
    rows.forEach((row) => {
      ns.bindRowEvents(row, fieldsMeta);
      ns.onFieldChange(row, fieldsMeta);
    });
    if (rowsContainer) ns.initSelect2IfAvailable(rowsContainer);
    ns.updateApplyVisibility();
  };

  // Formset helpers
  ns.getTotalForms = function (formsetEl) {
    const totalInput = formsetEl.querySelector('input[name="filters-TOTAL_FORMS"]');
    return totalInput ? parseInt(totalInput.value || '0', 10) : 0;
  };

  ns.setTotalForms = function (formsetEl, total) {
    const totalInput = formsetEl.querySelector('input[name="filters-TOTAL_FORMS"]');
    if (totalInput) totalInput.value = String(total);
  };

  // Toggle Apply button visibility based on rows count
  ns.updateApplyVisibility = function () {
    const rowsContainer = document.getElementById('filters-rows');
    const applyBtn = document.getElementById('filters-apply');
    if (!rowsContainer || !applyBtn) return;
    const count = rowsContainer.querySelectorAll('.filter-row').length;
    if (count > 0) {
      applyBtn.classList.remove('hidden');
      applyBtn.style.display = '';
    } else {
      applyBtn.classList.add('hidden');
      applyBtn.style.display = 'none';
    }
  };

  // Add a new row from the empty form template and initialize it
  ns.addNewRow = function (fieldsMeta) {
    const tmpl = document.getElementById('filters-empty-form-template');
    const rowsContainer = document.getElementById('filters-rows');
    const formEl = document.getElementById('filters-formset');
    if (!tmpl || !rowsContainer || !formEl) return;

    const currentTotal = ns.getTotalForms(formEl);
    const nextIndex = currentTotal;

    const html = tmpl.innerHTML.replace(/__prefix__/g, String(nextIndex));
    const wrap = document.createElement('div');
    wrap.innerHTML = html.trim();
    const rowEl = wrap.firstElementChild;
    rowsContainer.appendChild(rowEl);

    rowEl.setAttribute('data-index', String(nextIndex));

    ns.bindRowEvents(rowEl, fieldsMeta);
    // Default field to first non-empty option and trigger change
    const fieldSel = rowEl.querySelector('select[name^="filters-"][name$="-field"]');
    if (fieldSel) {
      const $ = window.jQuery || (window.django && window.django.jQuery);
      const firstVal = ns.firstNonEmptyOption(fieldSel);
      if ($) {
        $(fieldSel).val(firstVal).trigger('change');
      } else {
        fieldSel.value = firstVal;
        fieldSel.dispatchEvent(new Event('change'));
      }
    } else {
      ns.onFieldChange(rowEl, fieldsMeta);
    }

    ns.setTotalForms(formEl, nextIndex + 1);

    ns.initSelect2IfAvailable(rowEl);
    // Re-apply field change after Select2 init to ensure operator/value widgets match the selected field
    ns.onFieldChange(rowEl, fieldsMeta);
    ns.updateApplyVisibility();
  };

  // Clear all rows and reset the formset counter
  ns.resetFilters = function () {
    const rowsContainer = document.getElementById('filters-rows');
    const formEl = document.getElementById('filters-formset');
    if (!rowsContainer || !formEl) return;
    const $ = window.jQuery || (window.django && window.django.jQuery);
    if ($) {
      // Destroy any select2 instances within rows before removing
      $(rowsContainer).find('select').each(function () {
        const $sel = $(this);
        if ($sel.data('select2')) {
          try { $sel.select2('destroy'); } catch (e) { /* noop */ }
        }
      });
    }
    while (rowsContainer.firstChild) rowsContainer.removeChild(rowsContainer.firstChild);
    ns.setTotalForms(formEl, 0);
    console.debug('[Cell Info Filters][rows.state] resetFilters: cleared all rows');
    ns.updateApplyVisibility();
  };

  console.debug('[Cell Info Filters][rows.state] init');
})();
