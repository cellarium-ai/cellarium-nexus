(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Load precomputed suggestions from the template JSON script tag
  let PRECOMPUTED_SUGGESTIONS = {};
  let PRECOMPUTED_SUGGESTIONS_ALL = {};
  try {
    const tag = document.getElementById('precomputed-suggestions');
    if (tag && tag.textContent) {
      PRECOMPUTED_SUGGESTIONS = JSON.parse(tag.textContent) || {};
    }
    const tagAll = document.getElementById('precomputed-suggestions-all');
    if (tagAll && tagAll.textContent) {
      PRECOMPUTED_SUGGESTIONS_ALL = JSON.parse(tagAll.textContent) || {};
    }
  } catch (_) { /* noop */ }

  ns.initSelect2IfAvailable = function (container) {
    try {
      if (typeof window.unfoldInitSelect2 === 'function') {
        window.unfoldInitSelect2(container || document);
        // Sanitize field/operator selects: disable allowClear and placeholder if Unfold enabled them
        const $jq = window.jQuery || (window.django && window.django.jQuery);
        if ($jq) {
          $jq(container || document).find('select.admin-autocomplete').each(function () {
            const $el = $jq(this);
            const nameAttr = String($el.attr('name') || '');
            const isFieldOrOperator = /filters-\d+-(field|operator)$/.test(nameAttr);
            const isMultiple = !!$el.prop('multiple');
            if ($el.data('select2')) {
              try { $el.select2('destroy'); } catch (e) { /* noop */ }
            }
            if (isFieldOrOperator) {
              $el.select2({ width: 'style' });
            } else if (!isMultiple) {
              // value single-select (e.g., boolean) — keep placeholder, no allowClear
              let placeholder = '';
              const $first = $el.find('option').first();
              if ($first.length && ($first.val() === '' || $first.val() === null)) {
                placeholder = $first.text() || '';
              }
              $el.select2({ width: 'style', placeholder: placeholder || 'Select value', allowClear: false });
            } else {
              // categorical multi
              $el.select2({ width: 'style' });
            }
          });
        }
        return;
      }
      if (typeof window.initSelect2 === 'function') {
        window.initSelect2(container || document);
        return;
      }
    } catch (e) {
      // fall through
    }

    const $ = window.jQuery || (window.django && window.django.jQuery);
    if ($) {
      $(container || document).find('select.admin-autocomplete').each(function () {
        const $el = $(this);
        if ($el.data('select2')) return;
        const isMultiple = !!$el.prop('multiple');
        const nameAttr = String($el.attr('name') || '');
        const isFieldOrOperator = /filters-\d+-(field|operator)$/.test(nameAttr);
        // Derive placeholder from first empty option if present
        let placeholder = '';
        const $first = $el.find('option').first();
        if ($first.length && ($first.val() === '' || $first.val() === null)) {
          placeholder = $first.text() || '';
        }
        const baseOpts = { width: 'style' };
        const singleOpts = isMultiple
          ? {}
          : (isFieldOrOperator
              ? { width: 'style' }
              : { placeholder: placeholder || 'Select value', allowClear: false });
        $el.select2({ ...baseOpts, ...singleOpts });
      });
    }
  };

  ns.configureCategoricalSelect = function (selectEl, fieldMeta, fieldKey) {
    selectEl.setAttribute('multiple', 'multiple');
    selectEl.classList.add('unfold-admin-autocomplete', 'admin-autocomplete');
    selectEl.setAttribute('data-theme', 'admin-autocomplete');

    const $ = window.jQuery || (window.django && window.django.jQuery);
    if ($ && typeof $(selectEl).select2 === 'function') {
      // Resolve the suggestions map per current dataset if ALL is present
      function resolveDsMap() {
        const ds = (ns.getDataset && ns.getDataset()) || '';
        if (PRECOMPUTED_SUGGESTIONS_ALL && ds && PRECOMPUTED_SUGGESTIONS_ALL[ds]) {
          return PRECOMPUTED_SUGGESTIONS_ALL[ds] || {};
        }
        return PRECOMPUTED_SUGGESTIONS || {};
      }
      $(selectEl).select2({
        width: 'style',
        multiple: true,
        ajax: {
          transport: function (params, success /*, failure */) {
            const term = (params && params.data && params.data.term) ? String(params.data.term) : '';
            const DS_MAP = resolveDsMap();
            const arr = (DS_MAP && DS_MAP[fieldKey]) ? DS_MAP[fieldKey] : [];
            const t = term.toLowerCase();
            const filtered = t ? arr.filter(v => String(v).toLowerCase().includes(t)) : arr.slice();
            const limited = filtered.slice(0, 100);
            const results = limited.map(v => ({ id: v, text: v }));
            success({ results });
          },
          processResults: function (data) { return data; },
        },
        minimumInputLength: 0,
        allowClear: false,
        tags: false,
      });
    }
  };

  console.debug('[Cell Info Filters][select2] init');
})();

