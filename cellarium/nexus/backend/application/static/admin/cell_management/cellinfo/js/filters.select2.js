(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

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

    const minChars = typeof fieldMeta.suggest_min_chars === 'number' ? fieldMeta.suggest_min_chars : 1;
    const $ = window.jQuery || (window.django && window.django.jQuery);

    if ($ && typeof $(selectEl).select2 === 'function') {
      $(selectEl).select2({
        width: 'style',
        multiple: true,
        ajax: {
          transport: function (params, success, failure) {
            const term = params.data.term || '';
            if (term.length < minChars) {
              success({ results: [] });
              return;
            }
            const url = ns.URLS.suggest;
            const ds = ns.getDataset();
            const qs = new URLSearchParams({ field: fieldKey, q: term, dataset: ds });
            fetch(url + '?' + qs.toString(), { credentials: 'same-origin' })
              .then(r => r.json())
              .then(data => {
                const results = (data.suggestions || []).map(v => ({ id: v, text: v }));
                success({ results });
              })
              .catch(failure);
          },
          processResults: function (data) { return data; },
        },
        tags: false,
        allowClear: false,
      });
    }
  };

  console.debug('[Cell Info Filters][select2] init');
})();

