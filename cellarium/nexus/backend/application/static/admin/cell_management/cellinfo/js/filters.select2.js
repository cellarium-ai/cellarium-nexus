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
        // Keep dropdown open after each selection for faster multi-pick UX
        closeOnSelect: false,
        ajax: {
          cache: false,
          transport: function (params, success /*, failure */) {
            const term = (params && params.data && params.data.term) ? String(params.data.term) : '';
            const DS_MAP = resolveDsMap();
            const arr = (DS_MAP && DS_MAP[fieldKey]) ? DS_MAP[fieldKey] : [];
            // Exclude currently selected values to satisfy: hide once selected, reappear if deselected
            // Use jQuery value (more reliable with Select2's async state than reading <option selected>)
            const $jq = window.jQuery || (window.django && window.django.jQuery);
            const selectedVals = ($jq && typeof $jq.fn.val === 'function') ? ($jq(selectEl).val() || []) : Array.from(selectEl.options).filter(o => o.selected).map(o => o.value);
            // Also consult Select2's data immediately (handles timing where .val() hasn't updated yet)
            let selectedFromData = [];
            try {
              if ($jq && typeof $jq(selectEl).select2 === 'function') {
                selectedFromData = ($jq(selectEl).select2('data') || []).map(x => x && (x.id != null ? String(x.id) : String(x.text)));
              }
            } catch (_) { /* noop */ }
            const selected = new Set([...(selectedVals || []).map(v => String(v)), ...selectedFromData]);
            const t = term.toLowerCase();
            const base = t ? arr.filter(v => String(v).toLowerCase().includes(t)) : arr.slice();
            const filtered = base.filter(v => !selected.has(String(v)));
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

      // Force results refresh in the still-open dropdown after each select/unselect
      // Minimal and safe: simulate input in the search field to re-run the transport
      try {
        $(selectEl).on('select2:select select2:unselect', function () {
          const $el = $(this);
          const inst = $el.data('select2');
          // Try multiple locations for the search field depending on theme/mode
          let $search = null;
          if (inst && inst.dropdown && inst.dropdown.$search) {
            $search = inst.dropdown.$search;
          } else if (inst && inst.selection && inst.selection.$search) {
            $search = inst.selection.$search;
          } else if (inst && inst.$dropdown) {
            $search = inst.$dropdown.find('.select2-search__field');
          }
          // Nudge the search input to retrigger transport filtering
          if ($search && $search.length) {
            const current = $search.val();
            // Staggered triggers to avoid debounce swallowing updates
            setTimeout(function () {
              $search.val(String(current || '') + ' ').trigger('input').trigger('keyup');
            }, 0);
            setTimeout(function () {
              $search.val(current || '').trigger('input').trigger('keyup');
            }, 20);
          }
          // Also explicitly re-open (keeps open, but forces results render in some Select2 themes)
          setTimeout(function () {
            try { $el.select2('open'); } catch (_) { /* noop */ }
          }, 0);

          // As an immediate fallback, prune any options in the currently visible
          // dropdown that are now selected, so they disappear without needing a
          // transport re-query (covers themes/modes without a live search box).
          try {
            const selectedNow = new Set((($el.val() || [])).map(String));
            if (inst && inst.$dropdown) {
              inst.$dropdown.find('.select2-results__option').each(function () {
                const $opt = $(this);
                const data = $opt.data('data');
                const id = data && (data.id != null ? String(data.id) : String(data.text || ''));
                if (id && selectedNow.has(id)) {
                  $opt.remove();
                }
              });
              // Additionally remove any option flagged selected by Select2 ARIA state
              inst.$dropdown.find('.select2-results__option[aria-selected="true"]').remove();
            }
          } catch (_) { /* noop */ }
        });
      } catch (e) { /* noop */ }
    }
  };

  console.debug('[Cell Info Filters][select2] init');
})();

