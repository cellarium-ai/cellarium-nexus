// Bootstrap only: wire split modules and log init
(function () {
  function init() {
    const ns = window.CellInfoFilters;
    if (!ns) return;

    const block = document.getElementById('filters-block');
    if (!block) return;

    ns.readyFieldsMeta().then((fieldsMeta) => {
      ns.initExistingRows(fieldsMeta);
      // Start with zero filters on initial load
      ns.resetFilters();

      const addBtn = document.getElementById('filters-add');
      if (addBtn) {
        addBtn.addEventListener('click', () => ns.addNewRow(fieldsMeta));
      }

      const applyBtn = document.getElementById('filters-apply');
      if (applyBtn) {
        applyBtn.addEventListener('click', () => {
          const dataset = ns.getDataset();
          const filters = ns.collectFiltersPayload();
          const payload = { dataset, filters };
          const countEl = document.getElementById('cell-count');
          if (countEl) countEl.textContent = '…';
          applyBtn.disabled = true;
          ns.postCount(payload)
            .then((data) => {
              const count = (data && typeof data.count === 'number') ? data.count : 0;
              if (countEl) {
                try {
                  countEl.textContent = Number(count).toLocaleString(undefined);
                } catch (_) {
                  countEl.textContent = String(count);
                }
              }
            })
            .catch((err) => {
              console.error('[Cell Info Filters] Count request failed:', err);
              if (countEl) countEl.textContent = '—';
            })
            .finally(() => {
              applyBtn.disabled = false;
            });
        });
      }

      console.debug('[Cell Info Filters][bootstrap] init');
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();