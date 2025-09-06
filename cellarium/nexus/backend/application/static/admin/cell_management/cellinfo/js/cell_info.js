// Cell Info custom page interactions
// - Handles dataset selection toggling
// - Reflects selected dataset label (if present)
// - Updates mocked cell count from JSON provided via json_script("dataset-counts")

(function () {
  function fmt(n) {
    if (typeof n !== 'number') return '0';
    return n.toLocaleString(undefined);
  }

  function init() {
    const selector = document.getElementById('dataset-selector');
    if (!selector) return;

    const label = document.getElementById('selected-dataset-label');
    const descEl = document.getElementById('selected-dataset-description');
    const descToggle = document.getElementById('dataset-desc-toggle');
    const countEl = document.getElementById('cell-count');

    let counts = {};
    try {
      const dataTag = document.getElementById('dataset-counts');
      if (dataTag) {
        counts = JSON.parse(dataTag.textContent);
      }
    } catch (_) { /* noop */ }

    // Capture canonical className strings from server-rendered buttons
    const initialSelected = selector.querySelector('.dataset-card[aria-selected="true"]');
    const initialUnselected = selector.querySelector('.dataset-card[aria-selected="false"]');
    const selectedClassName = initialSelected ? initialSelected.className : '';
    const unselectedClassName = initialUnselected ? initialUnselected.className : '';

    function needsToggle(el) {
      if (!el) return false;
      if (!el.textContent || !el.textContent.trim()) return false;
      // If content overflows the element's current box, we need a toggle
      return el.scrollHeight > el.clientHeight + 1; // allow for sub-pixel rounding
    }

    function select(btn) {
      const cards = selector.querySelectorAll('.dataset-card');
      cards.forEach(function (card) {
        const isTarget = card === btn;
        card.setAttribute('aria-selected', isTarget ? 'true' : 'false');
        if (selectedClassName && unselectedClassName) {
          card.className = isTarget ? selectedClassName : unselectedClassName;
        }
      });
      const ds = btn && btn.dataset ? btn.dataset.dataset : null;
      const desc = btn && btn.dataset ? (btn.dataset.description || '') : '';
      if (label && ds) label.textContent = ds;
      if (descEl) {
        descEl.textContent = desc;
        descEl.classList.remove('expanded');
      }
      if (descToggle) {
        descToggle.textContent = 'Show more';
        // Defer measurement until next frame to ensure layout updated
        requestAnimationFrame(function () {
          const show = needsToggle(descEl);
          descToggle.style.visibility = show ? 'visible' : 'hidden';
        });
      }
      if (countEl && ds && counts && Object.prototype.hasOwnProperty.call(counts, ds)) {
        countEl.textContent = fmt(counts[ds]);
      }
      // Reset filters on dataset change
      try {
        if (window.CellInfoFilters && typeof window.CellInfoFilters.resetFilters === 'function') {
          window.CellInfoFilters.resetFilters();
        }
      } catch (_) { /* noop */ }
    }

    selector.addEventListener('click', function (e) {
      const btn = e.target.closest('.dataset-card');
      if (!btn) return;
      select(btn);
    });

    if (descToggle && descEl) {
      descToggle.addEventListener('click', function () {
        const isExpanded = descEl.classList.toggle('expanded');
        descToggle.textContent = isExpanded ? 'Show less' : 'Show more';
        // If collapsed again and text actually fits without truncation, hide the toggle
        if (!isExpanded) {
          requestAnimationFrame(function () {
            const show = descEl.scrollHeight > descEl.clientHeight + 1;
            descToggle.style.visibility = show ? 'visible' : 'hidden';
          });
        }
      });
    }

    if (initialSelected) select(initialSelected);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
