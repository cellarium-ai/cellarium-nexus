(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  // Create value widget from prototypes or fallback programmatically
  ns.cloneValuePrototype = function (type) {
    const protoRoot = document.getElementById('filters-prototypes');
    if (!protoRoot) {
      // Programmatic fallbacks
      if (type === 'number') {
        const inp = document.createElement('input');
        inp.type = 'number';
        inp.step = 'any';
        inp.className = 'vTextField w-full h-9 px-3 py-2 rounded border border-base-200 dark:border-base-700 bg-white dark:bg-base-900 text-base-700 dark:text-base-200 placeholder-base-400 focus:outline-none focus:ring focus:ring-primary-300 focus:border-primary-600 shadow-sm appearance-none';
        inp.style.width = '100%';
        return inp;
      }
      if (type === 'boolean') {
        const sel = document.createElement('select');
        const optBlank = document.createElement('option'); optBlank.value = ''; optBlank.textContent = 'Select value';
        const optT = document.createElement('option'); optT.value = 'true'; optT.textContent = 'True';
        const optF = document.createElement('option'); optF.value = 'false'; optF.textContent = 'False';
        sel.append(optBlank, optT, optF);
        sel.className = 'unfold-admin-autocomplete admin-autocomplete';
        sel.setAttribute('data-theme', 'admin-autocomplete');
        sel.style.width = '100%';
        return sel;
      }
      if (type === 'categorical') {
        const sel = document.createElement('select');
        sel.multiple = true;
        sel.className = 'unfold-admin-autocomplete admin-autocomplete';
        sel.setAttribute('data-theme', 'admin-autocomplete');
        sel.style.width = '100%';
        return sel;
      }
      const txt = document.createElement('input');
      txt.type = 'text';
      txt.className = 'vTextField w-full h-9 px-3 py-2 rounded border border-base-200 dark:border-base-700 bg-white dark:bg-base-900 text-base-700 dark:text-base-200 placeholder-base-400 focus:outline-none focus:ring focus:ring-primary-300 focus:border-primary-600 shadow-sm appearance-none';
      txt.style.width = '100%';
      return txt;
    }

    // From prototypes
    let id = null;
    if (type === 'number') id = 'filter-proto-number';
    else if (type === 'boolean') id = 'filter-proto-boolean';
    else if (type === 'categorical') id = 'filter-proto-categorical';
    else id = 'filter-proto-text';
    const proto = protoRoot.querySelector('#' + id);
    if (proto) return proto.cloneNode(true);

    // Fallback if prototype missing
    if (type === 'number') {
      const inp = document.createElement('input');
      inp.type = 'number';
      inp.step = 'any';
      inp.className = 'vTextField w-full h-9 px-3 py-2 rounded border border-base-200 dark:border-base-700 bg-white dark:bg-base-900 text-base-700 dark:text-base-200 placeholder-base-400 focus:outline-none focus:ring focus:ring-primary-300 focus:border-primary-600 shadow-sm appearance-none';
      inp.style.width = '100%';
      return inp;
    }
    if (type === 'boolean') {
      const sel = document.createElement('select');
      const optBlank = document.createElement('option'); optBlank.value = ''; optBlank.textContent = 'Select value';
      const optT = document.createElement('option'); optT.value = 'true'; optT.textContent = 'True';
      const optF = document.createElement('option'); optF.value = 'false'; optF.textContent = 'False';
      sel.append(optBlank, optT, optF);
      sel.className = 'unfold-admin-autocomplete admin-autocomplete';
      sel.setAttribute('data-theme', 'admin-autocomplete');
      sel.style.width = '100%';
      return sel;
    }
    if (type === 'categorical') {
      const sel = document.createElement('select');
      sel.multiple = true;
      sel.className = 'unfold-admin-autocomplete admin-autocomplete';
      sel.setAttribute('data-theme', 'admin-autocomplete');
      sel.style.width = '100%';
      return sel;
    }
    const txt = document.createElement('input');
    txt.type = 'text';
    txt.className = 'vTextField w-full h-9 px-3 py-2 rounded border border-base-200 dark:border-base-700 bg-white dark:bg-base-900 text-base-700 dark:text-base-200 placeholder-base-400 focus:outline-none focus:ring focus:ring-primary-300 focus:border-primary-600 shadow-sm appearance-none';
    txt.style.width = '100%';
    return txt;
  };

  // Replace value cell content with a new widget and initialize select2
  ns.replaceValueWidget = function (rowEl, type, index, fieldMeta, fieldKey) {
    const cell = rowEl.querySelector('.value-cell');
    if (!cell) return;

    // Destroy any select2 inside the value cell before replacing
    const $ = window.jQuery || (window.django && window.django.jQuery);
    if ($) {
      $(cell).find('select').each(function () {
        const $sel = $(this);
        if ($sel.data('select2')) {
          try { $sel.select2('destroy'); console.debug('[Cell Info Filters][widgets] destroyed previous Select2 in value cell'); } catch (e) { /* ignore */ }
        }
      });
    }

    while (cell.firstChild) cell.removeChild(cell.firstChild);

    const input = ns.cloneValuePrototype(type);
    if (!input) {
      const txt = document.createElement('input');
      txt.type = 'text';
      txt.style.width = '100%';
      cell.appendChild(txt);
      return;
    }

    const name = `filters-${index}-value`;
    const id = `id_filters-${index}-value`;
    input.name = name;
    input.id = id;

    cell.appendChild(input);
    console.debug('[Cell Info Filters][widgets] appended value widget', { index, type, fieldKey });

    if (type === 'categorical') {
      ns.configureCategoricalSelect(input, fieldMeta || {}, fieldKey);
    }

    ns.initSelect2IfAvailable(cell);
    console.debug('[Cell Info Filters][widgets] initialized select2 if available for value cell', { index, type });
  };

  console.debug('[Cell Info Filters][widgets] init');
})();
