(function () {
  window.CellInfoFilters = window.CellInfoFilters || {};
  const ns = window.CellInfoFilters;

  ns.URLS = {
    count: '/admin/cell_management/cellinfo/filters/count/',
  };

  ns.OP_LABELS = {
    eq: 'equals',
    not_eq: 'not equals',
    in: 'in',
    not_in: 'not in',
    gt: '>',
    gte: '≥',
    lt: '<',
    lte: '≤',
  };

  console.debug('[Cell Info Filters][constants] init');
})();
