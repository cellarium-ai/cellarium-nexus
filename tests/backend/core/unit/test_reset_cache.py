import types

import pytest

from cellarium.nexus.backend.core.utils import reset_cache as reset_cache_module


def test_reset_cache_calls_cell_management_and_returns_summary(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Assert reset_cache_and_repopulate() delegates to cell_management cache reset
    with repopulate=True and returns the expected summary list.
    """
    calls = types.SimpleNamespace(count=0, args=None, kwargs=None)

    def _fake_reset_cellinfo_filter_cache(*args, **kwargs):  # noqa: D401
        """Record call args/kwargs and return dummy stats."""
        calls.count += 1
        calls.args = args
        calls.kwargs = kwargs
        return {"deleted": ["countcache:*", "bqcache:*"]}

    monkeypatch.setattr(
        reset_cache_module.cm_cache_reset,
        "reset_cellinfo_filter_cache",
        _fake_reset_cellinfo_filter_cache,
        raising=True,
    )

    result = reset_cache_module.reset_cache_and_repopulate()

    assert calls.count == 1
    assert calls.kwargs == {"repopulate": True}
    assert result == ["countcache:*", "bqcache:*"]
