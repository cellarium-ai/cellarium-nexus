from __future__ import annotations

from io import StringIO

from django.core.management import call_command

from cellarium.nexus.backend.cell_management.utils import cache_reset
from cellarium.nexus.backend.core.utils import reset_cache


def test_list_urls_output_includes_custom_admin() -> None:
    """
    Invoke ``list_urls`` management command and ensure custom admin path exists.

    :return: None
    """

    stdout = StringIO()

    call_command("list_urls", stdout=stdout)

    output = stdout.getvalue()
    assert "admin/cell_management/cellinfo/" in output


def test_reset_cache_and_repopulate_delegates(monkeypatch) -> None:
    """
    Ensure ``reset_cache_and_repopulate`` delegates to cell management helper.

    :param monkeypatch: Pytest monkeypatch fixture
    """

    captured: dict[str, int] = {"deleted": 2, "repopulated": 5}

    def _fake_reset(*, repopulate: bool) -> dict[str, int]:
        assert repopulate is True
        return captured

    monkeypatch.setattr(cache_reset, "reset_cellinfo_filter_cache", _fake_reset)

    result = reset_cache.reset_cache_and_repopulate()

    assert result == ["countcache:*", "bqcache:*"]
