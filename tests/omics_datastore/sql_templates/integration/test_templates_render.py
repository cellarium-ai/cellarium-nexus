from __future__ import annotations

from pathlib import Path

from cellarium.nexus.omics_datastore.bq_ops.bq_sql import query as query_module
from cellarium.nexus.omics_datastore.bq_ops.bq_sql import template_data as td_module

TEMPLATES_ROOT = (
    Path(__file__).resolve().parents[3] / "cellarium" / "nexus" / "omics_datastore" / "bq_ops" / "sql_templates"
)


def minimal_kwargs_for_templates() -> dict[str, dict[str, object]]:
    """
    Provide minimal keyword arguments for known templates to render successfully.

    Keys are relative paths from the sql_templates root. Values are kwargs passed as other_kwargs
    to TemplateData so that Mako templates can render without missing variables.
    """
    default = {
        "table_name": "cell_info",
        "filter_statements": {},
        "column_name": "c.organism",
        "column_names": ["c.id"],
        "extract_table_prefix": "",
        "select_columns": ["c.id"],
        "start_bin": 0,
        "end_bin": 1,
        "limit": 10,
    }

    # You can override per-template here if any template needs special values
    overrides: dict[str, dict[str, object]] = {}

    # Build full mapping using defaults for all templates; apply overrides when present
    mapping: dict[str, dict[str, object]] = {}
    for p in sorted(TEMPLATES_ROOT.rglob("*.sql.mako")):
        rel = str(p.relative_to(TEMPLATES_ROOT))
        mapping[rel] = {**default, **overrides.get(rel, {})}

    return mapping


def test_all_sql_templates_render() -> None:
    """
    Render all SQL Mako templates with minimal TemplateData and ensure rendering succeeds.
    """
    manifest = minimal_kwargs_for_templates()

    for rel_path, kwargs in manifest.items():
        path = TEMPLATES_ROOT / rel_path
        td = td_module.TemplateData(project="dummy-project", dataset="dummy_dataset", **kwargs)
        rendered = query_module.render(
            template_path=str(path),
            template_data=td,
            sql_query_validator_on=False,
        )
        assert isinstance(rendered, str) and rendered.strip(), f"Rendered SQL empty for template: {rel_path}"
