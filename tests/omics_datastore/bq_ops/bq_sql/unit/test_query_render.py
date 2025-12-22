from __future__ import annotations

from pathlib import Path

from cellarium.nexus.omics_datastore.bq_ops.bq_sql import query as query_module
from cellarium.nexus.omics_datastore.bq_ops.bq_sql.validation.query_validators import base_query_validator as base_val


def write_tmp_template(tmp_path: Path, content: str) -> str:
    p = tmp_path / "t.sql.mako"
    p.write_text(content)
    return str(p)


def test_render_basic_lowercases_keywords(tmp_path: Path) -> None:
    """
    Verify that render() lowercases SQL keywords after Mako rendering.
    """
    template_path = write_tmp_template(tmp_path=tmp_path, content="SELECT 1")
    td = query_module.TemplateData(project="p")

    out = query_module.render(template_path=template_path, template_data=td)

    assert out.strip() == "select 1"


class FakeValidator(base_val.SQLSyntaxValidator):
    called = False

    @classmethod
    def validate_syntax(cls, sql_query: str) -> None:  # type: ignore[override]
        cls.called = True


def test_render_calls_validator_when_enabled(tmp_path: Path) -> None:
    """
    Ensure that the provided SQLSyntaxValidator is invoked when validation is enabled.
    """
    template_path = write_tmp_template(tmp_path=tmp_path, content="select 1")
    td = query_module.TemplateData(project="p")

    FakeValidator.called = False
    _ = query_module.render(
        template_path=template_path,
        template_data=td,
        sql_query_validator_class=FakeValidator,  # type: ignore[arg-type]
        sql_query_validator_on=True,
    )

    assert FakeValidator.called is True
