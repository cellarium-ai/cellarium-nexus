from django.contrib import admin
from django.core.files.base import ContentFile
from django.forms.models import BaseInlineFormSet
from django.http import HttpRequest
from django.utils.text import slugify
from django.utils.timezone import now
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin, StackedInline, TabularInline

from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.admin import forms
from cellarium.nexus.backend.ingest_management.utils import soma_schema_utils


class SomaVarSchemaInlineFormSet(BaseInlineFormSet):
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop("request", None)
        kwargs.pop("per_page", None)
        super().__init__(*args, **kwargs)

    def clean(self):
        super().clean()
        if any(self.errors):
            return

        if self.instance.pk:
            return

        valid_forms = [form for form in self.forms if not self._should_delete_form(form)]
        has_csv = any(form.cleaned_data.get("csv_file") for form in valid_forms)
        if not has_csv and valid_forms:
            valid_forms[0].add_error("csv_file", _("Upload a CSV file to create the var schema."))


class SomaVarSchemaInline(StackedInline):
    model = models.SomaVarSchema
    form = forms.SomaVarSchemaInlineForm
    formset = SomaVarSchemaInlineFormSet
    extra = 1
    max_num = 1
    can_delete = False
    fields = (
        "allow_subsets",
        "csv_file",
        "feature_count",
        "var_columns",
        "created_at",
        "updated_at",
    )
    readonly_fields = ("feature_count", "var_columns", "created_at", "updated_at")
    verbose_name = _("SOMA var schema")
    verbose_name_plural = _("SOMA var schemas")

    def get_fields(self, request: HttpRequest, obj=None):
        if obj is None:
            return ("allow_subsets", "csv_file")
        return super().get_fields(request, obj)


class SomaObsColumnSchemaInline(TabularInline):
    model = models.SomaObsColumnSchema
    extra = 1
    fields = ("name", "dtype", "nullable")
    can_delete = True
    verbose_name = _("SOMA obs column")
    verbose_name_plural = _("SOMA obs columns")


@admin.register(models.IngestSchema)
class IngestSchemaAdmin(ModelAdmin):
    list_display = ("name", "measurement_name", "x_validation_type", "obs_column_count", "updated_at")
    search_fields = ("name",)
    list_filter = ("x_validation_type",)
    readonly_fields = ("created_at", "updated_at")
    inlines = [SomaVarSchemaInline, SomaObsColumnSchemaInline]
    fieldsets = (
        (None, {"fields": ("name", "description", "measurement_name", "x_validation_type")}),
        (_("Timestamps"), {"fields": ("created_at", "updated_at"), "classes": ("collapse",)}),
    )

    def obs_column_count(self, obj: models.IngestSchema) -> int:
        return obj.obs_columns.count()

    obs_column_count.short_description = _("Obs columns")

    def save_formset(self, request: HttpRequest, form, formset, change) -> None:
        instances = formset.save(commit=False)
        parent = formset.instance

        for inline_form in formset.forms:
            instance = inline_form.instance

            if not isinstance(instance, models.SomaVarSchema):
                continue

            instance.ingest_schema = parent

            parsed_df = getattr(inline_form, "_parsed_df", None)
            if parsed_df is not None:
                instance.feature_count = len(parsed_df.index)
                instance.var_columns = list(parsed_df.columns)
                parquet_bytes = soma_schema_utils.dataframe_to_parquet_bytes(df=parsed_df)
                slug = slugify(parent.name) or "soma-var-schema"
                timestamp = now().strftime("%Y%m%d%H%M%S")
                filename = f"{slug}-{timestamp}.parquet"
                instance.var_parquet_file.save(filename, ContentFile(parquet_bytes), save=False)

        for obj in instances:
            obj.save()

        for obj in formset.deleted_objects:
            obj.delete()

        formset.save_m2m()

    def get_fieldsets(self, request: HttpRequest, obj=None):
        if obj is None:
            return ((None, {"fields": ("name", "description", "measurement_name", "x_validation_type")}),)
        return super().get_fieldsets(request, obj)


@admin.register(models.SomaVarSchema)
class SomaVarSchemaAdmin(ModelAdmin):
    list_display = ("ingest_schema", "allow_subsets", "feature_count", "updated_at")
    search_fields = ("ingest_schema__name",)
    list_filter = ("allow_subsets",)
    readonly_fields = ("var_parquet_file", "feature_count", "var_columns", "created_at", "updated_at")
    fields = (
        "ingest_schema",
        "allow_subsets",
        "var_parquet_file",
        "feature_count",
        "var_columns",
        "created_at",
        "updated_at",
    )

    def has_add_permission(self, request: HttpRequest) -> bool:
        return False
