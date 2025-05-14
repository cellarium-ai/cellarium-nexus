"""
Custom filters for Django admin.

This module provides custom filter implementations for Django admin interface,
specifically designed to work with Django Unfold.
"""

import logging
from typing import Any

from django import forms
from django.contrib import admin
from django.db.models import Field, QuerySet
from django.http import HttpRequest
from django.http.request import HttpRequest as WSGIRequest
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _
from unfold.widgets import UnfoldAdminTextInputWidget

from cellarium.nexus.backend.core.admin.filters import GenericDropdownFilter, GenericMultiDropdownFilter

logger = logging.getLogger(__name__)


class SearchWithExcludeForm(forms.Form):
    """
    A form that includes a text input field and an exclude checkbox.

    This form extends the standard Django form to include both a text input
    for search terms and a checkbox to indicate whether to exclude the matches.
    """

    class Media:
        css = {"all": ("cell_management/css/custom_filters.css",)}

    def __init__(self, name: str, label: str, *args, **kwargs) -> None:
        """
        Initialize the form with a text field and exclude checkbox.

        :param name: The name of the field
        :param label: The label for the field
        """
        super().__init__(*args, **kwargs)

        # Add the main search field
        self.fields[name] = forms.CharField(
            label=label,
            required=False,
            widget=UnfoldAdminTextInputWidget,
        )

        # Add the exclude checkbox with a better label that includes a tilde (~) icon
        self.fields[f"{name}_exclude"] = forms.BooleanField(
            label=mark_safe('<span class="exclude-checkbox-label"><span class="tilde-icon">~</span> Exclude</span>'),
            required=False,
            widget=forms.CheckboxInput(
                attrs={
                    "class": "exclude-checkbox",
                }
            ),
            help_text=_("Exclude these values instead of including them"),
        )


class OntologyTermFilter(admin.FieldListFilter):
    """
    A specialized filter for ontology term IDs with exclusion support.

    This filter is specifically designed for filtering ontology term IDs,
    with appropriate title and parameter formatting. It also supports
    excluding the specified values instead of including them.
    """

    title = _("Ontology Term ID")
    form_class = SearchWithExcludeForm
    template = "cell_management/filters/ontology_filter.html"

    # Define lookup parameters at class level to ensure they exist before expected_parameters is called
    lookup_kwarg = ""
    lookup_kwarg_exclude = ""

    def __init__(
        self,
        field: Field,
        request: WSGIRequest,
        params: dict,
        model: Any,
        model_admin: admin.ModelAdmin,
        field_path: str,
    ) -> None:
        """
        Initialize the filter.

        :param field: The field object
        :param request: The HTTP request
        :param params: The request parameters
        :param model: The model class
        :param model_admin: The model admin instance
        :param field_path: The field path
        """
        # Define the lookup parameters
        self.lookup_kwarg = field_path
        self.lookup_kwarg_exclude = f"{field_path}_exclude"

        # Get the values directly from the request GET parameters
        self.lookup_val = request.GET.get(self.lookup_kwarg, "")
        self.lookup_val_exclude = request.GET.get(self.lookup_kwarg_exclude, "")

        # Call parent init
        super().__init__(field, request, params, model, model_admin, field_path)

        # Set a more descriptive title based on the field path
        field_name = field_path.replace("_ontology_term_id", "")
        self.title = _(f"{field_name.replace('_', ' ').title()} Ontology Term ID")

    def has_output(self) -> bool:
        """
        Return whether this filter is active or has an output.

        :return: True if the filter has an output
        """
        return True

    def expected_parameters(self) -> list[str]:
        """
        Return the list of expected parameters for this filter.

        :return: List of parameter names
        """
        return [self.lookup_kwarg, self.lookup_kwarg_exclude]

    def choices(self, changelist):
        """
        Return the choices for the filter form.

        :param changelist: The changelist instance
        :return: A tuple containing a dictionary with the form
        """
        exclude_active = self.value_exclude()

        return (
            {
                "form": self.form_class(
                    label=_("By %(filter_title)s") % {"filter_title": self.title},
                    name=self.lookup_kwarg,  # Use the exact field path without __icontains
                    data={
                        self.lookup_kwarg: self.value(),  # Use the exact field path without __icontains
                        f"{self.lookup_kwarg}_exclude": "on" if exclude_active else "",
                    },
                ),
                "attrs": {"class": "filter-exclude-active" if exclude_active else ""},
            },
        )

    def value_exclude(self) -> bool:
        """
        Get the exclude value from the request parameters.

        :return: True if exclusion is enabled, False otherwise
        """
        return self.lookup_val_exclude == "on"

    def value(self) -> str:
        """
        Get the filter value from the request parameters.

        :return: The filter value as a string
        """
        return self.lookup_val or ""

    def queryset(self, request: HttpRequest, queryset: QuerySet) -> QuerySet:
        """
        Filter the queryset based on the filter value and exclusion setting.

        :param request: The HTTP request
        :param queryset: The queryset to filter

        :return: The filtered queryset
        """
        value = self.value()
        exclude = self.value_exclude()

        if not value:
            return queryset

        # Check if the value contains commas (indicating multiple values)
        if "," in value:
            # Split by comma and clean each value
            values = [v.strip() for v in value.split(",") if v.strip()]
            if not values:
                return queryset

            # Use __in lookup for multiple values
            filter_kwargs = {f"{self.field_path}__in": values}

            # Apply filter or exclude based on the exclude checkbox
            if exclude:
                filtered_qs = queryset.exclude(**filter_kwargs)
            else:
                filtered_qs = queryset.filter(**filter_kwargs)

            return filtered_qs

        # For single values
        filter_kwargs = {f"{self.field_path}__exact": value}

        # Apply filter or exclude based on the exclude checkbox
        if exclude:
            filtered_qs = queryset.exclude(**filter_kwargs)
        else:
            filtered_qs = queryset.filter(**filter_kwargs)

        return filtered_qs


class SexDropdownFilter(GenericDropdownFilter):
    title = _("sex")
    parameter_name = "sex"
    field_name = "sex"


class DonorDropdownFilter(GenericMultiDropdownFilter):
    title = _("donor id")
    parameter_name = "donor_id"
    field_name = "donor_id"


class SuspensionTypeDropdownFilter(GenericMultiDropdownFilter):
    title = _("suspension type")
    parameter_name = "suspension_type"
    field_name = "suspension_type"


class OrganismDropdownFilter(GenericMultiDropdownFilter):
    title = _("organism")
    parameter_name = "organism"
    field_name = "organism"


class TagDropdownFilter(GenericMultiDropdownFilter):
    title = _("tag")
    parameter_name = "tag"
    field_name = "tag"


class CellTypeDropdownFilter(GenericMultiDropdownFilter):
    title = _("cell type")
    parameter_name = "cell_type"
    field_name = "cell_type"


class DiseaseDropdownFilter(GenericMultiDropdownFilter):
    title = _("disease")
    parameter_name = "disease"
    field_name = "disease"
