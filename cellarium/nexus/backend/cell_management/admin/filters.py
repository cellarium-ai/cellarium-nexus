"""
Custom filters for Django admin.

This module provides custom filter implementations for Django admin interface,
specifically designed to work with Django Unfold.
"""

import logging
from typing import Any, List

from django import forms
from django.contrib import admin
from django.db.models import Field, QuerySet
from django.http import HttpRequest
from django.http.request import HttpRequest as WSGIRequest
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _
from unfold.contrib.filters.admin import FieldTextFilter
from unfold.widgets import UnfoldAdminTextInputWidget

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


class MultiValueTextFilter(FieldTextFilter):
    """
    A custom filter that supports comma-separated values for __in lookups.

    This filter extends Django Unfold's FieldTextFilter to support filtering
    by multiple values using comma-separated input. It automatically converts
    the input to an __in lookup when commas are detected.
    """

    def queryset(self, request: HttpRequest, queryset: QuerySet) -> QuerySet:
        """
        Filter the queryset based on the filter value.

        If the value contains commas, it will be split and used for an __in lookup.
        Otherwise, it will use the standard exact lookup.

        :param request: The HTTP request
        :param queryset: The queryset to filter

        :return: The filtered queryset
        """
        value = self.value()
        if not value:
            return queryset

        # Check if the value contains commas (indicating multiple values)
        if "," in value:
            # Split by comma and clean each value
            values = [v.strip() for v in value.split(",") if v.strip()]
            if values:
                # Use __in lookup for multiple values
                filter_kwargs = {f"{self.field_path}__in": values}

                # Debug logging to print the filter parameters and SQL query
                logger.info(f"MultiValueTextFilter: field_path={self.field_path}, values={values}")

                filtered_qs = queryset.filter(**filter_kwargs)

                # Log the SQL query (this will only work if DEBUG=True in settings)
                # Add a note about proper quoting in the actual query
                logger.info(f"Generated SQL: {filtered_qs.query}")
                logger.info(
                    "Note: The actual SQL will have proper quoting for string values, even if not shown in the log"
                )

                return filtered_qs
            return queryset

        # For single values, use __exact lookup explicitly
        # This ensures we match exactly what the user entered and handles quoting properly
        filter_kwargs = {f"{self.field_path}__exact": value}

        # Apply the filter
        filtered_qs = queryset.filter(**filter_kwargs)

        # Debug logging for single value filter
        logger.info(f"Exact filter: field_path={self.field_path}, value={value}")
        logger.info(f"Generated SQL: {filtered_qs.query}")
        logger.info("Note: The actual SQL will have proper quoting for string values, even if not shown in the log")

        # Ensure we never return None
        return filtered_qs


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

        # Log the raw request parameters for debugging
        logger.info(f"OntologyTermFilter init: field_path={field_path}")
        logger.info(f"OntologyTermFilter raw params: {request.GET}")

        # Get the values directly from the request GET parameters
        self.lookup_val = request.GET.get(self.lookup_kwarg, "")
        self.lookup_val_exclude = request.GET.get(self.lookup_kwarg_exclude, "")

        # Log the lookup values for debugging
        logger.info(
            f"OntologyTermFilter lookup values from request: lookup_kwarg={self.lookup_kwarg}, lookup_val={self.lookup_val}"
        )

        # Call parent init
        super().__init__(field, request, params, model, model_admin, field_path)

        # Set a more descriptive title based on the field path
        field_name = field_path.replace("_ontology_term_id", "")
        self.title = _(f"{field_name.replace('_', ' ').title()} Ontology Term ID")

        # Log the final lookup values for debugging
        logger.info(f"OntologyTermFilter final values: lookup_kwarg={self.lookup_kwarg}, lookup_val={self.lookup_val}")

    def has_output(self) -> bool:
        """
        Return whether this filter is active or has an output.

        :return: True if the filter has an output
        """
        return True

    def expected_parameters(self) -> List[str]:
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

        # Log the form data for debugging
        logger.info(
            f"OntologyTermFilter choices: lookup_kwarg={self.lookup_kwarg}, value={self.value()}, exclude={exclude_active}"
        )

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

        # Log the actual filter parameters being used
        logger.info(f"OntologyTermFilter queryset: field_path={self.field_path}, value={value}, exclude={exclude}")

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

            # Debug logging
            logger.info(f"OntologyTermFilter: field_path={self.field_path}, values={values}, exclude={exclude}")

            # Apply filter or exclude based on the exclude checkbox
            if exclude:
                filtered_qs = queryset.exclude(**filter_kwargs)
                logger.info(f"Excluding with: {filter_kwargs}")
            else:
                filtered_qs = queryset.filter(**filter_kwargs)
                logger.info(f"Filtering with: {filter_kwargs}")

            return filtered_qs

        # For single values
        filter_kwargs = {f"{self.field_path}__exact": value}

        # Apply filter or exclude based on the exclude checkbox
        if exclude:
            filtered_qs = queryset.exclude(**filter_kwargs)
            logger.info(f"Excluding with: {filter_kwargs}")
        else:
            filtered_qs = queryset.filter(**filter_kwargs)
            logger.info(f"Filtering with: {filter_kwargs}")

        return filtered_qs
