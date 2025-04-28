from django.core.management.base import BaseCommand
from django.urls import get_resolver


class Command(BaseCommand):
    """
    List all registered URLs in the Django project.

    This command prints all URL patterns registered in the project,
    including admin URLs, with their names and view functions.

    :raise: None

    :return: None
    """

    help = "List all registered URLs in the Django project"

    def handle(self, *args, **options):
        """
        Execute the command to list all URLs.

        :param args: Command arguments
        :param options: Command options

        :raise: None

        :return: None
        """
        resolver = get_resolver()
        self._print_url_patterns(resolver.url_patterns)

    def _print_url_patterns(self, url_patterns, prefix=""):
        """
        Recursively print URL patterns.

        :param url_patterns: List of URL patterns
        :param prefix: URL prefix for nested patterns

        :raise: None

        :return: None
        """
        for pattern in url_patterns:
            if hasattr(pattern, "url_patterns"):
                # This is an include - recurse into it
                self._print_url_patterns(pattern.url_patterns, prefix + str(pattern.pattern))
            else:
                # This is a URL pattern
                name = getattr(pattern, "name", None)
                view_name = self._get_view_name(pattern)
                pattern_str = prefix + str(pattern.pattern)
                self.stdout.write(f"{pattern_str} - Name: {name} - View: {view_name}")

    def _get_view_name(self, pattern):
        """
        Get the view name from a URL pattern.

        :param pattern: URL pattern

        :raise: None

        :return: String representation of the view
        """
        if hasattr(pattern, "callback"):
            if hasattr(pattern.callback, "__name__"):
                return f"{pattern.callback.__module__}.{pattern.callback.__name__}"
            return str(pattern.callback)
        return "No view"
