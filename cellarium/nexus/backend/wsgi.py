"""
WSGI config for backend project.

It exposes the WSGI callable as a module-level variable named ``backend``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application

os.environ.setdefault(key='DJANGO_SETTINGS_MODULE', value='nexus.backend.settings')

application = get_wsgi_application()
