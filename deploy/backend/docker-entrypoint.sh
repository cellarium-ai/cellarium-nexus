#!/bin/bash
set -e

# Start the gunicorn server
echo "Starting gunicorn server..."
exec gunicorn --bind 0.0.0.0:8080 \
     --workers 3 \
     --timeout 1200 \
     --log-level info \
     cellarium.nexus.backend.application.wsgi:application
