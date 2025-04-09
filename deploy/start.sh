#!/bin/bash
set -e

echo "Verifying environment..."
python -c "
import sys
print(f'Python version: {sys.version}')
print(f'Python path: {sys.path}')
print('Attempting imports...')
import cellarium
import cellarium.nexus
import cellarium.nexus.backend
import cellarium.nexus.backend.wsgi
print('All imports successful')
"

echo "Starting Gunicorn..."
exec gunicorn \
    --bind 0.0.0.0:8000 \
    --workers 3 \
    --timeout 120 \
    --log-level debug \
    --error-logfile - \
    --capture-output \
    --preload \
    cellarium.nexus.backend.wsgi:application
