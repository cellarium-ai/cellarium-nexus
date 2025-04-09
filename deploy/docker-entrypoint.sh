#!/bin/bash
set -e

# Initialize the Django application
echo "Running database migrations..."
python /app/cellarium/nexus/backend/manage.py migrate --noinput

echo "Collecting static files..."
python /app/cellarium/nexus/backend/manage.py collectstatic --noinput

# Function to check if there are any superusers
check_superuser_exists() {
    echo "Checking if superuser exists..."
    python -c "
import django
django.setup()
from django.contrib.auth import get_user_model
User = get_user_model()
superusers = User.objects.filter(is_superuser=True)
exists = superusers.exists()
print(f'Superuser exists: {exists}')
exit(0 if exists else 1)
"
}

# Create superuser if none exists
if ! check_superuser_exists; then
    echo "No superuser found. Creating superuser..."
    
    # Path to the .env file
    ENV_FILE="/app/conf/.env"
    
    # Check if .env file exists
    if [ -f "$ENV_FILE" ]; then
        # Read superuser credentials from .env file
        if grep -q "DJANGO_SUPERUSER_USERNAME=" "$ENV_FILE" && \
           grep -q "DJANGO_SUPERUSER_PASSWORD=" "$ENV_FILE" && \
           grep -q "DJANGO_SUPERUSER_EMAIL=" "$ENV_FILE"; then
            
            # Export variables from .env file
            export $(grep -E "^DJANGO_SUPERUSER_(USERNAME|PASSWORD|EMAIL)=" "$ENV_FILE" | xargs)
            
            # Create superuser
            python /app/cellarium/nexus/backend/manage.py createsuperuser --noinput
            echo "Superuser created successfully from .env file."
        else
            echo "Warning: Cannot create superuser. Required variables not found in .env file."
        fi
    else
        echo "Warning: .env file not found at $ENV_FILE. Cannot create superuser."
    fi
else
    echo "Superuser already exists. Skipping creation."
fi

# Execute the command passed to the script
echo "Starting server..."
exec "$@"
