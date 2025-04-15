#!/bin/bash
set -e

# Initialize the Django application
echo "Running database migrations..."
python /app/cellarium/nexus/backend/manage.py migrate --noinput

echo "Collecting static files..."
python /app/cellarium/nexus/backend/manage.py collectstatic --noinput

# Function to check if a superuser with specific username exists
check_superuser_exists() {
    local username=$1
    echo "Checking if superuser with username '$username' exists..."
    python -c "
import django
django.setup()
from django.contrib.auth import get_user_model
User = get_user_model()
user_exists = User.objects.filter(username='$username', is_superuser=True).exists()
print(f'Superuser with username $username exists: {user_exists}')
exit(0 if user_exists else 1)
"
}

# Path to the .env file
ENV_FILE="/app/conf/.env"

# Check if .env file exists and contains required variables
if [ -f "$ENV_FILE" ] && \
   grep -q "DJANGO_SUPERUSER_USERNAME=" "$ENV_FILE" && \
   grep -q "DJANGO_SUPERUSER_PASSWORD=" "$ENV_FILE" && \
   grep -q "DJANGO_SUPERUSER_EMAIL=" "$ENV_FILE"; then
    
    # Export variables from .env file
    export $(grep -E "^DJANGO_SUPERUSER_(USERNAME|PASSWORD|EMAIL)=" "$ENV_FILE" | xargs)
    
    # Check if superuser with this username already exists
    if ! check_superuser_exists "$DJANGO_SUPERUSER_USERNAME"; then
        echo "Creating superuser with username: $DJANGO_SUPERUSER_USERNAME"
        python /app/cellarium/nexus/backend/manage.py createsuperuser --noinput
        echo "Superuser created successfully."
    else
        echo "Superuser with username '$DJANGO_SUPERUSER_USERNAME' already exists. Skipping creation."
    fi
else
    echo "Warning: Cannot create superuser. Required variables not found in .env file or file not found at $ENV_FILE."
fi

# Execute the command passed to the script
echo "Starting server..."
exec "$@"
