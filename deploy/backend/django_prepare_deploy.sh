#!/bin/bash

set -e

# Color codes (enable only if output is a TTY)
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    NC='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; NC=''
fi

# Function to check command status
check_command() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Success${NC}"
    else
        echo -e "${RED}✗ Failed${NC}"
        exit 1
    fi
}

# Parse command line arguments
ENV_FILE=".env"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --env-file) ENV_FILE="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

## Use absolute path for manage.py instead of changing directory

echo -e "${GREEN}Running Django Deployment Preparation: Migrations, Static Files, and Superuser Creation${NC}\n"

# Run database migrations
echo -e "\n${YELLOW}Running database migrations...${NC}"
python /app/cellarium/nexus/backend/manage.py migrate --noinput
check_command
echo -e "${GREEN}Database migrations completed successfully!${NC}"

# Create cache table
echo -e "\n${YELLOW}Creating cache table...${NC}"
python /app/cellarium/nexus/backend/manage.py createcachetable
check_command
echo -e "${GREEN}Cache table created successfully!${NC}"

# Collect static files
echo -e "\n${YELLOW}Collecting static files...${NC}"
python /app/cellarium/nexus/backend/manage.py collectstatic --noinput
check_command
echo -e "${GREEN}Static files collected successfully!${NC}"

# Create superuser if it doesn't exist
echo -e "\n${YELLOW}Checking for superuser...${NC}"
if [ -f "$ENV_FILE" ] && grep -q "DJANGO_SUPERUSER_USERNAME=" "$ENV_FILE" && grep -q "DJANGO_SUPERUSER_PASSWORD=" "$ENV_FILE" && grep -q "DJANGO_SUPERUSER_EMAIL=" "$ENV_FILE"; then
    export $(grep -E "^DJANGO_SUPERUSER_(USERNAME|PASSWORD|EMAIL)=" "$ENV_FILE" | xargs)

    # Check if superuser exists
    python -c "
import django
django.setup()
from django.contrib.auth import get_user_model
User = get_user_model()
exists = User.objects.filter(username=\"$DJANGO_SUPERUSER_USERNAME\", is_superuser=True).exists()
print(f\"Superuser exists: {exists}\")
exit(0 if exists else 1)
" || python /app/cellarium/nexus/backend/manage.py createsuperuser --noinput

    check_command
    echo -e "${GREEN}Superuser check/creation completed successfully!${NC}"
else
    echo -e "${YELLOW}Warning: Cannot create superuser. Required variables not found in .env file or file not found.${NC}"
fi

echo -e "\n${GREEN}Django deployment preparation completed successfully!${NC}"
