#!/bin/bash

set -e

# Color codes for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to check command status
check_command() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Success${NC}"
    else
        echo -e "${RED}✗ Failed${NC}"
        exit 1
    fi
}

# Function to prompt for a required value
prompt_value() {
    local prompt_text=$1
    local var_name=$2
    local default_value=${3:-}
    
    if [ -n "$default_value" ]; then
        echo -e "${YELLOW}$prompt_text${NC} [default: $default_value]"
        echo -e "Leave blank for default"
    else
        echo -e "${YELLOW}$prompt_text${NC}"
    fi
    
    read -r value
    
    if [ -z "$value" ] && [ -n "$default_value" ]; then
        value=$default_value
        echo -e "${GREEN}Using default: $value${NC}"
    fi
    
    while [ -z "$value" ]; do
        echo -e "${RED}This value is required. Please enter a value:${NC}"
        read -r value
    done
    
    eval "$var_name=\"$value\""
}

echo -e "${GREEN}Welcome to Cellarium Nexus Cloud Run Deployment Script${NC}\n"

# Load environment variables from setup
ENV_FILE="./deploy.env"
if [ ! -f "$ENV_FILE" ]; then
    echo -e "${RED}Error: deploy.env file not found. Please run setup.sh first.${NC}"
    exit 1
fi

source "$ENV_FILE"

# Cloud Run settings
prompt_value "Enter Cloud Run service name" "SERVICE_NAME" "cellarium-nexus"
prompt_value "Enter minimum number of instances" "MIN_INSTANCES" "1"
prompt_value "Enter maximum number of instances" "MAX_INSTANCES" "10"
prompt_value "Enter CPU limit" "CPU" "2"
prompt_value "Enter memory limit (in GB)" "MEMORY" "4"
prompt_value "Enter request timeout (in seconds)" "TIMEOUT" "3600"

echo -e "\n${YELLOW}Deploying to Cloud Run...${NC}"

# Deploy to Cloud Run
gcloud run deploy "${SERVICE_NAME}" \
    --image="${IMAGE_PATH}" \
    --platform=managed \
    --region="${REPO_LOCATION}" \
    --project="${PROJECT_ID}" \
    --min-instances="${MIN_INSTANCES}" \
    --max-instances="${MAX_INSTANCES}" \
    --cpu="${CPU}" \
    --memory="${MEMORY}Gi" \
    --timeout="${TIMEOUT}s" \
    --port=8000 \
    --service-account="${BACKEND_SA_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --set-secrets="ENV_FILE=${ENV_SECRET_NAME}:latest" \
    --allow-unauthenticated
check_command

# Get the service URL
SERVICE_URL=$(gcloud run services describe "${SERVICE_NAME}" \
    --platform=managed \
    --region="${REPO_LOCATION}" \
    --project="${PROJECT_ID}" \
    --format='value(status.url)')

# Update local env file with SITE_URL
echo "export SITE_URL=\"${SERVICE_URL}\"" >> "$ENV_FILE"

# Create temporary env file for secret update
tmp_secret="/tmp/env_secret_$$"
while IFS= read -r line; do
    # Remove 'export' and keep the rest
    clean_line="${line#export }"
    echo "${clean_line}" >> "$tmp_secret"
done < "$ENV_FILE"

# Add SITE_URL to secret manager
echo "SITE_URL=${SERVICE_URL}" >> "$tmp_secret"

# Update secret with new environment variables
echo -e "\n${YELLOW}Updating environment secret with SITE_URL...${NC}"
gcloud secrets versions add "${ENV_SECRET_NAME}" \
    --data-file="${tmp_secret}" \
    --project="${PROJECT_ID}"
check_command

# Clean up temporary file
rm "$tmp_secret"

echo -e "\n${GREEN}Deployment completed successfully!${NC}"
echo -e "${YELLOW}Service URL:${NC} ${SERVICE_URL}"
echo -e "\n${YELLOW}Next steps:${NC}"
echo -e "1. Access your application at: ${SERVICE_URL}"
echo -e "2. Monitor the service in Cloud Console"
echo -e "3. Set up custom domain if needed"
echo -e "\nLocal environment file and Secret Manager have been updated with SITE_URL"
