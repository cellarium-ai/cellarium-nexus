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

# Function to print script usage
print_usage() {
    echo "Usage: $0 --image-path IMAGE_PATH --project-id PROJECT_ID --sql-connection-name SQL_CONNECTION_NAME [--repo-location REPO_LOCATION] [--backend-sa-name BACKEND_SA_NAME] [--env-secret-name ENV_SECRET_NAME] [--service-name SERVICE_NAME] [--min-instances MIN_INSTANCES] [--max-instances MAX_INSTANCES] [--cpu CPU] [--memory MEMORY] [--timeout TIMEOUT]"
}

# Default values
REPO_LOCATION="us-central1"
BACKEND_SA_NAME="cellarium-backend-sa"
ENV_SECRET_NAME="cellarium-nexus"
REMOTE_ENV_FILE="/app/conf/.env"
SERVICE_NAME="cellarium-nexus"
MIN_INSTANCES=1
MAX_INSTANCES=10
CPU=2
MEMORY=4
TIMEOUT=3600

# Parse input parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --image-path=*) IMAGE_PATH="${1#*=}" ;;
        --image-path) IMAGE_PATH="$2"; shift ;;
        --project-id=*) PROJECT_ID="${1#*=}" ;;
        --project-id) PROJECT_ID="$2"; shift ;;
        --repo-location=*) REPO_LOCATION="${1#*=}" ;;
        --repo-location) REPO_LOCATION="$2"; shift ;;
        --backend-sa-name=*) BACKEND_SA_NAME="${1#*=}" ;;
        --backend-sa-name) BACKEND_SA_NAME="$2"; shift ;;
        --env-secret-name=*) ENV_SECRET_NAME="${1#*=}" ;;
        --env-secret-name) ENV_SECRET_NAME="$2"; shift ;;
        --sql-connection-name=*) SQL_CONNECTION_NAME="${1#*=}" ;;
        --sql-connection-name) SQL_CONNECTION_NAME="$2"; shift ;;
        --service-name=*) SERVICE_NAME="${1#*=}" ;;
        --service-name) SERVICE_NAME="$2"; shift ;;
        --min-instances=*) MIN_INSTANCES="${1#*=}" ;;
        --min-instances) MIN_INSTANCES="$2"; shift ;;
        --max-instances=*) MAX_INSTANCES="${1#*=}" ;;
        --max-instances) MAX_INSTANCES="$2"; shift ;;
        --cpu=*) CPU="${1#*=}" ;;
        --cpu) CPU="$2"; shift ;;
        --memory=*) MEMORY="${1#*=}" ;;
        --memory) MEMORY="$2"; shift ;;
        --timeout=*) TIMEOUT="${1#*=}" ;;
        --timeout) TIMEOUT="$2"; shift ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo -e "${RED}Unknown parameter: $1${NC}"; print_usage; exit 1 ;;
    esac
    shift
done

# Check required parameters
if [ -z "$IMAGE_PATH" ] || [ -z "$PROJECT_ID" ] || [ -z "$SQL_CONNECTION_NAME" ]; then
    echo -e "${RED}Missing required parameters!${NC}"
    print_usage
    exit 1
fi

echo -e "${GREEN}Welcome to Cellarium Nexus Cloud Run Deployment Script${NC}\n"
echo -e "${YELLOW}This script will deploy Cellarium Nexus to Google Cloud Run.${NC}"
echo -e "${YELLOW}You will be prompted for necessary configuration values.${NC}\n"
echo -e "${YELLOW}Creating Cloud Run job for Django deployment preparation...${NC}"

# Define job name
JOB_NAME="django-prepare-deploy"

# Check if the Cloud Run job already exists
if gcloud run jobs describe "${JOB_NAME}" --region="${REPO_LOCATION}" --project="${PROJECT_ID}" &>/dev/null; then
    echo -e "${YELLOW}Updating existing Cloud Run job...${NC}"
    # Update the existing job
    gcloud run jobs update "${JOB_NAME}" \
        --image="${IMAGE_PATH}" \
        --region="${REPO_LOCATION}" \
        --project="${PROJECT_ID}" \
        --service-account="${BACKEND_SA_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
        --set-cloudsql-instances="${SQL_CONNECTION_NAME}" \
        --set-secrets="${REMOTE_ENV_FILE}=${ENV_SECRET_NAME}:latest" \
        --command="/bin/bash" \
        --args="-c,/app/deploy/backend/django_prepare_deploy.sh --env-file ${REMOTE_ENV_FILE}"
    check_command
else
    echo -e "${YELLOW}Creating new Cloud Run job...${NC}"
    # Create a new job
    gcloud run jobs create "${JOB_NAME}" \
        --image="${IMAGE_PATH}" \
        --region="${REPO_LOCATION}" \
        --project="${PROJECT_ID}" \
        --service-account="${BACKEND_SA_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
        --set-cloudsql-instances="${SQL_CONNECTION_NAME}" \
        --set-secrets="${REMOTE_ENV_FILE}=${ENV_SECRET_NAME}:latest" \
        --command="/bin/bash" \
        --args="-c,/app/deploy/backend/django_prepare_deploy.sh --env-file ${REMOTE_ENV_FILE}"
    check_command
fi

echo -e "\n${YELLOW}Executing Django deployment preparation job...${NC}"
gcloud run jobs execute "${JOB_NAME}" \
    --region="${REPO_LOCATION}" \
    --project="${PROJECT_ID}" \
    --wait
check_command

echo -e "\n${GREEN}Django deployment preparation completed successfully!${NC}"

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
    --port=8080 \
    --service-account="${BACKEND_SA_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --add-cloudsql-instances="${SQL_CONNECTION_NAME}" \
    --set-secrets="${REMOTE_ENV_FILE}=${ENV_SECRET_NAME}:latest" \
    --allow-unauthenticated
check_command

# Get the service URL
SERVICE_URL=$(gcloud run services describe "${SERVICE_NAME}" \
    --platform=managed \
    --region="${REPO_LOCATION}" \
    --project="${PROJECT_ID}" \
    --format='value(status.url)')

echo -e "\n${GREEN}Deployment completed successfully!${NC}"
echo -e "${YELLOW}Service URL:${NC} ${SERVICE_URL}"
echo -e "\n${YELLOW}Next steps:${NC}"
echo -e "1. Access your application at: ${SERVICE_URL}"
echo -e "2. Monitor the service in Cloud Console"
echo -e "3. Set up custom domain if needed"
echo -e "\nLocal env & secret updates are handled externally"
