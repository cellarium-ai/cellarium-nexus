#!/bin/bash

set -e

# Color codes for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Welcome to Cellarium Nexus Setup Script${NC}\n"

# Function to generate a random secret
generate_secret() {
    openssl rand -hex 32
}

# Function to prompt for a secret value with auto-generation option
prompt_secret() {
    local prompt_text=$1
    local var_name=$2
    read -rp "$(echo -e "${YELLOW}${prompt_text} (leave blank for random):${NC} ")" value
    if [ -z "$value" ]; then
        value=$(generate_secret)
        echo -e "${GREEN}Generated secret: $value${NC}"
    fi
    echo ""  # single blank line separation
    eval "$var_name=\"$value\""
}

# Function to prompt for a required value
prompt_value() {
    local prompt_text=$1
    local var_name=$2
    local default_value=${3:-}

    if [ -n "$default_value" ]; then
        read -rp "$(echo -e "${YELLOW}${prompt_text} [default: ${default_value}]:${NC} ")" value
        if [ -z "$value" ]; then
            value=$default_value
            echo -e "${GREEN}Using default: $value${NC}"
        fi
    else
        read -rp "$(echo -e "${YELLOW}${prompt_text}:${NC} ")" value
    fi

    while [ -z "$value" ]; do
        echo -e "${RED}This value is required. Please enter a value:${NC}"
        read -rp "$(echo -e "${YELLOW}${prompt_text}:${NC} ")" value
    done

    echo ""  # single blank line separation
    eval "$var_name=\"$value\""
}

# Docker image
prompt_value "Enter the Docker image URL for Cellarium Nexus" "IMAGE_PATH"

# Nexus pipeline base Docker image URL
prompt_value "Enter the Nexus pipeline base Docker image URL" "PIPELINE_BASE_IMAGE"

# Project settings
prompt_value "Enter your Google Cloud Project ID" "GCP_PROJECT_ID"

# Application billing label
prompt_value "Enter application billing label" "GCP_APPLICATION_BILLING_LABEL" "cellarium-nexus"

# Database instance settings
prompt_value "Enter database instance name" "DB_INSTANCE_NAME" "cellarium-nexus-instance"
prompt_value "Enter number of CPUs" "DB_CPU" "2"
prompt_value "Enter memory in GB" "DB_MEMORY" "8"
prompt_value "Enter storage size in GB" "DB_STORAGE" "10"

# Database settings
prompt_value "Enter the region for the database instance" "DB_REGION" "us-central1"
DB_INSTANCE_CONNECTION_NAME="${GCP_PROJECT_ID}:${DB_REGION}:${DB_INSTANCE_NAME}"

echo "Constructed DB Instance Connection Name: ${DB_INSTANCE_CONNECTION_NAME}"
echo 
prompt_value "Enter database name" "DB_NAME" "cellarium-nexus-db"
prompt_value "Enter database user" "DB_USER" "cellarium-nexus-db-user"
prompt_secret "Enter database password" "DB_PASSWORD"

# Secret Manager settings
prompt_value "Enter environment secret name" "ENV_SECRET_NAME" "cellarium-nexus"

# Service account settings
prompt_value "Enter Vertex AI Pipelines service account name" "PIPELINE_SA_NAME" "vertex-pipelines-sa"
prompt_value "Enter Vertex AI Pipelines service account display name" "PIPELINE_SA_DISPLAY_NAME" "Vertex AI Pipelines Service Account"
prompt_value "Enter Backend service account name" "BACKEND_SA_NAME" "cellarium-backend-sa"
prompt_value "Enter Backend service account display name" "BACKEND_SA_DISPLAY_NAME" "Cellarium Backend Service Account"

# GCS bucket settings
# Generate unique 6-character suffixes for buckets
PIPELINE_SUFFIX=$(uuidgen | tr -d '-' | tr '[:upper:]' '[:lower:]' | head -c6)
PRIVATE_SUFFIX=$(uuidgen | tr -d '-' | tr '[:upper:]' '[:lower:]' | head -c6)
PUBLIC_SUFFIX=$(uuidgen | tr -d '-' | tr '[:upper:]' '[:lower:]' | head -c6)

prompt_value "Enter pipeline bucket name" "PIPELINE_BUCKET" "${GCP_PROJECT_ID}-pipeline-root-${PIPELINE_SUFFIX}"
prompt_value "Enter private bucket name" "BUCKET_NAME_PRIVATE" "cellarium-nexus-file-system-${PRIVATE_SUFFIX}"
prompt_value "Enter public bucket name" "BUCKET_NAME_PUBLIC" "cellarium-nexus-public-files-${PUBLIC_SUFFIX}"

# Artifact Registry settings
prompt_value "Enter repository name" "REPO_NAME" "cellarium-images"
prompt_value "Enter repository location" "REPO_LOCATION" "us-central1"

# Cloud Run settings
prompt_value "Enter Cloud Run service name" "SERVICE_NAME" "nexus"
prompt_value "Enter minimum number of instances" "MIN_INSTANCES" "1"
prompt_value "Enter maximum number of instances" "MAX_INSTANCES" "10"
prompt_value "Enter CPU limit" "CPU" "2"
prompt_value "Enter memory limit (in GB)" "MEMORY" "4"
prompt_value "Enter request timeout (in seconds)" "TIMEOUT" "3600"

# Generate Django secret key
prompt_secret "Enter Django secret key" "SECRET_KEY"

# Django superuser settings
prompt_value "Enter Django superuser username" "DJANGO_SUPERUSER_USERNAME" "admin"
prompt_value "Enter Django superuser email" "DJANGO_SUPERUSER_EMAIL" "admin@example.com"

# Prompt for Django superuser password with 8-character random default
echo -e "${YELLOW}Enter Django superuser password (leave blank for random 8-char):${NC}"
read -r DJANGO_SUPERUSER_PASSWORD
if [ -z "$DJANGO_SUPERUSER_PASSWORD" ]; then
    DJANGO_SUPERUSER_PASSWORD=$(openssl rand -hex 4)
    echo -e "${GREEN}Generated superuser password: $DJANGO_SUPERUSER_PASSWORD${NC}"
fi

# Prompt for environment
echo -e "\n${YELLOW}Select environment (default: development):${NC}"
echo "1) development"
echo "2) local"
echo "3) production"
read -p "Enter choice [1-3] (default: 1): " env_choice

case $env_choice in
    2) ENVIRONMENT="local" ;;
    3) ENVIRONMENT="production" ;;
    ""|1) ENVIRONMENT="development" ;;
    *) 
        echo -e "${RED}Invalid choice. Please select 1, 2, or 3${NC}"
        exit 1
        ;;
esac

echo -e "\n${GREEN}Setting environment to: ${ENVIRONMENT}${NC}"

# Set initial URL values based on environment
if [ "$ENVIRONMENT" = "local" ]; then
    SITE_URL="http://localhost:8000"
    MAIN_HOST_ALLOWED="localhost"
    echo -e "${GREEN}Using local development URLs:${NC}"
else
    # Set temporary values for non-local environments
    # These will be updated after Cloud Run deployment
    SITE_URL="http://pending-cloud-run-url"
    MAIN_HOST_ALLOWED="pending-cloud-run-domain"
    echo -e "${YELLOW}Using temporary URLs (will be updated after Cloud Run deployment):${NC}"
fi

echo -e "SITE_URL: ${SITE_URL}"
echo -e "MAIN_HOST_ALLOWED: ${MAIN_HOST_ALLOWED}"

# Create output file
output_file=".env.generated"
echo "# Generated environment variables for Cellarium Nexus" > "$output_file"
echo "# Generated on: $(date)" >> "$output_file"
echo "" >> "$output_file"

# Construct full service account emails and pipeline paths
PIPELINE_SERVICE_ACCOUNT="${PIPELINE_SA_NAME}@${GCP_PROJECT_ID}.iam.gserviceaccount.com"
PIPELINE_ROOT_PATH="gs://${PIPELINE_BUCKET}/pipeline-root"

# Write all variables to file
declare -a vars=(
    "GCP_PROJECT_ID"
    "GCP_APPLICATION_BILLING_LABEL"
    "DB_INSTANCE_CONNECTION_NAME"
    "DB_NAME"
    "DB_USER"
    "DB_PASSWORD"
    "PIPELINE_SERVICE_ACCOUNT"
    "PIPELINE_ROOT_PATH"
    "PIPELINE_BASE_IMAGE"
    "BUCKET_NAME_PRIVATE"
    "BUCKET_NAME_PUBLIC"
    "SECRET_KEY"
    "ENVIRONMENT"
    "DJANGO_SUPERUSER_USERNAME"
    "DJANGO_SUPERUSER_PASSWORD"
    "DJANGO_SUPERUSER_EMAIL"
    "SITE_URL"
    "MAIN_HOST_ALLOWED"
)

# Create a temporary file for secret manager (without exports)
secret_file="/tmp/env_secret_$$"
for var in "${vars[@]}"; do
    # Write to local env file with exports (for shell sourcing)
    echo "export ${var}=\"${!var}\"" >> "$output_file"
    # Write to secret file without exports (for env file format)
    echo "${var}=${!var}" >> "$secret_file"
done

# Source the generated environment variables
source "$output_file"

echo -e "\n${GREEN}Starting deployment process...${NC}"

# Verify project access
echo -e "\n${YELLOW}Verifying Google Cloud project access...${NC}"
if ! gcloud projects describe "${GCP_PROJECT_ID}" --project="${GCP_PROJECT_ID}" > /dev/null 2>&1; then
    echo -e "${RED}Unable to access project ${GCP_PROJECT_ID}. Please ensure:${NC}"
    echo -e "1. The project ID is correct"
    echo -e "2. You have the necessary permissions"
    echo -e "3. The project exists and is active"
    exit 1
fi
echo -e "${GREEN}Project access verified successfully${NC}"

# Function to check command status
check_command() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Success${NC}"
        return 0
    else
        echo -e "${RED}✗ Failed${NC}"
        return 1
    fi
}

# Function to handle operation retry
retry_operation() {
    local operation_name=$1
    local cmd=$2

    while true; do
        echo -e "\n${YELLOW}Attempting: ${operation_name}${NC}"
        if eval "$cmd"; then
            echo -e "${GREEN}Operation successful${NC}"
            return 0
        else
            echo -e "${RED}Operation failed${NC}"
            echo -e "\nChoose an option:"
            echo "1) Skip this operation (resource may already exist)"
            echo "2) Retry the operation"
            read -r choice

            case $choice in
                1)
                    echo -e "${YELLOW}Skipping operation - assuming resource exists${NC}"
                    return 0  # Return success when skipping
                    ;;
                2)
                    echo -e "${YELLOW}Retrying operation...${NC}"
                    continue
                    ;;
                *)
                    echo -e "${RED}Invalid choice. Retrying...${NC}"
                    ;;
            esac
        fi
    done
}

# Function to create service account
create_service_account() {
    local sa_name=$1
    local sa_display_name=$2

    retry_operation "Creating service account ${sa_name}" \
        "gcloud iam service-accounts create \"${sa_name}\" \
        --display-name=\"${sa_display_name}\" \
        --project=\"${GCP_PROJECT_ID}\""
}

echo -e "\n${YELLOW}Enabling required Google Cloud APIs...${NC}"
gcloud services enable \
    secretmanager.googleapis.com \
    aiplatform.googleapis.com \
    storage.googleapis.com \
    artifactregistry.googleapis.com \
    run.googleapis.com \
    cloudbuild.googleapis.com \
    sqladmin.googleapis.com \
    bigquery.googleapis.com \
    --project="${GCP_PROJECT_ID}"
check_command

echo -e "\n${YELLOW}Creating GCS buckets...${NC}"
# Create private bucket
echo "Creating private bucket: ${BUCKET_NAME_PRIVATE}"
retry_operation "Creating private bucket ${BUCKET_NAME_PRIVATE}" \
    "gcloud storage buckets create \"gs://${BUCKET_NAME_PRIVATE}\" \
    --project=\"${GCP_PROJECT_ID}\" \
    --location=us-central1 \
    --uniform-bucket-level-access \
    --labels=application=${GCP_APPLICATION_BILLING_LABEL}"

# Create public bucket
echo "Creating public bucket: ${BUCKET_NAME_PUBLIC}"
retry_operation "Creating public bucket ${BUCKET_NAME_PUBLIC}" \
    "gcloud storage buckets create \"gs://${BUCKET_NAME_PUBLIC}\" \
    --project=\"${GCP_PROJECT_ID}\" \
    --location=us-central1 \
    --uniform-bucket-level-access \
    --labels=application=${GCP_APPLICATION_BILLING_LABEL}"

echo -e "\n${YELLOW}Configuring CORS for public bucket...${NC}"
gcloud storage buckets update "gs://${BUCKET_NAME_PUBLIC}" \
    --cors-file=deploy/setup/cors.json \
    --project="${GCP_PROJECT_ID}"
check_command

echo -e "\n${YELLOW}Configuring public bucket IAM policy...${NC}"
# Generate bucket policy with project-specific values
cat > /tmp/bucket-policy.yaml << EOF
bindings:
- members:
  - projectEditor:${GCP_PROJECT_ID}
  - projectOwner:${GCP_PROJECT_ID}
  role: roles/storage.legacyBucketOwner
- members:
  - projectViewer:${GCP_PROJECT_ID}
  role: roles/storage.legacyBucketReader
- members:
  - projectEditor:${GCP_PROJECT_ID}
  - projectOwner:${GCP_PROJECT_ID}
  role: roles/storage.legacyObjectOwner
- members:
  - allUsers
  - projectViewer:${GCP_PROJECT_ID}
  role: roles/storage.objectViewer
EOF

gcloud storage buckets set-iam-policy "gs://${BUCKET_NAME_PUBLIC}" /tmp/bucket-policy.yaml
check_command
rm /tmp/bucket-policy.yaml

echo -e "\n${YELLOW}Creating service accounts...${NC}"

# Function to grant roles to a service account
grant_roles() {
    local sa_name=$1
    shift
    local roles=($@)

    for role in "${roles[@]}"; do
        retry_operation "Granting ${role} to ${sa_name}" \
            "gcloud projects add-iam-policy-binding \"${GCP_PROJECT_ID}\" \
            --member=\"serviceAccount:${sa_name}@${GCP_PROJECT_ID}.iam.gserviceaccount.com\" \
            --role=\"${role}\""
    done
}

# Create and configure Vertex AI Pipelines service account
echo -e "\n${YELLOW}Setting up Vertex AI Pipelines service account...${NC}"
create_service_account "${PIPELINE_SA_NAME}" "${PIPELINE_SA_DISPLAY_NAME}"

# Create and configure Backend service account
echo -e "\n${YELLOW}Setting up Backend service account...${NC}"
create_service_account "${BACKEND_SA_NAME}" "${BACKEND_SA_DISPLAY_NAME}"

echo -e "\n${YELLOW}Granting IAM roles...${NC}"
# Grant roles to Vertex AI Pipelines service account
declare -a pipeline_roles=(
    "roles/aiplatform.user"
    "roles/bigquery.admin"
    "roles/bigquery.jobUser"
    "roles/storage.admin"
)

# Grant roles to Backend service account
declare -a backend_roles=(
    "roles/aiplatform.user"
    "roles/bigquery.dataEditor"
    "roles/bigquery.jobUser"
    "roles/storage.admin"
    "roles/secretmanager.secretAccessor"
)

# Grant roles (will work with both existing and new service accounts)
grant_roles "${PIPELINE_SA_NAME}" "${pipeline_roles[@]}"
grant_roles "${BACKEND_SA_NAME}" "${backend_roles[@]}"

# Grant backend service account permission to impersonate pipeline service account
echo -e "\n${YELLOW}Granting pipeline service account impersonation to backend service account...${NC}"
retry_operation "Granting serviceAccountUser role to backend service account" \
    "gcloud iam service-accounts add-iam-policy-binding \
    ${PIPELINE_SA_NAME}@${GCP_PROJECT_ID}.iam.gserviceaccount.com \
    --member=serviceAccount:${BACKEND_SA_NAME}@${GCP_PROJECT_ID}.iam.gserviceaccount.com \
    --role=roles/iam.serviceAccountUser"

echo -e "\n${YELLOW}Creating Cloud SQL instance...${NC}"
retry_operation "Creating Cloud SQL instance ${DB_INSTANCE_NAME}" \
    "gcloud sql instances create \"${DB_INSTANCE_NAME}\" \
    --database-version=POSTGRES_16 \
    --edition=ENTERPRISE \
    --cpu=\"${DB_CPU}\" \
    --memory=\"${DB_MEMORY}GB\" \
    --storage-size=\"${DB_STORAGE}\" \
    --storage-auto-increase \
    --region=${DB_REGION} \
    --require-ssl \
    --root-password=\"${DB_PASSWORD}\" \
    --labels=application=${GCP_APPLICATION_BILLING_LABEL} \
    --project=\"${GCP_PROJECT_ID}\""

echo -e "\n${YELLOW}Creating database and user...${NC}"
# Create database
retry_operation "Creating database ${DB_NAME}" \
    "gcloud sql databases create \"${DB_NAME}\" \
    --instance=\"${DB_INSTANCE_NAME}\" \
    --project=\"${GCP_PROJECT_ID}\""

# Create user
retry_operation "Creating user ${DB_USER}" \
    "gcloud sql users create \"${DB_USER}\" \
    --instance=\"${DB_INSTANCE_NAME}\" \
    --password=\"${DB_PASSWORD}\" \
    --project=\"${GCP_PROJECT_ID}\""

echo -e "\n${YELLOW}Creating Secret Manager secrets...${NC}"
# Create environment variables secret
retry_operation "Creating secret ${ENV_SECRET_NAME}" \
    "gcloud secrets create \"${ENV_SECRET_NAME}\" \
    --replication-policy=\"automatic\" \
    --labels=application=${GCP_APPLICATION_BILLING_LABEL} \
    --project=\"${GCP_PROJECT_ID}\""

# Update secret with environment variables
gcloud secrets versions add "${ENV_SECRET_NAME}" \
    --data-file="${secret_file}" \
    --project="${GCP_PROJECT_ID}"
check_command

# Clean up temporary file
rm "${secret_file}"

echo -e "\n${YELLOW}Creating Vertex AI Pipeline bucket...${NC}"
retry_operation "Creating bucket ${PIPELINE_BUCKET}" \
    "gcloud storage buckets create \"gs://${PIPELINE_BUCKET}\" \
    --project=\"${GCP_PROJECT_ID}\" \
    --location=us-central1 \
    --uniform-bucket-level-access \
    --labels=application=${GCP_APPLICATION_BILLING_LABEL}"

echo -e "\n${YELLOW}Deploying to Cloud Run...${NC}"

# Add Cloud SQL connection permissions to service account
retry_operation "Granting Cloud SQL Client role" \
    "gcloud projects add-iam-policy-binding \"${GCP_PROJECT_ID}\" \
    --member=\"serviceAccount:${BACKEND_SA_NAME}@${GCP_PROJECT_ID}.iam.gserviceaccount.com\" \
    --role=\"roles/cloudsql.client\""

# Get the SQL instance connection name
SQL_CONNECTION_NAME=$(gcloud sql instances describe "${DB_INSTANCE_NAME}" \
    --project="${GCP_PROJECT_ID}" \
    --format="value(connectionName)")

bash "$(dirname "$0")/deploy-cloudrun.sh" \
    --image-path "${IMAGE_PATH}" \
    --project-id "${GCP_PROJECT_ID}" \
    --sql-connection-name "${SQL_CONNECTION_NAME}" \
    --repo-location "${REPO_LOCATION}" \
    --backend-sa-name "${BACKEND_SA_NAME}" \
    --env-secret-name "${ENV_SECRET_NAME}" \
    --service-name "${SERVICE_NAME}" \
    --min-instances "${MIN_INSTANCES}" \
    --max-instances "${MAX_INSTANCES}" \
    --cpu "${CPU}" \
    --memory "${MEMORY}" \
    --timeout "${TIMEOUT}" \
    --application-label "${GCP_APPLICATION_BILLING_LABEL}"
check_command

# Get the service URL
echo -e "\n${YELLOW}Getting service URL...${NC}"
SERVICE_URL=$(gcloud run services describe "${SERVICE_NAME}" \
    --platform=managed \
    --region="${REPO_LOCATION}" \
    --project="${GCP_PROJECT_ID}" \
    --format='value(status.address.url)')

if [ -z "$SERVICE_URL" ]; then
    echo -e "${RED}Failed to get service URL, falling back to status.url...${NC}"
    SERVICE_URL=$(gcloud run services describe "${SERVICE_NAME}" \
        --platform=managed \
        --region="${REPO_LOCATION}" \
        --project="${GCP_PROJECT_ID}" \
        --format='value(status.url)')
fi

echo -e "${GREEN}Service URL: ${SERVICE_URL}${NC}"

# Create a temporary file for the environment variables
tmp_secret=$(mktemp)
# Add header to temp secret to ensure comments at top
echo "# Generated environment variables for Cellarium Nexus" > "$tmp_secret"
echo "# Generated on: $(date)" >> "$tmp_secret"
echo "" >> "$tmp_secret"

# First, get current secret content if it exists
if gcloud secrets versions access latest --secret="${ENV_SECRET_NAME}" --project="${GCP_PROJECT_ID}" >> "$tmp_secret" 2>/dev/null; then
    echo -e "${GREEN}Retrieved existing secret content${NC}"
else
    echo -e "${YELLOW}No existing secret found, creating new one${NC}"
    echo "ENVIRONMENT=\"${ENVIRONMENT}\"" >> "$tmp_secret"
    echo "SECRET_KEY=\"${SECRET_KEY}\"" >> "$tmp_secret"
    echo "GCP_PROJECT_ID=\"${GCP_PROJECT_ID}\"" >> "$tmp_secret"
fi

# Add environment variables from the output file
while IFS= read -r line; do
    # Skip comment and blank lines
    [[ "$line" =~ ^# ]] && continue
    [[ -z "$line" ]] && continue
    # Remove 'export ' prefix and keep the rest
    clean_line="${line#export }"
    # Extract variable name before the equals sign
    var_name="${clean_line%%=*}"
    # Only add if not already in the secret file
    if ! grep -q "^${var_name}=" "$tmp_secret"; then
        echo "${clean_line}" >> "$tmp_secret"
    fi
done < "$output_file"

if [ "$ENVIRONMENT" != "local" ]; then
    # Update URLs with actual Cloud Run values
    SITE_URL="${SERVICE_URL}"
    MAIN_HOST_ALLOWED=$(echo "${SERVICE_URL}" | sed 's|^https://||')
    echo -e "${GREEN}Updating URLs with Cloud Run values:${NC}"
    echo -e "SITE_URL: ${SITE_URL}"
    echo -e "MAIN_HOST_ALLOWED: ${MAIN_HOST_ALLOWED}"

    # Update the secret with new URLs
    # Use temp files for sed on macOS
    sed "/^SITE_URL=/d" "$tmp_secret" > "${tmp_secret}.tmp" && mv "${tmp_secret}.tmp" "$tmp_secret"
    sed "/^MAIN_HOST_ALLOWED=/d" "$tmp_secret" > "${tmp_secret}.tmp" && mv "${tmp_secret}.tmp" "$tmp_secret"
    echo "SITE_URL=${SITE_URL}" >> "$tmp_secret"
    echo "MAIN_HOST_ALLOWED=${MAIN_HOST_ALLOWED}" >> "$tmp_secret"
fi

echo -e "${GREEN}Added SITE_URL and MAIN_HOST_ALLOWED to secret file${NC}"

# Update secret with new environment variables
echo -e "\n${YELLOW}Updating environment secret with SITE_URL and MAIN_HOST_ALLOWED...${NC}"
gcloud secrets versions add "${ENV_SECRET_NAME}" \
    --data-file="${tmp_secret}" \
    --project="${GCP_PROJECT_ID}"
check_command

# Clean up temporary file
rm "$tmp_secret"

# Generate Cloud Console URL
CONSOLE_URL="https://console.cloud.google.com/run/detail/${REPO_LOCATION}/${SERVICE_NAME}/metrics?project=${GCP_PROJECT_ID}"

echo -e "\n${GREEN}Deployment completed successfully!${NC}"
echo -e "${YELLOW}Service URL:${NC} ${SERVICE_URL}"
echo -e "\n${YELLOW}Next steps:${NC}"
echo -e "1. Access your application at: ${SERVICE_URL}"
echo -e "2. Monitor the service in Cloud Console: ${CONSOLE_URL}"
echo -e "3. Set up custom domain if needed"
echo -e "\nConfiguration file has been saved to: ${output_file}"
