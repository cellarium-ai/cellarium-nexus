Cellarium Nexus
===============
This is a monorepository of the codebase responsible for managing and administering omics data and the pipelines
related to it.

Configuration Variables
--------------------

Set these variables according to your environment:

.. code-block:: bash

    # Project settings
    export PROJECT_ID="your-project-id"

    # Database instance settings
    export DB_INSTANCE_NAME="cellarium-nexus-instance"
    export DB_VERSION="POSTGRES_16"
    export DB_REGION="us-central1"

    # Database credentials
    export DB_NAME="cellarium-nexus-db"
    export DB_USER="cellarium-nexus-db-user"
    export DB_PASSWORD="your-strong-password"

    # Secret Manager settings
    export ENV_SECRET_NAME="cellarium-nexus"

    # Service account settings
    export PIPELINE_SA_NAME="vertex-pipelines-sa"
    export PIPELINE_SA_DISPLAY_NAME="Vertex AI Pipelines Service Account"
    export BACKEND_SA_NAME="cellarium-backend-sa"
    export BACKEND_SA_DISPLAY_NAME="Cellarium Backend Service Account"

    # Vertex AI Pipelines settings
    export PIPELINE_BUCKET="${PROJECT_ID}-pipeline-root"

    # Artifact Registry settings
    export REPO_NAME="cellarium-images"
    export REPO_LOCATION="us-central1"
    export IMAGE_PATH="$REPO_LOCATION-docker.pkg.dev/$PROJECT_ID/$REPO_NAME/cellarium-nexus"

    # Cloud Run settings
    export SERVICE_NAME="cellarium-nexus"

    # GCS bucket settings
    export BUCKET_NAME_PRIVATE="cellarium-nx-file-system"
    export BUCKET_NAME_PUBLIC="cellarium-nx-public-files"

Enable Required APIs
-------------------

Enable all necessary Google Cloud APIs for the project:

.. code-block:: bash

    gcloud services enable \
        secretmanager.googleapis.com \
        aiplatform.googleapis.com \
        storage.googleapis.com

Google Cloud Storage Configuration
-----------------------------

1. Create the required buckets:

.. code-block:: bash

    # Create private bucket with uniform bucket-level access
    gcloud storage buckets create gs://${BUCKET_NAME_PRIVATE} \
        --project=${PROJECT_ID} \
        --location=us-central1 \
        --uniform-bucket-level-access

    # Create public bucket with uniform bucket-level access
    gcloud storage buckets create gs://${BUCKET_NAME_PUBLIC} \
        --project=${PROJECT_ID} \
        --location=us-central1 \
        --uniform-bucket-level-access

2. Configure the public bucket for static files:

   a. Set up CORS policy:

   .. code-block:: bash

       # Apply CORS configuration
       gsutil cors set deploy/cors.json gs://cellarium-nx-public-files

   b. Set up IAM policy for public access:

   .. code-block:: bash

       # Apply bucket IAM policy
       gcloud storage buckets set-iam-policy gs://cellarium-nx-public-files deploy/bucket-policy.yaml

        artifactregistry.googleapis.com \
        run.googleapis.com \
        cloudbuild.googleapis.com \
        sqladmin.googleapis.com \
        bigquery.googleapis.com

Service Accounts Setup
--------------------

Create and configure service accounts for different components of the system:

1. Create service account for Vertex AI Pipelines:

   .. code-block:: bash

       gcloud iam service-accounts create $PIPELINE_SA_NAME \
           --display-name="$PIPELINE_SA_DISPLAY_NAME"

       # Get the full service account email
       export PIPELINE_SA_EMAIL="$PIPELINE_SA_NAME@$PROJECT_ID.iam.gserviceaccount.com"

       # Grant BigQuery Admin role for data operations
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$PIPELINE_SA_EMAIL" \
           --role="roles/bigquery.admin"

       # Grant Storage Admin role for GCS operations
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$PIPELINE_SA_EMAIL" \
           --role="roles/storage.admin"

2. Create service account for backend application:

   .. code-block:: bash

       gcloud iam service-accounts create $BACKEND_SA_NAME \
           --display-name="$BACKEND_SA_DISPLAY_NAME"

       # Get the full service account email
       export BACKEND_SA_EMAIL="$BACKEND_SA_NAME@$PROJECT_ID.iam.gserviceaccount.com"

       # Grant Secret Manager Secret Accessor role
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$BACKEND_SA_EMAIL" \
           --role="roles/secretmanager.secretAccessor"

       # Grant Cloud SQL Client role
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$BACKEND_SA_EMAIL" \
           --role="roles/cloudsql.client"

       # Grant Vertex AI User role for pipeline submission
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$BACKEND_SA_EMAIL" \
           --role="roles/aiplatform.user"

       # Grant BigQuery Admin role for dataset creation and data access
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$BACKEND_SA_EMAIL" \
           --role="roles/bigquery.admin"

       # Grant Storage Admin role for GCS operations
       gcloud projects add-iam-policy-binding $PROJECT_ID \
           --member="serviceAccount:$BACKEND_SA_EMAIL" \
           --role="roles/storage.admin"

Cloud SQL Setup
-------------

This document describes how to set up and configure a Cloud SQL database for Cellarium Nexus.

Create Database Instance
~~~~~~~~~~~~~~~~~~~~~

1. Create the Cloud SQL instance:

   .. code-block:: bash

       gcloud sql instances create $DB_INSTANCE_NAME \
           --database-version=$DB_VERSION \
           --region=$DB_REGION \
           --storage-type=SSD

Database Configuration
--------------------

1. Create database:

   .. code-block:: bash

       gcloud sql databases create $DB_NAME \
           --instance=$DB_INSTANCE_NAME

2. Create user:

   .. code-block:: bash

       gcloud sql users create $DB_USER \
           --instance=$DB_INSTANCE_NAME \
           --password=$DB_PASSWORD

Create Secret in Google Secret Manager
----------------------------------

Create and configure secrets in Google Secret Manager:

1. Create the secret:

   .. code-block:: bash

       # Create secret
       gcloud secrets create $ENV_SECRET_NAME \
           --replication-policy="automatic"

Vertex AI Pipelines Setup
-----------------------

Configure Vertex AI Pipelines for running data ingestion and extraction workflows:

1. Create a pipeline root bucket for storing artifacts:

   .. code-block:: bash

       gsutil mb -p $PROJECT_ID gs://$PIPELINE_BUCKET



After completing these steps, the system will be configured to run data ingestion and extraction pipelines using Vertex AI. The pipelines can be monitored and managed through the Vertex AI Pipelines UI in the Google Cloud Console.

Artifact Registry Setup
--------------------

Set up Artifact Registry to store Docker images:

1. Create a Docker repository:

   .. code-block:: bash

       gcloud artifacts repositories create $REPO_NAME \
           --repository-format=docker \
           --location=$REPO_LOCATION \
           --description="Docker repository for Cellarium Nexus images"

2. Configure Docker to use the repository:

   .. code-block:: bash

       gcloud auth configure-docker $REPO_LOCATION-docker.pkg.dev

Cloud Run Setup
-------------

Deploy the backend application to Cloud Run:

1. Build and push the Docker image:

   .. code-block:: bash

       # Build the image
       docker build -t $IMAGE_PATH:latest .

       # Push to Artifact Registry
       docker push $IMAGE_PATH:latest

2. Deploy to Cloud Run with environment variables from Secret Manager:

   .. code-block:: bash

       # Deploy and capture the service URL
       SERVICE_URL=$(gcloud run deploy $SERVICE_NAME \
           --image=$IMAGE_PATH:latest \
           --region=$REPO_LOCATION \
           --platform=managed \
           --service-account=$BACKEND_SA_EMAIL \
           --set-secrets="/cellarium/nexus/.env=$ENV_SECRET_NAME:latest" \
           --allow-unauthenticated \
           --format="get(status.url)")

   This command:
   - Deploys the application to Cloud Run
   - Uses the service account created earlier
   - Mounts the secret at `/cellarium/nexus/.env` in the container
   - Makes the service publicly accessible (remove `--allow-unauthenticated` for private access)
   - Captures the service URL in the SERVICE_URL variable

Local Environment Setup
---------------------

After deploying to Cloud Run, generate a local .env file for development:

.. code-block:: bash

    # Generate random secret key
    SECRET_KEY=$(openssl rand -hex 32)

    # Create .env file with Cloud Run service URL
    cat > .env << EOL
    SECRET_KEY=${SECRET_KEY}
    ENVIRONMENT=development
    SITE_URL=${SERVICE_URL}
    DB_HOST=localhost
    DB_PORT=5432
    DB_NAME=$DB_NAME
    DB_USER=$DB_USER
    DB_PASSWORD=$DB_PASSWORD
    GCP_PROJECT_ID=$PROJECT_ID
    BUCKET_NAME_PRIVATE=$BUCKET_NAME_PRIVATE
    BUCKET_NAME_PUBLIC=$BUCKET_NAME_PUBLIC
    PIPELINE_ROOT=gs://${PIPELINE_BUCKET}
    PIPELINE_SERVICE_ACCOUNT=$PIPELINE_SA_EMAIL
    EOL

Connection Configuration
~~~~~~~~~~~~~~~~~~~~~

Set up Cloud SQL Proxy for local development:

.. code-block:: bash

    # Download and install Cloud SQL Proxy
    curl -o cloud_sql_proxy https://dl.google.com/cloudsql/cloud_sql_proxy.darwin.amd64
    chmod +x cloud_sql_proxy

    # Start proxy
    ./cloud_sql_proxy \
        --instances=$PROJECT_ID:$DB_REGION:$DB_INSTANCE_NAME=tcp:5432

Best Practices
~~~~~~~~~~~~

* Use strong passwords and store them securely
* Enable automatic backups
* Use Cloud SQL Proxy for local development
* Monitor and optimize slow queries
* Implement proper database indexing strategy
