Deployment Setup
================

This directory contains deployment configurations and scripts for Cellarium Nexus.

Directory Structure
-------------------

- ``backend/``: Backend service deployment files
- ``workflows/``: Workflow service deployment files
- ``setup/``: Setup scripts and configuration files

Requirements
------------

Before running the setup script, ensure you have:

1. **Google Cloud SDK (gcloud)** installed and configured
2. **Billing** enabled on the Google Cloud project
3. The following **IAM roles** assigned to your user account:

   - ``roles/serviceusage.serviceUsageAdmin`` - To enable required Google Cloud APIs
   - ``roles/iam.serviceAccountAdmin`` - To create service accounts for pipelines and backend
   - ``roles/resourcemanager.projectIamAdmin`` - To grant IAM roles to service accounts
   - ``roles/secretmanager.admin`` - To create and manage secrets for environment variables
   - ``roles/cloudsql.admin`` - To create and configure the PostgreSQL database instance
   - ``roles/storage.admin`` - To create GCS buckets
   - ``roles/run.admin`` - To deploy and manage Cloud Run services
   - ``roles/aiplatform.admin`` - For Vertex AI Pipeline configurations
   - ``roles/billing.user`` - To associate resources with billing accounts


Setup Instructions
------------------

The setup process can be run from a local computer with the requirements above:

.. code-block:: bash

    ./deploy/setup/setup.sh

The script will:

1. Create necessary service accounts with appropriate permissions
2. Set up Cloud SQL PostgreSQL database
3. Create Secret Manager secrets for environment variables
4. Create GCS buckets for pipeline data and file storage
5. Deploy the application to Cloud Run