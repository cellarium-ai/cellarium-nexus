======================
Cellarium Nexus Deploy
======================

This directory contains deployment configurations and scripts for Cellarium Nexus.

Directory Structure
------------------

- ``backend/``: Backend service deployment files
- ``workflows/``: Workflow service deployment files
- ``setup/``: Setup scripts and configuration files

Setup Instructions
-----------------

The setup process can be run from a local computer, but requires:

1. An authorized gcloud utility logged in as a project owner
2. Access to the Google Cloud project where Nexus will be deployed

To set up Cellarium Nexus:

.. code-block:: bash

    cd deploy/setup
    ./setup.sh

Django Migrations and Static Files
---------------------------------

Cellarium Nexus uses a combined Cloud Run job approach for handling Django database migrations and static file collection.

Benefits:
~~~~~~~~

- Simplifies the deployment process
- Reduces container startup time
- Ensures operations are performed in the correct order
- Works well with CI/CD pipelines
- Keeps containers stateless and quick to start

Implementation:
~~~~~~~~~~~~~

The ``django_prepare_deploy.sh`` script creates a Cloud Run job that handles:

1. Database migrations
2. Static file collection
3. Superuser creation (if needed)

Execute this job after deploying your application:

.. code-block:: bash

    gcloud run jobs execute django-prepare-deploy --region=${REPO_LOCATION} --project=${PROJECT_ID}

Alternative Approaches:
~~~~~~~~~~~~~~~~~~~~~

- **For very small applications**: Running operations at container startup might be acceptable
- **For complex deployments**: Separate jobs for migrations, static files, and superuser creation might be preferred

Remember that Cloud Run containers should be stateless and quick to start. Moving operations like migrations and static file collection to dedicated Cloud Run jobs helps achieve this goal.