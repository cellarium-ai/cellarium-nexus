Cellarium Nexus
===============

A comprehensive platform for managing, processing, and administering omics data using Google Cloud Platform services and Kubeflow pipelines.

Project Overview
---------------

Cellarium Nexus is a monorepository containing the codebase responsible for:

1. Data ingestion into BigQuery
2. Data extraction from BigQuery into AnnData files
3. Processing and validation of omics datasets
4. Orchestration of data workflows through Kubeflow pipelines

The system leverages Google Cloud Platform services including:

- BigQuery for scalable data storage and querying
- Cloud Storage for file management
- Vertex AI Pipelines for workflow orchestration
- Cloud SQL for application data storage
- Secret Manager for secure configuration

System Architecture
-------------------

Cellarium Nexus is structured around these key components:

Backend
~~~~~~~

- Django-based API server with admin dashboard for system interaction
- Manages metadata, users, and job coordination
- Provides intuitive interface for monitoring and controlling pipelines
- Modules for cell management, ingest management, and curriculum management

Coordinator
~~~~~~~~~~~

- Bridges backend services with data operations
- Handles data ingestion and extraction workflows
- Manages validation of AnnData files

Omics Datastore
~~~~~~~~~~~~~~

- Specialized storage for omics datasets
- Schema definitions for biological data types
- Query interfaces for cellular data

Workflows
~~~~~~~~~

- Kubeflow pipelines for scalable data processing
- Components for creating ingest files, data ingestion, and extraction
- Parallelized execution for high-throughput data operations

Key Workflows
-------------

Data Ingestion
~~~~~~~~~~~~~

1. ``ingest_data_pipeline``: Orchestrates the ingestion process

   - Creates ingest files from input data sources in parallel
   - Ingests prepared data into BigQuery tables
   - Handles multiple input files concurrently

Data Extraction
~~~~~~~~~~~~~~

1. ``extract_data_pipeline``: Manages the extraction process

   - Prepares extract tables with specified features and filters
   - Extracts data to AnnData files in parallel
   - Marks curriculum as finished upon completion

Development
-----------

Project Structure
~~~~~~~~~~~~~~~~~

.. code-block:: text

    cellarium/
    ├── nexus/
    │   ├── backend/              # Django-based API server
    │   ├── client/               # Client library for communicating with backend
    │   ├── coordinator/          # Data processing coordination
    │   ├── omics_datastore/      # Specialized data storage for omics data (powered by BigQuery)
    │   ├── shared/               # Shared utilities and schemas
    │   └── workflows/            # Kubeflow pipeline definitions
    │       └── kubeflow/
    │           ├── components.py # Pipeline component definitions
    │           ├── pipelines.py  # Pipeline orchestration
    │           └── utils/        # Utilities for pipeline operations

Working with Kubeflow Pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The system uses Kubeflow pipelines to orchestrate data workflows:

1. Components (``components.py``): Define individual processing steps

   - ``create_ingest_files_job``: Prepares data for ingestion
   - ``ingest_data_to_bigquery_job``: Loads data into BigQuery
   - ``prepare_extract_job``: Sets up extraction parameters
   - ``extract_job``: Extracts data to AnnData files

2. Pipelines (``pipelines.py``): Orchestrate components into workflows

   - ``ingest_data_pipeline``: Manages the complete ingestion process
   - ``extract_data_pipeline``: Handles the full extraction workflow

3. Configuration (``component_configs.py``): Defines parameters for each component

   - Configurations are passed as YAML files via GCS

Environment and Tooling (Poetry)
--------------------------------

This repository uses Poetry for dependency management, virtual environments, tasks, and packaging.

Install tooling
~~~~~~~~~~~~~~~

Using pip (in a virtual environment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   python3 -m venv .venv-poetry
   source .venv-poetry/bin/activate
   pip install poetry
   poetry self add "poetry-dynamic-versioning[plugin]"
   poetry self add poetry-plugin-export
   poetry --version

Create the environment and install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   poetry install --with dev,test,backend

   # Optional: keep the venv inside the repo
   poetry config virtualenvs.in-project true

Run tests
~~~~~~~~~

.. code-block:: bash

   # All tests
   poetry run poe test

   # Subsets (pytest markers defined under tests/)
   poetry run poe unit
   poetry run poe integration

Lint and format
~~~~~~~~~~~~~~~

.. code-block:: bash

   # Lint check (Ruff + Black --check)
   poetry run poe lint

   # Auto-format (Ruff fixes + Black)
   poetry run poe format

Export pinned requirements for Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   poetry run poe export-backend-reqs
   poetry run poe export-workflows-reqs

Dockerfiles consume the exported files under ``deploy/requirements/`` and install the package with ``pip install .`` for reproducible builds.

Contributing
------------

When contributing to this repository, please follow these guidelines:

1. Use built-in type annotations for all function signatures
2. Write docstrings in imperative mood and reST format
3. Include proper error documentation with ``:raise:`` sections
4. Use absolute imports throughout the codebase
