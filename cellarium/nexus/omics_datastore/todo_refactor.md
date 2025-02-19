# TODO Refactor Plan

## 1. **Create a New Controller in `bq_ops`**

### **File:** `repo/cellarium/nexus/omics_datastore/bq_ops/bq_datastore_controller.py`

- [x] **Create** the file `bq_datastore_controller.py` inside `repo/cellarium/nexus/omics_datastore/bq_ops/`.
- [x] **Define** a new class, for example `BQDatastoreController`, which will handle:
  - Creating/managing BigQuery datasets/tables.
  - Preparing extract tables.
  - Extracting metadata.
  - Ingesting data to BigQuery.
- [x] **Initialize** the class with:
  ```python
  def __init__(
      self,
      *,
      client: bigquery.Client,
      project: str,
      dataset: str,
  ) -> None:
      ...
  ```
- [x] **Import** and wrap logic from existing modules:
  - `metadata_extractor.py`
  - `prepare_extract.py`
  - `ingest_data_to_bigquery.py` (where appropriate)
  - `create_bq_tables.py`  

---

## 2. **Refactor BQ-Specific Logic**

### **File:** `repo/cellarium/nexus/omics_datastore/controller.py`

- [x] **Identify** methods doing direct BigQuery tasks (e.g. `create_bigquery_dataset`, `ingest_data_to_bigquery`).
- [x] **Move** or **delegate** those methods to the new `BQDatastoreController`.
- [x] **Retain** only the Nexus-backend interactions in `NexusDataController` that are not strictly BigQuery-related.
- [x] **Replace** BigQuery logic calls with calls to `BQDatastoreController`.

---

## 3. **Combine "prepare_extract" Logic**

### **Files to review in `bq_ops`:**
- `repo/cellarium/nexus/omics_datastore/bq_ops/extract/prepare_extract.py`
- `repo/cellarium/nexus/omics_datastore/bq_ops/extract/metadata_extractor.py`  
 
- [x] **Choose** one version of `prepare_extract_tables` to keep.
- [x] **Migrate** or **wrap** its logic into a method (e.g. `prepare_extract_tables`) in `BQDatastoreController`.
- [x] **Decide** if `ExtractTablePreparer` remains as a separate class or merges fully into `BQDatastoreController`.
- [x] **Ensure** consistent docstrings and naming inside your new code to align with your docstring style (reST, imperative mood, etc).

---

## 4. **Integrate `metadata_extractor.py`**

### **File:** `repo/cellarium/nexus/omics_datastore/bq_ops/extract/metadata_extractor.py`

- [x] **Migrate** or **embed** (`MetadataExtractor`) inside `BQDatastoreController`.
- [x] **Wrap** calls to methods like `get_extract_metadata`, `save_metadata` in new higher-level methods (e.g. `extract_and_save_metadata`).
- [x] **Remove** duplication of code if `MetadataExtractor` references appear elsewhere.

---

## 5. **Remove or Unify Duplicate Code**

- [x] **Check** for overlapping logic between `prepare_extract.py` and `prepare_extract_refactored.py`.
- [x] **Consolidate** calls to `ingest_data_to_bigquery` in one place (within the new controller).
- [x] **Ensure** only a single "source of truth" for each operation (e.g., dataset creation, metadata extraction, etc.).

---

## 6. **Test Thoroughly**

- [ ] **Run** unit tests to ensure the new `BQDatastoreController` works as intended.
- [ ] **Verify** integration tests that rely on table creation, ingestion, or metadata extraction still pass.
- [ ] **Fix** or **update** any test code pointing to old modules you've changed/removed.

---

## 7. **Clean Up Old References**

- [ ] **Keep** `prepare_extract.py` for now since it's still being used directly.
- [ ] **Plan** future migration to remove direct usage of `prepare_extract.py`:
  - Update imports in other modules using prepare_extract.py
  - Once all usages are migrated, remove the file
- [ ] **Remove** or adapt old references to "controller.py" if they're replaced by the new, consolidated approach.
- [ ] **Eliminate** duplicated code across `metadata_extractor.py`, `prepare_extract.py`, etc.

---

## 8. **Finalize & Document**

- [ ] **Update** your repository's primary README or relevant docfiles to clarify usage of `bq_datastore_controller.py` as the new single interface for:
  - Creating datasets & tables  
  - Ingesting data  
  - Preparing extracts  
  - Managing metadata  
- [ ] **Provide** usage examples in docstrings or separate docs so others can follow the new structure.
- [ ] **Add** migration guide for moving from old modules to BQDatastoreController.

---

## Future Tasks

### **Update Django Admin Interface**
- [ ] **Update** admin.py to use BQDatastoreController:
  - Update BigQueryDatasetAdmin to use BQDatastoreController for dataset creation
  - Update CellInfoAdmin to use BQDatastoreController for preparing extract tables
  - Update IngestFileInfoAdmin to use BQDatastoreController for data ingestion
- [ ] **Fix** linter errors in admin.py:
  - Fix attribute access for `obs_mappings` and `var_mappings` in ColumnMapping model
  - Fix type conversion for `object_id` in `get_object` function
  - Fix return type for None in HttpResponse
  - Fix attribute access for `obs_mappings` and `var_mappings` in ColumnMappingAdmin
