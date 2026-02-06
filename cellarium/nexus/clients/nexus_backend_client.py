from datetime import datetime
from enum import Enum
from typing import Any, Literal

from cellarium.nexus.clients.api_schemas import (
    BackendResetAPISchema,
    CurriculumAPISchema,
    IngestInfoAPISchema,
    ValidationReportAPISchema,
    ValidationReportItemAPISchema,
)
from cellarium.nexus.clients.base import BaseAPIHTTPClient


class ApiEndpoints:
    CREATE_INGEST: str = "api/ingest-management/ingest/create"
    UPDATE_INGEST: str = "api/ingest-management/ingest/{id}"
    CREATE_INGEST_FILE: str = "api/ingest-management/ingest-file/create"
    UPDATE_INGEST_FILE_INFO: str = "api/ingest-management/ingest-file/{id}"
    CELL_INFO_RESERVE_INDEXES: str = "api/ingest-management/reserve-indexes/cell-info/"
    FEATURE_INFO_RESERVE_INDEXES: str = "api/ingest-management/reserve-indexes/feature-info/"
    REGISTER_CURRICULUM: str = "api/curriculum/curriculums/"
    UPDATE_CURRICULUM: str = "api/curriculum/curriculums/{name}/"
    CREATE_VALIDATION_REPORT: str = "api/ingest-management/validation-report/create/"
    CREATE_VALIDATION_REPORT_ITEM: str = "api/ingest-management/validation-report-item/create/"
    RESET_CACHE: str = "api/ingest-management/reset-cache/"


class ReserveIndexesModelType(Enum):
    CELL_INFO: str = "cell_info"
    FEATURE_INFO: str = "feature_info"


class NexusBackendAPIClient(BaseAPIHTTPClient):
    def create_ingest_file_info(
        self,
        *,
        omics_dataset: str | None = None,
        ingest_id: int | None = None,
        gcs_file_path: str | None = None,
        tag: str | None = None,
    ) -> IngestInfoAPISchema:
        """
        Create a new ingest file info record.

        :param omics_dataset: Name of the Omics dataset
        :param ingest_id: Optional parent ingest ID
        :param gcs_file_path: Optional GCS file path
        :param tag: Optional tag value

        :raise HTTPError: if the request fails
        :raise ValueError: if the response is invalid

        :return: Created ingest file info record
        """
        data: dict[str, Any] = {}
        if omics_dataset is not None:
            data["omics_dataset"] = omics_dataset
        if ingest_id is not None:
            data["ingest_id"] = ingest_id
        if gcs_file_path is not None:
            data["gcs_file_path"] = gcs_file_path
        if tag is not None:
            data["tag"] = tag

        api_out = self.post_json(endpoint=ApiEndpoints.CREATE_INGEST_FILE, data=data)
        return IngestInfoAPISchema(**api_out)

    def update_ingest_metadata_extra(self, ingest_id: int, new_metadata_extra: dict[str, Any]) -> IngestInfoAPISchema:
        """
        Update the metadata_extra field of an ingest file info record.

        :param ingest_id: ID of the ingest record to update
        :param new_metadata_extra: New metadata dictionary to set

        :raise HTTPError: if the request fails
        :raise ValueError: if the response is invalid

        :return: Updated ingest file info record
        """
        endpoint = ApiEndpoints.UPDATE_INGEST_FILE_INFO.format(id=ingest_id)
        data = {"metadata_extra": new_metadata_extra}
        api_out = self.patch_json(endpoint=endpoint, data=data)
        return IngestInfoAPISchema(**api_out)

    def update_ingest_status(
        self,
        ingest_id: int,
        new_status: Literal["STARTED", "SUCCEEDED", "FAILED"],
        ingest_finish_timestamp: datetime | None = None,
    ) -> IngestInfoAPISchema:
        """
        Update the status of an ingest file info record.

        :param ingest_id: ID of the ingest record to update
        :param new_status: New status to set (either "SUCCEEDED" or "FAILED")
        :param ingest_finish_timestamp: Optional timestamp for when the ingest finished

        :raise HTTPError: if the request fails
        :raise ValueError: if the response is invalid

        :return: Updated ingest file info record
        """
        endpoint = ApiEndpoints.UPDATE_INGEST_FILE_INFO.format(id=ingest_id)
        data = {"status": new_status}

        if ingest_finish_timestamp is not None:
            data["ingest_finish_timestamp"] = ingest_finish_timestamp.isoformat()

        api_out = self.patch_json(endpoint=endpoint, data=data)
        return IngestInfoAPISchema(**api_out)

    def __reserve_indexes_for_model(
        self, *, bigquery_dataset: str, batch_size: int, reserve_model_type: ReserveIndexesModelType
    ) -> tuple[int, int]:
        """
        Reserve a batch of indexes for a specific model type and dataset.

        :param bigquery_dataset: Name of the BigQuery dataset to scope the reservation
        :param batch_size: Number of indexes to reserve
        :param reserve_model_type: Type of model to reserve indexes for

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid

        :return: Tuple of (start_index, end_index)
        """
        match reserve_model_type:
            case ReserveIndexesModelType.CELL_INFO:
                endpoint = ApiEndpoints.CELL_INFO_RESERVE_INDEXES
            case ReserveIndexesModelType.FEATURE_INFO:
                endpoint = ApiEndpoints.FEATURE_INFO_RESERVE_INDEXES
            case _:
                raise ValueError("Not supported type")

        api_out = self.post_json(
            endpoint=endpoint,
            data={
                "omics_dataset": bigquery_dataset,
                "batch_size": batch_size,
            },
        )
        return api_out["index_start"], api_out["index_end"]

    def reserve_indexes_cell_info(self, *, bigquery_dataset: str, batch_size: int) -> tuple[int, int]:
        """
        Reserve a batch of indexes for cell info records.

        :param bigquery_dataset: Name of the BigQuery dataset to scope the reservation
        :param batch_size: Number of indexes to reserve

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid

        :return: Tuple of (start_index, end_index)
        """
        return self.__reserve_indexes_for_model(
            bigquery_dataset=bigquery_dataset,
            batch_size=batch_size,
            reserve_model_type=ReserveIndexesModelType.CELL_INFO,
        )

    def reserve_indexes_feature_info(self, *, bigquery_dataset: str, batch_size: int) -> tuple[int, int]:
        """
        Reserve a batch of indexes for feature info records.

        :param bigquery_dataset: Name of the BigQuery dataset to scope the reservation
        :param batch_size: Number of indexes to reserve

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid

        :return: Tuple of (start_index, end_index)
        """
        return self.__reserve_indexes_for_model(
            bigquery_dataset=bigquery_dataset,
            batch_size=batch_size,
            reserve_model_type=ReserveIndexesModelType.FEATURE_INFO,
        )

    def register_curriculum(
        self,
        *,
        name: str,
        creator_id: int,
        extract_bin_size: int,
        filters_json: dict[str, Any] | None = None,
    ) -> CurriculumAPISchema:
        """
        Register a new curriculum.

        :param name: Name of the curriculum
        :param creator_id: ID of the creator
        :param extract_bin_size: Size of extraction bins
        :param filters_json: Optional filters in JSON format

        :raise: HTTPError, ValueError

        :return: Registered curriculum object
        """
        data = {
            "name": name,
            "creator_id": creator_id,
            "extract_bin_size": extract_bin_size,
        }

        if filters_json is not None:
            data["filters_json"] = filters_json

        api_out = self.post_json(endpoint=ApiEndpoints.REGISTER_CURRICULUM, data=data)
        return CurriculumAPISchema(**api_out)

    def update_curriculum(
        self,
        *,
        name: str,
        status: Literal["EXTRACTING", "SUCCEEDED", "FAILED"] | None = None,
        cell_count: int | None = None,
        extract_bin_count: int | None = None,
        extract_files_path: str | None = None,
        metadata_file_path: str | None = None,
    ) -> CurriculumAPISchema:
        """
        Update curriculum information including status and extract metadata.

        :param name: Name of the curriculum to update
        :param status: New status to set (either "EXTRACTING", "SUCCEEDED", or "FAILED")
        :param cell_count: Total number of cells in the extract
        :param extract_bin_count: Number of extract bins
        :param extract_files_path: Directory containing the extract files
        :param metadata_file_path: Path to the metadata file

        :raise: HTTPError, ValueError

        :return: Updated curriculum object
        """
        data = {}
        if status is not None:
            data["status"] = status
        if cell_count is not None:
            data["cell_count"] = cell_count
        if extract_bin_count is not None:
            data["extract_bin_count"] = extract_bin_count
        if extract_files_path is not None:
            data["extract_files_path"] = extract_files_path
        if metadata_file_path is not None:
            data["metadata_file_path"] = metadata_file_path

        api_out = self.patch_json(endpoint=ApiEndpoints.UPDATE_CURRICULUM.format(name=name), data=data)
        return CurriculumAPISchema(**api_out)

    def create_validation_report(self) -> ValidationReportAPISchema:
        """
        Create a new validation report.

        The creator is set from the authenticated request context on the backend.

        :raise: HTTPError if the request fails
        :raise: ValueError if the response is invalid

        :return: Created validation report
        """
        api_out = self.post_json(endpoint=ApiEndpoints.CREATE_VALIDATION_REPORT, data={})
        return ValidationReportAPISchema(**api_out)

    def create_validation_report_item(
        self,
        *,
        report_id: int,
        input_file_path: str,
        validator_name: str,
        is_valid: bool,
        message: str | None = None,
        sanitized_file_path: str | None = None,
    ) -> ValidationReportItemAPISchema:
        """
        Create a validation report item for an existing validation report.

        :param report_id: ID of the existing validation report
        :param input_file_path: Path to the input file that was validated
        :param validator_name: Name of the validator that performed the validation
        :param is_valid: Whether the validation passed or failed
        :param message: Optional message with validation details
        :param sanitized_file_path: Optional path to the sanitized output file

        :raise: HTTPError if the request fails
        :raise: ValueError if the response is invalid

        :return: Created validation report item
        """
        data = {
            "report_id": report_id,
            "input_file_path": input_file_path,
            "validator_name": validator_name,
            "is_valid": is_valid,
        }

        if message is not None:
            data["message"] = message

        if sanitized_file_path is not None:
            data["sanitized_file_path"] = sanitized_file_path

        api_out = self.post_json(endpoint=ApiEndpoints.CREATE_VALIDATION_REPORT_ITEM, data=data)
        return ValidationReportItemAPISchema(**api_out)

    def reset_backend_cache(self) -> BackendResetAPISchema:
        """
        Reset and repopulate the cache for all cached filters.

        This method triggers a complete cache reset and repopulation of all
        cached filters used in the admin interface. It requires admin privileges
        to execute.

        :raise: HTTPError if the request fails or if the user doesn't have admin privileges

        :return: Dictionary containing status, message, repopulated keys, and count
        """
        api_out = self.post_json(endpoint=ApiEndpoints.RESET_CACHE)
        return BackendResetAPISchema(**api_out)
