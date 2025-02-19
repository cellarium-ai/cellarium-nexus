from datetime import datetime
from enum import Enum
from typing import Any, Literal
from uuid import UUID

from nexus.clients.api_schemas import CellInfoAPISchema, FeatureInfoAPISchema, IngestInfoAPISchema
from nexus.clients.base import BaseAPIHTTPClient
from pydantic import BaseModel


class ApiEndpoints:
    CREATE_INGEST_FILE: str = "api/cell_management/ingest/create"
    UPDATE_INGEST_FILE_INFO: str = "api/cell_management/ingest/{id}"
    CREATE_CELL_INFO_BULK: str = "api/cell_management/cell-info/create-bulk"
    CREATE_FEATURE_INFO_BULK: str = "api/cell_management/feature-info/create-bulk"
    CELL_INFO_RESERVE_INDEXES: str = "api/cell_management/cell-info/reserve-indexes"
    FEATURE_INFO_RESERVE_INDEXES: str = "api/cell_management/feature-info/reserve-indexes"
    INGEST_FROM_AVRO: str = "api/cell_management/ingest-from-avro/"


class ReserveIndexesModelType(Enum):
    CELL_INFO: str = "cell_info"
    FEATURE_INFO: str = "feature_info"


class NexusBackendAPIClient(BaseAPIHTTPClient):
    def create_ingest_file_info(self, bigquery_dataset: str) -> IngestInfoAPISchema:
        """
        Create a new ingest file info record.

        :param bigquery_dataset: Name of the BigQuery dataset

        :raise HTTPError: if the request fails
        :raise ValueError: if the response is invalid

        :return: Created ingest file info record
        """
        data = {"bigquery_dataset": bigquery_dataset}
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
        new_status: Literal["SUCCEEDED", "FAILED"],
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

    def create_cell_info_bulk(self, cell_info_api_schema_objects: list[CellInfoAPISchema]) -> None:
        """
        Create multiple cell info records in bulk.

        :param cell_info_api_schema_objects: List of cell info objects to create

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid
        """
        data = [obj.model_dump() for obj in cell_info_api_schema_objects]
        _ = self.post_json(endpoint=ApiEndpoints.CREATE_CELL_INFO_BULK, data=data)

    def create_feature_info_bulk(self, feature_info_api_schema_objects: list[FeatureInfoAPISchema]) -> None:
        """
        Create multiple feature info records in bulk.

        :param feature_info_api_schema_objects: List of feature info objects to create

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid
        """
        data = [obj.model_dump() for obj in feature_info_api_schema_objects]
        _ = self.post_json(endpoint=ApiEndpoints.CREATE_CELL_INFO_BULK, data=data)

    def __reserve_indexes_for_model(
        self, batch_size: int, reserve_model_type: ReserveIndexesModelType
    ) -> tuple[int, int]:
        """
        Reserve a batch of indexes for a specific model type.

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

        api_out = self.post_json(endpoint=endpoint, data={"batch_size": batch_size})
        return api_out["start_index"], api_out["end_index"]

    def reserve_indexes_cell_info(self, batch_size: int) -> tuple[int, int]:
        """
        Reserve a batch of indexes for cell info records.

        :param batch_size: Number of indexes to reserve

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid

        :return: Tuple of (start_index, end_index)
        """
        return self.__reserve_indexes_for_model(
            batch_size=batch_size, reserve_model_type=ReserveIndexesModelType.CELL_INFO
        )

    def reserve_indexes_feature_info(self, batch_size: int) -> tuple[int, int]:
        """
        Reserve a batch of indexes for feature info records.

        :param batch_size: Number of indexes to reserve

        :raise HTTPError: if the request fails
        :raise ValueError: if the request data is invalid

        :return: Tuple of (start_index, end_index)
        """
        return self.__reserve_indexes_for_model(
            batch_size=batch_size, reserve_model_type=ReserveIndexesModelType.FEATURE_INFO
        )

    def ingest_from_avro(self, stage_dir: str, ingest_nexus_uuid: str) -> dict[str, int]:
        """
        Trigger ingestion of CellInfo and FeatureInfo from Avro files.

        :param stage_dir: Base staging directory path where the Avro files are located
        :param ingest_nexus_uuid: UUID of the ingest process

        :raise HTTPError: if the request fails
        :raise ValueError: if the response is invalid

        :return: Dictionary containing counts of ingested records
        """
        data = {"stage_dir": stage_dir, "ingest_nexus_uuid": ingest_nexus_uuid}

        api_out = self.post_json(endpoint=ApiEndpoints.INGEST_FROM_AVRO, data=data)

        return {"cell_info_count": api_out["cell_info_count"], "feature_info_count": api_out["feature_info_count"]}
