from storages.backends.gcloud import GoogleCloudStorage


class CustomGoogleCloudStorage(GoogleCloudStorage):
    def gs_path(self, name: str) -> str:
        """
        Returns the full `gs://` path for a given file name.
        """
        return f"gs://{self.bucket.name}/{self._normalize_name(name)}"
