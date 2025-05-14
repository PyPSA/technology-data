"""Classes for source management and processing, for pre-packaged and user-provided data sources."""

import logging
import pathlib
import subprocess
from collections.abc import Iterable
from datetime import datetime
from typing import Any

import frictionless as ftl
import pandas as pd
import requests
import savepagenow

logger = logging.getLogger(__name__)

DATASOURCES_PATH = pathlib.Path(__file__).parent / "datasources"
SPECIFICATIONS_PATH = DATASOURCES_PATH / "specification"


# Generate a list of all available sources currently available in the package's datasources folder.
def _get_available_sources() -> dict[str, pathlib.Path]:
    """
    Determine all available sources based on the folders in datasources.

    Returns
    -------
    dict[str, Path]
        A dictionary with the source name as the key and the path to the source folder as the value.

    """
    # Valid sources are folders and that contain at least one file named sources.csv and
    # one additional other file (a feature) ending in `.csv`
    sources = [
        folder
        for folder in DATASOURCES_PATH.iterdir()
        if folder.glob("sources.csv") and len(list(folder.glob("*.csv"))) > 1
    ]

    return {folder.stem: folder for folder in sources}


# Dictionary of available sources packaged and available for the user
AVAILABLE_SOURCES = _get_available_sources()


class Source:
    """
    A source of data that can provide one or more data features.

    Attributes
    ----------
    name : str
        A unique name for distinguishing the source, either the name of a pre-packaged source or a user provided name.
    path : Path
        Path to the source's folder.
    details : pd.DataFrame
        Details on the source, containing author, title, URL and other information.
    available_features : list[str]
        List of features provided by the source (e.g., 'Technologies').

    """

    def __init__(self, name: str, path: pathlib.Path | str | None = None) -> None:
        """
        Create a source of data that can provide one or more data features.

        Some data sources are already pre-packaged. In addition, data sources can be loaded from local paths.
        Possible data features are features that are understood by this package, e.g. 'Technologies'.
        Data sources are usually included in the package pre-processed, in the case of automatically extractd data
        they also come with a processing script that allows to re-process and thereby reproduce the data.

        Parameters
        ----------
        name : str
            A unique name for distinguishing the source, either indicating a pre-packaged source or a new source from a local path.
        path : Path | str
            Path to the source's folder. Not required for loading a pre-packaged source.

        Examples
        --------
        >>> from technologydata import Source
        >>> source = Source("example01")
        >>> source = Source("example01", "./some/local/path/") # folder contains sources.csv and other data files

        """
        self.name = name

        if name in AVAILABLE_SOURCES and path is None:
            # Use prepackaged source
            path = AVAILABLE_SOURCES[name]
        elif name not in AVAILABLE_SOURCES and path is None:
            raise ValueError(
                f"No pre-packaged source with the name {name} found. "
                f"Check the available sources in {DATASOURCES_PATH} or provide a path to a local source."
            )

        # Ensure path is a Path object
        if isinstance(path, str):
            path = pathlib.Path(path)
        self.path = path

        # pd.DataFrame: Details on the source, containing author, title, URL and other information, loaded from the folder
        self.details = self._load_details()
        # list[str]: List of features provided by the source (e.g., 'Technologies')
        self.available_features = self._detect_features()

    def _load_details(self) -> pd.DataFrame | None:
        """Load the sources.csv file into a DataFrame."""
        if self.path is None:
            logger.error("The path attribute is not set.")
            return None

        source_file = self.path / "sources.csv"
        schema_file = SPECIFICATIONS_PATH / "sources.schema.json"

        if not source_file.exists():
            raise FileNotFoundError(f"Missing sources.csv in {self.path}")

        resource = ftl.Resource(path=str(source_file), schema=str(schema_file))
        report = resource.validate()

        if not report.valid:
            raise ValueError(f"source.csv in {self.path} is invalid: {report}")

        details = resource.to_pandas()

        # Add the source name as first column
        details.insert(0, "source_name", self.name)

        return details

    def _detect_features(self) -> list[str]:
        """Detect which features are available in the source folder."""
        # Find all CSV files in the source folder

        if self.path is None:
            logger.error("The path attribute is not set.")
            return []  # return an empty list

        files = {file.stem for file in self.path.glob("*.csv")}
        files = files - {"sources"}  # Exclude source.csv

        # TODO can check for supported features
        # TODO make more flexible instead of hardcoding
        nicer_names = {
            "technologies": "Technologies",
        }
        features = list({nicer_names[f] for f in files})
        features = sorted(features)

        return features

    def process(self, trusted_execution: bool = False) -> bool:
        """
        Process the source data.

        If a `process.py` script is present in the source folder, this method is available and run the script.
        Otherwise, a message will be logged indicating that the process method is unavailable.

        For security reasons, `process.py` is only enabled by default for `process.py` scripts inside the package,
        for sources with a folder outside of the package, `trusted_execution` must be set to True.

        Parameters
        ----------
        trusted_execution : bool
            By default, only `process.py` scripts inside the package are trusted and executed,
            outside data source processing scripts are not executed. To allow their execution,
            set this parameter explicitly to True.

        Returns
        -------
        bool
            True if the process.py script was executed successfully, False otherwise.

        """
        raise NotImplementedError(
            "This functionality is not implemented yet and should probably be redesigned. "
            "For now enter new data manually. Automatic processing can still be used, but is"
            "not yet stitched into the package."
        )
        process_script = self.path / "process.py"

        # Check if the process.py script exists or raise an error that the source doesn't support automatic processing
        if not process_script.exists():
            logger.info(
                f"No process.py script found for source: {self.path.stem} - Source does not support automatic processing."
            )
            return False

        # Check if the source is part of the package, in which case the script is trusted and can be executed
        if self.path.parent == DATASOURCES_PATH:
            trusted_execution = True
        elif not trusted_execution:
            logger.warning(
                f"Untrusted execution of process.py script for source: {self.path.stem}. Set trusted_execution=True to execute."
            )
            return False

        logger.debug(f"Executing process.py for source: {self.path.stem}")
        p = subprocess.popen(
            ["python", str(process_script)],
            cwd=self.path,  # Run the script in the source folder
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        (stdoutdata, stderrdata) = p.communicate()

        if p.returncode != 0:
            logger.debug(
                f"Process.py execution failed for source: {self.path.stem} with error code {p.returncode}.\n\n"
                f"Error message: {stderrdata.decode()}\n"
                f"{stdoutdata.decode()}"
            )
            return False
        else:
            logger.debug(
                f"Process.py executed successfully for source: {self.path.stem}"
            )
            return True

    def ensure_snapshot(self) -> None:
        """
        Ensure that the Source object has the url_archived and url_date fields populated.
        If not, check if the URL already has a snapshot stored. If not, store it and populate
        url_archived and url_date.

        Returns
        -------
        None

        """
        if self.details is None:
            logger.error(f"The details attribute of the source {self.name} is not set.")
            return None

        if self.details.url is None:
            logger.error(f"The url attribute of the source {self.name} is not set.")
            return None

        # Check if the url_date is populated. If not, store a snapshot
        if (
            self.details["url_date"].isna().all()
            or not self.details["url_archived"].isna().all()
        ):
            archived_info = self.store_snapshot_on_wayback(
                self.details["url"].to_numpy()[0]
            )
            if archived_info is not None:
                archived_url, new_capture_flag, timestamp = archived_info
                if new_capture_flag:
                    logger.info(
                        f"A new snapshot has been stored for the url {self.details.url} with timestamp {timestamp} and Archive.org url {archived_url}."
                    )
                else:
                    logger.info(
                        f"There is already a snapshot for the url {self.details.url}."
                    )
                self.details.url_date = timestamp
                self.details.url_archived = archived_url

                # Update the existing .csv file with the new attributes
                self.details.to_csv(self.path, index=False)
        else:
            logger.info(
                f"Both url_date and url_archived are present for the source {self.name}"
            )

    @staticmethod
    def store_snapshot_on_wayback(
        url_to_archive: str,
    ) -> tuple[Any, bool | None, str | None] | None:
        """
        Store a snapshot of the given URL on the Wayback Machine and extract the timestamp.
        This method captures the specified URL using the Wayback Machine and retrieves the
        corresponding archive URL along with a formatted timestamp. The timestamp is extracted
        from the archive URL and converted to a more readable format.

        Parameters
        ----------
        url_to_archive : str
            The URL that you want to archive on the Wayback Machine.

        Returns
        -------
        tuple[str, bool, str] | None
            A tuple containing the archive URL, a boolean indicating if a new capture was conducted (if the boolean is True,
            archive.org conducted a new capture. If it is False, archive.org has returned a recently cached capture
            instead, likely taken in the previous minutes) and the formatted timestamp (with format YYYY-MM-DD hh:mm:ss)
            if the operation is successful. Returns None if the timestamp cannot be extracted due to a ValueError (e.g.,
            if the expected substrings are not found in the archive URL).

        """
        archive_url = savepagenow.capture_or_cache(url_to_archive)
        try:
            # The timestamp is between "web/" and the next "/" afterward
            # Find the starting index of "web/"
            start_index = archive_url[0].index("web/") + len("web/")
            # Find the ending index of the timestamp by locating the next "/" after the start_index
            end_index = archive_url[0].index("/", start_index)
            # Extract the timestamp substring
            timestamp = archive_url[0][start_index:end_index]
            output_timestamp = Source.change_datetime_format(
                timestamp,
                "%Y%m%d%H%M%S",
                "%Y-%m-%d %H:%M:%S",
            )
            return archive_url[0], archive_url[1], output_timestamp
        except ValueError:
            # If "web/" or next "/" not found, return empty string
            return None

    def download_file_from_wayback(self) -> pathlib.Path | None:
        """
        Download a file from the Wayback Machine and save it to a specified path.

        This method retrieves an archived file from the Wayback Machine using the URL
        stored in the `details` attribute of the instance. The file is saved in the
        specified format based on its Content-Type field in the Response Header.
        Supported formats include:
        - Plain text (.txt)
        - PDF (.pdf)
        - Excel (.xls and .xlsx)
        - Parquet (.parquet)

        Returns
        -------
        pathlib.Path | None
            The specified path where the file is stored, or None if an error occurs.

        Raises
        ------
        requests.exceptions.RequestException
            If there is an issue with the HTTP request.

        Notes
        -----
        - The `details` attribute must contain a key "url_archived" with a valid URL.
        - The `path` and `name` attributes must be defined in the instance for saving the file.

        """
        if self.details is None:
            logger.error(f"The details attribute of the source {self.name} is not set.")
            return None

        if (
            "url_archived" not in self.details
            or self.details["url_archived"].isna().all()
        ):
            logger.error(
                f"The url_archived attribute of source {self.name} is not set."
            )
            return None

        if self.path is None:
            logger.error(f"The path attribute of the source {self.name} is not set.")
            return None

        url_archived = self.details["url_archived"].to_numpy()[0]
        save_path = self._get_save_path(url_archived)

        if save_path is None:
            return None

        return self._download_file(url_archived, save_path)

    def _get_save_path(self, url_archived: str) -> pathlib.Path | None:
        """
        Determine the save path based on the content type.
        This method retrieves the content type of the archived URL and determines the appropriate
        file extension based on the content type. It constructs the full save path using the
        instance's `path` and `name` attributes.

        Parameters
        ----------
        url_archived : str
            The URL of the archived file from which the content type will be determined.

        Returns
        -------
        pathlib.Path | None
            The full path where the file should be saved, including the appropriate file extension,
            or None if the content type is unsupported or an error occurs.

        Notes
        -----
            - Supported content types and their corresponding file extensions:
                - "text/plain" -> ".txt"
                - "application/pdf" -> ".pdf"
                - "application/vnd.ms-excel" -> ".xls"
                - "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet" -> ".xlsx"
                - "application/parquet" -> ".parquet"
            - If the content type is unsupported, a warning is logged.

        """
        content_type = self._get_content_type(url_archived)
        if content_type is None:
            return None

        extension_map = {
            "text/plain": ".txt",
            "application/pdf": ".pdf",
            "application/vnd.ms-excel": ".xls",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet": ".xlsx",
            "application/parquet": ".parquet",
        }

        extension = extension_map.get(content_type)
        if extension is None:
            logger.warning(f"Unsupported content type: {content_type}")
            return None

        if self.path is not None and self.name is not None:
            return pathlib.Path(self.path, self.name + extension)
        else:
            return None

    @staticmethod
    def _get_content_type(url_archived: str) -> Any:
        """
        Fetch the content type of the archived URL.

        This method sends a HEAD request to the specified archived URL to retrieve the
        Content-Type from the response headers. It returns the content type as a string
        if the request is successful; otherwise, it logs an error and returns None.

        Parameters
        ----------
        url_archived : str
            The URL of the archived resource for which the content type is to be fetched.

        Returns
        -------
        str | None
            The Content-Type of the archived URL if the request is successful, or None
            if an error occurs during the request.

        Raises
        ------
        requests.exceptions.RequestException
            If there is an issue with the HTTP request, an error is logged, and None is returned.

        """
        try:
            response = requests.head(url_archived)
            response.raise_for_status()
            return response.headers.get("Content-Type")
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to retrieve content type: {e}")
            return None

    @staticmethod
    def _download_file(
        url_archived: str, save_path: pathlib.Path
    ) -> pathlib.Path | None:
        """
        Download the file and save it to the specified path.

        This method retrieves the content from the specified archived URL and saves it
        to the provided file path. It handles HTTP errors and logs appropriate messages
        based on the outcome of the download operation.

        Parameters
        ----------
        url_archived : str
            The URL of the archived file to be downloaded.

        save_path : pathlib.Path
            The path where the downloaded file will be saved, including the file name.

        Returns
        -------
        pathlib.Path | None
            The path where the file has been saved if the download is successful, or None
            if an error occurs during the download process.

        Raises
        ------
        requests.exceptions.RequestException
            If there is an issue with the HTTP request, an error is logged, and None is returned.

        """
        try:
            response = requests.get(url_archived)
            response.raise_for_status()  # Check for HTTP errors

            with open(save_path, "wb") as file:
                file.write(response.content)

            logger.info(f"File downloaded successfully and saved to {save_path}")
            return save_path
        except requests.exceptions.RequestException as e:
            logger.error(f"An error occurred during file download: {e}")
            return None

    @staticmethod
    def change_datetime_format(
        input_datetime_string: str,
        input_datetime_format: str,
        output_datetime_format: str,
    ) -> str | None:
        """
        Change the format of a given datetime string from one format to another. This method takes a
        datetime string and its current format, then converts it to a specified output format.
        If the input string does not match the provided input format, it logs an error and returns None.

        Parameters
        ----------
        input_datetime_string : str
            datetime string that needs to be reformatted

        input_datetime_format : str
            format of the input datetime string, following the strftime format codes

        output_datetime_format : str
            desired format for the output datetime string, following the strftime format codes.

        Returns
        -------
           str | None
               reformatted datetime string if successful, otherwise None

        Raises
        ------
        ValueError
            If the input datetime string does not match the input format.

        """
        try:
            dt = datetime.strptime(input_datetime_string, input_datetime_format)
            logger.info(
                f"The datetime string follows the format {input_datetime_format}"
            )
            output_datetime_string = dt.strftime(output_datetime_format)
            logger.info(f"The format is now changed to {output_datetime_format}")
            return output_datetime_string
        except ValueError as e:
            logger.info(f"Error during datetime formatting: {e}")
            return None

    @staticmethod
    def is_wayback_snapshot_available(
        input_url: str, input_timestamp: str
    ) -> str | None | Any:
        """
        The method queries the Internet Archive's Wayback Machine to check for the availability
        of archived snapshots of a given URL. It constructs an API request to the Wayback Machine
        and processes the response to extract the closest archived snapshot information.

        Parameters
        ----------
        input_url : str
            URL for which to retrieve the Wayback Machine snapshot
        input_timestamp : str
            timestamp for which to retrieve the Wayback Machine snapshot

        Returns
        -------
        Tuple[str, str, str] | None
            The tuple contains:
                -) archived_url : str URL of the archived snapshot
                -) timestamp : str timestamp of the archived snapshot with format YYYY-MM-DD hh:mm:ss
                -) status : str status of the archived snapshot
            Returns None if no archived snapshot is found or if an error occurs during the API request.

        Raises
        ------
        requests.RequestException
            If there is an issue with the HTTP request to the Wayback Machine API

        """
        api_timestamp = Source.change_datetime_format(
            input_timestamp,
            "%Y-%m-%d %H:%M:%S",
            "%Y%m%d%H%M%S",
        )
        api_url = f"http://archive.org/wayback/available?url={input_url}&timestamp={api_timestamp}"
        try:
            response = requests.get(api_url)
            # Raise an error for bad HTTP status codes
            response.raise_for_status()
            data = response.json()

            # Extract the closest archived snapshot information
            closest = data.get("archived_snapshots", {}).get("closest", {})
            if closest:
                available = closest.get("available", False)
                archived_url = closest.get("url", "")
                timestamp = closest.get("timestamp", "")
                status = closest.get("status", "")

                logger.info(f"Available: {available}")
                logger.info(f"Archived URL: {archived_url}")
                logger.info(f"Timestamp: {timestamp}")
                logger.info(f"Status: {status}")

                output_timestamp = Source.change_datetime_format(
                    timestamp,
                    "%Y%m%d%H%M%S",
                    "%Y-%m-%d %H:%M:%S",
                )

                return archived_url, output_timestamp, status
            else:
                logger.info("No archived snapshot found.")
                return None

        except requests.RequestException as e:
            logger.info(f"Error during API request: {e}")
            return None


class Sources:
    """
    A collection of data sources that can be loaded and used.

    Attributes
    ----------
    sources : list[Source]
        A list of Source objects.
    details : pd.DataFrame
        A DataFrame containing the details of all the sources, like author, title, URL, etc. .
    available_features : pd.DataFrame
        A DataFrame containing all source names and their available features.

    """

    schema_name = "sources"
    # Load the datapackage schema to be able to validate against it
    schema = ftl.Schema(str(SPECIFICATIONS_PATH / (schema_name + ".schema.json")))

    def __init__(
        self, sources: str | Source | list[str | Source] | dict[str, pathlib.Path]
    ) -> None:
        """
        Create a collection of data sources.

        Parameters
        ----------
        sources : str | Source | list[Source | str] | dict[str:Path]
            A list of Source objects representing the data sources, names of pre-packaged sources,
            or a dictionary with source names as keys and paths as values.

        """
        self.sources = []

        if isinstance(sources, str):
            # Single source name specified
            if sources in AVAILABLE_SOURCES:
                self.sources = [Source(sources)]
            else:
                raise ValueError(
                    f"Source {sources} not found. Check the available sources in {DATASOURCES_PATH}"
                )
        elif isinstance(sources, Source):
            self.sources = [sources]
        elif isinstance(sources, dict):
            # All sources need to be loaded when specified this way
            self.sources = [
                Source(
                    name, path if isinstance(path, pathlib.Path) else pathlib.Path(path)
                )
                for name, path in sources.items()
            ]
        elif isinstance(sources, Iterable):
            # Directly add already loaded sources to the object
            loaded_sources = [
                source for source in sources if isinstance(source, Source)
            ]
            # Sources specified by name need to be loaded first
            unloaded_sources = [
                Source(source) for source in sources if isinstance(source, str)
            ]
            # Contains all the loaded sources
            self.sources = loaded_sources + unloaded_sources
        else:
            raise TypeError(f"Unsupported type for sources: {type(sources)}")

        # A pd.DataFrame, containing the details of all the sources
        self.details = pd.concat(
            [source.details for source in self.sources], ignore_index=True
        ).sort_values(by="source_name", ascending=True)

        # Ad pd.DataFrame, containing the sources and their available features
        self.available_features = pd.DataFrame(
            [
                {
                    "source_name": source.name,
                    "available_features": source.available_features,
                }
                for source in self.sources
            ]
        )

    def to_csv(self, path: str | pathlib.Path) -> None:
        """
        Save the details of the sources to a CSV file.

        Parameters
        ----------
        path : str | Path
            The path to save the CSV file.

        """
        self.details.to_csv(path, index=False)

    def to_datapackage(self, path: str | pathlib.Path, overwrite: bool = False) -> None:
        """
        Export the data to a folder following the datapackage specification.

        Params
        -------
        path : str|Path
            The path to save the data including sources and schema to.
            Must be a non-existing or empty folder.
            To overwrite existing files in a folder, use the `overwrite` parameter.
        overwrite : bool
            Existing files with the same name in the target path will be overwritten, default is False.
        """
        path = pathlib.Path(path)

        # Check if the path exists and is empty
        if path.exists():
            if not path.is_dir():
                raise ValueError(f"Path {path} is not a directory.")
            if len(list(path.iterdir())) > 0 and not overwrite:
                raise ValueError(
                    f"Path {path} is not empty. To overwrite existing files, set `overwrite=True`."
                )

        # Safe to write beyond this point
        # Write the object data and its schema
        data = self.details.drop(
            columns=["source_name"]
        )  # not part of the schema, drop manually here
        ftl.Resource(data).write(path=str(path / f"{self.schema_name}.csv"))
        self.schema.to_json(path=str(path / f"{self.schema_name}.schema.json"))
