"""Classes for source management and processing, for pre-packaged and user-provided data sources."""

import logging
import subprocess
from collections.abc import Iterable
from datetime import datetime
from pathlib import Path
from typing import Any

import frictionless as ftl
import pandas as pd
import requests

logger = logging.getLogger(__name__)

DATASOURCES_PATH = Path(__file__).parent / "datasources"
SPECIFICATIONS_PATH = DATASOURCES_PATH / "specification"


# Generate a list of all available sources currently available in the package's datasources folder.
def _get_available_sources() -> dict[str, Path]:
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

    def __init__(self, name: str, path: Path | str | None = None) -> None:
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
            path = Path(path)
        self.path = path

        # pd.DataFrame: Details on the source, containing author, title, URL and other information, loaded from the folder
        self.details = self._load_details()
        # list[str]: List of features provided by the source (e.g., 'Technologies')
        self.available_features = self._detect_features()

    def _load_details(self) -> pd.DataFrame:
        """Load the source.csv file into a DataFrame."""
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
        features = []
        # Find all CSV files in the source folder
        files = {file.stem for file in self.path.glob("*.csv")}
        files = files - {"sources"}  # Exclude source.csv

        # TODO can check for supported features
        # TODO make more flexible instead of hardcoding
        nicer_names = {
            "technologies": "Technologies",
        }
        features = {nicer_names[f] for f in files}
        features = sorted(features)

        return features

    def process(self, trusted_execution: bool = False) -> None:
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

    @staticmethod
    def archive_file_on_internet_archive(url) -> str | None | Any:
        """
        Ask Internet Archive to save the URL via the Save Page Now API.
        Returns the archived URL or None if failed.
        """
        save_api = f"https://web.archive.org/save/{url}"  # TODO: move it to some config so that it is not hardcoded
        logger.info(f"Requesting archiving for: {url}")
        try:
            response = requests.get(save_api, timeout=30)
            if response.status_code == 200 or response.status_code == 201:
                # The API may redirect to the archived page URL in 'Content-Location' header
                archived_path = response.headers.get("Content-Location")

                if archived_path:
                    archived_url = "https://web.archive.org" + archived_path
                    logger.info(f"Successfully archived URL: {archived_url}")
                    return archived_url
                else:
                    logger.info(
                        "Archive request succeeded but no Content-Location header found."
                    )
                    return None
            else:
                logger.info(
                    f"Archive request failed with status code: {response.status_code}"
                )
                return None
        except Exception as e:
            logger.info(f"Exception during archiving request: {e}")
            return None

    @staticmethod
    def change_datetime_format(
        input_datetime_string, input_datetime_format, output_datetime_format
    ) -> str | None:
        """
        Change the format of a given datetime string from one format to another. This function takes a
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
            None

    @staticmethod
    def get_wayback_snapshot(url) -> str | None | Any:
        """
        The function queries the Internet Archive's Wayback Machine to check for the availability
        of archived snapshots of a given URL. It constructs an API request to the Wayback Machine
        and processes the response to extract the closest archived snapshot information.

        Parameters
        ----------
        url : str
            URL for which to retrieve the Wayback Machine snapshot

        Returns
        -------
        Tuple[str, str, str] | None
            The tuple contains:
                -) archived_url : str URL of the archived snapshot
                -) timestamp : str timestamp of the archived snapshot
                -) status : str status of the archived snapshot
            Returns None if no archived snapshot is found or if an error occurs during the API request.

        Raises
        ------
        requests.RequestException
            If there is an issue with the HTTP request to the Wayback Machine API

        """
        api_url = f"http://archive.org/wayback/available?url={url}"  # TODO: better to put it in the configs
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
                # reformatted_timestamp = timestamp
                reformatted_timestamp = Source.change_datetime_format(
                    timestamp,
                    "%Y%m%d%H%M%S",
                    "%Y-%m-%d %H:%M:%S",
                )
                return archived_url, reformatted_timestamp, status
            else:
                logger.info("No archived snapshot found.")
                var = None

        except requests.RequestException as e:
            logger.info(f"Error during API request: {e}")
            var = None

    @staticmethod
    def download_file(url, local_path=None) -> str | None:
        """

        Download a file from the given URL.
        Saves it to local_path if given, otherwise saves to temporary file with name from URL.
        Returns the local file path or None if failed.
        """
        if local_path is None:
            local_path = url.split("/")[-1].split("?")[0]
            if not local_path:
                local_path = "downloaded_file"

        logger.info(f"Downloading file from {url} to {local_path} ...")
        try:
            with requests.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(local_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                logger.info("Download successful.")
                return local_path
        except Exception as e:
            logger.info(f"Error downloading file: {e}")
            return None

    @staticmethod
    def retrieve_from_internet_archive(archive_url):
        # Check if the file exists in the Internet Archive
        # If it does, download it; if not, return an error message
        pass

    @staticmethod
    def compare_versions(original_file, archived_file):
        # Compare the two files and return a warning if they differ
        pass

    @staticmethod
    def manage_data_source(url):
        # Step 1: Archive the current version of the file
        archive_url = Source.archive_file_on_internet_archive(url)

        # Step 2: Check if the file is still available
        if Source.is_file_available(url):
            # Step 3: Retrieve the file from the original source
            return Source.download_file(url)
        else:
            # Step 4: Retrieve the file from the Internet Archive
            return Source.retrieve_from_internet_archive(archive_url)


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
        self, sources: str | Source | list[str | Source] | dict[str, Path]
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
            self.sources = [Source(name, path) for name, path in sources.items()]
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

    def to_csv(self, path: str | Path) -> None:
        """
        Save the details of the sources to a CSV file.

        Parameters
        ----------
        path : str | Path
            The path to save the CSV file.

        """
        self.details.to_csv(path, index=False)

    def to_datapackage(self, path: str | Path, overwrite: bool = False) -> None:
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
        path = Path(path)

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
