# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""
Source class for representing bibliographic and web sources, with archiving support.

Examples
--------
>>> src = Source(title="Example Source", authors="The Authors")
>>> src.store_in_wayback()
>>> src.retrieve_from_wayback()

"""

import logging
import pathlib
import typing

import pydantic
import requests
import savepagenow

import technologydata

logger = logging.getLogger(__name__)


class Source(pydantic.BaseModel):  # type: ignore
    """
    Represent a data source, including bibliographic and web information.

    Attributes
    ----------
    title : str
        Title of the source.
    authors : str
        Authors of the source.
    url : Optional[str]
        URL of the source.
    url_archive : Optional[str]
        Archived URL.
    url_date : Optional[str]
        Date the URL was accessed.
    url_date_archive : Optional[str]
        Date the URL was archived.

    """

    title: typing.Annotated[str, pydantic.Field(description="Title of the source.")]
    authors: typing.Annotated[str, pydantic.Field(description="Authors of the source.")]
    url: typing.Annotated[
        str | None, pydantic.Field(description="URL of the source.")
    ] = None
    url_archive: typing.Annotated[
        str | None, pydantic.Field(description="Archived URL.")
    ] = None
    url_date: typing.Annotated[
        str | None, pydantic.Field(description="Date the URL was accessed.")
    ] = None
    url_date_archive: typing.Annotated[
        str | None, pydantic.Field(description="Date the URL was archived.")
    ] = None

    def __eq__(self, other: object) -> bool:
        """
        Check for equality with another Source object based on non-None attributes.

        Compares all attributes of the current instance with those of the other object.
        Only compares attributes that are not None in both instances.

        Parameters
        ----------
        other : object
            The object to compare with. Expected to be an instance of Source.

        Returns
        -------
        bool
            True if all non-None attributes are equal between self and other, False otherwise.
            Returns False if other is not a Source instance.

        Notes
        -----
        This method considers only attributes that are not None in both objects.
        If an attribute is None in either object, it is ignored in the comparison.

        """
        if not isinstance(other, Source):
            logger.error("The object is not a Source instance.")
            return False

        for field in self.__class__.model_fields.keys():
            value_self = getattr(self, field)
            value_other = getattr(other, field)
            if value_self != value_other:
                return False
        return True

    def __hash__(self) -> int:
        """
        Return a hash value for the Source instance based on all attributes.

        This method computes a combined hash of the instance's attributes to
        uniquely identify the object in hash-based collections such as sets and dictionaries.

        Returns
        -------
        int
            The hash value of the Source instance.

        """
        # Retrieve all attribute values dynamically
        attribute_values = tuple(
            getattr(self, field) for field in self.__class__.model_fields.keys()
        )
        return hash(attribute_values)

    def __str__(self) -> str:
        """
        Return a string representation of the Source, including all available attributes.

        Returns
        -------
        str
            A string detailing the source's information.

        """
        parts = [f"'{self.authors}': '{self.title}'"]
        if self.url:
            parts.append(f"from url '{self.url}'")
        if self.url_date:
            parts.append(f"last accessed on '{self.url_date}'")
        if self.url_archive:
            parts.append(f"archived at '{self.url_archive}'")
        if self.url_date_archive:
            parts.append(f"on '{self.url_date_archive}'.")
        return ", ".join(parts)

    def ensure_in_wayback(self) -> None:
        """
        Ensure that the source URL is archived in the Wayback Machine.

        This method checks if the source's `url` attribute is set and whether
        an archived URL or archive date is already present. If neither is available, it attempts to archive the
        URL using the Wayback Machine and updates the corresponding attributes.

        Parameters
        ----------
        None

        Returns
        -------
        None
            This method updates the Source object's `url_archive` and `url_date_archive` attributes in place.

        Raises
        ------
        ValueError
            If the `url` attribute is not set (None or NaN).

        Examples
        --------
        >>> from technologydata import Source
        >>> source = Source(url="http://example.com", title="Example Site", authors="The Authors")
        >>> source.ensure_in_wayback()
        A new snapshot has been stored for the url http://example.com with timestamp 2023-10-01T12:00:00Z and Archive.org url http://web.archive.org/web/20231001120000/http://example.com.
        >>> source.url_archive
        'http://web.archive.org/web/20231001120000/http://example.com'
        >>> source.url_date_archive
        '2023-10-01T12:00:00Z'

        """
        if self.url is None:
            raise ValueError(
                f"The url attribute of the source {self.title} is not set or contains a NaN value."
            )

        if self.url_archive is None and self.url_date_archive is None:
            archived_info = self.store_in_wayback(self.url)
            if archived_info is not None:
                archived_url, new_capture_flag, timestamp = archived_info
                if new_capture_flag:
                    logger.info(
                        f"A new snapshot has been stored for the url {self.url} with timestamp {timestamp} and Archive.org url {archived_url}."
                    )
                else:
                    logger.info(
                        f"There is already a snapshot for the url {self.url} with timestamp {timestamp} and Archive.org url {archived_url}."
                    )
                self.url_date_archive = timestamp
                self.url_archive = archived_url

    @staticmethod
    def store_in_wayback(
        url_to_archive: str,
    ) -> tuple[typing.Any, bool | None, str | None] | None:
        """
        Store a snapshot of the given URL on the Wayback Machine and extract the timestamp.

        The method captures the specified URL using the Wayback Machine and retrieves the
        corresponding archive URL along with a formatted timestamp. The timestamp is extracted
        from the archive URL and converted to a more readable format.

        Parameters
        ----------
        url_to_archive : str
            The URL that you want to archive on the Wayback Machine.

        Returns
        -------
        tuple[str, bool, str] | None
            A tuple containing the archive URL, a boolean indicating if a new capture was conducted (if the boolean is
            True, archive.org conducted a new capture. If it is False, archive.org has returned a recently cached capture
            instead, likely taken in the previous minutes) and the formatted timestamp (with format YYYY-MM-DD hh:mm:ss)
            if the operation is successful. Returns None if the timestamp cannot be extracted due to a ValueError (e.g.,
            if the expected substrings are not found in the archive URL).

        Examples
        --------
        >>> from technologydata import Source
        >>> some_url = "some_url"
        >>> archived_info = Source.store_in_wayback(some_url)

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
            output_timestamp = technologydata.Commons.change_datetime_format(
                timestamp,
                technologydata.DateFormatEnum.SOURCES_CSV,
            )
            return archive_url[0], archive_url[1], output_timestamp
        except ValueError:
            # If "web/" or next "/" not found, return empty string
            return None

    def retrieve_from_wayback(
        self, download_directory: pathlib.Path
    ) -> pathlib.Path | None:
        """
        Download a file from the Wayback Machine and save it to a specified path.

        The method retrieves an archived file from the Wayback Machine using the URL
        from the url_archive attribute of the instance. The file is saved in the
        specified format based on its Content-Type field in the Response Header or the extension
        that can be extracted from the URL.

        Parameters
        ----------
        download_directory : pathlib.Path
            The base path where the file will be saved.


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
        - The attribute "url_archived" should contain a valid URL.

        Examples
        --------
        >>> from technologydata import Source
        >>> source = Source(title="example01", authors="The Authors")
        >>> output_path = source.retrieve_from_wayback(pathlib.Path("base_path"))

        """
        if self.url_archive is None:
            logger.error(
                f"The url_archive attribute of source {self.title} is not set."
            )
            return None
        if download_directory is None:
            logger.error(f"The base path of the source {self.title} is not set.")
            return None

        source_title = technologydata.Commons.replace_special_characters(self.title)
        save_path = self._get_save_path(
            self.url_archive, download_directory, source_title
        )

        if save_path is None:
            logger.warning(
                f"It was not possible to determine a file path for the source {source_title}."
            )
            return None

        if save_path.is_file():
            logger.warning(
                f"There is already a file stored at the path {save_path}. Not downloading or overwriting this file."
            )
            return None

        storage_path = self._download_file(self.url_archive, save_path)
        return storage_path

    @staticmethod
    def _get_save_path(
        url_archived: str, source_path: pathlib.Path, source_title: str
    ) -> pathlib.Path | None:
        """
        Determine the save path based on the content type or archived URL.

        The method retrieves the content type of the archived URL and determines the appropriate
        file extension based on the content type or based on the archived URL. It constructs the full save path using
        the provided source path and source title.

        Parameters
        ----------
        url_archived : str
            The URL of the archived file from which the content type will be determined.
        source_path : pathlib.Path
            The base path where the file will be saved.
        source_title : str
            The title of the given source from sources.csv, used as the filename.

        Returns
        -------
        pathlib.Path | None
            The full path where the file should be saved, including the appropriate file extension,
            or None if the content type is unsupported or an error occurs.

        Raises
        ------
        ValueError
            If the extension is not among the supported ones.

        """
        content_type = Source._get_content_type(url_archived)
        if content_type is None:
            return None

        extension = technologydata.FileExtensionEnum.get_extension(
            content_type
        ) or technologydata.FileExtensionEnum.search_file_extension_in_url(url_archived)
        if extension is None:
            raise ValueError(
                f"Unable to infer file extension from content type: {content_type} or URL: {url_archived}"
            )

        if source_path is not None and source_title is not None:
            return pathlib.Path(source_path, source_title + extension)
        else:
            return None

    @staticmethod
    def _get_content_type(url_archived: str) -> typing.Any:
        """
        Fetch the content type of the archived URL.

        The method sends a HEAD request to the specified archived URL to retrieve the
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
            raise requests.exceptions.RequestException(
                f"Failed to retrieve content type: {e}"
            )

    @staticmethod
    def _download_file(
        url_archived: str, save_path: pathlib.Path
    ) -> pathlib.Path | None:
        """
        Download the file and save it to the specified path.

        The method retrieves the content from the specified archived URL and saves it
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
            requests.exceptions.RequestException(
                f"An error occurred during file download: {e}"
            )
            return None
