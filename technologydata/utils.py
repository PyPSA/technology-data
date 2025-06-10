"""Classes for utils methods."""

import logging
import re
from enum import Enum
from typing import Any

from dateutil import parser

logger = logging.getLogger(__name__)


class DateFormatEnum(str, Enum):
    """
    Enum for date formats used in different sources.

    Attributes
    ----------
    SOURCES_CSV : str
        Date format for CSV sources, e.g., "2023-10-01 12:00:00".
    WAYBACK : str
        Date format for Wayback Machine, e.g., "20231001120000".
    NONE : str
        Represents an empty date format.

    """

    SOURCES_CSV = "%Y-%m-%d %H:%M:%S"
    WAYBACK = "%Y%m%d%H%M%S"
    NONE = ""


class FileExtensionEnum(Enum):
    """
    An enumeration that maps various file extensions to their corresponding MIME types.

    This Enum provides a structured way to associate common file extensions with their respective
    MIME types, facilitating easy retrieval of file extensions based on content types. Each member
    of the enumeration is a tuple containing the file extension and its associated MIME type.

    Members
    --------
    TEXT_PLAIN : tuple
        Represents the MIME type "text/plain" with the file extension ".txt".
    TEXT_HTML : tuple
        Represents the MIME type "text/html" with the file extension ".html".
    TEXT_CSV : tuple
        Represents the MIME type "text/csv" with the file extension ".csv".
    TEXT_XML : tuple
        Represents the MIME type "text/xml" with the file extension ".xml".
    APPLICATION_MS_EXCEL : tuple
        Represents the MIME type "application/vnd.ms-excel" with the file extension ".xls".
    APPLICATION_ODS : tuple
        Represents the MIME type "application/vnd.oasis.opendocument.spreadsheet" with the file extension ".ods".
    APPLICATION_OPENXML_EXCEL : tuple
        Represents the MIME type "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        with the file extension ".xlsx".
    APPLICATION_JSON : tuple
        Represents the MIME type "application/json" with the file extension ".json".
    APPLICATION_XML : tuple
        Represents the MIME type "application/xml" with the file extension ".xml".
    APPLICATION_PDF : tuple
        Represents the MIME type "application/pdf" with the file extension ".pdf".
    APPLICATION_PARQUET : tuple
        Represents the MIME type "application/parquet" with the file extension ".parquet".
    APPLICATION_VDN_PARQUET : tuple
        Represents the MIME type "application/vdn.apache.parquet" with the file extension ".parquet".
    APPLICATION_RAR_WINDOWS : tuple
        Represents the MIME type "application/x-rar-compressed" with the file extension ".rar".
    APPLICATION_RAR : tuple
        Represents the MIME type "application/vnd.rar" with the file extension ".rar".
    APPLICATION_ZIP : tuple
        Represents the MIME type "application/zip" with the file extension ".zip".
    APPLICATION_ZIP_WINDOWS : tuple
        Represents the MIME type "application/x-zip-compressed" with the file extension ".zip".
    """

    TEXT_PLAIN = (".txt", "text/plain")
    TEXT_HTML = (".html", "text/html")
    TEXT_CSV = (".csv", "text/csv")
    TEXT_XML = (".xml", "text/xml")
    APPLICATION_MS_EXCEL = (".xls", "application/vnd.ms-excel")
    APPLICATION_ODS = (".ods", "application/vnd.oasis.opendocument.spreadsheet")
    APPLICATION_OPENXML_EXCEL = (
        ".xlsx",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    APPLICATION_JSON = (".json", "application/json")
    APPLICATION_XML = (".xml", "application/xml")
    APPLICATION_PDF = (".pdf", "application/pdf")
    APPLICATION_PARQUET = (".parquet", "application/parquet")
    APPLICATION_VDN_PARQUET = (".parquet", "application/vdn.apache.parquet")
    APPLICATION_RAR_WINDOWS = (".rar", "application/x-rar-compressed")
    APPLICATION_RAR = (".rar", "application/vnd.rar")
    APPLICATION_ZIP = (".zip", "application/zip")
    APPLICATION_ZIP_WINDOWS = (".zip", "application/x-zip-compressed")

    @classmethod
    def get_extension(cls, content_type: str) -> str | None:
        """
        Retrieve the file extension associated with a given MIME type.

        Parameters
        ----------
        content_type : str
            The MIME type for which the corresponding file extension is to be retrieved.

        Returns
        -------
        str | None
            The file extension associated with the given MIME type, or None if the
            MIME type is not supported.

        Examples
        --------
        >>> FileExtensionEnum.get_extension("application/pdf")
        >>> '.pdf'
        >>> FileExtensionEnum.get_extension("application/unknown")
        >>> None

        """
        for member in cls:
            if member.value[1] == content_type:
                return member.value[0]
        return None

    @classmethod
    def search_file_extension_in_url(cls, url: str) -> str | None:
        """
        Search for the file extension in a given URL.

        Parameters
        ----------
        url : str
            The URL to search for the file extension.

        Returns
        -------
        str | None
            The file extension, or None if no match is found.

        Examples
        --------
        >>> FileExtensionEnum.search_file_extension_in_url("https://example.com/file.pdf")
        '.pdf'
        >>> FileExtensionEnum.search_file_extension_in_url("https://example.com/file.unknown")
        None

        """
        for member in cls:
            if re.search(r"\b" + re.escape(member.value[0]) + r"\b", url):
                return member.value[0]
        return None


class Utils:
    """
    A utility class for various helper functions.

    The class contains static methods that provide utility functions for
    common tasks, such as changing the format of datetime strings. The methods
    in this class are designed to be stateless and can be called without
    instantiating the class.

    Methods
    -------
    change_datetime_format(input_datetime_string: str, output_datetime_format: DateFormatEnum) -> str | None:
        Change the format of a given datetime string to a specified output format.

    """

    @staticmethod
    def change_datetime_format(
        input_datetime_string: str,
        output_datetime_format: DateFormatEnum,
    ) -> str | Any:
        """
        Change the format of a given datetime string to a specified output format.

        The method takes a datetime string and automatically detects its format, then converts it to the specified output format.
        If the input string cannot be parsed, it logs an error and returns None.

        Parameters
        ----------
        input_datetime_string : str
            datetime string that needs to be reformatted

        output_datetime_format : DateFormatEnum
            desired format for the output datetime string, following the strftime format codes.

        Returns
        -------
           str | None
               reformatted datetime string if successful, otherwise None

        Raises
        ------
        ValueError
            If the input datetime string cannot be parsed.

        """
        try:
            # Automatically detect the format of the input datetime string
            dt = parser.parse(input_datetime_string)
            logger.debug(f"The datetime string has been parsed successfully: {dt}")
            output_datetime_string = dt.strftime(output_datetime_format.value)
            logger.debug(f"The format is now changed to {output_datetime_format.value}")
            return output_datetime_string
        except ValueError as e:
            raise ValueError(f"Error during datetime formatting: {e}")

    @staticmethod
    def replace_special_characters(input_string: str) -> str:
        """
        Replace special characters and spaces in a string.

        The method replaces special characters and spaces in a string with underscores,
        collapsing multiple consecutive underscores into a single underscore. Finally, it lowercases all characters of the string and removes leading or
        trailing underscores.

        Parameters
        ----------
        input_string : str
            The input string from which special characters and spaces will be replaced.

        Returns
        -------
        str
            A new string with all special characters and spaces replaced by a single underscore
            where consecutive underscores occur.

        Examples
        --------
        >>> replace_special_characters("Hello, World! Welcome to Python @ 2023.")
        'hello_world_welcome_to_python_2023'

        >>> replace_special_characters("Special#Characters$Are%Fun!")
        'special_characters_are_fun'

        """
        # Replace any character that is not a word character or whitespace with underscore
        replaced = re.sub(r"[^\w\s]", "_", input_string)
        # Replace whitespace with underscore
        replaced = replaced.replace(" ", "_")
        # Collapse multiple consecutive underscores into a single underscore
        replaced = re.sub(r"_+", "_", replaced)
        # Remove leading and trailing underscores
        replaced = replaced.strip("_")
        # Lower case the string
        replaced = replaced.casefold()
        return replaced
