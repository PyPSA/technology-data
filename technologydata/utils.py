"""Classes for utils methods."""

import logging
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
    An enumeration that maps file extensions to their corresponding MIME
    (https://en.wikipedia.org/wiki/MIME) types. This Enum provides a way to associate common
    file extensions with their respective MIME types, allowing for easy retrieval of file extensions
    based on content types.

    Members
    --------
    TEXT_PLAIN : tuple
        Represents the MIME type "text/plain" with the file extension ".txt".
    APPLICATION_PDF : tuple
        Represents the MIME type "application/pdf" with the file extension ".pdf".
    MS_EXCEL : tuple
        Represents the MIME type "application/vnd.ms-excel" with the file extension ".xls".
    OPENXML_EXCEL : tuple
        Represents the MIME type "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        with the file extension ".xlsx".
    PARQUET : tuple
        Represents the MIME type "application/parquet" with the file extension ".parquet".
    """

    TEXT_PLAIN = (".txt", "text/plain")
    APPLICATION_PDF = (".pdf", "application/pdf")
    MS_EXCEL = (".xls", "application/vnd.ms-excel")
    OPENXML_EXCEL = (
        ".xlsx",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    PARQUET = (".parquet", "application/parquet")

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
        '.pdf'
        >>> FileExtensionEnum.get_extension("application/unknown")
        None

        """
        for member in cls:
            if member.value[1] == content_type:
                return member.value[0]
        return None


class Utils:
    """
    A utility class for various helper functions.
    This class contains static methods that provide utility functions for
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
        Change the format of a given datetime string to a specified output format. This method takes a
        datetime string and automatically detects its format, then converts it to the specified output format.
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
