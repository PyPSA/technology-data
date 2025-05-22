"""Classes for utils methods."""
from enum import Enum
from typing import Any, Tuple, Dict, Type
from dateutil import parser
import logging

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


class Utils:

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