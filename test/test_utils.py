"""Test the utility methods."""

import pathlib
import sys
from typing import Any

import pytest

import technologydata as td

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "input_datetime_string, date_format, expected_date",
    [
        (
            "2025-05-20 14:45:00",
            td.DateFormatEnum.SOURCES_CSV,
            "2025-05-20 14:45:00",
        ),
        (
            "20250520144500",
            td.DateFormatEnum.SOURCES_CSV,
            "2025-05-20 14:45:00",
        ),
        ("2025-05-20 14:45:00", td.DateFormatEnum.WAYBACK, "20250520144500"),
        ("20250520144500", td.DateFormatEnum.WAYBACK, "20250520144500"),
        ("2025-05-20 14:45:00", td.DateFormatEnum.NONE, ""),
        ("invalid-date-string", td.DateFormatEnum.SOURCES_CSV, ValueError),
        ("2025/13/01", td.DateFormatEnum.SOURCES_CSV, ValueError),
    ],
)  # type: ignore
def test_change_datetime_format(
    input_datetime_string: str,
    date_format: td.DateFormatEnum,
    expected_date: str | Any,
) -> None:
    """Check if the datetime is correctly transformed to a new format."""
    if expected_date is ValueError:
        with pytest.raises(ValueError, match="Error during datetime formatting"):
            td.Utils.change_datetime_format(input_datetime_string, date_format)
    else:
        result = td.Utils.change_datetime_format(input_datetime_string, date_format)
        assert result == expected_date


@pytest.mark.parametrize(
    "input_string, expected_string",
    [
        (
            "Hello, World! Welcome to Python @ 2023.",
            "hello_world_welcome_to_python_2023",
        ),
        (
            "  Special#Characters$Are%Fun!  ",
            "special_characters_are_fun",
        ),
        (
            "!!!LeadingAndTrailing!!!",
            "leadingandtrailing",
        ),
    ],
)  # type: ignore
def test_replace_special_characters(
    input_string: str,
    expected_string: str,
) -> None:
    """Check if the special characters are removed from a string and the string is lowercased."""
    assert td.Utils.replace_special_characters(input_string) == expected_string


@pytest.mark.parametrize(
    "input_string, expected_string",
    [
        ("text/plain", ".txt"),
        ("text/html", ".html"),
        ("text/csv", ".csv"),
        ("text/xml", ".xml"),
        ("application/vnd.ms-excel", ".xls"),
        ("application/vnd.oasis.opendocument.spreadsheet", ".ods"),
        ("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", ".xlsx"),
        ("application/json", ".json"),
        ("application/xml", ".xml"),
        ("application/pdf", ".pdf"),
        ("application/parquet", ".parquet"),
        ("application/vdn.apache.parquet", ".parquet"),
        ("application/x-rar-compressed", ".rar"),
        ("application/vnd.rar", ".rar"),
        ("application/zip", ".zip"),
        ("application/x-zip-compressed", ".zip"),
    ],
)  # type: ignore
def test_get_extension(
    input_string: str,
    expected_string: str,
) -> None:
    """Check if the correct file extension is associated to a given MIME type."""
    assert td.FileExtensionEnum.get_extension(input_string) == expected_string


@pytest.mark.parametrize(
    "input_string, expected_string",
    [
        ("https://example.com/file.txt", ".txt"),
        ("https://example.com/file.html", ".html"),
        ("https://example.com/file.csv", ".csv"),
        ("https://example.com/file.xml", ".xml"),
        ("https://example.com/file.xls", ".xls"),
        ("https://example.com/file.ods", ".ods"),
        ("https://example.com/file.xlsx", ".xlsx"),
        ("https://example.com/file.json", ".json"),
        ("https://example.com/file.pdf", ".pdf"),
        ("https://example.com/file.parquet", ".parquet"),
        ("https://example.com/file.rar", ".rar"),
        ("https://example.com/file.zip", ".zip"),
        ("https://example.com/file.unknown", None),
    ],
)  # type: ignore
def test_search_file_extension_in_url(
    input_string: str,
    expected_string: str,
) -> None:
    """Check if the correct file extension is found in a given url."""
    assert (
        td.FileExtensionEnum.search_file_extension_in_url(input_string)
        == expected_string
    )
