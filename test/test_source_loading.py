"""Test different ways of loading and initializing the Source and Sources objects."""

import pathlib
import sys
from datetime import datetime
from typing import Any

import pandas as pd
import pytest

import technologydata as td

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


def test_source_loading() -> None:
    """Simple approach for custom user sources."""
    name = "example01"
    fp = pathlib.Path(path_cwd, "technologydata", "datasources", "example01")

    # Check if the source can be loaded from a string or a path
    assert td.Source(name, fp)
    assert td.Source(name, str(fp))
    assert td.Source(name, fp.absolute())
    assert td.Source(name, str(fp.absolute()))


@pytest.mark.parametrize(
    "example_source, expected_path, expected_features",
    [
        (
            {
                "source_name": "example01",
                "source_path": pathlib.Path(
                    "technologydata", "datasources", "example01"
                ),
            },
            pathlib.Path(path_cwd, "technologydata", "datasources", "example01"),
            ["Technologies"],
        ),
        (
            {
                "source_name": "example02",
                "source_path": pathlib.Path(
                    "technologydata", "datasources", "example02"
                ),
            },
            pathlib.Path(path_cwd, "technologydata", "datasources", "example02"),
            ["Technologies"],
        ),
    ],
    indirect=["example_source"],
)  # type: ignore
def test_source_initialization(
    example_source: td.Source, expected_path: pathlib.Path, expected_features: list[str]
) -> None:
    """Test the initialization of the Source class."""
    # Path should be the provided path
    assert example_source.path == expected_path
    # Check datatype
    assert isinstance(example_source.details, pd.DataFrame)
    # Ensure the DataFrame is not empty
    assert not example_source.details.empty
    # Check for expected features
    assert example_source.available_features == expected_features
    # Check whether packaged sources can also be loaded
    loaded_source = td.Source(example_source.name)
    assert (
        loaded_source.name == example_source.name
        and loaded_source.path is not None  # Ensure path is not None
        and loaded_source.path.absolute() == example_source.path.absolute()
    )


@pytest.mark.parametrize(
    "example_source",
    [
        {
            "source_name": "example01",
            "source_path": pathlib.Path("technologydata", "datasources", "example01"),
        },
        {
            "source_name": "example02",
            "source_path": pathlib.Path("technologydata", "datasources", "example02"),
        },
    ],
    indirect=True,
)  # type: ignore
def test_sources_initialization(example_source: td.Source) -> None:
    """Test different ways of initializing a Sources object."""
    # See if we can directly load all packaged sources directly / load from dict
    assert td.Sources(td.AVAILABLE_SOURCES)
    # Load a single prepackaged source by name
    assert td.Sources([example_source.name])
    # Load single source not as list
    assert td.Sources(example_source.name)
    # Convert a single source into a Sources object
    assert td.Sources(td.Source(example_source.name))
    # Combine multiple loaded sources into one
    assert td.Sources([example_source, example_source])
    # Mixed loading of named and loaded sources
    assert td.Sources([example_source.name, example_source])


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
    "url, timestamp, expected",
    [
        (
            "https://ens.dk/media/3273/download",
            "2025-05-06 16:02:04",
            (
                "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                "2025-05-06 16:02:04",
                "200",
            ),
        ),
        ("https://ens.dk/media/1/download", "2025-05-06 16:02:04", None),
    ],
)  # type: ignore
def test_is_wayback_snapshot_available(
    url: str, timestamp: str, expected: tuple[str, str, str] | None
) -> None:
    """Check if the example source is available on the Internet Archive Wayback Machine."""
    output = td.Source.is_wayback_snapshot_available(url, timestamp)
    assert output == expected


@pytest.mark.parametrize(
    "example_source",
    [
        {
            "source_name": "example03",
            "source_path": pathlib.Path("technologydata", "datasources", "example03"),
        }
    ],
    indirect=True,
)  # type: ignore
def test_download_file_from_wayback(example_source: td.Source) -> None:
    """Check if the example source is downloaded from the Internet Archive Wayback Machine."""
    storage_paths = example_source.download_file_from_wayback()
    # Check if storage_paths is not None and is a list
    assert storage_paths is not None, (
        "Expected a valid storage path list, but got None."
    )
    assert isinstance(storage_paths, list), "Expected storage_paths to be a list."
    for storage_path in storage_paths:
        # Check if each storage_path is not None
        assert storage_path is not None, "Expected a valid storage path, but got None."
        assert storage_path.is_file(), (
            f"Expected {storage_path} to be a file, but it does not exist."
        )
        # Delete the downloaded file
        storage_path.unlink(missing_ok=True)


def test_store_snapshot_on_wayback() -> None:
    """Check if a given url is correctly stored as a snapshot on Internet Archive Wayback Machine."""
    url_to_archive = "https://openenergytransition.org/outputs.html"
    archived_info = td.Source.store_snapshot_on_wayback(url_to_archive)

    # Check if archived_info is None
    assert archived_info is not None, "archived_info should not be None"

    archived_url, new_capture, output_timestamp = archived_info

    assert "https://web.archive.org/web/" in archived_url
    assert url_to_archive in archived_url
    assert isinstance(new_capture, bool)

    assert output_timestamp is not None, "output_timestamp should not be None"
    try:
        datetime.strptime(output_timestamp, td.DateFormatEnum.SOURCES_CSV)
    except ValueError:
        pytest.fail("Valid date-time string did not match the format")


@pytest.mark.parametrize(
    "example_source",
    [
        {
            "source_name": "example03",
            "source_path": pathlib.Path("technologydata", "datasources", "example03"),
        },
        {
            "source_name": "example04",
            "source_path": pathlib.Path("technologydata", "datasources", "example04"),
        },
    ],
    indirect=["example_source"],
)  # type: ignore
def test_ensure_snapshot(example_source: td.Source) -> None:
    """Check if a given sources.csv file contains the fields url_archive_date and url_archived."""
    # Construct the file path
    file_name = "sources_modified.csv"

    assert example_source.path is not None, "path should not be None"
    file_path = pathlib.Path(path_cwd, example_source.path, file_name)

    # Ensure the snapshot is created
    example_source.ensure_snapshot(file_name)

    # Read the output DataFrame
    try:
        output_df = pd.read_csv(file_path)
    except Exception as e:
        pytest.fail(f"Failed to read CSV file: {e}")

    # Assert that the archived URL is now filled
    assert not output_df["url_archived"].isna().all(), (
        "Some values in url_archived are null"
    )

    # Assert that the timestamp is now filled
    assert not output_df["url_archive_date"].isna().all(), (
        "Some values in url_archive_date are null"
    )

    # Clean up by removing the file if it exists
    file_path.unlink(missing_ok=True)
