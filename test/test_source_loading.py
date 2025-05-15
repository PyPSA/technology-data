"""Test different ways of loading and initializing the Source and Sources objects."""

import pathlib
import sys

import pandas as pd
import pytest

from technologydata import AVAILABLE_SOURCES, Source, Sources

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


def test_source_loading() -> None:
    """Simple approach for custom user sources."""
    name = "example01"
    fp = pathlib.Path(path_cwd, "technologydata", "datasources", "example01")

    # Check if the source can be loaded from a string or a path
    assert Source(name, fp)
    assert Source(name, str(fp))
    assert Source(name, fp.absolute())
    assert Source(name, str(fp.absolute()))


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
    example_source: Source, expected_path: pathlib.Path, expected_features: list[str]
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
    loaded_source = Source(example_source.name)
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
def test_sources_initialization(example_source: Source) -> None:
    """Test different ways of initializing a Sources object."""
    # See if we can directly load all packaged sources directly / load from dict
    assert Sources(AVAILABLE_SOURCES)
    # Load a single prepackaged source by name
    assert Sources([example_source.name])
    # Load single source not as list
    assert Sources(example_source.name)
    # Convert a single source into a Sources object
    assert Sources(Source(example_source.name))
    # Combine multiple loaded sources into one
    assert Sources([example_source, example_source])
    # Mixed loading of named and loaded sources
    assert Sources([example_source.name, example_source])


@pytest.mark.parametrize(
    "input_date, input_format, output_format, expected_date",
    [
        (
            "20250507105201",
            "%Y%m%d%H%M%S",
            "%Y-%m-%d %H:%M:%S",
            "2025-05-07 10:52:01",
        ),
        (
            "2025-05-07 10:52:01",
            "%Y-%m-%d %H:%M:%S",
            "%Y%m%d%H%M%S",
            "20250507105201",
        ),
        (
            "2025-05-07 10:52:01",
            "%Y-%m-%d",
            "%Y%m%d%H%M%S",
            None,
        ),
    ],
)  # type: ignore
def test_change_datetime_format(
    input_date: str, input_format: str, output_format: str, expected_date: str
) -> None:
    """Check if the datetime is correctly transformed to a new format."""
    output_date = Source.change_datetime_format(input_date, input_format, output_format)
    assert output_date == expected_date


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
    output = Source.is_wayback_snapshot_available(url, timestamp)
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
def test_download_file_from_wayback(example_source: Source) -> None:
    """Check if the example source is downloaded from the Internet Archive Wayback Machine."""
    storage_path = example_source.download_file_from_wayback()

    # Check if storage_path is not None
    assert storage_path is not None, "Expected a valid storage path, but got None."

    assert storage_path.is_file()

    # Delete the downloaded file
    pathlib.Path(storage_path).unlink(missing_ok=True)


def test_store_snapshot_on_wayback() -> None:
    """Check if a given url is correctly stored as a snapshot on Internet Archive Wayback Machine."""
    url_to_archive = "https://openenergytransition.org/outputs.html"
    archived_info = Source.store_snapshot_on_wayback(url_to_archive)

    # Check if archived_info is None
    assert archived_info is not None, "archived_info should not be None"

    archived_url, new_capture, output_timestamp = archived_info
    assert (
        archived_url
        == "https://web.archive.org/web/20250513133237/https://openenergytransition.org/outputs.html"
    )
    assert output_timestamp == "2025-05-13 13:32:37"


@pytest.mark.parametrize(
    "example_source",
    [
        {
            "source_name": "example03",
            "source_path": pathlib.Path(
                "technologydata", "datasources", "example03"
            ),
        },
        {
            "source_name": "example04",
            "source_path": pathlib.Path(
                "technologydata", "datasources", "example04"
            ),
        },
    ],
    indirect=["example_source"],
)  # type: ignore
def test_ensure_snapshot(
    example_source: Source
) -> None:
    """Check if a given sources.csv file contains the fields url_date and url_archived."""
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
        f"Some values in url_archived are null"
    )

    # Assert that the timestamp is now filled
    assert not output_df["url_date"].isna().all(), (
        f"Some values in url_date are null"
    )

    # Clean up by removing the file if it exists
    # file_path.unlink(missing_ok=True)
