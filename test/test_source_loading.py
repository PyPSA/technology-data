"""Test different ways of loading and initializing the Source and Sources objects."""

import pathlib
import sys

import pandas as pd
import pytest

import technologydata as td
from technologydata import Source

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


def test_source_loading():
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
)
def test_source_initialization(example_source, expected_path, expected_features):
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
    assert (
        td.Source(example_source.name).name == example_source.name
        and td.Source(example_source.name).path.absolute()
        == example_source.path.absolute()
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
)
def test_sources_initialization(example_source):
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
)
def test_change_datetime_format(
    input_date, input_format, output_format, expected_date
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
)
def test_is_wayback_snapshot_available(url, timestamp, expected) -> None:
    """Check if the example source is available on the Internet Archive Wayback Machine."""
    output = Source.is_wayback_snapshot_available(url, timestamp)
    print(output)
    assert output == expected
