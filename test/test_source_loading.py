"""Test different ways of loading and initializing the Source and Sources objects."""

from pathlib import Path

import pandas as pd
import pytest

import technologydata as td


@pytest.fixture()
def example_source():
    """Fixture to create an example source."""
    return td.Source(
        "example01",
        "technologydata/datasources/example01",
    )


def test_source_loading():
    """Simple approach for custom user sources."""
    name = "example01"
    fp = Path("technologydata/datasources/example01")

    # Check if the source can be loaded from a string or a path
    assert td.Source(name, fp)
    assert td.Source(name, str(fp))
    assert td.Source(name, fp.absolute())
    assert td.Source(name, str(fp.absolute()))


def test_source_initialization(example_source):
    """Test the initialization of the Source class."""
    assert example_source.name == "example01"
    # Path should be the provided path
    assert example_source.path == Path("technologydata/datasources/example01")
    # Check datatype
    assert isinstance(example_source.details, pd.DataFrame)
    # Ensure the DataFrame is not empty
    assert not example_source.details.empty
    # Check for expected features
    assert example_source.available_features == ["Technologies"]
    # Check whether packaged sources can also be loaded
    assert (
        td.Source("example01").name == example_source.name
        and td.Source("example01").path.absolute() == example_source.path.absolute()
    )


def test_sources_initialization(example_source):
    """Test different ways of initializing a Sources object."""
    # See if we can directly load all packaged sources directly / load from dict
    assert td.Sources(td.AVAILABLE_SOURCES)
    # Load a single prepackaged source by name
    assert td.Sources(["example01"])
    # Load single source not as list
    assert td.Sources("example01")
    # Convert a single source into a Sources object
    assert td.Sources(td.Source("example01"))
    # Combine multiple loaded sources into one
    assert td.Sources([example_source, example_source])
    # Mixed loading of named and loaded sources
    assert td.Sources(["example01", example_source])


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
    "url, expected",
    [
        (
            "https://ens.dk/media/3273/download",
            (
                "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                "2025-05-06 16:02:04",
                "200",
            ),
        ),
        ("https://ens.dk/media/1/download", None),
    ],
)
def test_is_wayback_snapshot_available(url, expected) -> None:
    """Check if the example source is available on the Internet Archive Wayback Machine."""
    output = Source.is_wayback_snapshot_available(url)
    assert output == expected
