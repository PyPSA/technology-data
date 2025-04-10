import technologydata as td
import pandas as pd
import pytest
from pathlib import Path


@pytest.fixture()
def example_source():
    """Fixture to create an example source."""
    return td.Source(
        "example01",
        "technologydata/datasources/example01",
    )


def test_source_loading():
    name = "example01"
    fp = Path("technologydata/datasources/example01")

    # Check if the source can be loaded from a string or a path
    assert td.Source(name, fp)
    assert td.Source(name, str(fp))


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
    # Combine multiple loaded sources into one
    assert td.Sources([example_source, example_source])
    # Mixed loading of named and loaded sources
    assert td.Sources(["example01", example_source])
