import pathlib
import sys

import pytest

import technologydata as td

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


@pytest.fixture()
def example_source():
    """Fixture to create an example source."""
    return td.Source(
        "example01",
        pathlib.Path(path_cwd, "technologydata", "datasources", "example01"),
    )
