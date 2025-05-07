import pathlib
import sys

import pytest

import technologydata as td

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


@pytest.fixture(scope="function")
def example_source(request):
    """Fixture to create an example source."""
    source_name = request.param.get("source_name", "example01")
    source_path = request.param.get(
        "source_path", pathlib.Path("technologydata", "datasources", "example01")
    )

    def load_example_source():
        """Inner function to create the source."""
        return td.Source(
            source_name,
            pathlib.Path(path_cwd, source_path),
        )

    # Call the inner function and return the result
    return load_example_source()
