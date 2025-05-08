"""
The module is used to define pytest fixtures that can be shared across
multiple test files in the project.
Fixtures are a way to provide a fixed baseline upon which tests can
rely. They allow for setup code to be reused and can help manage
resources needed for tests, such as database connections, test data,
or configuration settings.
By placing fixtures in this file, they become accessible to all test
modules in the same directory and subdirectories, promoting code
reusability and organization.
"""

import pathlib
import sys

import pytest

import technologydata as td

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


@pytest.fixture(scope="function")
def example_source(request) -> td.Source:
    """Fixture to create an example source."""
    source_name = request.param.get("source_name", "example01")
    source_path = request.param.get(
        "source_path", pathlib.Path("technologydata", "datasources", "example01")
    )

    def load_example_source() -> td.Source:
        """Inner function to create the source."""
        return td.Source(
            source_name,
            pathlib.Path(path_cwd, source_path),
        )

    # Call the inner function and return the result
    return load_example_source()
