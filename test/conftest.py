"""
The module is used to define sharable pytest fixtures.

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


@pytest.fixture(scope="function")  # type: ignore
def example_source(request: pytest.FixtureRequest) -> td.Source:
    """Fixture to create an example source."""
    # Fetch the necessary values from the request object
    source_name = request.param.get("source_name", "example01")
    source_path = request.param.get(
        "source_path", pathlib.Path("technologydata", "datasources", "example01")
    )

    def load_example_source() -> td.Source:
        """Inner function to create the source object."""
        return td.Source(
            source_name,
            pathlib.Path(path_cwd, source_path),
        )

    # Call the inner function and return the result
    return load_example_source()


@pytest.fixture(scope="function")  # type: ignore
def example_technologies(request: pytest.FixtureRequest) -> td.Technologies:
    """Fixture to provide the example technologies."""
    # Fetch the necessary values from the request object
    technologies_name = request.param.get("technologies_name", "example01")
    technologies_path = request.param.get(
        "technologies_path", pathlib.Path("technologydata", "datasources", "example01")
    )

    def load_example_technologies() -> td.Technologies:
        """Inner function to create the source object."""
        return td.Technologies({technologies_name: technologies_path})

    # Call the inner function and return the result
    return load_example_technologies()
