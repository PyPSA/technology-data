# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

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

import technologydata

sys.path.append("./technology-data")
path_cwd = pathlib.Path.cwd()


def pytest_addoption(parser: pytest.Parser) -> None:
    """
    Add custom command-line options to pytest.

    This function adds a custom option `--run_webarchive` to the pytest command line.
    When this option is specified, it allows the execution of tests marked with the `webarchive` marker.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser object used to add custom options.

    """
    parser.addoption(
        "--run_webarchive",
        action="store_true",
        default=False,
        help="run the webarchive tests",
    )


def pytest_configure(config: pytest.Config) -> None:
    """
    Configure pytest with custom markers.

    This function adds a custom marker `webarchive` to pytest. Tests marked with this marker
    will only be run if the `--run_webarchive` option is specified.

    Parameters
    ----------
    config : pytest.Config
        The pytest configuration object.

    """
    config.addinivalue_line("markers", "webarchive: mark test as webarchive to run")


def pytest_collection_modifyitems(config: pytest.Config, items: pytest.Item) -> None:
    """
    Modify the test items collection based on command-line options.

    This function modifies the collection of test items. If the `--run_webarchive` option is not specified,
    it skips tests marked with the `webarchive` marker.
    By default, tests marked with `webarchive` will be skipped unless the option `--run_webarchive` is provided.
    This is because the CI was failing due to rate limits on anonymous capture webarchive requests using savepagenow.
    For further details, see https://github.com/open-energy-transition/technology-data/pull/41.

    Parameters
    ----------
    config : pytest.Config
        The pytest configuration object.
    items : list
        The list of test items collected by pytest.

    """
    if config.getoption("--run_webarchive"):
        return
    skip_webarchive = pytest.mark.skip(reason="need --run_webarchive option to run")
    for item in items:
        if "webarchive" in item.keywords:
            item.add_marker(skip_webarchive)


def create_source_from_params(params: dict[str, str]) -> technologydata.Source:
    """
    Create a Source object from a parameter dictionary with validation.

    This function takes a dictionary of parameters and validates that the required fields are present.
    If any required fields are missing, a ValueError is raised. If all required fields are present,
    a new Source object is created and returned.

    Parameters
    ----------
    params : dict[str, str]
        A dictionary containing the parameters for creating a Source object.
        Must include the keys "source_title" and "source_authors".
        Other keys are optional.

    Returns
    -------
    technologydata.Source
        A Source object initialized with the provided parameters.

    Raises
    ------
    ValueError
        If any of the required fields ("source_title", "source_authors") are missing from the params.

    """
    return technologydata.Source(
        title=params["source_title"],
        authors=params["source_authors"],
        url=params.get("source_url"),
        url_archive=params.get("source_url_archive"),
        url_date=params.get("source_url_date"),
        url_date_archive=params.get("source_url_date_archive"),
    )


@pytest.fixture(scope="function")  # type: ignore
def example_source(request: pytest.FixtureRequest) -> technologydata.Source:
    """Fixture to create an example source."""
    return create_source_from_params(request.param)


@pytest.fixture(scope="function")  # type: ignore
def example_source_collection(
    request: pytest.FixtureRequest,
) -> technologydata.SourceCollection:
    """
    Fixture to create an example SourceCollection from a list of parameter dicts.

    Each dict in the list must contain the following keys:
        - source_title
        - source_authors

    This fixture is compatible with pytest parametrize.
    """
    sources_params: list[dict[str, str]] = request.param
    sources = [create_source_from_params(params) for params in sources_params]
    return technologydata.SourceCollection(sources=sources)


@pytest.fixture(scope="function")  # type: ignore
def example_parameter(request: pytest.FixtureRequest) -> technologydata.Parameter:
    """Fixture to create an example parameter."""
    source_list = request.param.pop("parameter_sources", [])
    return technologydata.Parameter(
        magnitude=request.param.get("parameter_magnitude"),
        units=request.param.get("parameter_units"),
        carrier=request.param.get("parameter_carrier"),
        heating_value=request.param.get("parameter_heating_value"),
        provenance=request.param.get("parameter_provenance"),
        note=request.param.get("parameter_note"),
        sources=technologydata.SourceCollection(sources=source_list),
    )
