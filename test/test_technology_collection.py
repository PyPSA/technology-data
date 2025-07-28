# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the TechnologyCollection class."""

import pathlib

import pandas
import pytest

import technologydata

path_cwd = pathlib.Path.cwd()


def test_to_csv() -> None:
    """Check if the example technology collection is exported to CSV."""
    input_file = pathlib.Path(
        path_cwd,
        "test",
        "test_data",
        "solar_photovoltaics_example_03",
        "technologies.json",
    )
    technology_collection = technologydata.TechnologyCollection.from_json(input_file)
    output_file = pathlib.Path(path_cwd, "technologies.csv")
    technology_collection.to_csv(path_or_buf=output_file)
    assert output_file.is_file()
    output_file.unlink(missing_ok=True)


def test_to_dataframe() -> None:
    """Check if the example technology collection is exported to pandas dataframe."""
    input_file = pathlib.Path(
        path_cwd,
        "test",
        "test_data",
        "solar_photovoltaics_example_03",
        "technologies.json",
    )
    technology_collection = technologydata.TechnologyCollection.from_json(input_file)
    assert isinstance(technology_collection.to_dataframe(), pandas.DataFrame)


def test_to_json() -> None:
    """Check if the example technology collection is exported to JSON."""
    input_file = pathlib.Path(
        path_cwd,
        "test",
        "test_data",
        "solar_photovoltaics_example_03",
        "technologies.json",
    )
    technology_collection = technologydata.TechnologyCollection.from_json(input_file)
    output_file = pathlib.Path(path_cwd, "technologies.json")
    schema_file = pathlib.Path(path_cwd, "technologies.schema.json")
    technology_collection.to_json(pathlib.Path(output_file))
    assert output_file.is_file()
    assert schema_file.is_file()
    output_file.unlink(missing_ok=True)
    schema_file.unlink(missing_ok=True)


def test_from_json() -> None:
    """Check if the example technology collection is imported from JSON."""
    input_file = pathlib.Path(
        path_cwd,
        "test",
        "test_data",
        "solar_photovoltaics_example_03",
        "technologies.json",
    )
    technology_collection = technologydata.TechnologyCollection.from_json(input_file)
    assert isinstance(technology_collection, technologydata.TechnologyCollection)
    assert len(technology_collection) == 2


@pytest.mark.parametrize(
    "name, region, year, case, detailed_technology",
    [
        ["Solar photovoltaics", "DEU", 2022, "example-scenario", "Si-HC"],
        ["Solar photovoltaics", "DEU", 2022, "example-project", "Si-HC"],
    ],
)  # type: ignore
def test_get(
    name: str, region: str, year: int, case: str, detailed_technology: str
) -> None:
    """Check if the example technology collection is filtered as requested."""
    input_file = pathlib.Path(
        path_cwd,
        "test",
        "test_data",
        "solar_photovoltaics_example_03",
        "technologies.json",
    )
    technologies_collection = technologydata.TechnologyCollection.from_json(input_file)
    result = technologies_collection.get(
        name=name,
        region=region,
        year=year,
        case=case,
        detailed_technology=detailed_technology,
    )
    assert isinstance(result, technologydata.TechnologyCollection)
    assert len(result.technologies) == 1
