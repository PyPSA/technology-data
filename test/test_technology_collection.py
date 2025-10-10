# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the TechnologyCollection class."""

import pathlib

import pandas
import pytest

import technologydata

path_cwd = pathlib.Path.cwd()


class TestTechnologyCollection:
    """Test suite for the TechnologyCollection class in the technologydata module."""

    def test_to_csv(self) -> None:
        """Check if the example technology collection is exported to CSV."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technology_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        output_file = pathlib.Path(path_cwd, "technologies.csv")
        technology_collection.to_csv(path_or_buf=output_file)
        assert output_file.is_file()
        output_file.unlink(missing_ok=True)

    def test_to_dataframe(self) -> None:
        """Check if the example technology collection is exported to pandas dataframe."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technology_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        assert isinstance(technology_collection.to_dataframe(), pandas.DataFrame)

    def test_to_json(self) -> None:
        """Check if to_json works as expected."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technology_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        output_file = pathlib.Path(path_cwd, "technologies.json")
        schema_file = pathlib.Path(path_cwd, "technologies.schema.json")
        technology_collection.to_json(pathlib.Path(output_file))
        assert output_file.is_file()
        assert schema_file.is_file()
        output_file.unlink(missing_ok=True)
        schema_file.unlink(missing_ok=True)

    def test_from_json(self) -> None:
        """Check if from_json works as expected."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technology_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        assert isinstance(technology_collection, technologydata.TechnologyCollection)
        assert len(technology_collection) == 2

    def test_from_json_to_json(self) -> None:
        """Check whether reading with from_json and exporting with to_json yields the same file."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technology_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        output_file = pathlib.Path("to_json_test.json")
        schema_file = pathlib.Path(path_cwd, "to_json_test.schema.json")
        technology_collection.to_json(output_file)

        # Read files and strip trailing whitespace/newlines before comparing
        with open(input_file) as f1, open(output_file) as f2:
            assert f1.read().rstrip() == f2.read().rstrip(), "Files are not identical"
        assert output_file.is_file()
        assert schema_file.is_file()
        output_file.unlink(missing_ok=True)
        schema_file.unlink(missing_ok=True)

    @pytest.mark.parametrize(
        "name, region, year, case, detailed_technology",
        [
            ["Solar photovoltaics", "DEU", 2022, "example-scenario", "Si-HC"],
            ["Solar photovoltaics", "DEU", 2022, "example-project", "Si-HC"],
        ],
    )  # type: ignore
    def test_get(
        self, name: str, region: str, year: int, case: str, detailed_technology: str
    ) -> None:
        """Check if the example technology collection is filtered as requested."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        technologies_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        result = technologies_collection.get(
            name=name,
            region=region,
            year=year,
            case=case,
            detailed_technology=detailed_technology,
        )
        assert isinstance(result, technologydata.TechnologyCollection)
        assert len(result.technologies) == 1
