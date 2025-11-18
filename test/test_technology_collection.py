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

    @pytest.mark.parametrize(
        "output_schema",
        [True, False],
    )  # type: ignore
    def test_to_json(self, output_schema: bool) -> None:
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
        technology_collection.to_json(
            pathlib.Path(output_file), output_schema=output_schema
        )
        assert output_file.is_file()
        if output_schema:
            assert schema_file.is_file()
        else:
            assert not schema_file.is_file()
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
        technology_collection.to_json(output_file, output_schema=True)

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

    def test_fit_linear_growth(self) -> None:
        """Test TechnologyCollection.fit with LinearGrowth model."""
        tech = technologydata.Technology(
            name="Amazing technology",
            detailed_technology="",
            region="",
            case="",
            year=2020,
            parameters={
                "total units": technologydata.Parameter(magnitude=2020),
            },
        )

        tc = technologydata.TechnologyCollection(
            technologies=[
                tech,
                tech.model_copy(
                    deep=True,
                    update={
                        "year": 2030,
                        "parameters": {
                            "total units": technologydata.Parameter(magnitude=2030),
                        },
                    },
                ),
                tech.model_copy(
                    deep=True,
                    update={
                        "year": 2040,
                        "parameters": {
                            "total units": technologydata.Parameter(magnitude=2040),
                        },
                    },
                ),
            ]
        )

        # Fit 'total units' parameter with LinearGrowth
        from technologydata.technologies.growth_models import LinearGrowth

        model = LinearGrowth()
        fitted = tc.fit("total units", model, p0={"x0": 0, "m": 1, "A": 0})
        assert isinstance(fitted, LinearGrowth)
        print(fitted)
        assert pytest.approx(fitted.x0) == 0
        assert pytest.approx(fitted.m) == 1
        assert pytest.approx(fitted.A) == 0

    def test_project_linear_growth(self) -> None:
        """Test TechnologyCollection.project with LinearGrowth model."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )
        tc = technologydata.TechnologyCollection.from_json(input_file)
        from technologydata.technologies.growth_models import LinearGrowth

        projected_tc = tc.project(
            to_years=[2030],
            parameters={"capacity": LinearGrowth(x0=2020)},
        )
        assert isinstance(projected_tc, technologydata.TechnologyCollection)
        assert projected_tc.technologies[0].year == 2030
        assert "capacity" in projected_tc.technologies[0].parameters
        assert isinstance(
            projected_tc.technologies[0].parameters["capacity"].magnitude, float
        )

        # non-projected parameters should not be present
        assert "investment" not in projected_tc.technologies[0].parameters

    def test_project_other_parameter_options(self) -> None:
        """Test projection of parameters using 'mean', 'closest', and 'NaN' options."""
        tech = technologydata.Technology(
            name="Amazing technology",
            detailed_technology="",
            region="",
            case="",
            year=2020,
            parameters={
                "total units": technologydata.Parameter(magnitude=2000),
            },
        )

        tc = technologydata.TechnologyCollection(
            technologies=[
                tech,
                tech.model_copy(
                    deep=True,
                    update={
                        "year": 2030,
                        "parameters": {
                            "total units": technologydata.Parameter(magnitude=3000),
                        },
                    },
                ),
            ]
        )

        ptc = tc.project(
            to_years=[2025],
            parameters={
                "total units": "mean",
            },
        )

        assert (
            pytest.approx(
                (
                    tc.technologies[0].parameters["total units"].magnitude
                    + tc.technologies[1].parameters["total units"].magnitude
                )
                / 2,
            )
            == ptc.technologies[0].parameters["total units"].magnitude
        )

        ptc = tc.project(
            to_years=[2025],
            parameters={
                "total units": "NaN",
            },
        )

        assert pandas.isna(ptc.technologies[0].parameters["total units"].magnitude)

        # "closest" currently raises NotImplementedError
        with pytest.raises(NotImplementedError):
            _ = tc.project(
                to_years=[2025],
                parameters={
                    "total units": "closest",
                },
            )

    def test_to_currency(self) -> None:
        """Test currency conversion for all technologies in the collection."""
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

        # Convert to USD_2020
        converted = technology_collection.to_currency("USD_2020")

        # Check that a new TechnologyCollection object was returned
        assert isinstance(converted, technologydata.TechnologyCollection)
        assert converted is not technology_collection

        # Check that the number of technologies is the same
        assert len(converted) == len(technology_collection)

        # Check that all technologies have been converted
        for tech in converted:
            for param_name, param in tech.parameters.items():
                if param.units and "EUR_2022" in param.units:
                    # Should not happen - all EUR_2022 should be converted
                    assert False, f"Parameter {param_name} still has EUR_2022"

    def test_to_currency_with_overwrite_country(self) -> None:
        """Test currency conversion with a different country for inflation adjustment."""
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

        # Convert using a different country (USA) for inflation adjustment
        converted = technology_collection.to_currency(
            "USD_2023", overwrite_country="USA"
        )

        assert isinstance(converted, technologydata.TechnologyCollection)
        assert len(converted) == len(technology_collection)

        # Check that currencies have been converted
        for tech in converted:
            for param in tech.parameters.values():
                if param.units and "USD" in param.units:
                    assert "USD_2023" in param.units

    def test_to_currency_with_source(self) -> None:
        """Test currency conversion with different inflation data sources."""
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

        # Convert using worldbank source
        converted_wb = technology_collection.to_currency(
            "USD_2022", source="international_monetary_fund"
        )
        assert isinstance(converted_wb, technologydata.TechnologyCollection)
        assert len(converted_wb) == len(technology_collection)

        # Convert using IMF source
        converted_imf = technology_collection.to_currency("USD_2022", source="imf")
        assert isinstance(converted_imf, technologydata.TechnologyCollection)
        assert len(converted_imf) == len(technology_collection)

    def test_to_currency_preserves_technology_attributes(self) -> None:
        """Test that currency conversion preserves technology attributes."""
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

        converted = technology_collection.to_currency("USD_2020")

        # Check that technology attributes are preserved
        for orig_tech, conv_tech in zip(
            technology_collection.technologies, converted.technologies
        ):
            assert orig_tech.name == conv_tech.name
            assert orig_tech.detailed_technology == conv_tech.detailed_technology
            assert orig_tech.case == conv_tech.case
            assert orig_tech.region == conv_tech.region
            assert orig_tech.year == conv_tech.year
