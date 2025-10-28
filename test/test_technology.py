# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the Technology class."""

import pathlib

import technologydata

path_cwd = pathlib.Path.cwd()


class TestTechnology:
    """Test suite for the Technology class in the technologydata module."""

    def test_technology_creation(self) -> None:
        """Test the creation of a Technology instance."""
        tech = technologydata.Technology(
            name="Solar photovoltaics",
            detailed_technology="Si-HC",
            case="example-scenario",
            region="DEU",
            year=2022,
            parameters={
                "investment": technologydata.Parameter(
                    magnitude=500,
                    units="EUR_2022/kW",
                    provenance="Industry report",
                    note="Average overnight cost",
                ),
                "lifetime": technologydata.Parameter(
                    magnitude=25,
                    units="years",
                    provenance="Literature",
                ),
            },
        )
        assert tech.name == "Solar photovoltaics"
        assert tech.detailed_technology == "Si-HC"
        assert tech.case == "example-scenario"
        assert tech.region == "DEU"
        assert tech.year == 2022
        assert len(tech.parameters) == 2
        assert "investment" in tech.parameters
        assert "lifetime" in tech.parameters

    def test_technology_to_currency(self) -> None:
        """Test currency conversion for all parameters in a Technology."""
        tech = technologydata.Technology(
            name="Solar photovoltaics",
            detailed_technology="Si-HC",
            case="example-scenario",
            region="DEU",
            year=2022,
            parameters={
                "investment": technologydata.Parameter(
                    magnitude=500,
                    units="EUR_2022/kW",
                    provenance="Industry report",
                    note="Average overnight cost",
                ),
                "fixed_om": technologydata.Parameter(
                    magnitude=10,
                    units="EUR_2022/kW/year",
                    provenance="Literature",
                ),
                "lifetime": technologydata.Parameter(
                    magnitude=25,
                    units="years",
                    provenance="Literature",
                ),
            },
        )

        # Convert to USD_2020
        converted = tech.to_currency("USD_2020")

        # Check that a new Technology object was returned
        assert isinstance(converted, technologydata.Technology)
        assert converted is not tech  # Should be a different object

        # Check that basic attributes are preserved
        assert converted.name == tech.name
        assert converted.detailed_technology == tech.detailed_technology
        assert converted.case == tech.case
        assert converted.region == tech.region
        assert converted.year == tech.year

        # Check that parameters with currency units are converted
        assert "USD_2020" in str(converted.parameters["investment"].units)
        assert "EUR_2022" not in str(converted.parameters["investment"].units)
        assert "USD_2020" in str(converted.parameters["fixed_om"].units)
        assert "EUR_2022" not in str(converted.parameters["fixed_om"].units)

        # Check that parameters without currency units remain unchanged
        assert (
            converted.parameters["lifetime"].units == tech.parameters["lifetime"].units
        )
        assert (
            converted.parameters["lifetime"].magnitude
            == tech.parameters["lifetime"].magnitude
        )

        # Check that magnitude changed for currency parameters
        assert (
            converted.parameters["investment"].magnitude
            != tech.parameters["investment"].magnitude
        )
        assert (
            converted.parameters["fixed_om"].magnitude
            != tech.parameters["fixed_om"].magnitude
        )

    def test_technology_to_currency_with_overwrite_country(self) -> None:
        """Test currency conversion with a different country for inflation adjustment."""
        tech = technologydata.Technology(
            name="Wind turbine",
            detailed_technology="Onshore",
            case="base-case",
            region="DEU",
            year=2022,
            parameters={
                "investment": technologydata.Parameter(
                    magnitude=1200,
                    units="EUR_2020/kW",
                    provenance="Industry data",
                ),
            },
        )

        # Convert using a different country (USA) for inflation adjustment
        converted = tech.to_currency("USD_2023", overwrite_country="USA")

        assert isinstance(converted, technologydata.Technology)
        assert "USD_2023" in str(converted.parameters["investment"].units)
        assert (
            converted.parameters["investment"].magnitude
            != tech.parameters["investment"].magnitude
        )

    def test_technology_to_currency_with_source(self) -> None:
        """Test currency conversion with different inflation data sources."""
        tech = technologydata.Technology(
            name="Battery storage",
            detailed_technology="Li-ion",
            case="scenario-1",
            region="USA",
            year=2022,
            parameters={
                "investment": technologydata.Parameter(
                    magnitude=300,
                    units="USD_2020/kWh",
                ),
            },
        )

        # Convert using worldbank source
        converted_wb = tech.to_currency("EUR_2022", source="worldbank")
        assert isinstance(converted_wb, technologydata.Technology)
        assert "EUR_2022" in str(converted_wb.parameters["investment"].units)

        # Convert using IMF source
        converted_imf = tech.to_currency("EUR_2022", source="imf")
        assert isinstance(converted_imf, technologydata.Technology)
        assert "EUR_2022" in str(converted_imf.parameters["investment"].units)

    def test_technology_to_currency_preserves_other_attributes(self) -> None:
        """Test that currency conversion preserves other parameter attributes."""
        tech = technologydata.Technology(
            name="Solar photovoltaics",
            detailed_technology="Si-HC",
            case="example-scenario",
            region="DEU",
            year=2022,
            parameters={
                "investment": technologydata.Parameter(
                    magnitude=500,
                    units="EUR_2022/kW",
                    carrier="electricity",
                    provenance="Industry report",
                    note="Average overnight cost",
                    sources=technologydata.SourceCollection(
                        sources=[
                            technologydata.Source(
                                title="Test Source",
                                authors="Test Authors",
                            )
                        ]
                    ),
                ),
            },
        )

        converted = tech.to_currency("USD_2020")

        # Check that non-currency attributes are preserved
        assert (
            converted.parameters["investment"].carrier
            == tech.parameters["investment"].carrier
        )
        assert (
            converted.parameters["investment"].provenance
            == tech.parameters["investment"].provenance
        )
        assert (
            converted.parameters["investment"].note
            == tech.parameters["investment"].note
        )
        assert len(converted.parameters["investment"].sources.sources) == len(
            tech.parameters["investment"].sources.sources
        )

    def test_technology_to_currency_empty_parameters(self) -> None:
        """Test currency conversion on a Technology with no parameters."""
        tech = technologydata.Technology(
            name="Solar photovoltaics",
            detailed_technology="Si-HC",
            case="example-scenario",
            region="DEU",
            year=2022,
            parameters={},
        )

        # Should handle empty parameters gracefully
        converted = tech.to_currency("USD_2020")

        assert isinstance(converted, technologydata.Technology)
        assert len(converted.parameters) == 0

    def test_technology_to_currency_mixed_currencies(self) -> None:
        """Test currency conversion when parameters have different currencies."""
        tech = technologydata.Technology(
            name="Hybrid system",
            detailed_technology="Solar+Wind",
            case="scenario-1",
            region="DEU",
            year=2022,
            parameters={
                "investment_solar": technologydata.Parameter(
                    magnitude=500,
                    units="EUR_2022/kW",
                ),
                "investment_wind": technologydata.Parameter(
                    magnitude=1000,
                    units="USD_2020/kW",
                ),
                "capacity": technologydata.Parameter(
                    magnitude=100,
                    units="MW",
                ),
            },
        )

        # Convert all to GBP_2021
        converted = tech.to_currency("GBP_2021")

        # All currency parameters should be converted to GBP_2021
        assert "GBP_2021" in str(converted.parameters["investment_solar"].units)
        assert "EUR_2022" not in str(converted.parameters["investment_solar"].units)
        assert "GBP_2021" in str(converted.parameters["investment_wind"].units)
        assert "USD_2020" not in str(converted.parameters["investment_wind"].units)

        # Non-currency parameter should remain unchanged
        assert (
            converted.parameters["capacity"].units == tech.parameters["capacity"].units
        )
        assert (
            converted.parameters["capacity"].magnitude
            == tech.parameters["capacity"].magnitude
        )
