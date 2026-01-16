# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the functions of the DEA energy storage parser."""

import pathlib
import typing

import pandas
import pytest

from technologydata.package_data.dea_energy_storage.dea_energy_storage import (
    build_technology_collection,
    clean_est_string,
    clean_parameter_string,
    clean_technology_string,
    drop_invalid_rows,
    extract_year,
    filter_parameters,
    format_val_number,
    standardize_units,
)


class TestDEAEnergyStorage:
    """Test suite for the functions of the DEA energy storage parser."""

    @pytest.mark.parametrize(
        "input_string, expected_string",
        [
            ("- Discharge efficiency [%]", "discharge efficiency"),
            ("Construction time [years]", "construction time"),
            ("- Another example [with brackets]", "another example"),
            ("No brackets or hyphen here", "no brackets or hyphen here"),
            ("- Another example [with brackets)", "another example"),
        ],
    )  # type: ignore
    def test_clean_parameter_string(
        self, input_string: str, expected_string: str
    ) -> None:
        """Check if the clean_parameter_string works as expected."""
        assert clean_parameter_string(input_string) == expected_string

    def test_drop_invalid_rows(self) -> None:
        """Check if the drop_invalid_rows works as expected."""
        input_dataframe = pandas.DataFrame(
            {
                "Technology": ["AI", "", "ML", "new_tech", "Tech", "Tech_1", "Tech_2"],
                "par": ["p1", "p2", "", "p4", None, "p5", "p6"],
                "val": ["<1", "   ", None, "abc", "456", "456,1", "1,1x10^3"],
                "year": [
                    "almost 2020",
                    "2021",
                    "2020",
                    "2022",
                    "",
                    "nearly 2020",
                    "2024",
                ],
            }
        )
        expected_dataframe = pandas.DataFrame(
            {
                "Technology": ["Tech_1", "Tech_2"],
                "par": ["p5", "p6"],
                "val": ["456,1", "1,1x10^3"],
                "year": ["nearly 2020", "2024"],
            }
        )
        output_dataframe = drop_invalid_rows(input_dataframe).reset_index(drop=True)
        comparison_df = output_dataframe.compare(expected_dataframe)
        assert comparison_df.empty

    @pytest.mark.parametrize(
        "input_string, expected_string",
        [
            ("192 string", "string"),
            (" STRING 123 ", "string 123"),
            (" string 123a ", "string 123a"),
            (" 123b string ", "string"),
            (" 123c ", ""),
        ],
    )  # type: ignore
    def test_clean_technology_string(
        self, input_string: str, expected_string: str
    ) -> None:
        """Check if clean_technology_string works as expected."""
        result = clean_technology_string(input_string)
        assert result == expected_string
        assert isinstance(result, str)

    @pytest.mark.parametrize(
        "input_year, expected_year",
        [
            ("some 1999", 1999),
            ("maybe 1", 1),
            ("1", 1),
        ],
    )  # type: ignore
    def test_extract_year(self, input_year: str, expected_year: int) -> None:
        """Check if extract_year works as expected."""
        result = extract_year(input_year)
        assert result == expected_year
        assert isinstance(result, int)

    @pytest.mark.parametrize(
        "input_number, expected_number",
        [
            ("1,1", 1.1),
            ("1", 1.0),
            ("1.31×10-23", 1.31e-23),
        ],
    )  # type: ignore
    def test_format_val_number(
        self, input_number: str, expected_number: int | None | typing.Any
    ) -> None:
        """Check if format_val_number works as expected, including exception handling."""
        result = format_val_number(input_number, 2)
        assert isinstance(result, float)
        assert result == expected_number

    @pytest.mark.parametrize(
        "input_string, expected_string",
        [
            ("ctrl", "control"),
            ("Upper", "upper"),
            (" Yeah ", "yeah"),
        ],
    )  # type: ignore
    def test_clean_est_string(self, input_string: str, expected_string: str) -> None:
        """Check if clean_est_string works as expected."""
        result = clean_est_string(input_string)
        assert isinstance(result, str)
        assert result == expected_string

    def test_standardize_units(self) -> None:
        """Check if standardize_units works as expected."""
        input_dataframe = pandas.DataFrame(
            {
                "par": [
                    "energy storage capacity for one unit",
                    "typical temperature difference in storage",
                    "fixed o&m",
                    "lifetime in total number of cycles",
                    "cycle life",
                    "p1",
                    "p2",
                    "p3",
                    "p4",
                    "p5",
                ],
                "unit": [
                    "",
                    "",
                    "",
                    None,
                    None,
                    "unit*pct./period",
                    "unit/⁰C",
                    "pct./30sec",
                    "m3",
                    "EUR_2020/tank/year",
                ],
            }
        )
        expected_dataframe = pandas.DataFrame(
            {
                "par": [
                    "energy storage capacity for one unit",
                    "typical temperature difference in storage",
                    "fixed o&m",
                    "lifetime in total number of cycles",
                    "cycle life",
                    "p1",
                    "p2",
                    "p3",
                    "p4",
                    "p5",
                ],
                "unit": [
                    "MWh",
                    "K",
                    "percent/year",
                    "cycles",
                    "cycles",
                    "unit*percent",
                    "unit/C",
                    "percent",
                    "meter**3",
                    "EUR_2020/year",
                ],
            }
        )
        input_dataframe[["par", "unit"]] = (
            input_dataframe[["par", "unit"]]
            .apply(standardize_units, axis=1)
            .reset_index(drop=True)
        )
        comparison_df = input_dataframe.compare(expected_dataframe)
        assert comparison_df.empty

    def test_filter_parameters_true_false(self) -> None:
        """Check if filter_parameters works as expected."""
        dataframe = pandas.DataFrame(
            {"par": ["technical lifetime", "other"], "Technology": ["A", "B"]}
        )
        filtered = filter_parameters(dataframe, True)
        assert all(filtered["par"] == "technical lifetime")
        filtered2 = filter_parameters(dataframe, False)
        assert len(filtered2) == 2

    @pytest.mark.webarchive  # type: ignore
    def test_build_technology_collection(self, tmp_path: pathlib.Path) -> None:
        """Check if build_technology_collection works as expected."""
        dataframe = pandas.DataFrame(
            {
                "est": ["base"],
                "year": [2020],
                "ws": ["Tech1"],
                "Technology": ["Solar PV"],
                "par": ["specific_investment"],
                "val": [1000],
                "unit": ["EUR_2020/kW"],
            }
        )
        sources_path = pathlib.Path(tmp_path, "sources.json")
        tech_collection = build_technology_collection(dataframe, sources_path)
        assert len(tech_collection.technologies) == 1
        assert "specific_investment" in tech_collection.technologies[0].parameters
