# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.append("./scripts")

from compile_cost_assumptions_usa import (
    calculate_fom_percentage,
    duplicate_fuel_cost,
    filter_atb_input_file,
    get_conversion_dictionary,
    get_query_string,
    pre_process_atb_input_file,
    pre_process_cost_input_file,
    pre_process_manual_input_usa,
    query_cost_dataframe,
    replace_value_name,
)

path_cwd = pathlib.Path.cwd()
additional_occ_query_string = "atb_year == @x.atb_year & core_metric_case == @x.core_metric_case & core_metric_parameter.str.casefold() == 'additional occ' & core_metric_variable == @x.core_metric_variable & display_name == @x.display_name & scenario == @x.scenario & technology == @x.technology & technology_alias == @x.technology_alias"
capex_query_string = "atb_year == @x.atb_year & core_metric_case == @x.core_metric_case & core_metric_parameter.str.casefold() == 'capex' & core_metric_variable == @x.core_metric_variable & display_name == @x.display_name & scenario == @x.scenario & technology == @x.technology & technology_alias == @x.technology_alias"


@pytest.mark.parametrize(
    "file_year, year, expected",
    [
        (
            2019,
            2020,
            "atb_e_2019 - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet",
        ),
        (2022, 2020, (2960, 10)),
        (2024, 2025, (3036, 10)),
        (2024, 2030, (3222, 10)),
        (2024, 2035, (3246, 10)),
        (2024, 2040, (3246, 10)),
        (2024, 2045, (3246, 10)),
        (2024, 2050, (3246, 10)),
    ],
)
def test_filter_atb_input_file(config, file_year, year, expected):
    """
    The test verifies what is returned by filter_atb_input_file.
    """
    list_columns_to_keep = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    list_core_metric_parameter_to_keep = config["nrel_atb"][
        "nrel_atb_core_metric_parameter_to_keep"
    ]
    nrel_atb_technology_to_remove = config["nrel_atb"]["nrel_atb_technology_to_remove"]
    input_file_path = pathlib.Path(
        path_cwd, "inputs", "US", f"atb_e_{file_year}.parquet"
    )
    if file_year in [2022, 2024]:
        input_file = filter_atb_input_file(
            input_file_path,
            year,
            list_columns_to_keep,
            list_core_metric_parameter_to_keep,
            nrel_atb_technology_to_remove,
        )
        assert (
            len(
                set(input_file["core_metric_parameter"].unique()).difference(
                    set(list_core_metric_parameter_to_keep)
                )
            )
            == 0
        )
        assert int(input_file["core_metric_variable"].unique().item()) == year
        assert input_file.drop_duplicates(keep="first").shape == expected
        assert "aeo" not in input_file["technology"].astype(str).str.casefold().unique()
        assert (
            "coal-ccs-95% -> transformational tech"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "coal-max-ccs -> transformational tech"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "coal-new -> transformational tech"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "ng combined cycle 95% ccs (f-frame basis -> transformational tech)"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "ng combined cycle 95% ccs (h-frame basis -> transformational tech)"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "ng combined cycle max ccs (f-frame basis -> transformational tech)"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
        assert (
            "ng combined cycle max ccs (h-frame basis -> transformational tech)"
            not in input_file["technology"].astype(str).str.casefold().unique()
        )
    else:
        with pytest.raises(Exception) as excinfo:
            input_file = filter_atb_input_file(
                input_file_path,
                year,
                list_columns_to_keep,
                list_core_metric_parameter_to_keep,
                nrel_atb_technology_to_remove,
            )
        assert str(excinfo.value) == expected


@pytest.mark.parametrize(
    "parameter_value, columns_to_exclude, expected",
    [
        ("additional occ", ["units", "value"], additional_occ_query_string),
        ("capex", ["units", "value"], capex_query_string),
        (
            "fail_test",
            ["random_column", "value"],
            "The following columns ['random_column'] are not included in the original list",
        ),
    ],
)
def test_get_query_string(config, parameter_value, columns_to_exclude, expected):
    """
    The test verifies what is returned by get_query_string.
    """
    columns_list = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    if parameter_value == "fail_test":
        with pytest.raises(Exception) as excinfo:
            output_string = get_query_string(
                columns_list, columns_to_exclude, parameter_value
            )
            print(str(excinfo.value))
        assert str(excinfo.value) == expected
    else:
        output_string = get_query_string(
            columns_list, columns_to_exclude, parameter_value
        )
        assert output_string == expected


@pytest.mark.parametrize(
    "display_name, expected",
    [
        ("Coal-new", 2.13),
        ("Coal-95%-CCS", 2.06),
        ("Coal-99%-CCS", 2.05),
        ("Coal-IGCC", 2.38),
        ("Coal-IGCC-90%-CCS", 2.37),
        ("Coal integrated retrofit 90%-CCS", 7.37),
        ("Coal integrated retrofit 95%-CCS", 7.22),
    ],
)
def test_calculate_fom_percentage(config, display_name, expected):
    """
    The test verifies what is returned by calculate_fom_percentage.
    """
    columns_list = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    test_df = pd.read_csv(pathlib.Path(path_cwd, "test", "test_data", "coal_test.csv"))
    test_df["value"] = test_df.apply(
        lambda x: calculate_fom_percentage(x, test_df, columns_list), axis=1
    )
    assert (
        test_df.loc[
            (test_df["display_name"] == display_name)
            & (test_df["core_metric_parameter"] == "Fixed O&M")
        ]["value"].item()
        == expected
    )


def test_replace_value_name():
    """
    The test verifies what is returned by replace_value_name.
    """
    test_df = pd.DataFrame(
        {
            "Name": ["Tom", "Paul", "John", "Sarah"],
            "Age": [31, 42, 12, 56],
            "Country": ["US", "DE", "UK", "IT"],
        }
    )
    reference_df = pd.DataFrame(
        {
            "Name": ["Tom", "Paul", "John", "Sarah"],
            "Age": [31, 42, 12, 56],
            "Country": ["United States", "Germany", "United Kingdom", "IT"],
        }
    )
    conversion_dict = {
        "US": "United States",
        "DE": "Germany",
        "UK": "United Kingdom",
        "ES": "Spain",
    }
    output_df = replace_value_name(test_df, conversion_dict, "Country")
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


def test_pre_process_cost_input_file(tmpdir, cost_dataframe):
    """
    The test verifies what is returned by pre_process_cost_input_file.
    """
    reference_df = pd.DataFrame(
        {
            "technology": ["coal", "coal", "another_tech"],
            "parameter": ["co2 intensity", "lifetime", "investment"],
            "value": [1.0, 1.0, 3.0],
            "unit": ["unit", "unit", "unit"],
            "source": ["source", "source", "source"],
            "further description": ["g", "h", "i"],
            "currency_year": [2020, 2020, 2020],
            "financial_case": [np.nan, np.nan, np.nan],
            "scenario": [np.nan, np.nan, np.nan],
        }
    )
    input_file_path = pathlib.Path(tmpdir, "tmp_costs.csv")
    cost_dataframe.to_csv(input_file_path, index=False)
    output_df = pre_process_cost_input_file(
        input_file_path, ["financial_case", "scenario"]
    )
    comparison_df = output_df.compare(reference_df)
    pathlib.Path(input_file_path).unlink(missing_ok=True)
    assert comparison_df.empty


@pytest.mark.parametrize(
    "input_file_year, year, expected",
    [(2022, 2020, (2960, 9)), (2024, 2050, (3246, 9))],
)
def test_pre_process_atb_input_file(config, input_file_year, year, expected):
    """
    The test verifies what is returned by pre_process_atb_input_file.
    """
    input_file_path = pathlib.Path(
        path_cwd, "inputs", "US", f"atb_e_{input_file_year}.parquet"
    )
    nrel_atb_columns_to_keep = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = config["nrel_atb"][
        "nrel_atb_core_metric_parameter_to_keep"
    ]
    nrel_atb_technology_to_remove = config["nrel_atb"]["nrel_atb_technology_to_remove"]
    nrel_atb_source_link = config["nrel_atb"]["nrel_atb_source_link"]
    nrel_atb_further_description = config["nrel_atb"]["nrel_atb_further_description"]
    output_df = pre_process_atb_input_file(
        input_file_path,
        nrel_atb_source_link,
        nrel_atb_further_description,
        year,
        nrel_atb_columns_to_keep,
        nrel_atb_core_metric_parameter_to_keep,
        nrel_atb_technology_to_remove,
    )
    reference_parameter_list = sorted(["investment", "CF", "FOM", "VOM", "fuel"])
    output_parameter_list = sorted(list(output_df["parameter"].unique()))
    assert output_df.shape == expected
    assert len(output_parameter_list) == len(reference_parameter_list)
    assert all(
        [x == y for x, y in zip(reference_parameter_list, output_parameter_list)]
    )
    units_df = output_df[["parameter", "unit"]].drop_duplicates(keep="first")
    assert (
        units_df.loc[units_df["parameter"].astype(str).str.casefold() == "investment"][
            "unit"
        ].item()
        == "USD/kW"
    )
    assert (
        units_df.loc[units_df["parameter"].astype(str).str.casefold() == "cf"][
            "unit"
        ].item()
        == "per unit"
    )
    assert (
        units_df.loc[units_df["parameter"].astype(str).str.casefold() == "fom"][
            "unit"
        ].item()
        == "%/year"
    )
    assert (
        units_df.loc[units_df["parameter"].astype(str).str.casefold() == "vom"][
            "unit"
        ].item()
        == "USD/MWh"
    )
    assert (
        units_df.loc[units_df["parameter"].astype(str).str.casefold() == "fuel"][
            "unit"
        ].item()
        == "USD/MWh"
    )


def test_query_cost_dataframe(cost_dataframe):
    """
    The test verifies what is returned by query_cost_dataframe.
    """
    reference_df = pd.DataFrame(
        {
            "technology": ["coal", "coal", "another_tech"],
            "parameter": ["co2 intensity", "lifetime", "investment"],
            "value": [1.0, 1.0, 3.0],
            "unit": ["unit", "unit", "unit"],
            "source": ["source", "source", "source"],
            "further description": ["g", "h", "i"],
            "currency_year": [2020, 2020, 2020],
        }
    )
    pypsa_technology_dictionary = get_conversion_dictionary("pypsa_technology_name")
    parameter_dictionary = get_conversion_dictionary("parameter")
    output_df = query_cost_dataframe(
        cost_dataframe, pypsa_technology_dictionary, parameter_dictionary
    )
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


def test_duplicate_fuel_cost(config):
    """
    The test verifies what is returned by duplicate_fuel_cost.
    """
    input_file_path = pathlib.Path(path_cwd, "inputs", "US", "fuel_costs_usa.csv")
    output_df = duplicate_fuel_cost(input_file_path, config["years"])
    assert output_df.shape == (21, 10)

    # The row corresponding to the coal technology for 2025 is replicated for any later year
    assert (
        output_df.loc[
            (output_df["technology"] == "coal")
            & (output_df["year"].isin([2025, 2030, 2035, 2040, 2045, 2050]))
        ]["value"]
        .unique()
        .item()
        == 8.12
    )
    assert (
        output_df.loc[
            (output_df["technology"] == "coal")
            & (output_df["year"].isin([2025, 2030, 2035, 2040, 2045, 2050]))
        ]["further description"]
        .unique()
        .item()
        == "46.97 USD/short ton of bituminous coal with energy content = 19.73 million BTU/short ton"
    )

    # The row corresponding to the gas technology for 2030 is replicated for any later year
    assert (
        output_df.loc[
            (output_df["technology"] == "gas")
            & (output_df["year"].isin([2030, 2035, 2040, 2045, 2050]))
        ]["value"]
        .unique()
        .item()
        == 14.05
    )

    # The row corresponding to the oil technology for 2030 is replicated for any later year
    assert (
        output_df.loc[
            (output_df["technology"] == "oil")
            & (output_df["year"].isin([2030, 2035, 2040, 2045, 2050]))
        ]["value"]
        .unique()
        .item()
        == 44.22
    )


@pytest.mark.parametrize(
    "year, expected",
    [
        (2020, (130, 9)),
        (2025, (130, 9)),
        (2030, (130, 9)),
        (2035, (130, 9)),
        (2040, (130, 9)),
        (2045, (130, 9)),
        (2050, (130, 9)),
    ],
)
def test_pre_process_manual_input_usa(config, year, expected):
    """
    The test verifies what is returned by pre_process_manual_input_usa.
    """
    list_of_years = config["years"]
    manual_input_usa_file_path = pathlib.Path(
        path_cwd, "inputs", "US", "manual_input_usa.csv"
    )
    year = 2020
    output_dataframe = pre_process_manual_input_usa(
        manual_input_usa_file_path,
        list_of_years,
        year,
    )
    assert output_dataframe.shape == expected


def test_final_output(tmpdir, cost_dataframe, atb_cost_dataframe):
    """
    The test verifies what is returned by the concatenation of the existing cost file and NREL/ATB.
    """
    reference_df = pd.DataFrame(
        {
            "technology": [
                "coal",
                "coal",
                "another_tech",
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
            ],
            "parameter": [
                "co2 intensity",
                "lifetime",
                "investment",
                "investment",
                "FOM",
                "VOM",
                "fuel",
                "investment",
                "discount rate",
            ],
            "value": [1.0, 1.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            "unit": [
                "unit",
                "unit",
                "unit",
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
            ],
            "source": [
                "source",
                "source",
                "source",
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
            ],
            "further description": ["g", "h", "i", "a", "b", "c", "d", "e", "f"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020],
            "financial_case": [
                np.nan,
                np.nan,
                np.nan,
                "R&D",
                "R&D",
                "R&D",
                "R&D",
                "R&D",
                "R&D",
            ],
            "scenario": [
                np.nan,
                np.nan,
                np.nan,
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
            ],
        }
    )
    input_cost_path = pathlib.Path(tmpdir, "tmp_costs.csv")
    cost_dataframe.to_csv(input_cost_path, index=False)

    cost_df = pre_process_cost_input_file(
        input_cost_path, ["financial_case", "scenario"]
    )

    output_df = pd.concat([cost_df, atb_cost_dataframe]).reset_index(drop=True)

    comparison_df = output_df.compare(reference_df)
    pathlib.Path(input_cost_path).unlink(missing_ok=True)
    assert comparison_df.empty
