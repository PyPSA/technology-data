#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys
import pandas as pd
import numpy as np

sys.path.append("./scripts")

from compile_cost_assumptions_nrel import calculate_fom_percentage, filter_input_file, get_convertion_dictionary, get_query_string, pre_process_input_file, replace_value_name, update_cost_values

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "file_year, year, expected",
    [(2019, 2020, "atb_e_2019 - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet"), (2022, 2020, (2960, 10)), (2024, 2025, (3036, 10)), (2024, 2030, (3222, 10)), (2024, 2035, (3246, 10)), (2024, 2040, (3246, 10)), (2024, 2045, (3246, 10)), (2024, 2050, (3246, 10))],
)
def test_filter_input_file(config, file_year, year, expected):
    """
    The test verifies what is returned by filter_input_file.
    """
    list_columns_to_keep = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    list_core_metric_parameter_to_keep = config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_technology_to_remove = config["nrel_atb"]["nrel_atb_technology_to_remove"]
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(file_year))
    if file_year in [2022, 2024]:
        input_file = filter_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_atb_technology_to_remove)
        assert input_file.shape == expected
        assert "aeo" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "coal-ccs-95% -> transformational tech" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "coal-max-ccs -> transformational tech" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "coal-new -> transformational tech" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "ng combined cycle 95% ccs (f-frame basis -> transformational tech)" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "ng combined cycle 95% ccs (h-frame basis -> transformational tech)" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "ng combined cycle max ccs (f-frame basis -> transformational tech)" not in input_file["technology"].astype(str).str.casefold().unique()
        assert "ng combined cycle max ccs (h-frame basis -> transformational tech)" not in input_file["technology"].astype(str).str.casefold().unique()
    else:
        with pytest.raises(Exception) as excinfo:
            input_file = filter_input_file(input_file_path, year, list_columns_to_keep,
                                           list_core_metric_parameter_to_keep, nrel_atb_technology_to_remove)
        assert str(excinfo.value) == expected


def test_replace_value_name():
    """
    The test verifies what is returned by replace_value_name.
    """
    test_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["US", "DE", "UK", "IT"]})
    reference_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["United States", "Germany", "United Kingdom", "IT"]})
    conversion_dict = {"US": "United States", "DE": "Germany", "UK": "United Kingdom", "ES": "Spain"}
    output_df = replace_value_name(test_df, conversion_dict, "Country")
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


@pytest.mark.parametrize(
    "display_name, expected",
    [("Coal-new", 2.13), ("Coal-95%-CCS", 2.06), ("Coal-99%-CCS", 2.05), ("Coal-IGCC", 2.38), ("Coal-IGCC-90%-CCS", 2.37), ("Coal integrated retrofit 90%-CCS", 7.37), ("Coal integrated retrofit 95%-CCS", 7.22)],
)
def test_calculate_fom_percentage(config, display_name, expected):
    """
    The test verifies what is returned by calculate_fom_percentage.
    """
    columns_list = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    test_df = pd.read_csv(pathlib.Path(path_cwd, "test", "test_data", "coal_test.csv"))
    test_df["value"] = test_df.apply(lambda x: calculate_fom_percentage(x, test_df, columns_list), axis=1)
    assert test_df.loc[(test_df["display_name"] == display_name) & (test_df["core_metric_parameter"] == "Fixed O&M")]["value"].item() == expected


@pytest.mark.parametrize(
    "input_file_year, year, expected", [(2022, 2020, (2960, 9)), (2024, 2050, (3246, 9))],
)
def test_pre_process_input_file(config, input_file_year, year, expected):
    """
    The test verifies what is returned by pre_process_input_file.
    """
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(input_file_year))
    nrel_atb_columns_to_keep = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_technology_to_remove = config["nrel_atb"]["nrel_atb_technology_to_remove"]
    nrel_atb_source_link = config["nrel_atb"]["nrel_atb_source_link"]
    output_df = pre_process_input_file(input_file_path, year, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link, nrel_atb_technology_to_remove)
    reference_parameter_list = sorted(["investment", "CF", "FOM", "VOM", "fuel"])
    output_parameter_list = sorted(list(output_df["parameter"].unique()))
    assert output_df.shape == expected
    assert len(output_parameter_list) == len(reference_parameter_list)
    assert all([x == y for x, y in zip(reference_parameter_list, output_parameter_list)])


def test_update_cost_values(cost_dataframe, atb_cost_dataframe):
    """
    The test verifies what is returned by update_cost_values.
    """
    reference_df = pd.DataFrame(
        {
            "technology": ["coal", "coal", "coal", "coal", "coal", "coal", "coal", "coal"],
            "parameter": ["co2 intensity", "lifetime", "investment", "FOM", "VOM", "fuel", "investment", "discount rate"],
            "value": [1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            "unit": ["unit", "unit", "unit_atb", "unit_atb", "unit_atb", "unit_atb", "unit_atb", "unit_atb"],
            "source": ["source", "source", "source_atb", "source_atb", "source_atb", "source_atb", "source_atb", "source_atb"],
            "further description": ["g", "h", "a", "b", "c", "d", "e", "f"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020],
            "financial_case": [np.nan, np.nan, "R&D", "R&D", "R&D", "R&D", "R&D", "R&D"],
            "scenario": [np.nan, np.nan, "Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate"]
        }
    )
    pypsa_technology_dictionary = get_convertion_dictionary("pypsa_technology_name")
    parameter_dictionary = get_convertion_dictionary("parameter")
    columns_to_add_list = ["financial_case", "scenario"]
    output_df = update_cost_values(cost_dataframe, atb_cost_dataframe, pypsa_technology_dictionary, parameter_dictionary, columns_to_add_list)
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


@pytest.mark.parametrize(
        "parameter_value, columns_to_exclude, expected", [("additional occ", ["units", "value"], "atb_year == @x.atb_year & core_metric_case == @x.core_metric_case & core_metric_parameter.str.casefold() == 'additional occ' & core_metric_variable == @x.core_metric_variable & display_name == @x.display_name & scenario == @x.scenario & technology == @x.technology & technology_alias == @x.technology_alias"), ("capex", ["units", "value"], "atb_year == @x.atb_year & core_metric_case == @x.core_metric_case & core_metric_parameter.str.casefold() == 'capex' & core_metric_variable == @x.core_metric_variable & display_name == @x.display_name & scenario == @x.scenario & technology == @x.technology & technology_alias == @x.technology_alias"), ("fail_test", ["random_column", "value"], "The following columns ['random_column'] are not included in the original list")],
    )
def test_get_query_string(config, parameter_value, columns_to_exclude, expected):
    """
    The test verifies what is returned by get_query_string.
    """
    columns_list = config["nrel_atb"]["nrel_atb_columns_to_keep"]
    if parameter_value == "fail_test":
        with pytest.raises(Exception) as excinfo:
            output_string = get_query_string(columns_list, columns_to_exclude, parameter_value)
            print(str(excinfo.value))
        assert str(excinfo.value) == expected
    else:
        output_string = get_query_string(columns_list, columns_to_exclude, parameter_value)
        assert output_string == expected


# def test_test(config):
#     list_columns_to_keep = config["nrel_atb"]["nrel_atb_columns_to_keep"]
#     list_core_metric_parameter_to_keep = config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
#     nrel_atb_technology_to_remove = config["nrel_atb"]["nrel_atb_technology_to_remove"]
#     input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(2022))
#     input_file = filter_input_file(input_file_path, 2022, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_atb_technology_to_remove)
#     assert False