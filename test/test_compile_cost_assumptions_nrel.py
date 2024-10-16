#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys
import pandas as pd

sys.path.append("./scripts")

from test.conftest import get_config_dict
from compile_cost_assumptions_nrel import calculate_fom_percentage, concatenate_columns, filter_input_file, pre_process_input_file, repeat_values, replace_value_name

path_cwd = pathlib.Path.cwd()

@pytest.mark.parametrize(
    "file_year, year, expected",
    [(2019, 2020, "atb_e_2019 - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet"), (2022, 2020, (3092, 12)), (2024, 2025, (3213, 12)), (2024, 2030, (3405, 12)), (2024, 2035, (3429, 12)), (2024, 2040, (3429, 12)), (2024, 2045, (3429, 12)), (2024, 2050, (3429, 12))],
)
def test_filter_input_file(get_config_dict, file_year, year, expected):
    """
    Verify what returned by filter_input_file.
    """
    config_dict = get_config_dict
    list_columns_to_keep = config_dict["nrel_atb"]["nrel_atb_columns_to_keep"]
    list_core_metric_parameter_to_keep = config_dict["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(file_year))
    if file_year in [2022, 2024]:
        input_file = filter_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep)
        assert input_file.shape == expected
    else:
        with pytest.raises(Exception) as excinfo:
            input_file = filter_input_file(input_file_path, year, list_columns_to_keep,
                                           list_core_metric_parameter_to_keep)
        assert str(excinfo.value) == expected


def test_replace_value_name():
    """
    Verify what returned by replace_value_name.
    """
    test_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["US", "DE", "UK", "IT"]})
    reference_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["United States", "Germany", "United Kingdom", "IT"]})
    conversion_dict = {"US": "United States", "DE": "Germany", "UK": "United Kingdom", "ES": "Spain"}
    output_df = replace_value_name(test_df, conversion_dict, "Country")
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


def test_concatenate_columns():
    """
    Verify what returned by concatenate_columns.
    """
    test_df = pd.DataFrame({"Name": ["Tom", "Paul", "Sarah"], "Age": [31, 42, 56]})
    reference_df = pd.DataFrame({"Name": ["Tom", "Paul", "Sarah"], "Age": [31, 42, 56], "Name_Age": ["Tom_31", "Paul_42", "Sarah_56"]})
    output_df = concatenate_columns(test_df, "Name_Age", ["Name", "Age"])
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty


def test_repeat_values():
    """
    Verify what returned by repeat_values.
    """
    test_df = pd.DataFrame({"Name": ["Tom", "Paul", "Sarah", "Sabrina"], "Age": [31, 31, 31, 12],
                            "City": ["Rome", "Rome", "Rome", "Berlin"]})

    # Numeric values
    numeric_reference_df = pd.DataFrame(
        {"Name": ["Sabrina", "Tom", "Tom", "Tom", "Paul", "Paul", "Paul", "Sarah", "Sarah", "Sarah"],
         "Age": [12, 40, 41, 42, 40, 41, 42, 40, 41, 42],
         "City": ["Berlin", "Rome", "Rome", "Rome", "Rome", "Rome", "Rome", "Rome", "Rome", "Rome"]
         })
    numeric_output_df = repeat_values(test_df, [40, 41, 42], "Age", 31, 3)
    numeric_comparison_df = numeric_output_df.compare(numeric_reference_df)

    # String values
    string_reference_df = pd.DataFrame(
        {"Name": ["Sabrina", "Tom", "Tom", "Tom", "Paul", "Paul", "Paul", "Sarah", "Sarah", "Sarah"],
         "Age": [12, 31, 31, 31, 31, 31, 31, 31, 31, 31],
         "City": ["Berlin", "Munich", "Paris", "Mainz", "Munich", "Paris", "Mainz", "Munich", "Paris", "Mainz"]
         })
    string_output_df = repeat_values(test_df, ["Munich", "Paris", "Mainz"], "City", "rome", 3)
    string_comparison_df = string_output_df.compare(string_reference_df)

    assert numeric_comparison_df.empty
    assert string_comparison_df.empty

@pytest.mark.parametrize(
    "display_name, expected",
    [("Coal-new", 2.13), ("Coal-95%-CCS", 2.06), ("Coal-99%-CCS", 2.05), ("Coal-IGCC", 2.38), ("Coal-IGCC-90%-CCS", 2.37), ("Coal integrated retrofit 90%-CCS", 7.37), ("Coal integrated retrofit 95%-CCS", 7.22)],
)
def test_calculate_fom_percentage(display_name, expected):
    """
    Verify what returned by calculate_fom_percentage.
    """
    test_df = pd.read_csv(pathlib.Path(path_cwd, "test", "test_data", "coal_test.csv"))
    test_df["value"] = test_df.apply(lambda x: calculate_fom_percentage(x, test_df), axis=1)
    assert test_df.loc[(test_df["display_name"] == display_name) & (test_df["core_metric_parameter"] == "Fixed O&M")]["value"].item() == expected

@pytest.mark.parametrize(
    "input_file_year, year, expected", [(2022, 2020, (3098, 11)), (2024, 2050, (3435, 11))],
)
def test_pre_process_input_file(get_config_dict, input_file_year, year, expected):
    config_dict = get_config_dict
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(input_file_year))
    nrel_atb_columns_to_keep = config_dict["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = config_dict["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_source_link = config_dict["nrel_atb"]["nrel_atb_source_link"]
    output_df = pre_process_input_file(input_file_path, year, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)
    reference_parameter_list = sorted(["investment", "CF", "FOM", "VOM", "fuel", "discount rate"])
    output_parameter_list = sorted(list(output_df["parameter"].unique()))
    print(output_parameter_list, reference_parameter_list)
    assert output_df.shape == expected
    assert len(output_parameter_list) == len(reference_parameter_list)
    assert all([x == y for x, y in zip(reference_parameter_list, output_parameter_list)])
