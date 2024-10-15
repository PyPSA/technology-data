#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys
import pandas as pd

sys.path.append("./scripts")

from test.conftest import get_config_dict
from compile_cost_assumptions_nrel import calculate_fom_percentage, concatenate_columns, filter_input_file, repeat_values, replace_value_name

path_cwd = pathlib.Path.cwd()

@pytest.mark.parametrize(
    "year, expected",
    [(2019, "atb_e_2019 - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet"), (2022, (3092, 12)), (2024, (20334, 12))],
)
def test_filter_input_file(get_config_dict, year, expected):
    """
    Verify what returned by filter_input_file.
    """
    config_dict = get_config_dict
    list_years = config_dict["years"]
    list_columns_to_keep = config_dict["nrel_atb"]["nrel_atb_columns_to_keep"]
    list_core_metric_parameter_to_keep = config_dict["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(year))
    if year in [2022, 2024]:
        input_file = filter_input_file(input_file_path, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep)
        assert input_file.shape == expected
    else:
        with pytest.raises(Exception) as excinfo:
            input_file = filter_input_file(input_file_path, list_years, list_columns_to_keep,
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
    dataframe = pd.read_csv("/Users/fabriziofinozzi/Desktop/OpenEnergyTransition/repo/technology-data/intermediate_atp_2022.csv")
    repeat_values(dataframe, "technology_alias_detail")
    assert False

# def test_calculate_fom_percentage():
#     test_df = pd.DataFrame({"Name": ["Tom", "Paul", "Sarah"], "Age": [31, 42, 56]})
#     calculate_fom_percentage(row, test_df)


