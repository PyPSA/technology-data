#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys
import pandas as pd
import numpy as np

sys.path.append("./scripts")

from test.conftest import get_config_dict
from compile_cost_assumptions_nrel import calculate_fom_percentage, filter_input_file, get_convertion_dictionary, pre_process_input_file, replace_value_name, update_cost_values

path_cwd = pathlib.Path.cwd()

@pytest.mark.parametrize(
    "file_year, year, expected",
    [(2019, 2020, "atb_e_2019 - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet"), (2022, 2020, (3002, 11)), (2024, 2025, (3126, 11)), (2024, 2030, (3312, 11)), (2024, 2035, (3336, 11)), (2024, 2040, (3336, 11)), (2024, 2045, (3336, 11)), (2024, 2050, (3336, 11))],
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
        assert "aeo" not in input_file["technology"].astype(str).str.casefold().unique()

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
    "input_file_year, year, expected", [(2022, 2020, (3002, 11)), (2024, 2050, (3336, 11))],
)
def test_pre_process_input_file(get_config_dict, input_file_year, year, expected):
    config_dict = get_config_dict
    input_file_path = pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(input_file_year))
    nrel_atb_columns_to_keep = config_dict["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = config_dict["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_source_link = config_dict["nrel_atb"]["nrel_atb_source_link"]
    output_df = pre_process_input_file(input_file_path, year, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)
    reference_parameter_list = sorted(["investment", "CF", "FOM", "VOM", "fuel"])
    output_parameter_list = sorted(list(output_df["parameter"].unique()))
    assert output_df.shape == expected
    assert len(output_parameter_list) == len(reference_parameter_list)
    assert all([x == y for x, y in zip(reference_parameter_list, output_parameter_list)])

def test_update_cost_values():

    test_atb_df = pd.DataFrame({
        "technology": ["coal", "CCGT", "hydro", "ror", "offwind", "onwind", "Offshore Wind - Class 3", "Offshore Wind - Class 10"],
        "financial_case": ["Market", "R&D", "Market", "Market", "Market", "R&D", "Market", "R&D"],
        "scenario": ["Advanced", "Conservative", "Moderate", "Advanced", "Conservative", "Moderate", "Conservative", "Moderate"],
        "tax_credit_case": ["ITC", "ITC", "ITC", "ITC", "ITC", "ITC", "ITC", "ITC"]
    })

    test_cost_df = pd.DataFrame()
    test_cost_df["technology"] = ["coal", "CCGT", "hydro", "ror", "offwind", "onwind", "BEV Bus City", "Ammonia cracker"]

    reference_df = pd.DataFrame({
        "technology": ["BEV Bus City", "Ammonia cracker", "coal", "CCGT", "hydro", "ror", "offwind", "onwind", "Offshore Wind - Class 3", "Offshore Wind - Class 10"],
        "financial_case": [np.nan, np.nan, "Market", "R&D", "Market", "Market", "Market", "R&D", "Market", "R&D"],
        "scenario": [np.nan, np.nan, "Advanced", "Conservative", "Moderate", "Advanced", "Conservative", "Moderate", "Conservative", "Moderate"],
        "tax_credit_case": [np.nan, np.nan, "ITC", "ITC", "ITC", "ITC", "ITC", "ITC", "ITC", "ITC"]
    })

    technology_conversion_dictionary = get_convertion_dictionary("technology")
    output_df = update_cost_values(test_cost_df, test_atb_df, technology_conversion_dictionary, ["financial_case", "scenario", "tax_credit_case"])
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty
