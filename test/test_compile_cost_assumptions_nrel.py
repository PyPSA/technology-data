#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys
import pandas as pd

sys.path.append("./scripts")

from test.conftest import get_config_dict
from compile_cost_assumptions_nrel import filter_input_file, replace_value_name

path_cwd = pathlib.Path.cwd()

@pytest.mark.parametrize(
    "year, expected",
    [(2022, (3098, 12)), (2024, (20334, 12))],
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
    input_file = filter_input_file(input_file_path, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep)
    assert input_file.shape == expected


def test_replace_value_name():
    """
    Verify what returned by replace_value_name.
    """
    test_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["US", "DE", "UK", "IT"]})
    reference_df = pd.DataFrame({"Name": ["Tom", "Paul", "John", "Sarah"], "Age": [31, 42, 12, 56], "Country": ["United States", "Germany", "United Kingdom", "IT"]})
    conversion_dict = {"US": "United States", "DE": "Germany", "UK": "United Kingdom"}
    output_df = replace_value_name(test_df, conversion_dict, "Country")
    comparison_df = output_df.compare(reference_df)
    assert comparison_df.empty
