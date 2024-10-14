#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import pytest
import sys

sys.path.append("./scripts")


from test.conftest import get_config_dict
from compile_cost_assumptions_nrel import read_atb_input_file

path_cwd = pathlib.Path.cwd()

@pytest.mark.parametrize(
    "year,expected",
    [(2022, (286998, 14)), (2024, (572232, 19))],
)
def test_read_atb_input_file(year, expected):
    list_inputs = [pathlib.Path(path_cwd, "inputs", "atb_e_{}.parquet".format(year))]
    input_file = read_atb_input_file(list_inputs)
    assert input_file[0].shape == expected
