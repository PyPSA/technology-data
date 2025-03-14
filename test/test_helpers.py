# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import copy
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.append("./scripts")

from _helpers import prepare_inflation_rate

path_cwd = pathlib.Path.cwd()


def test_prepare_inflation_rate():
    """
    The test verifies what is returned by prepare_inflation_rate.
    """

    output_series = prepare_inflation_rate("/Users/fabriziofinozzi/Desktop/OpenEnergyTransition/repo/technology-data/inputs/Eurostat_inflation_rates.xlsx", "usd").round(decimals=3)
    print(output_series)
    assert False
#    reference_output_series = pd.Series(
#        [0.02, 0.015, 0.025, 0.018],
#        index=[2001, 2002, 2003, 2004],
#        name="European Union - 27 countries (from 2020)",
#    )
#    comparison_series = output_series.compare(reference_output_series)
#    assert comparison_series.size == 0