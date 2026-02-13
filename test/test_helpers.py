# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import pathlib
import sys

import numpy as np
import pytest

sys.path.append("./scripts")

from _helpers import prepare_inflation_rate

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "currency_to_use, expected_series_name, expected_index, expected_values",
    [
        (
            "eur",
            "European Economic Area (EEA18-1995, EEA28-2004, EEA30-2007, EEA31-2013, EEA30-2020)",
            [2017, 2022, 2023],
            [0.017, 0.092, 0.064],
        ),
        ("usd", "United States", [2017, 2022, 2023], [0.018, 0.087, 0.03]),
        (
            "EuR",
            "European Economic Area (EEA18-1995, EEA28-2004, EEA30-2007, EEA31-2013, EEA30-2020)",
            [2017, 2022, 2023],
            [0.017, 0.092, 0.064],
        ),
        ("USD", "United States", [2017, 2022, 2023], [0.018, 0.087, 0.03]),
    ],
)
def test_prepare_inflation_rate(
    currency_to_use, expected_series_name, expected_index, expected_values
):
    """
    The test verifies what is returned by prepare_inflation_rate.
    """
    inflation_rate_input_file_path = pathlib.Path(
        path_cwd, "inputs", "prc_hicp_aind__custom_20097956_spreadsheet.xlsx"
    )
    output_series = prepare_inflation_rate(
        inflation_rate_input_file_path, currency_to_use
    ).round(decimals=3)
    assert np.array_equal(
        output_series.loc[expected_index].values, np.array(expected_values)
    )
    assert output_series.name == expected_series_name
