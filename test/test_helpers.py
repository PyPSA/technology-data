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
            "European Union - 27 countries (from 2020)",
            [2017, 2022, 2023],
            [0.016, 0.092, 0.064],
        ),
        ("usd", "United States", [2017, 2022, 2023], [0.018, 0.087, 0.03]),
        (
            "EuR",
            "European Union - 27 countries (from 2020)",
            [2017, 2022, 2023],
            [0.016, 0.092, 0.064],
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
        path_cwd, "inputs", "Eurostat_inflation_rates.xlsx"
    )
    output_series = prepare_inflation_rate(
        inflation_rate_input_file_path, currency_to_use
    ).round(decimals=3)
    assert np.array_equal(
        output_series.loc[expected_index].values, np.array(expected_values)
    )
    assert output_series.name == expected_series_name
