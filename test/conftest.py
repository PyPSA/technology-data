# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import pathlib

import pandas as pd
import pytest
import yaml


@pytest.fixture(scope="session")
def config():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config.yaml")
    with open(path_config) as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def cost_dataframe():
    return pd.DataFrame(
        {
            "technology": [
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
                "coal",
                "another_tech",
            ],
            "parameter": [
                "investment",
                "FOM",
                "VOM",
                "fuel",
                "investment",
                "discount rate",
                "co2 intensity",
                "lifetime",
                "investment",
            ],
            "value": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0],
            "unit": [
                "unit",
                "unit",
                "unit",
                "unit",
                "unit",
                "unit",
                "unit",
                "unit",
                "unit",
            ],
            "source": [
                "source",
                "source",
                "source",
                "source",
                "source",
                "source",
                "source",
                "source",
                "source",
            ],
            "further description": ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020],
        },
        index=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    )


@pytest.fixture(scope="function")
def atb_cost_dataframe():
    return pd.DataFrame(
        {
            "technology": ["coal", "coal", "coal", "coal", "coal", "coal"],
            "parameter": [
                "investment",
                "FOM",
                "VOM",
                "fuel",
                "investment",
                "discount rate",
            ],
            "value": [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            "unit": [
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
                "unit_atb",
            ],
            "source": [
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
                "source_atb",
            ],
            "further description": ["a", "b", "c", "d", "e", "f"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020],
            "financial_case": ["R&D", "R&D", "R&D", "R&D", "R&D", "R&D"],
            "scenario": [
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
                "Moderate",
            ],
        },
        index=[0, 1, 2, 3, 4, 5],
    )
