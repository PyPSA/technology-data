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
    """
    Fixture to load configuration from config.yaml file.

    Returns
    -------
    dict
        configuration dictionary
    """
    path_config = pathlib.Path(pathlib.Path.cwd(), "config.yaml")
    try:
        with open(path_config) as file:
            config_dict = yaml.safe_load(file)
    except FileNotFoundError:
        pytest.fail(f"Configuration file {path_config} not found.")
    return config_dict


@pytest.fixture(scope="function")
def cost_dataframe():
    """
    Fixture to provide a sample cost dataframe.

    Returns
    -------
    pandas.DataFrame
        sample data for cost
    """
    return pd.DataFrame(
        {
            "technology": ["coal"] * 8 + ["another_tech"],
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
            "value": [1.0] * 8 + [3.0],
            "unit": ["unit"] * 9,
            "source": ["source"] * 9,
            "further description": list("abcdefghi"),
            "currency_year": [2020] * 9,
        }
    )


@pytest.fixture(scope="function")
def atb_cost_dataframe():
    """
    Fixture to provide a sample ATB cost dataframe.

    Returns
    -------
    pandas.DataFrame
        sample data for ATB cost
    """
    return pd.DataFrame(
        {
            "technology": ["coal"] * 6,
            "parameter": [
                "investment",
                "FOM",
                "VOM",
                "fuel",
                "investment",
                "discount rate",
            ],
            "value": [2.0] * 6,
            "unit": ["unit_atb"] * 6,
            "source": ["source_atb"] * 6,
            "further description": list("abcdef"),
            "currency_year": [2020] * 6,
            "financial_case": ["R&D"] * 6,
            "scenario": ["Moderate"] * 6,
        }
    )


@pytest.fixture(scope="function")
def mock_input_data():
    """
    Fixture to provide a sample mock input dataframe.

    Returns
    -------
    pandas.DataFrame
        sample mock input data
    """
    return pd.DataFrame(
        {
            "technology": ["random_tech"] * 50
            + [
                "central air-sourced heat pump",
                "central geothermal-sourced heat pump",
                "central gas boiler",
                "central resistive heater",
                "decentral air-sourced heat pump",
                "decentral gas boiler",
                "decentral ground-sourced heat pump",
            ]
            + ["fuel cell"] * 4,
            "unit": [
                "kW",
                "€",
                " per ",
                " / ",
                " /",
                "J/s",
                "$",
                "₤",
                "MEUR",
                "mio EUR",
                "mill. EUR",
                "1000EUR",
                "k EUR",
                "r/kW",
                "r/GWh",
                "r/GJ",
                " a year",
                "2015EUR",
                "2015-EUR",
                "2020-EUR",
                "EUR2015",
                "EUR-2015",
                "MWe",
                "EUR/MW of total input_e",
                "MWth",
                "MWheat",
                "MWhth",
                "MWhheat",
                "MWH Liquids",
                "MW Liquids",
                "MW Methanol",
                "MW/year FT Liquids/year",
                "MW/year Methanol",
                "MWh FT Liquids/year",
                "MWh methanol",
                "MW/year SNG",
                "MWh SNG",
                "MW SNG",
                "EUR/MWh of total input",
                "EUR/MWeh",
                "% -points of heat loss",
                "FT Liquids Output, MWh/MWh Total Input",
                "MW Ammonia output",
                "MW Ammonia",
                "MWh Ammonia",
                "EUR/MW/y",
                "EUR/MW",
                "EUR/MW/year",
                "EUR/MWh",
                "MW",
                "EUR/MW",
                "EUR/MW/year",
                "EUR/MWh",
                "MW",
                "EUR/MW",
                "EUR/MW/year",
                "EUR/MWh",
                "EUR/MW",
                "EUR/MW/year",
                "EUR/MWh",
                "MW",
            ],
            "value": [1000.0]
            + [1.0] * 7
            + [0.000001] * 3
            + [0.001] * 3
            + [1000.0]
            + [1.0] * 46,
        },
    ).set_index(["technology"])


@pytest.fixture(scope="function")
def mock_output_data():
    """
    Fixture to provide a sample mock output dataframe.

    Returns
    -------
    Callable
        a function to generate mock output data based on the source
    """

    def mock_output(source):
        if source == "dea":
            return pd.DataFrame(
                {
                    "technology": ["random_tech"] * 50
                    + [
                        "central air-sourced heat pump",
                        "central geothermal-sourced heat pump",
                        "central gas boiler",
                        "central resistive heater",
                        "decentral air-sourced heat pump",
                        "decentral gas boiler",
                        "decentral ground-sourced heat pump",
                    ]
                    + ["fuel cell"] * 4,
                    "unit": [
                        "MW",
                        "EUR",
                        "/",
                        "/",
                        "/",
                        "W",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "r/MW",
                        "r/MWh",
                        "r/MWh",
                        "/year",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "MW_e",
                        "EUR/MW_e",
                        "MW_th",
                        "MW_th",
                        "MWh_th",
                        "MWh_th",
                        "MWh_FT",
                        "MW_FT",
                        "MW_MeOH",
                        "MW_FT/year",
                        "MW_MeOH/year",
                        "MWh_FT",
                        "MWh_MeOH",
                        "MW_CH4/year",
                        "MWh_CH4",
                        "MW_CH4",
                        "EUR/MWh_e",
                        "EUR/MW_eh",
                        "MWh_th/MWh_el",
                        "MWh_FT/MWh_H2",
                        "MW_NH3",
                        "MW_NH3",
                        "MWh_NH3",
                        "EUR/MW/year",
                        "EUR/MW",
                        "EUR/MW/year",
                        "EUR/MWh",
                        "MW",
                        "EUR/MW_th",
                        "EUR/MW_th/year",
                        "EUR/MWh_th",
                        "MW_th",
                        "EUR/MW_th",
                        "EUR/MW_th/year",
                        "EUR/MWh_th",
                        "EUR/MW_e",
                        "EUR/MW_e/year",
                        "EUR/MWh_e",
                        "MW_e",
                    ],
                    "value": [1.0] * 6
                    + [
                        0.8917822267802202,
                        1.177107611177814,
                    ]
                    + [1.0] * 7
                    + [3.6]
                    + [1.0] * 45,
                },
            ).set_index(["technology"])
        else:
            return pd.DataFrame(
                {
                    "technology": ["random_tech"] * 50
                    + [
                        "central air-sourced heat pump",
                        "central geothermal-sourced heat pump",
                        "central gas boiler",
                        "central resistive heater",
                        "decentral air-sourced heat pump",
                        "decentral gas boiler",
                        "decentral ground-sourced heat pump",
                    ]
                    + ["fuel cell"] * 4,
                    "unit": [
                        "MW",
                        "EUR",
                        "/",
                        "/",
                        "/",
                        "W",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "r/MW",
                        "r/MWh",
                        "r/MWh",
                        "/year",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "EUR",
                        "MW_e",
                        "EUR/MW_e",
                        "MW_th",
                        "MW_th",
                        "MWh_th",
                        "MWh_th",
                        "MWh_FT",
                        "MW_FT",
                        "MW_MeOH",
                        "MW_FT/year",
                        "MW_MeOH/year",
                        "MWh_FT",
                        "MWh_MeOH",
                        "MW_CH4/year",
                        "MWh_CH4",
                        "MW_CH4",
                        "EUR/MWh_e",
                        "EUR/MW_eh",
                        "MWh_th/MWh_el",
                        "MWh_FT/MWh_H2",
                        "MW_NH3",
                        "MW_NH3",
                        "MWh_NH3",
                        "EUR/MW/year",
                        "EUR/MW",
                        "EUR/MW/year",
                        "EUR/MWh",
                        "MW",
                        "EUR/MW",
                        "EUR/MW/year",
                        "EUR/MWh",
                        "MW",
                        "EUR/MW",
                        "EUR/MW/year",
                        "EUR/MWh",
                        "EUR/MW",
                        "EUR/MW/year",
                        "EUR/MWh",
                        "MW",
                    ],
                    "value": [1.0] * 6
                    + [
                        0.8917822267802202,
                        1.177107611177814,
                    ]
                    + [1.0] * 7
                    + [3.6]
                    + [1.0] * 45,
                },
            ).set_index(["technology"])

    return mock_output
