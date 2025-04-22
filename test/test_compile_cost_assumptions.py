# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import copy
import pathlib

import numpy as np
import pandas as pd
import pytest

from scripts.compile_cost_assumptions import (
    add_carbon_capture,
    add_description,
    add_gas_storage,
    annuity,
    clean_up_units,
    convert_units,
    dea_sheet_names,
    geometric_series,
    get_data_from_DEA,
    get_dea_vehicle_data,
    get_excel_sheets,
    get_sheet_location,
    rename_pypsa_old,
    set_round_trip_efficiency,
    set_specify_assumptions,
)

path_cwd = pathlib.Path.cwd()

snakemake_input_dictionary = {
    "inflation_rate": "inputs/Eurostat_inflation_rate.xlsx",
    "pypsa_costs": "inputs/costs_PyPSA.csv",
    "fraunhofer_costs": "inputs/Fraunhofer_ISE_costs.csv",
    "fraunhofer_energy_prices": "inputs/Fraunhofer_ISE_energy_prices.csv",
    "fraunhofer_vehicles_costs": "inputs/Fraunhofer_ISE_vehicles_costs.csv",
    "EWG_costs": "inputs/EWG_costs.csv",
    "dea_transport": "inputs/energy_transport_data_sheet_dec_2017.xlsx",
    "dea_vehicles": "inputs/data_sheets_for_commercial_freight_and_passenger_transport_0.xlsx",
    "dea_renewable_fuels": "inputs/data_sheets_for_renewable_fuels.xlsx",
    "dea_storage": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
    "dea_generation": "inputs/technology_data_for_el_and_dh.xlsx",
    "dea_heating": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
    "dea_industrial": "inputs/technology_data_for_industrial_process_heat.xlsx",
    "dea_ship": "inputs/data_sheets_for_maritime_commercial_freight_and_passenger_transport.xlsx",
    "dea_ccts": "inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
    "pnnl_energy_storage": "inputs/pnnl-energy-storage-database.xlsx",
    "manual_input": "inputs/manual_input.csv",
}


@pytest.mark.parametrize("source", ["", "dea"])
def test_clean_up_units(mock_input_data, mock_output_data, source):
    """
    The test verifies what is returned by clean_up_units.
    """
    expected_df = mock_output_data(source)
    output_df = clean_up_units(
        mock_input_data.copy(deep=True), value_column="value", source=source
    )
    comparison_df = output_df.compare(expected_df)
    print(comparison_df)
    assert comparison_df.empty


def test_get_excel_sheets():
    """
    The test verifies what is returned by get_excel_sheets.
    """
    reference_output_dictionary = {
        "inputs/energy_transport_data_sheet_dec_2017.xlsx": 16,
        "inputs/data_sheets_for_commercial_freight_and_passenger_transport_0.xlsx": 19,
        "inputs/data_sheets_for_renewable_fuels.xlsx": 45,
        "inputs/technology_data_catalogue_for_energy_storage.xlsx": 15,
        "inputs/technology_data_for_el_and_dh.xlsx": 72,
        "inputs/technologydatafor_heating_installations_marts_2018.xlsx": 29,
        "inputs/technology_data_for_industrial_process_heat.xlsx": 32,
        "inputs/data_sheets_for_maritime_commercial_freight_and_passenger_transport.xlsx": 22,
        "inputs/technology_data_for_carbon_capture_transport_storage.xlsx": 31,
    }
    excel_files = [
        v for k, v in snakemake_input_dictionary.items() if "dea" in k.casefold()
    ]
    output_dict = get_excel_sheets(excel_files)
    comparison_dictionary = {}
    for key, value in output_dict.items():
        comparison_dictionary[key] = len(value)
    assert reference_output_dictionary == comparison_dictionary


def test_get_sheet_location():
    """
    The test verifies what is returned by get_sheet_location.
    """
    reference_output_dictionary = {
        "onwind": "inputs/technology_data_for_el_and_dh.xlsx",
        "offwind": "inputs/technology_data_for_el_and_dh.xlsx",
        "solar-utility": "inputs/technology_data_for_el_and_dh.xlsx",
        "solar-utility single-axis tracking": "inputs/technology_data_for_el_and_dh.xlsx",
        "solar-rooftop residential": "inputs/technology_data_for_el_and_dh.xlsx",
        "solar-rooftop commercial": "inputs/technology_data_for_el_and_dh.xlsx",
        "OCGT": "inputs/technology_data_for_el_and_dh.xlsx",
        "CCGT": "inputs/technology_data_for_el_and_dh.xlsx",
        "oil": "inputs/technology_data_for_el_and_dh.xlsx",
        "biomass CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "biomass EOP": "inputs/technology_data_for_el_and_dh.xlsx",
        "biomass HOP": "inputs/technology_data_for_el_and_dh.xlsx",
        "central coal CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "central gas CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "central gas CHP CC": "inputs/technology_data_for_el_and_dh.xlsx",
        "central solid biomass CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "central solid biomass CHP CC": "inputs/technology_data_for_el_and_dh.xlsx",
        "central solid biomass CHP powerboost CC": "inputs/technology_data_for_el_and_dh.xlsx",
        "central air-sourced heat pump": "inputs/technology_data_for_el_and_dh.xlsx",
        "central geothermal heat source": "inputs/technology_data_for_el_and_dh.xlsx",
        "central excess-heat-sourced heat pump": "inputs/technology_data_for_el_and_dh.xlsx",
        "central ground-sourced heat pump": "inputs/technology_data_for_el_and_dh.xlsx",
        "central resistive heater": "inputs/technology_data_for_el_and_dh.xlsx",
        "central gas boiler": "inputs/technology_data_for_el_and_dh.xlsx",
        "decentral gas boiler": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        "direct firing gas": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "direct firing gas CC": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "direct firing solid fuels": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "direct firing solid fuels CC": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "decentral ground-sourced heat pump": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        "decentral air-sourced heat pump": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        "central water pit storage": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "central water tank storage": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "decentral water tank storage": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "fuel cell": "inputs/technology_data_for_el_and_dh.xlsx",
        "hydrogen storage underground": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "hydrogen storage tank type 1 including compressor": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "micro CHP": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        "biogas": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "biogas CC": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "biogas upgrading": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "battery": "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        "industrial heat pump medium temperature": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "industrial heat pump high temperature": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "electric boiler steam": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "gas boiler steam": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "solid biomass boiler steam": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "solid biomass boiler steam CC": "inputs/technology_data_for_industrial_process_heat.xlsx",
        "biomass boiler": "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        "electrolysis": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "direct air capture": "inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
        "biomass CHP capture": "inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
        "cement capture": "inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
        "BioSNG": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "BtL": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "biomass-to-methanol": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "biogas plus hydrogen": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "methanolisation": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "Fischer-Tropsch": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "central hydrogen CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "Haber-Bosch": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "air separation unit": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "waste CHP": "inputs/technology_data_for_el_and_dh.xlsx",
        "waste CHP CC": "inputs/technology_data_for_el_and_dh.xlsx",
        "biochar pyrolysis": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "electrolysis small": "inputs/data_sheets_for_renewable_fuels.xlsx",
        "random tech": "Sheet not found",
    }
    dea_sheet_names_dict = copy.deepcopy(dea_sheet_names)
    dea_sheet_names_dict["random tech"] = "random sheet"
    excel_files = [
        v for k, v in snakemake_input_dictionary.items() if "dea" in k.casefold()
    ]
    output_dict = get_excel_sheets(excel_files)
    sheet_location_dictionary = {}
    for tech, dea_tech in dea_sheet_names_dict.items():
        technology_location = get_sheet_location(
            tech, dea_sheet_names_dict, output_dict
        )
        sheet_location_dictionary[tech] = technology_location
    assert sheet_location_dictionary == reference_output_dictionary


def test_get_data_from_dea(config):
    """
    The test verifies what is returned by get_data_from_DEA.
    """
    reference_output_dictionary = {
        "onwind": (4, 9),
        "offwind": (4, 9),
        "solar-utility": (4, 9),
        "solar-utility single-axis tracking": (4, 9),
        "solar-rooftop residential": (4, 9),
        "solar-rooftop commercial": (4, 9),
        "OCGT": (8, 9),
        "CCGT": (8, 9),
        "oil": (8, 9),
        "biomass CHP": (13, 9),
        "biomass EOP": (13, 9),
        "biomass HOP": (9, 9),
        "central coal CHP": (7, 9),
        "central gas CHP": (8, 9),
        "central gas CHP CC": (8, 9),
        "central solid biomass CHP": (13, 9),
        "central solid biomass CHP CC": (13, 9),
        "central solid biomass CHP powerboost CC": (13, 9),
        "central air-sourced heat pump": (6, 9),
        "central geothermal heat source": (11, 9),
        "central excess-heat-sourced heat pump": (6, 9),
        "central ground-sourced heat pump": (5, 9),
        "central resistive heater": (7, 9),
        "central gas boiler": (6, 9),
        "decentral gas boiler": (8, 9),
        "direct firing gas": (7, 9),
        "direct firing gas CC": (7, 9),
        "direct firing solid fuels": (7, 9),
        "direct firing solid fuels CC": (7, 9),
        "decentral ground-sourced heat pump": (9, 9),
        "decentral air-sourced heat pump": (9, 9),
        "central water pit storage": (10, 9),
        "central water tank storage": (10, 9),
        "decentral water tank storage": (9, 9),
        "fuel cell": (8, 9),
        "hydrogen storage underground": (10, 9),
        "hydrogen storage tank type 1 including compressor": (10, 9),
        "micro CHP": (8, 9),
        "biogas": (4, 9),
        "biogas CC": (4, 9),
        "biogas upgrading": (4, 9),
        "battery": (13, 9),
        "industrial heat pump medium temperature": (6, 9),
        "industrial heat pump high temperature": (6, 9),
        "electric boiler steam": (7, 9),
        "gas boiler steam": (7, 9),
        "solid biomass boiler steam": (7, 9),
        "solid biomass boiler steam CC": (7, 9),
        "biomass boiler": (6, 9),
        "electrolysis": (7, 9),
        "direct air capture": (8, 9),
        "biomass CHP capture": (10, 9),
        "cement capture": (10, 9),
        "BioSNG": (6, 9),
        "BtL": (6, 9),
        "biomass-to-methanol": (8, 9),
        "biogas plus hydrogen": (6, 9),
        "methanolisation": (7, 9),
        "Fischer-Tropsch": (6, 9),
        "central hydrogen CHP": (8, 9),
        "Haber-Bosch": (10, 9),
        "air separation unit": (7, 9),
        "waste CHP": (16, 9),
        "waste CHP CC": (16, 9),
        "biochar pyrolysis": (7, 9),
        "electrolysis small": (7, 9),
        "random tech": (0, 0),
    }
    excel_files = [
        v for k, v in snakemake_input_dictionary.items() if "dea" in k.casefold()
    ]
    input_dea_files_dict = get_excel_sheets(excel_files)
    dea_sheet_names_dict = copy.deepcopy(dea_sheet_names)
    dea_sheet_names_dict["random tech"] = "random sheet"
    output_dictionary = get_data_from_DEA(
        config["years"],
        dea_sheet_names_dict,
        input_dea_files_dict,
        expectation=config["expectation"],
    )
    comparison_dictionary = {}
    for key, value in output_dictionary.items():
        comparison_dictionary[key] = value.shape
    assert comparison_dictionary == reference_output_dictionary


def test_set_specify_assumptions():
    """
    The test verifies what is returned by set_specify_assumptions.
    """
    input_df = pd.DataFrame(
        {
            "technology": [
                "central resistive heater",
                "decentral gas boiler",
                "decentral gas boiler",
                "decentral gas boiler",
                "biogas upgrading",
                "biogas upgrading",
                "solar-rooftop",
                "heat pump",
            ],
            "parameter": [
                "Nominal investment, 400/690 V; 1-5 MW",
                "Heat efficiency, annual average, net",
                "Possible additional specific investment",
                "Technical lifetime",
                "investment",
                "Technical lifetime",
                "PV module conversion efficiency [p.u.]",
                "Heat efficiency, annual average, net, radiators",
            ],
            "2020": [1.0] * 8,
            "source": ["source"] * 8,
            "unit": ["unit"] * 8,
        }
    ).set_index(["technology", "parameter"])

    reference_output_df = pd.DataFrame(
        {
            "technology": [
                "biogas upgrading",
                "biogas upgrading",
                "decentral gas boiler",
                "decentral gas boiler connection",
                "decentral gas boiler connection",
                "heat pump",
            ],
            "parameter": [
                "Technical lifetime",
                "investment (upgrading, methane redution and grid injection)",
                "Technical lifetime",
                "Possible additional specific investment",
                "Technical lifetime",
                "Heat efficiency, annual average, net, radiators, existing one family house",
            ],
            "2020": [1.0, 1.0, 1.0, 1.0, 50.0, 1.0],
            "source": ["source"] * 6,
            "unit": ["unit"] * 6,
        }
    )
    list_of_years = ["2020"]
    output_df = set_specify_assumptions(list_of_years, input_df)
    output_df = output_df.reset_index(drop=False).rename(
        columns={"level_0": "technology", "level_1": "parameter"}
    )
    comparison_df = output_df.compare(reference_output_df)
    assert comparison_df.empty


def test_set_round_trip_efficiency():
    """
    The test verifies what is returned by set_round_trip_efficiency.
    """
    input_df = pd.DataFrame(
        {
            "technology": [
                "hydrogen storage underground",
                "hydrogen storage tank type 1 including compressor",
            ]
            + ["battery"] * 5,
            "parameter": [
                "Round trip efficiency",
                "Round trip efficiency",
                "Round trip efficiency DC",
                "Output capacity expansion cost",
                "Technical lifetime",
                "Fixed O&M",
                "Energy storage expansion cost",
            ],
            "2020": [1.0] * 7,
            "source": ["source"] * 7,
            "unit": ["unit"] * 7,
        }
    ).set_index(["technology", "parameter"])

    reference_output_df = pd.DataFrame(
        {
            "technology": ["battery inverter"] * 4
            + ["battery storage"] * 2
            + [
                "hydrogen storage tank type 1 including compressor",
                "hydrogen storage underground",
            ],
            "parameter": [
                "Fixed O&M",
                "Output capacity expansion cost investment",
                "Round trip efficiency DC",
                "Technical lifetime",
                "Energy storage expansion cost investment",
                "Technical lifetime",
                "Round trip efficiency",
                "Round trip efficiency",
            ],
            "2020": [
                1.0,
                1.0,
                1.0,
                10.0,
                1.0,
                1.0,
                100.0,
                100.0,
            ],
            "source": ["source"] * 3 + ["source, Note K."] + ["source"] * 4,
            "unit": ["unit"] * 8,
        }
    )
    list_of_years = ["2020"]
    output_df = set_round_trip_efficiency(list_of_years, input_df)
    output_df = output_df.reset_index(drop=False).rename(
        columns={"level_0": "technology", "level_1": "parameter"}
    )
    comparison_df = output_df.compare(reference_output_df)
    assert comparison_df.empty


def test_add_description():
    """
    The test verifies what is returned by add_description.
    """
    technology_series = pd.Series(
        [
            "offwind",
            "technology_name_2",
        ]
    )

    parameter_series = pd.Series(
        [
            "investment",
            "parameter_name_2",
        ]
    )
    input_df = pd.DataFrame(
        {
            "2020": [1.0, 1.0],
            "source": [
                "source",
                "source",
            ],
            "unit": [
                "unit",
                "unit",
            ],
            "further description": [
                "text",
                "text",
            ],
        }
    ).set_index([technology_series, parameter_series])

    reference_output_df = pd.DataFrame(
        {
            "technology": ["offwind", "technology_name_2"],
            "parameter": ["investment", "parameter_name_2"],
            "2020": [1.0, 1.0],
            "unit": [
                "unit",
                "unit",
            ],
            "source": [
                "source",
                "source",
            ],
            "further description": [
                "21 Offshore turbines:  text grid connection costs subtracted from investment costs",
                ":  text",
            ],
        }
    )
    list_of_years = ["2020"]
    output_df = add_description(list_of_years, input_df).reset_index()
    comparison_df = output_df.compare(reference_output_df)
    assert comparison_df.empty


def test_convert_units(config):
    """
    The test verifies what is returned by convert_units.
    """
    technology_series = pd.Series(
        [
            "technology_name_1",
            "technology_name_2",
        ]
    )

    parameter_series = pd.Series(
        [
            "efficiency",
            "efficiency-heat",
        ]
    )
    input_df = pd.DataFrame(
        {
            "2020": [100.0, 100.0],
            "2025": [100.0, 100.0],
            "2030": [100.0, 100.0],
            "2035": [100.0, 100.0],
            "2040": [100.0, 100.0],
            "2045": [100.0, 100.0],
            "2050": [100.0, 100.0],
            "source": [
                "source",
                "source",
            ],
            "unit": [
                "unit",
                "unit",
            ],
            "further description": [
                "text",
                "text",
            ],
        }
    ).set_index([technology_series, parameter_series])

    reference_output_df = pd.DataFrame(
        {
            "technology": ["technology_name_1", "technology_name_2"],
            "parameter": ["efficiency", "efficiency-heat"],
            "2020": [1.0, 1.0],
            "2025": [1.0, 1.0],
            "2030": [1.0, 1.0],
            "2035": [1.0, 1.0],
            "2040": [1.0, 1.0],
            "2045": [1.0, 1.0],
            "2050": [1.0, 1.0],
            "source": [
                "source",
                "source",
            ],
            "unit": [
                "per unit",
                "per unit",
            ],
            "further description": [
                "text",
                "text",
            ],
        }
    )
    list_of_years = [str(x) for x in config["years"]]
    output_df = (
        convert_units(list_of_years, input_df)
        .reset_index()
        .rename(columns={"level_0": "technology", "level_1": "parameter"})
    )
    comparison_df = output_df.compare(reference_output_df)
    assert comparison_df.empty


@pytest.mark.parametrize(
    "discount_rate_value, expected_annuity", [(0.07, 1.07), (-0.1, 1.0)]
)
def test_annuity(discount_rate_value, expected_annuity):
    """
    The test verifies what is returned by annuity.
    """
    assert annuity(n=1.0, r=discount_rate_value) == expected_annuity


@pytest.mark.parametrize(
    "nom_val, den_val, n_terms, start_val, expected_val",
    [(1.0, 6.5, 3, 0, 1.18), (1.0, 2.0, 3, 0, 1.75)],
)
def test_geometric_series(nom_val, den_val, n_terms, start_val, expected_val):
    """
    The test verifies what is returned by annuity.
    """
    assert (
        np.round(geometric_series(nom_val, den_val, n_terms, start_val), 2)
        == expected_val
    )


def test_add_gas_storage(config):
    """
    The test verifies what is returned by add_gas_storage.
    """
    input_file_path = pathlib.Path(
        path_cwd, "inputs", "technology_data_catalogue_for_energy_storage.xlsx"
    )
    list_of_years = ["2020"]
    technology_series = pd.Series(
        [
            "gas_storage",
            "gas_storage",
            "gas storage charger",
            "gas storage discharger",
            "gas_storage",
        ]
    )
    parameter_series = pd.Series(
        [
            "investment",
            "lifetime",
            "investment",
            "investment",
            "FOM",
        ]
    )
    technology_dataframe = pd.DataFrame(
        {
            "2020": [np.nan] * 5,
        }
    ).set_index([technology_series, parameter_series])
    output_df = add_gas_storage(input_file_path, list_of_years, technology_dataframe)
    cleaned_df = (
        output_df.dropna(how="all")
        .reset_index(drop=False)
        .rename(columns={"level_0": "technology", "level_1": "parameter"})
    )

    reference_output_df = pd.DataFrame(
        {
            "technology": [
                "gas storage charger",
                "gas storage discharger",
                "gas storage",
                "gas storage",
                "gas storage",
            ],
            "parameter": ["investment"] * 3 + ["lifetime", "FOM"],
            "2020": [
                14.33885157711495,
                4.77961719237165,
                0.03290254335105841,
                100.0,
                3.5918748566291763,
            ],
            "source": ["Danish Energy Agency"] * 3
            + ["TODO no source"]
            + ["Danish Energy Agency"],
            "further description": [
                "150 Underground Storage of Gas, Process equipment (units converted)"
            ]
            * 2
            + [
                "150 Underground Storage of Gas, Establishment of one cavern (units converted)"
            ]
            + [
                "estimation: most underground storage are already build, they do have a long lifetime"
            ]
            + [
                "150 Underground Storage of Gas, Operation and Maintenance, salt cavern (units converted)"
            ],
            "unit": ["EUR/kW"] * 2 + ["EUR/kWh"] + ["years", "%"],
            "currency_year": [2015.0, np.nan, 2015.0, np.nan, np.nan],
        }
    )
    comparison_df = cleaned_df.compare(reference_output_df)
    assert comparison_df.empty


def test_add_carbon_capture(config):
    """
    The test verifies what is returned by add_carbon_capture.
    """
    list_of_years = ["2020"]
    technology_series = pd.Series(
        ["cement capture"] * 9
        + ["biomass CHP capture"] * 9
        + ["direct air capture"] * 9,
    )
    input_parameter_series = pd.Series(
        [
            "Ax) CO2 capture rate, net",
            "Specific investment",
            "Fixed O&M",
            "C2) Eletricity input ",
            "C1) Heat  input ",
            "C1) Heat out ",
            "CO₂ compression and dehydration - Electricity input",
            "CO₂ compression and dehydration - Heat out",
            "lifetime",
        ]
        * 3,
    )
    technology_dataframe = pd.DataFrame(
        {
            "2020": [50, 100, 10, 40, 90, 9, 30, 80, 10] * 3,
            "source": ["source"] * 27,
        }
    ).set_index([technology_series, input_parameter_series])

    output_parameter_series = pd.Series(
        [
            "capture_rate",
            "investment",
            "FOM",
            "electricity-input",
            "heat-input",
            "heat-output",
            "compression-electricity-input",
            "compression-heat-output",
            "lifetime",
        ]
        * 3,
    )
    new_technology_dataframe = pd.DataFrame(
        {
            "2020": [np.nan] * 27,
            "source": ["source"] * 27,
        }
    ).set_index([technology_series, output_parameter_series])

    output_df = add_carbon_capture(
        list_of_years, dea_sheet_names, new_technology_dataframe, technology_dataframe
    )

    assert output_df.loc["cement capture", "capture_rate"].equals(
        pd.Series(
            [0.5, "source", "per unit", "401.c Post comb - Cement kiln"],
            name="(cement capture, capture_rate)",
            index=["2020", "source", "unit", "further description"],
        )
    )
    assert output_df.loc["cement capture", "FOM"].equals(
        pd.Series(
            [10.0, "source", "%/year", "401.c Post comb - Cement kiln"],
            name="(cement capture, FOM)",
            index=["2020", "source", "unit", "further description"],
        )
    )
    assert output_df.loc["cement capture", "compression-heat-output"].equals(
        pd.Series(
            [80.0, "source", "MWh/tCO2", "401.c Post comb - Cement kiln"],
            name="(cement capture, compression-heat-output)",
            index=["2020", "source", "unit", "further description"],
        )
    )
    assert output_df.loc["cement capture", "lifetime"].equals(
        pd.Series(
            [np.nan, "source", np.nan, "401.c Post comb - Cement kiln"],
            name="(cement capture, lifetime)",
            index=["2020", "source", "unit", "further description"],
        )
    )
    assert output_df.loc["biomass CHP capture", "electricity-input"].equals(
        pd.Series(
            [40.0, "source", "MWh/tCO2", "401.a Post comb - small CHP"],
            name="(cement capture, lifetime)",
            index=["2020", "source", "unit", "further description"],
        )
    )
    assert output_df.loc["direct air capture", "investment"].equals(
        pd.Series(
            [100000000.0, "source", "EUR/(tCO2/h)", "403.a Direct air capture"],
            name="(cement capture, lifetime)",
            index=["2020", "source", "unit", "further description"],
        )
    )


def test_get_dea_vehicle_data(config):
    """
    The test verifies what is returned by get_dea_vehicle_data.
    """
    df = get_dea_vehicle_data(
        snakemake_input_dictionary["dea_vehicles"], ["2020"], pd.DataFrame()
    )
    assert df.shape == (90, 5)
    assert sorted(list(df.columns)) == sorted(
        ["2020", "unit", "source", "currency_year", "further description"]
    )


def test_rename_pypsa_old():
    """
    The test verifies what is returned by rename_pypsa_old.
    """
    data = {
        ("decentral water tank storage", "investment"): [46.8],
        ("central CHP", "investment"): [3000],
        ("hydrogen underground storage", "investment"): [2000],
        ("retrofitting I", "investment"): [1000],
        ("retrofitting II", "investment"): [1500],
    }
    index = pd.MultiIndex.from_tuples(data.keys())
    cost_dataframe_pypsa = pd.DataFrame(data.values(), index=index, columns=["value"])
    expected_data = {
        ("decentral water tank storage", "investment"): [1.0],
        ("central gas CHP", "investment"): [3000],
        ("hydrogen storage underground", "investment"): [2000],
    }
    expected_index = pd.MultiIndex.from_tuples(expected_data.keys())
    expected_df = pd.DataFrame(
        expected_data.values(), index=expected_index, columns=["value"]
    )
    expected_df.loc[("decentral water tank storage", "investment"), "unit"] = "EUR/kWh"
    output_df = rename_pypsa_old(cost_dataframe_pypsa)
    comparison_df = output_df.compare(expected_df)
    assert comparison_df.empty
