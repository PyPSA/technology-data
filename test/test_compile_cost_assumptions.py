# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import pathlib
import sys

sys.path.append("./scripts")

from compile_cost_assumptions import get_excel_sheets, get_data_from_DEA

path_cwd = pathlib.Path.cwd()

snakemake_input_dictionary = {
    "inflation_rate": "inputs/prc_hicp_aind__custom_9928419_spreadsheet.xlsx",
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


def test_get_excel_sheets():
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

def test_get_data_from_DEA(config):
    excel_files = [
        v for k, v in snakemake_input_dictionary.items() if "dea" in k.casefold()
    ]
    input_dea_files_dict = get_excel_sheets(excel_files)
    output_df = get_data_from_DEA(input_dea_files_dict, expectation=config["expectation"])
    print(output_df.shape)
    #print(output_df)
    assert False
