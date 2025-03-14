# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

import pathlib
from shutil import rmtree


configfile: "config.yaml"


rule compile_cost_assumptions:
    input:
        inflation_rate="inputs/Eurostat_inflation_rates.xlsx",
        pypsa_costs="inputs/costs_PyPSA.csv",
        fraunhofer_costs="inputs/Fraunhofer_ISE_costs.csv",
        fraunhofer_energy_prices="inputs/Fraunhofer_ISE_energy_prices.csv",
        fraunhofer_vehicles_costs="inputs/Fraunhofer_ISE_vehicles_costs.csv",
        EWG_costs="inputs/EWG_costs.csv",
        dea_transport="inputs/energy_transport_data_sheet_dec_2017.xlsx",
        dea_vehicles="inputs/data_sheets_for_commercial_freight_and_passenger_transport_0.xlsx",
        dea_renewable_fuels="inputs/data_sheets_for_renewable_fuels.xlsx",
        dea_storage="inputs/technology_data_catalogue_for_energy_storage.xlsx",
        dea_generation="inputs/technology_data_for_el_and_dh.xlsx",
        dea_heating="inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        dea_industrial="inputs/technology_data_for_industrial_process_heat.xlsx",
        dea_ship="inputs/data_sheets_for_maritime_commercial_freight_and_passenger_transport.xlsx",
        dea_ccts="inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
        pnnl_energy_storage="inputs/pnnl-energy-storage-database.xlsx",
        manual_input="inputs/manual_input.csv",
    output:
        expand("outputs/costs_{year}.csv", year=config["years"]),
    threads: 1
    resources:
        mem=500,
    conda:
        "environment.yaml"
    log:
        pathlib.Path("logs", "compile_cost_assumptions.log"),
    script:
        "scripts/compile_cost_assumptions.py"


rule compile_cost_assumptions_usa:
    input:
        cost_files_to_modify=expand("outputs/costs_{year}.csv", year=config["years"]),
        nrel_atb_input_files=expand(
            "inputs/US/atb_e_{year}.parquet",
            year=config["nrel_atb"]["nrel_atb_input_years"],
        ),
        nrel_atb_manual_input_usa="inputs/US/manual_input_usa.csv",
        eur_inflation_rate="inputs/Eurostat_inflation_rates.xlsx",
        nrel_atb_input_discount_rate="inputs/US/discount_rates_usa.csv",
        nrel_atb_input_fuel_costs="inputs/US/fuel_costs_usa.csv",
    output:
        expand("outputs/US/costs_{year}.csv", year=config["years"]),
    threads: 1
    resources:
        mem=500,
    conda:
        "environment.yaml"
    log:
        pathlib.Path("logs", "compile_cost_assumptions_usa.log"),
    script:
        "scripts/compile_cost_assumptions_usa.py"


# rule convert_fraunhofer:
#     input:
#         fraunhofer = "docu/Anhang-Studie-Wege-zu-einem-klimaneutralen-Energiesystem.pdf"
#     output:
#         costs = "inputs/Fraunhofer_ISE_costs.csv",
#         energy_prices = "inputs/Fraunhofer_ISE_energy_prices.csv"
#     threads: 1
#     resources: mem=500
#     conda: "environment.yaml"
#     script: "scripts/convert_pdf_fraunhofer_to_dataframe.py"


rule convert_EWG:
    input:
        EWG="docu/EWG_LUT_100RE_All_Sectors_Global_Report_2019.pdf",
    output:
        costs="inputs/EWG_costs.csv",
    threads: 1
    resources:
        mem=500,
    conda:
        "environment.yaml"
    script:
        "scripts/convert_pdf_EWG_to_dataframe.py"


rule all:
    input:
        rules.compile_cost_assumptions.output,
        rules.compile_cost_assumptions_usa.output,
    default_target: True


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated outputs? [y/N] "
        )
        if do_purge == "y":
            rmtree("outputs/", ignore_errors=True)
            print("Purging generated outputs.")
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")
