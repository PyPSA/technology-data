
configfile: "config.yaml"


rule compile_cost_assumptions:
    input:
        inflation_rate = "inputs/prc_hicp_aind__custom_9928419_spreadsheet.xlsx",
        pypsa_costs = "inputs/costs_PyPSA.csv",
        fraunhofer_costs = "inputs/Fraunhofer_ISE_costs.csv",
        fraunhofer_energy_prices = "inputs/Fraunhofer_ISE_energy_prices.csv",
	    fraunhofer_vehicles_costs = "inputs/Fraunhofer_ISE_vehicles_costs.csv",
        EWG_costs = "inputs/EWG_costs.csv",
        dea_transport = "inputs/energy_transport_data_sheet_dec_2017.xlsx",
        dea_vehicles = "inputs/data_sheets_for_commercial_freight_and_passenger_transport_0.xlsx",
        dea_renewable_fuels = "inputs/data_sheets_for_renewable_fuels.xlsx",
        dea_storage = "inputs/technology_data_catalogue_for_energy_storage.xlsx",
        dea_generation = "inputs/technology_data_for_el_and_dh.xlsx",
        dea_heating = "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
        dea_industrial = "inputs/technology_data_for_industrial_process_heat.xlsx",
        dea_ship = "inputs/data_sheets_for_maritime_commercial_freight_and_passenger_transport.xlsx",
        dea_ccts = "inputs/technology_data_for_carbon_capture_transport_storage.xlsx",
        pnnl_energy_storage = "inputs/pnnl-energy-storage-database.xlsx",
        manual_input = "inputs/manual_input.csv"
    output:
        expand("outputs/costs_{year}.csv", year = config["years"])
    threads: 1
    resources: mem=500
    conda: "environment.yaml"
    script: "scripts/compile_cost_assumptions.py"


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
        EWG = "docu/EWG_LUT_100RE_All_Sectors_Global_Report_2019.pdf"
    output:
        costs = "inputs/EWG_costs.csv",
    threads: 1
    resources: mem=500
    conda: "environment.yaml"
    script: "scripts/convert_pdf_EWG_to_dataframe.py"
