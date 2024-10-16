#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import pathlib
import numpy as np
from _helpers import mock_snakemake


def filter_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep):

    atb_input_df = pd.read_parquet(input_file_path)
    list_core_metric_parameter_to_keep = [str(x).casefold() for x in list_core_metric_parameter_to_keep]
    year_string = str(year).casefold()

    # --> select columns
    for column_name in list_columns_to_keep:
        if column_name not in atb_input_df:
            print("missing column", column_name)
            atb_input_df[column_name] = "not_available"
        else:
            pass
    atb_input_df = atb_input_df.loc[:, list_columns_to_keep]

    # --> select rows based on core_metric_parameter
    # Note: we do not apply any selection on core_metric_case and scenario.
    # Such selection can be done in the model config (e.g. PyPSA-Earth)
    atb_input_df = atb_input_df.loc[atb_input_df["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]

    # --> select rows based on core_metric_variable
    if input_file_path.name == "atb_e_2022.parquet":
        # Note: 2020 data are fetched from the input file atb_e_2022.
        atb_input_df = atb_input_df.loc[atb_input_df["core_metric_variable"].astype(str).str.casefold() == year_string]
    elif input_file_path.name == "atb_e_2024.parquet":
        # Note: 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
        atb_input_df = atb_input_df.loc[atb_input_df["core_metric_variable"].astype(str).str.casefold() == year_string]
    else:
        raise Exception("{} - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet".format(input_file_path.stem))

    # --> drop duplicated rows
    atb_input_df = atb_input_df.drop_duplicates(keep="first")

    # --> remove technology AEO
    atb_input_df = atb_input_df.query("technology.str.casefold() != 'aeo'")

    return atb_input_df


def calculate_fom_percentage(x, dataframe):

    # Note: for technologies as Coal Retrofit or Natural Gas Retrofit,
    # the Fixed O&M is normalized by Additional OCC. Else, the Fixed O&M is
    # normalized by the CAPEX

    if x["core_metric_parameter"].casefold() == "fixed o&m":
        if "retrofit" in x["technology"].casefold():
            fom_perc_value = x.value / dataframe.query("atb_year == @x.atb_year & core_metric_parameter.str.casefold() == 'additional occ' & core_metric_case == @x.core_metric_case & core_metric_variable == @x.core_metric_variable & technology == @x.technology & technology_alias == @x.technology_alias & display_name == @x.display_name & scenario == @x.scenario")["value"]*100.0
        else:
            fom_perc_value = x.value / dataframe.query("atb_year == @x.atb_year & core_metric_parameter.str.casefold() == 'capex' & core_metric_case == @x.core_metric_case & core_metric_variable == @x.core_metric_variable & technology == @x.technology & technology_alias == @x.technology_alias & display_name == @x.display_name & scenario == @x.scenario")["value"]*100.0
        return round(fom_perc_value.values[0], 2)
    else:
        return x.value


def replace_value_name(dataframe, conversion_dict, column_name):
    dataframe[column_name] = dataframe[column_name].replace(conversion_dict)
    return dataframe


def concatenate_columns(dataframe, column_name, column_name_list):
    dataframe[column_name] = dataframe[column_name_list].astype(str).agg('_'.join, axis=1)
    dataframe[column_name] = dataframe[column_name].str.replace(" ", "_")
    return dataframe


def repeat_values(dataframe, values_list, column_name, value_to_filter, n_repeat):
    dataframe_to_repeat = dataframe.loc[dataframe[column_name].astype(str).str.casefold() == str(value_to_filter).casefold()]
    dataframe_to_keep = dataframe.loc[dataframe[column_name].astype(str).str.casefold() != str(value_to_filter).casefold()]
    dataframe_repeated = pd.DataFrame(np.repeat(dataframe_to_repeat.values, n_repeat, axis=0), columns=dataframe_to_repeat.columns)
    for i_rep in range(n_repeat):
        dataframe_repeated.loc[i_rep::n_repeat, column_name] = values_list[i_rep]
    return pd.concat([dataframe_to_keep, dataframe_repeated], ignore_index=True)


def pre_process_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_source):

    # read inputs and filter relevant columns and relevant core_metric_variables (i.e. the years to consider)

    atb_input_df = filter_input_file(pathlib.Path(input_file_path), year, list_columns_to_keep, list_core_metric_parameter_to_keep)

    # Normalize Fixed O&M by CAPEX (or Additional OCC for retrofit technologies)
    atb_input_df["value"] = atb_input_df.apply(lambda x: calculate_fom_percentage(x, atb_input_df), axis=1)

    # Modify the unit of the normalized Fixed O&M to %-yr
    atb_input_df["units"] = atb_input_df.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)

    # Create a new column technology_alias_detail from the concatenation of the technology_alias and techdetail columns
    concatenate_columns(atb_input_df, "technology_alias_detail", ["technology_alias", "techdetail"])

    # For WACC_Real, Hydropower provides just one row. This one row has to be associated with all types of hydropower
    # The one WACC_Real value is replicated such that each type of hydropower has its own WACC_Real entry
    atb_input_df = repeat_values(atb_input_df, ["Hydropower_NPD1", "Hydropower_NSD1"], "technology_alias_detail", "Hydropower_*", 2)

    # Replace technology_alias_detail with PyPSA technology names
    technology_conversion_dict_atb = {
        "Coal_Coal-new": "coal",
        "Coal_newAvgCF2ndGen": "coal",
        "Coal_*": "coal",
        "Natural_Gas_NG_Combustion_Turbine_(F-Frame)": "CCGT",
        "Natural_Gas_CTAvgCF": "CCGT",
        "Natural_Gas_*": "CCGT",
        "Hydropower_NPD1": "hydro",
        "Hydropower_NSD1": "ror",
        "Pumped_Storage_Hydropower_NatlClass1": "PHS",
        "Pumped_Storage_Hydropower_*": "PHS",
        "Nuclear_Large": "nuclear",
        "Nuclear_Nuclear": "nuclear",
        "Nuclear_*": "nuclear",
        "Geothermal_HydroFlash": "geothermal",
        "Geothermal_*": "geothermal",
        "Land-Based_Wind_Class1": "onwind",
        "Land-Based_Wind_*": "onwind",
        "Offshore_Wind_Class1": "offwind",
        "Offshore_Wind_*": "offwind",
        "Utility_PV_Class1": "solar-utility",
        "Utility_PV_*": "solar-utility",
        "Commercial_PV_Class1": "solar-rooftop",
        "Commercial_PV_*": "solar-rooftop",
        "Utility-Scale_Battery_Storage_4Hr_Battery_Storage": "battery storage",
        "Biopower_Dedicated": "biomass",
        "Biopower_*": "biomass",
        "CSP_Class2": "csp-tower",
        "CSP_*": "csp-tower",
    }
    replace_value_name(atb_input_df, technology_conversion_dict_atb, "technology_alias_detail")

    # Add source column
    atb_input_df["source"] = nrel_source

    # Rename columns and consider just columns used in PyPSA
    column_rename_dict = {
        "technology_alias_detail": "technology",
        "core_metric_parameter": "parameter",
        "units": "unit",
        "display_name": "further description",
        "atb_year": "currency_year",
        "core_metric_case": "financial_case"
    }

    tuple_output_columns_to_keep = ("technology_alias_detail", "core_metric_parameter", "value", "units", "source", "display_name", "atb_year", "scenario", "core_metric_case", "core_metric_variable", "tax_credit_case")
    atb_input_df = atb_input_df.loc[:, tuple_output_columns_to_keep].rename(columns=column_rename_dict)

    # Replace parameter with PyPSA cost parameter names
    parameter_conversion_dict = {
        "CAPEX": "investment",
        "Fixed O&M": "FOM",
        "Variable O&M": "VOM",
        "Fuel": "fuel",
        "Additional OCC": "investment",
        "WACC Real": "discount rate"
    }
    replace_value_name(atb_input_df, parameter_conversion_dict, "parameter")

    return atb_input_df

    #new_atb_input_df_2022 = atb_input_df_2022.set_index(["technology", "parameter"])
    #new_atb_input_df_2024 = atb_input_df_2024.set_index(["technology", "parameter"])
    #return new_atb_input_df_2022, new_atb_input_df_2024


def update_cost_values(cost_dataframe, atb_dataframe):
    pass


if __name__ == "__main__":
    if 'snakemake' not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']
    input_file_list_atb = snakemake.input.nrel_atb_input_files
    nrel_atb_columns_to_keep = snakemake.config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = snakemake.config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_source_link = snakemake.config["nrel_atb"]["nrel_atb_source_link"]

    cost_atb_2022, cost_atb_2024 = pre_process_input_file(input_file_list_atb, year_list, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)

    for i, input_file_name in enumerate(snakemake.input.cost_files_to_modify):
        input_file_path = pathlib.Path(input_file_name)
        cost_df = pd.read_csv(input_file_path).reset_index()
        if input_file_path.name == "costs_2020.csv":
            update_cost_values(cost_df, cost_atb_2022)
        elif input_file_path.name in ["costs_2025.csv", "costs_2030.csv", "costs_2035.csv", "costs_2040.csv", "costs_2045.csv", "costs_2050.csv"]:
            #filtered_cost_atb_2024 =
            update_cost_values(cost_df, filtered_cost_atb_2024)
        else:
            raise Exception("{} is not among the costs files to consider".format(input_file_path.name))


        cost_df.to_csv(snakemake.output[i])

    print(cost_atb_2022.columns)
    cost_atb_2022.to_csv("2022_data.csv")
    print(cost_atb_2022["parameter"].unique())
    print(cost_atb_2022.shape)
    print(cost_atb_2024.columns)
    cost_atb_2024.to_csv("2024_data.csv")
    print(cost_atb_2024["parameter"].unique())
    print(cost_atb_2024.shape)
