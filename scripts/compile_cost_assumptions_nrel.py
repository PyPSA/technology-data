#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from _helpers import mock_snakemake


def read_atb_input_file(input_file_list):
    dataframe_list = []
    for file in input_file_list:
        dataframe_list.append(pd.read_parquet(file))
    return dataframe_list


def filter_input_file(input_file_list, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep):
    atb_input_df_list = read_atb_input_file(input_file_list)
    list_core_metric_parameter_to_keep = [str(x).casefold() for x in list_core_metric_parameter_to_keep]
    list_years = [str(x).casefold() for x in list_years]

    # --> select columns
    atb_input_df_2022 = atb_input_df_list[0].loc[:, list_columns_to_keep]
    atb_input_df_2024 = atb_input_df_list[1].loc[:, list_columns_to_keep]

    # --> select rows based on core_metric_parameter
    # Note: we do not apply any selection on core_metric_case and scenario.
    # Such selection can be done in the model config (e.g. PyPSA-Earth)
    atb_input_df_2022 = atb_input_df_2022.loc[
        atb_input_df_2022["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]
    atb_input_df_2024 = atb_input_df_2024.loc[
        atb_input_df_2024["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]

    # --> select rows based on core_metric_variable
    # Note: 2020 data are fetched from the input file atb_e_2022.*
    atb_input_df_2022 = atb_input_df_2022.loc[
        atb_input_df_2022["core_metric_variable"].astype(str).str.casefold() == list_years[0]]

    # Note: 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
    atb_input_df_2024 = atb_input_df_2024.loc[
        atb_input_df_2024["core_metric_variable"].astype(str).str.casefold().isin(list_years[1:])]

    # --> drop duplicated rows
    atb_input_df_2022 = atb_input_df_2022.drop_duplicates(keep="first")
    atb_input_df_2024 = atb_input_df_2024.drop_duplicates(keep="first")

    return atb_input_df_2022, atb_input_df_2024


def calculate_fom_percentage(x, dataframe):

    # Note: for technologies as Coal Retrofit or Natural Gas Retrofit,
    # the Fixed O&M is normalized by Additional OCC. Else, the Fixed O&M is
    # normalized by the CAPEX

    if x["core_metric_parameter"].casefold() == "fixed o&m":
        if "retrofit" in x["technology"].casefold():
            fom_perc_value = x.value / dataframe.query(
                "core_metric_parameter.str.casefold() == 'additional occ' & "
                "core_metric_case == @x.core_metric_case & "
                "core_metric_variable == @x.core_metric_variable & "
                "technology == @x.technology & "
                "technology_alias == @x.technology_alias & "
                "techdetail == @x.techdetail & "
                "display_name == @x.display_name & "
                "scenario == @x.scenario")["value"]*100.0
        else:
            fom_perc_value = x.value / dataframe.query(
                "core_metric_parameter.str.casefold() == 'capex' & "
                "core_metric_case == @x.core_metric_case & "
                "core_metric_variable == @x.core_metric_variable & "
                "technology == @x.technology & "
                "technology_alias == @x.technology_alias & "
                "techdetail == @x.techdetail & "
                "display_name == @x.display_name & "
                "scenario == @x.scenario")["value"]*100.0
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


def pre_process_input_file(input_file_list, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_source):

    atb_input_df_2022, atb_input_df_2024 = filter_input_file(input_file_list, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep)

    # normalize Fixed O&M by CAPEX (or Additional OCC for retrofit technologies)
    atb_input_df_2022["value"] = atb_input_df_2022.apply(lambda x: calculate_fom_percentage(x, atb_input_df_2022), axis=1)
    atb_input_df_2024["value"] = atb_input_df_2024.apply(lambda x: calculate_fom_percentage(x, atb_input_df_2024), axis=1)

    # modify the unit of the normalized Fixed O&M to %-yr
    atb_input_df_2022["units"] = atb_input_df_2022.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)
    atb_input_df_2024["units"] = atb_input_df_2024.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)

    # create a new column technology_alias_detail from the concatenation of the technology_alias and techdetail columns
    concatenate_columns(atb_input_df_2022, "technology_alias_detail", ["technology_alias", "techdetail"])
    concatenate_columns(atb_input_df_2024, "technology_alias_detail", ["technology_alias", "techdetail"])

    technology_conversion_dict = {
        "Coal_newAvgCF2ndGen": "coal",
        "Natural_Gas_CTAvgCF": "CCGT",
        "Hydropower_NPD1": "hydro",
        "Hydropower_NSD1": "ror",
        "Pumped_Storage_Hydropower_NatlClass1": "PHS",
        "Nuclear_Nuclear": "nuclear",
        "Geothermal_HydroFlash": "geothermal",
        "Land-Based_Wind_Class1": "onwind",
        "Offshore_Wind_Class1": "offwind",
        "Utility_PV_Class1": "solar-utility",
        "Commercial_PV_Class1": "solar-rooftop",
        "Utility-Scale_Battery_Storage_4Hr_Battery_Storage": "battery storage",
        "Biopower_Dedicated": "biomass",
    }

    # replace technology_alias_detail with PyPSA technology names
    replace_value_name(atb_input_df_2022, technology_conversion_dict, "technology_alias_detail")
    replace_value_name(atb_input_df_2024, technology_conversion_dict, "technology_alias_detail")

    # add source column
    atb_input_df_2022["source"] = nrel_source
    atb_input_df_2024["source"] = nrel_source

    column_rename_dict = {
        "technology_alias_detail": "technology",
        "core_metric_parameter": "parameter",
        "units": "unit",
        "display_name": "further description",
        "atb_year": "currency_year",
        "core_metric_case": "financial_case"
    }


    # here you need to keep the core_metric_variable because this is what you need to split the atb --> cost_year.csv
    # ("technology_alias_detail", "core_metric_parameter", "value", "units", "source", "display_name", "atb_year", "scenario", "core_metric_case") = nrel_atb_columns_to_keep - le righe che non mi servono
    atb_input_df_2022 = atb_input_df_2022.loc[:, ("technology_alias_detail", "core_metric_parameter", "value", "units", "source", "display_name", "atb_year", "scenario", "core_metric_case")].rename(columns=column_rename_dict)
    atb_input_df_2024 = atb_input_df_2024.loc[:, ("technology_alias_detail", "core_metric_parameter", "value", "units", "source", "display_name", "atb_year", "scenario", "core_metric_case")].rename(columns=column_rename_dict)

    parameter_conversion_dict = {
        "CAPEX": "investment",
        "Fixed O&M": "FOM",
        "Variable O&M": "VOM",
        "Fuel": "fuel",
        "Additional OCC": "investment",
        "WACC Real": "discount rate"
    }

    # replace parameter with PyPSA cost parameter names
    replace_value_name(atb_input_df_2022, parameter_conversion_dict, "parameter")
    replace_value_name(atb_input_df_2024, parameter_conversion_dict, "parameter")

    # Filter the technolies in atb_input_* with those in PyPSA
    # These rows should be replacing existing values in the next costs_....csv files
    # The rows which do not correspond to pypsa technologies are to be appended to the costs_....csv files

    # Note: core_metric_variable is used such that costs_(core_metric_variable).csv
    return atb_input_df_2022, atb_input_df_2024


if __name__ == "__main__":
    if 'snakemake' not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']
    input_file_list_atb = snakemake.input.nrel_atb_input_files
    nrel_atb_columns_to_keep = snakemake.config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = snakemake.config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_source_link = snakemake.config["nrel_atb"]["nrel_atb_source_link"]

    for i, input_file in enumerate(snakemake.input.cost_files_to_modify):
        cost_df = pd.read_csv(input_file).set_index("technology")
        cost_df.to_csv(snakemake.output[i])

    df_2022, df_2024 = pre_process_input_file(input_file_list_atb, year_list, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)

    print(df_2022.columns)
    df_2022.to_csv("2022_data.csv")
    print(df_2022["parameter"].unique())
    print(df_2022.shape)
    print(df_2024.columns)
    df_2024.to_csv("2024_data.csv")
    print(df_2024["parameter"].unique())
    print(df_2024.shape)
