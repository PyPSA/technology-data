#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import pathlib
import numpy as np
from _helpers import mock_snakemake


def filter_input_file(input_file_path, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep):

    atb_input_df = pd.read_parquet(input_file_path)
    list_core_metric_parameter_to_keep = [str(x).casefold() for x in list_core_metric_parameter_to_keep]
    list_years = [str(x).casefold() for x in list_years]

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
        atb_input_df = atb_input_df.loc[atb_input_df["core_metric_variable"].astype(str).str.casefold() == list_years[0]]
    elif input_file_path.name == "atb_e_2024.parquet":
        # Note: 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
        atb_input_df = atb_input_df.loc[atb_input_df["core_metric_variable"].astype(str).str.casefold().isin(list_years[1:])]
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
            fom_perc_value = x.value / dataframe.query("atb_year == @x.atb_year & core_metric_parameter.str.casefold() == 'additional occ' & core_metric_case == @x.core_metric_case & core_metric_variable == @x.core_metric_variable & technology == @x.technology & technology_alias == @x.technology_alias & techdetail == @x.techdetail & display_name == @x.display_name & scenario == @x.scenario")["value"]*100.0
        else:
            fom_perc_value = x.value / dataframe.query("atb_year == @x.atb_year & core_metric_parameter.str.casefold() == 'capex' & core_metric_case == @x.core_metric_case & core_metric_variable == @x.core_metric_variable & technology == @x.technology & technology_alias == @x.technology_alias & techdetail == @x.techdetail & display_name == @x.display_name & scenario == @x.scenario")["value"]*100.0
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


def repeat_values(dataframe, column_name):
    dataframe_to_repeat = dataframe.query("technology_alias_detail.str.casefold() == 'hydropower_*'")
    dataframe_to_keep = dataframe.query("technology_alias_detail.str.casefold() != 'hydropower_*'")
    dataframe_repeated = pd.DataFrame(np.repeat(dataframe_to_repeat.values, 2, axis=0), columns=dataframe_to_repeat.columns)
    dataframe_repeated.loc[0::2, column_name] = "Hydropower_NPD1"
    dataframe_repeated.loc[1::2, column_name] = "Hydropower_NSD1"
    return pd.concat([dataframe_to_keep, dataframe_repeated])


def pre_process_input_file(input_file_list, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_source):

    # read inputs and filter relevant columns and relevant core_metric_variables (i.e. the years to consider)
    dictionary_list = []
    for file_path in input_file_list:
        file_path_val = pathlib.Path(file_path)
        dictionary_list.append(filter_input_file(file_path_val, list_years, list_columns_to_keep, list_core_metric_parameter_to_keep))

    atb_input_df_2022, atb_input_df_2024 = dictionary_list[0], dictionary_list[1]

    # normalize Fixed O&M by CAPEX (or Additional OCC for retrofit technologies)
    atb_input_df_2022["value"] = atb_input_df_2022.apply(lambda x: calculate_fom_percentage(x, atb_input_df_2022), axis=1)
    atb_input_df_2024["value"] = atb_input_df_2024.apply(lambda x: calculate_fom_percentage(x, atb_input_df_2024), axis=1)

    # modify the unit of the normalized Fixed O&M to %-yr
    atb_input_df_2022["units"] = atb_input_df_2022.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)
    atb_input_df_2024["units"] = atb_input_df_2024.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)

    # create a new column technology_alias_detail from the concatenation of the technology_alias and techdetail columns
    concatenate_columns(atb_input_df_2022, "technology_alias_detail", ["technology_alias", "techdetail"])
    concatenate_columns(atb_input_df_2024, "technology_alias_detail", ["technology_alias", "techdetail"])

    atb_input_df_2022 = repeat_values(atb_input_df_2022, "technology_alias_detail")
    atb_input_df_2024 = repeat_values(atb_input_df_2024, "technology_alias_detail")

    # replace technology_alias_detail with PyPSA technology names
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
    replace_value_name(atb_input_df_2022, technology_conversion_dict_atb, "technology_alias_detail")
    replace_value_name(atb_input_df_2024, technology_conversion_dict_atb, "technology_alias_detail")

    # add source column
    atb_input_df_2022["source"] = nrel_source
    atb_input_df_2024["source"] = nrel_source

    # rename columns and consider just columns to keep for output (**)
    column_rename_dict = {
        "technology_alias_detail": "technology",
        "core_metric_parameter": "parameter",
        "units": "unit",
        "display_name": "further description",
        "atb_year": "currency_year",
        "core_metric_case": "financial_case"
    }

    tuple_output_columns_to_keep = ("technology_alias_detail", "core_metric_parameter", "value", "units", "source", "display_name", "atb_year", "scenario", "core_metric_case", "core_metric_variable", "tax_credit_case")
    atb_input_df_2022 = atb_input_df_2022.loc[:, tuple_output_columns_to_keep].rename(columns=column_rename_dict)
    atb_input_df_2024 = atb_input_df_2024.loc[:, tuple_output_columns_to_keep].rename(columns=column_rename_dict)

    # replace parameter with PyPSA cost parameter names
    parameter_conversion_dict = {
        "CAPEX": "investment",
        "Fixed O&M": "FOM",
        "Variable O&M": "VOM",
        "Fuel": "fuel",
        "Additional OCC": "investment",
        "WACC Real": "discount rate"
    }
    replace_value_name(atb_input_df_2022, parameter_conversion_dict, "parameter")
    replace_value_name(atb_input_df_2024, parameter_conversion_dict, "parameter")

    # Filter the technologies in atb_input_* with those in PyPSA
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
