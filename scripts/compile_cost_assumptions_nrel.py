#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import pathlib
from _helpers import mock_snakemake


def get_convertion_dictionary(flag):
    if flag.casefold() == "parameter":
        return {
            "CAPEX": "investment",
            "Fixed O&M": "FOM",
            "Variable O&M": "VOM",
            "Fuel": "fuel",
            "Additional OCC": "investment",
            "WACC Real": "discount rate"
        }
    elif flag.casefold() == "technology":
        return {
            "Coal-new -> 2nd Gen Tech": "coal",
            "Coal-new": "coal",
            "NG F-Frame CT": "CCGT",
            "NG Combustion Turbine (F-Frame)": "CCGT",
            "Hydropower - NPD 1": "hydro",
            "Hydropower - NSD 1": "ror",
            "Pumped Storage Hydropower - National Class 1": "PHS",
            "Nuclear - Large": "nuclear",
            "Nuclear - AP100": "nuclear",
            "Geothermal - Hydro / Flash": "geothermal",
            "Land-Based Wind - Class 1": "onwind",
            "Land-Based Wind - Class 1 - Technology 1": "onwind",
            "Offshore Wind - Class 1": "offwind",
            "Utility PV - Class 1": "solar-utility",
            "Commercial PV - Class 1": "solar-rooftop",
            "Utility-Scale Battery Storage - 6Hr": "battery storage",
            "Biopower - Dedicated": "biomass",
            "CSP - Class 2": "csp-tower",
        }
    elif flag.casefold() == "output_column":
        return {
            "display_name": "technology",
            "core_metric_parameter": "parameter",
            "units": "unit",
            "atb_year": "currency_year",
            "core_metric_case": "financial_case"
        }
    else:
        raise Exception("{} is not among the allowed choices: parameter, technology, output_column")


def filter_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep):

    atb_file_df = pd.read_parquet(input_file_path)
    list_core_metric_parameter_to_keep = [str(x).casefold() for x in list_core_metric_parameter_to_keep]
    year_string = str(year).casefold()

    # --> select columns
    for column_name in list_columns_to_keep:
        if column_name not in atb_file_df:
            print("missing column", column_name)
            atb_file_df[column_name] = pd.Series(dtype="str")
        else:
            pass
    atb_file_df = atb_file_df.loc[:, list_columns_to_keep]

    # --> select rows based on core_metric_parameter
    # Note: we do not apply any selection on core_metric_case and scenario.
    # Such selection can be done in the model config (e.g. PyPSA-Earth)
    atb_file_df = atb_file_df.loc[atb_file_df["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]

    # --> select rows based on core_metric_variable
    if input_file_path.name == "atb_e_2022.parquet":
        # Note: 2020 data are fetched from the input file atb_e_2022.
        atb_file_df = atb_file_df.loc[atb_file_df["core_metric_variable"].astype(str).str.casefold() == year_string]
    elif input_file_path.name == "atb_e_2024.parquet":
        # Note: 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
        atb_file_df = atb_file_df.loc[atb_file_df["core_metric_variable"].astype(str).str.casefold() == year_string]
    else:
        raise Exception("{} - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet".format(input_file_path.stem))

    # --> drop duplicated rows
    atb_file_df = atb_file_df.drop_duplicates(keep="first")

    # --> remove technology AEO
    atb_file_df = atb_file_df.query("technology.str.casefold() != 'aeo'")

    return atb_file_df


def get_query_string(column_list, column_to_exclude, parameter_value):
    if set(column_to_exclude).issubset(set(column_list)):
        column_to_use_list = list(set(column_list)-set(column_to_exclude))
        query_list = sorted(["{}.str.casefold() == '{}'".format(column_name, parameter_value) if column_name == "core_metric_parameter" else "{} == @x.{}".format(
            column_name, column_name) for column_name in column_to_use_list])
        query_string = " & ".join(query_list)
    else:
        exception_message = "The following columns {} are not included in the original list".format(list(set(column_to_exclude).difference(set(column_list))))
        raise Exception(exception_message)
    return query_string.strip()


def calculate_fom_percentage(x, dataframe, columns_list):

    # Note: for technologies as Coal Retrofit or Natural Gas Retrofit,
    # the Fixed O&M is normalized by Additional OCC. Else, the Fixed O&M is
    # normalized by the CAPEX

    if x["core_metric_parameter"].casefold() == "fixed o&m":
        if "retrofit" in x["technology"].casefold():
            query_string = get_query_string(columns_list, ["units", "value", "tax_credit_case"], "additional occ")
        else:
            query_string = get_query_string(columns_list, ["units", "value", "tax_credit_case"], "capex")
        fom_perc_value = x.value / dataframe.query(query_string)["value"]*100.0
        return round(fom_perc_value.values[0], 2)
    else:
        return x.value


def replace_value_name(dataframe, conversion_dict, column_name):
    dataframe[column_name] = dataframe[column_name].replace(conversion_dict)
    return dataframe


def pre_process_input_file(input_file_path, year, list_columns_to_keep, list_core_metric_parameter_to_keep, nrel_source):

    # Read inputs and filter relevant columns and relevant core_metric_variables (i.e. the years to consider)
    atb_input_df = filter_input_file(pathlib.Path(input_file_path), year, list_columns_to_keep, list_core_metric_parameter_to_keep)

    # Normalize Fixed O&M by CAPEX (or Additional OCC for retrofit technologies)
    atb_input_df["value"] = atb_input_df.apply(lambda x: calculate_fom_percentage(x, atb_input_df, list_columns_to_keep), axis=1)

    # Modify the unit of the normalized Fixed O&M to %-yr
    atb_input_df["units"] = atb_input_df.apply(lambda x: "%-yr" if x["core_metric_parameter"].casefold() == "fixed o&m" else x["units"], axis=1)

    # Replace the display_name column values with PyPSA technology names
    technology_conversion_dict_atb = get_convertion_dictionary("technology")
    atb_input_df = replace_value_name(atb_input_df, technology_conversion_dict_atb, "display_name")

    # Add source column
    atb_input_df["source"] = nrel_source

    # Add further description column
    atb_input_df["further description"] = pd.Series(dtype="str")

    # Rename columns and consider just columns used in PyPSA
    column_rename_dict = get_convertion_dictionary("output_column")
    tuple_output_columns_to_keep = ("display_name", "core_metric_parameter", "value", "units", "source", "further description", "atb_year", "scenario", "core_metric_case", "tax_credit_case")
    atb_input_df = atb_input_df.loc[:, tuple_output_columns_to_keep].rename(columns=column_rename_dict)

    # Replace parameter with PyPSA cost parameter names
    parameter_conversion_dict = get_convertion_dictionary("parameter")
    atb_input_df = replace_value_name(atb_input_df, parameter_conversion_dict, "parameter")

    # Cast currency_year from int to float
    atb_input_df["currency_year"] = atb_input_df["currency_year"].astype(float)

    return atb_input_df.reset_index(drop=True)


def update_cost_values(cost_dataframe, atb_dataframe, conversion_dictionary, columns_to_add_list):
    for column_name in columns_to_add_list:
        cost_dataframe[column_name] = pd.Series(dtype="str")
    values_to_keep = [str(x).casefold() for x in set(conversion_dictionary.values())]

    # query to keep the technologies currently NOT included in PyPSA
    query_string_part_one = "~technology.str.casefold().isin(@values_to_keep)"

    # query to keep the technologies currently included in PyPSA for which the columns ["financial_case", "scenario"] are valid values
    query_string_part_two = "technology.str.casefold().isin(@values_to_keep) & ~(scenario.isna() & financial_case.isna())"

    # query to replace the entries corresponding to technologies currently included in PyPSA with valus from NREL and keep value other technologies
    query_string = "{} | {}".format(query_string_part_one, query_string_part_two)

    new_cost_dataframe = pd.concat([cost_dataframe, atb_dataframe]).query(query_string).reset_index(drop=True)

    return new_cost_dataframe


def add_discount_rate_value(cost_dataframe, discount_rate_dataframe):
    discount_rate_dataframe["further description"] = pd.Series(dtype=str)
    return pd.concat([cost_dataframe, discount_rate_dataframe]).reset_index(drop=True)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']
    input_file_list_atb = snakemake.input.nrel_atb_input_files
    input_file_discount_rate = snakemake.input.nrel_atb_input_discount_rate
    cost_file_list = snakemake.input.cost_files_to_modify
    nrel_atb_columns_to_keep = snakemake.config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = snakemake.config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]
    nrel_atb_source_link = snakemake.config["nrel_atb"]["nrel_atb_source_link"]

    if len(year_list) != len(cost_file_list):
        raise Exception("The cost files {} are more than the considered years {}".format(year_list, cost_file_list))

    discount_rate_df = pd.read_csv(input_file_discount_rate)

    for year_val in year_list:

        # get the cost file to modify
        input_cost_path_list = [path for path in snakemake.input.cost_files_to_modify if str(year_val) in path]
        if len(input_cost_path_list) == 1:
            input_cost_path = input_cost_path_list[0]
            cost_df = pd.read_csv(input_cost_path).reset_index(drop=True)
        else:
            raise Exception("Please verify the list of cost files. It may contain duplicates.")

        # get the atb values for a given year
        if year_val == 2020:
            # choose atb_e_2022
            input_file_atb = input_file_list_atb[0]
            atb_e_df = pre_process_input_file(input_file_atb, year_val, nrel_atb_columns_to_keep,
                                   nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)
        elif year_val in year_list[1:]:
            # choose atb_e_2024
            input_file_atb = input_file_list_atb[1]
            atb_e_df = pre_process_input_file(input_file_atb, year_val, nrel_atb_columns_to_keep,
                                   nrel_atb_core_metric_parameter_to_keep, nrel_atb_source_link)
        else:
            raise Exception("{} is not a considered year".format(year_val))

        # get the discount rate file
        discount_rate_year_df = discount_rate_df.loc[discount_rate_df["year"] == year_val]
        discount_rate_year_df = discount_rate_year_df.loc[:, ("technology", "parameter", "value", "unit", "source", "further description", "currency_year", "financial_case", "scenario", "tax_credit_case")].reset_index(drop=True)

        # update the cost file
        updated_cost_df = update_cost_values(cost_df, atb_e_df, get_convertion_dictionary("technology"), ["financial_case", "scenario", "tax_credit_case"])

        # add discount rate
        updated_cost_df = add_discount_rate_value(updated_cost_df, discount_rate_year_df)

        # output the modified cost file
        output_cost_path_list = [path for path in snakemake.output if str(year_val) in path]
        if len(output_cost_path_list) == 1:
            output_cost_path = output_cost_path_list[0]
            updated_cost_df.to_csv(output_cost_path)
        else:
            raise Exception("Please verify the list of cost files. It may contain duplicates.")



