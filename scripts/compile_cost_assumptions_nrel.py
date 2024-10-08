#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from _helpers import mock_snakemake


def read_atb_input_file(input_file_list):
    dataframe_list = []
    for file in input_file_list:
        dataframe_list.append(pd.read_parquet(file))
    return dataframe_list


def pre_process_input_file(input_file_list, list_columns_to_keep, list_core_metric_parameter_to_keep):
    atb_input_df_list = read_atb_input_file(input_file_list)
    list_core_metric_parameter_to_keep = [str(x).casefold() for x in list_core_metric_parameter_to_keep]

    # 2020 data are fetched from the input file atb_e_2022.*
    atb_input_df_2022 = atb_input_df_list[0].loc[:, list_columns_to_keep]

    # 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
    atb_input_df_2024 = atb_input_df_list[1].loc[:, list_columns_to_keep]

    atb_input_df_2022 = atb_input_df_2022.loc[atb_input_df_2022["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]
    atb_input_df_2024 = atb_input_df_2024.loc[
        atb_input_df_2024["core_metric_parameter"].str.casefold().isin(list_core_metric_parameter_to_keep)]

    return atb_input_df_2022, atb_input_df_2024



if __name__ == "__main__":
    if 'snakemake' not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']
    input_file_list_atb = snakemake.input.nrel_atb_input_files
    nrel_atb_columns_to_keep = snakemake.config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = snakemake.config["nrel_atb"]["nrel_atb_core_metric_parameter_to_keep"]

    for i, input_file in enumerate(snakemake.input.cost_files_to_modify):
        cost_df = pd.read_csv(input_file).set_index("technology")
        cost_df.to_csv(snakemake.output[i])

    df_2022, df_2024 = pre_process_input_file(input_file_list_atb, nrel_atb_columns_to_keep, nrel_atb_core_metric_parameter_to_keep)

    print(df_2022.columns)
    print(df_2022["core_metric_parameter"].unique())
    print(df_2024.columns)
    print(df_2024["core_metric_parameter"].unique())
