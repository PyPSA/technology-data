#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from _helpers import mock_snakemake


def read_atb_input_file(input_file_list):
    dataframe_list = []
    for file in input_file_list:
        dataframe_list.append(pd.read_parquet(file))
    return dataframe_list


if __name__ == "__main__":
    if 'snakemake' not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']
    input_file_list_atb = snakemake.input.atb_input_files

    for i, input_file in enumerate(snakemake.input.cost_files_to_modify):
        cost_df = pd.read_csv(input_file).set_index("technology")
        print(snakemake.output[i])
        cost_df.to_csv(snakemake.output[i])

    atb_input_df = read_atb_input_file(input_file_list_atb)
    print(atb_input_df)
