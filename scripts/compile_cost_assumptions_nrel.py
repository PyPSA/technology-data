#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake("compile_cost_assumptions_nrel")

    year_list = snakemake.config['years']

    for i, input_file in enumerate(snakemake.input):
        cost_df = pd.read_csv(input_file).set_index("technology")
        print(snakemake.output[i])
        cost_df.to_csv(snakemake.output[i])

