# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8
"""
script to convert the technology data assumptions of the Fraunhofer ISE Study
"Wege zu einem klimaneutralen Energiesystem" (see pdf in folder docu) into a
.csv format
"""

import numpy as np
import pandas as pd
from tabula import read_pdf

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("convert_pdf_fraunhofer_to_dataframe")

    df_list = read_pdf(snakemake.input.fraunhofer, pages="3-15", multiple_tables=True)
    print(df_list)
    clean_df = []
    j = 0
    for i in range(len(df_list)):
        print(i)
        if len(df_list[i]) == 1:
            print("table is dropped ", i)
        if "Komponente" in df_list[i].iloc[0].unique():
            print("table is added ", i)
            clean_df.append(df_list[i])
            j += 1
            print(j)
            continue
        else:
            clean_df[j - 1] = clean_df[j - 1].append(df_list[i])
            print("table is appended ", i)

    for i in range(len(clean_df)):
        clean_df[i].columns = clean_df[i].iloc[0, :]
        clean_df[i] = clean_df[i].iloc[1:, :]
        clean_df[i].reset_index(drop=True, inplace=True)
        clean_df[i].dropna(axis=1, how="all", inplace=True)

    columns = [
        "Komponente",
        "Größe",
        "Einheit",
        2020,
        2025,
        2030,
        2035,
        2040,
        2045,
        2050,
    ]
    for table in range(len(clean_df)):
        print(table)
        test = clean_df[table]["Komponente"]

        counter = len(test) - 2
        for i in range(counter, -1, -1):
            if (isinstance(test.iloc[i], str)) and (isinstance(test.iloc[i - 1], str)):
                test.iloc[i - 1] = test.iloc[i - 1] + " " + test.iloc[i]
                test.iloc[i] = np.nan
        test.fillna(method="ffill", inplace=True)
        clean_df[table]["Komponente"] = test

        new = {}
        for i in range(len(test)):
            a = clean_df[table].loc[i].dropna()
            if len(a) == len(columns):
                new[i] = a
                new[i].index = columns
            else:
                print(a)

        clean_df[table] = pd.concat(new, axis=1).T
        clean_df[table].set_index(["Komponente", "Größe"], inplace=True)
        clean_df[table].dropna(how="all", inplace=True)

    total = pd.concat(clean_df)

    total.Einheit = total.Einheit.str.replace("€", "EUR")
    total.to_csv(snakemake.output.costs, encoding="iso-8859-15")

    energiepreise = read_pdf(snakemake.input.fraunhofer, pages="15")
    energiepreise.dropna(axis=1, how="all", inplace=True)
    energiepreise.dropna(axis=0, how="all", inplace=True)
    energiepreise = energiepreise.rename(columns={"Unnamed: 1": "Fuel"}).set_index(
        "Fuel"
    )
    energiepreise["unit"] = "Eur/MWh"
    energiepreise.to_csv(snakemake.output.energy_prices, encoding="iso-8859-15")
