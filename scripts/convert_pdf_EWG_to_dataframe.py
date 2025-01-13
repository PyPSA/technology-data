# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8
"""
script to convert the technology data assumptions of the Study
"Global Energy System based on 100% Renewable Energy" of Energywatchgroup/LTU University
http://energywatchgroup.org/wp-content/uploads/EWG_LUT_100RE_All_Sectors_Global_Report_2019.pdf
(see also pdf in folder docu) into a .csv format
"""

import numpy as np
import pandas as pd
from tabula import read_pdf

# Detect running outside of snakemake and mock snakemake for testing
if "snakemake" not in globals():
    from vresutils.snakemake import MockSnakemake

    snakemake = MockSnakemake()
    snakemake.input = dict(EWG="docu/EWG_LUT_100RE_All_Sectors_Global_Report_2019.pdf")
    snakemake.output = dict(costs="inputs/EWG_costs.csv")

df_list = read_pdf(snakemake.input["EWG"], pages="305-309", multiple_tables=True)
# %%
# wished columns
wished_columns = [
    "Technologies",
    "Type",
    "Units",
    "2015",
    "2020",
    "2025",
    "2030",
    "2035",
    "2040",
    "2045",
    "2050",
    "Ref",
]
# clean data frame
split_units = df_list[0]["Units 2015"].fillna(" ").str.split(" ", expand=True)
# check where split is too long
to_be_merged = split_units[split_units[2].apply(lambda x: x is not None)].index
split_units.loc[to_be_merged, 0] = (
    split_units.loc[to_be_merged, 0] + " " + split_units.loc[to_be_merged, 1]
)
split_units.loc[to_be_merged, 1] = split_units.loc[to_be_merged, 2]
df_list[0] = pd.concat(
    [df_list[0].drop("Units 2015", axis=1), split_units.iloc[:, 0:2]], axis=1
)

# renmae columns
df_list[0] = (
    df_list[0]
    .rename(columns={"Unnamed: 0": "Type", 0: "Units", 1: "2015"})
    .reindex(wished_columns, axis=1)
)
for page in range(1, len(df_list)):
    df_list[page] = df_list[page].T.reset_index().T
    if len(df_list[page].columns) != len(wished_columns):
        df_list[page][11] = df_list[page].loc[:, 11:].fillna("").agg(" ".join, axis=1)
        df_list[page] = df_list[page].iloc[:, :12]
    df_list[page].columns = wished_columns

df_list[4] = pd.concat(
    [df_list[4].iloc[:, :2], df_list[4].iloc[:, 2:].shift(axis=1)], axis=1
)
for sign in [" €", " kWh"]:
    split_units = df_list[4]["Type"].fillna(" ").str.split(sign, expand=True)
    # check where split is too long
    to_be_merged = split_units[split_units[1].apply(lambda x: x is not None)].index
    df_list[4].loc[to_be_merged, "Type"] = split_units.loc[to_be_merged, 0]
    df_list[4].loc[to_be_merged, "Units"] = sign + split_units.loc[to_be_merged, 1]

clean_df = pd.concat([df_list[page] for page in range(len(df_list))])
clean_df.reset_index(drop=True, inplace=True)
# %%
# fill missing rows with tech name
for row in range(len(clean_df)):
    if not str(clean_df.loc[row, "Technologies"]) == "nan":
        for end in range(row + 1, row + 5):
            if not (
                any(clean_df.loc[row:end, "Technologies"].isna())
                or any(
                    [
                        exclude in str(clean_df.loc[end, "Technologies"])
                        for exclude in ["Residential", "Battery", "DH"]
                    ]
                )
            ):
                clean_df.loc[row, "Technologies"] += (
                    " " + clean_df.loc[end, "Technologies"]
                )
            else:
                if any(
                    [
                        exclude in str(clean_df.loc[end, "Technologies"])
                        for exclude in ["Residential", "Battery", "DH"]
                    ]
                ):
                    end -= 1
                clean_df.loc[row + 1 : end, "Technologies"] = np.nan
                break

# convert to float
years = ["2015", "2020", "2025", "2030", "2035", "2040", "2045", "2050"]
clean_df[years] = clean_df[years].applymap(lambda x: str(x).replace(",", "."))
clean_df[years] = clean_df[years].apply(lambda x: pd.to_numeric(x, errors="coerce"))

clean_df = clean_df.loc[clean_df[years].dropna(how="all", axis=0).index]
clean_df.Technologies.fillna(method="ffill", inplace=True)
clean_df.Type.fillna(method="ffill", inplace=True)
clean_df.set_index(["Technologies", "Type"], inplace=True)


clean_df[years] = clean_df[years].fillna(axis=1, method="ffill")

rename_types = {
    "Lifetime": "lifetime",
    "Lifetime years": "lifetime",
    "Opex var": "VOM",
    "[120Opex var": "VOM",
    "[123,Opex fix": "FOM",
    "[125,Opex fix": "FOM",
    "[137Opex fix": "FOM",
    "Opex fix": "FOM",
    "Capex": "investment",
}
clean_df.rename(index=rename_types, level=1, inplace=True)
aggregate = {year: sum for year in years}
aggregate["Units"] = "first"
aggregate["Ref"] = "first"
final_df = clean_df.groupby(clean_df.index).agg(aggregate)
final_df.index = pd.MultiIndex.from_tuples(final_df.index)

fom = (
    final_df[years].xs("FOM", level=1).div(final_df[years].xs("investment", level=1))
    * 100
).dropna(how="all", axis=0)
fom.index = pd.MultiIndex.from_product([fom.index, ["FOM"]])
final_df.loc[fom.index, years] = fom
final_df.loc[fom.index, "Units"] = "%/year"
final_df.index.rename(["technology", "parameter"], inplace=True)
final_df.rename(columns={"Units": "unit", "Ref": "source"}, inplace=True)
final_df["unit"] = final_df["unit"].apply(lambda x: str(x).replace("€", "EUR"))

round(final_df, ndigits=2).to_csv(snakemake.output["costs"], encoding="iso-8859-15")
