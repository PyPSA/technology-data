#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""_summary_





"""

import numpy as np
import pandas as pd
from compile_cost_assumptions import (
    clean_up_units,
    get_excel_sheets,
    get_sheet_location,
    new_format,
    sheet_names,
    source_dict,
    uncrtnty_lookup,
)

# ---------- sources -------------------------------------------------------


dea_indonesia = {
    "onwind": "Wind Onshore",
    "offwind": "Wind Offshore",
    "solar-utility": "Ground-mounted PV",
    "solar-rooftop residential": "Rooftop PV",
    "solar-rooftop commercial": "Industrial PV",
    "OCGT": "SCGT",
    "CCGT": "CCGT",
    "oil": "Diesel power plant",
    "biomass CHP": "Biomass power plant (small)",
    "coal": "Coal Subcritical",
    "biogas": "Biogas power plant (small)",
    "geothermal": "Geothermal - large",
    "hydro": "Large hydro power plant",
    "Pumped-Storage-Hydro-store": "Hydro pumped storage",
}


dea_vietnam = {
    "onwind": "8 Wind onshore",
    "offwind": "8 Wind offshore Nearshore",
    "solar-utility": "7 Ground-mounted PV",
    "solar-rooftop residential": "7 Rooftop PV",
    "OCGT": "3 SCGT",
    "CCGT": "3 CCGT",
    "coal": "1 Coal supercritical",
    "oil": "14 ICE - Diesel",
    "biogas": "13 Biogas power plant (small)",
    "geothermal": "15 Geothermal - large",
    "hydro": "6 Large hydro",
    "Pumped-Storage-Hydro-store": "16 Hydro pumped storage",
}


def get_data_from_DEA(data_in, sheet_names, expectation=None):
    """
    saves technology data from DEA in dictionary d_by_tech
    """
    d_by_tech = {}

    for tech, dea_tech in sheet_names.items():
        print(f"{tech} in PyPSA corresponds to {dea_tech} in DEA database.")
        df = get_data_DEA(
            tech=tech, data_in=data_in, sheet_names=sheet_names, expectation=expectation
        ).fillna(0)
        d_by_tech[tech] = df

    return d_by_tech


#
def get_data_DEA(tech, data_in, sheet_names, expectation=None):
    """
    interpolate cost for a given technology from DEA database sheet

    uncertainty can be "optimist", "pessimist" or None|""
    """
    excel_file = get_sheet_location(tech, sheet_names, data_in)
    if excel_file is None:
        print("excel file not found for tech ", tech)
        return None

    if tech == "battery":
        usecols = "B:J"
    elif tech in ["direct air capture", "cement capture", "biomass CHP capture"]:
        usecols = "A:F"
    elif tech in [
        "industrial heat pump medium temperature",
        "industrial heat pump high temperature",
        "electric boiler steam",
        "gas boiler steam",
        "solid biomass boiler steam",
        "solid biomass boiler steam CC",
        "direct firing gas",
        "direct firing gas CC",
        "direct firing solid fuels",
        "direct firing solid fuels CC",
    ]:
        usecols = "A:E"
    elif tech in ["Fischer-Tropsch", "Haber-Bosch", "air separation unit"]:
        usecols = "B:F"
    else:
        usecols = "B:G"

    usecols += f",{uncrtnty_lookup[tech]}"

    if (tech in new_format) and (tech != "electrolysis"):
        skiprows = [0]
    else:
        skiprows = [0, 1]

    excel = pd.read_excel(
        excel_file,
        sheet_name=sheet_names[tech],
        index_col=0,
        usecols=usecols,
        skiprows=skiprows,
        na_values="N.A",
    )
    # print(excel)

    excel.dropna(axis=1, how="all", inplace=True)

    excel.index = excel.index.fillna(" ")
    excel.index = excel.index.astype(str)
    excel.dropna(axis=0, how="all", inplace=True)
    # print(excel)

    if 2020 not in excel.columns:
        selection = excel[excel.isin([2020])].dropna(how="all").index
        excel.columns = excel.loc[selection].iloc[0, :].fillna("Technology", limit=1)
        excel.drop(selection, inplace=True)

    uncertainty_columns = ["2050-optimist", "2050-pessimist"]
    if uncrtnty_lookup[tech]:
        # hydrogen storage sheets have reverse order of lower/upper estimates
        if tech in [
            "hydrogen storage tank type 1 including compressor",
            "hydrogen storage cavern",
        ]:
            uncertainty_columns.reverse()
        excel.rename(
            columns={
                excel.columns[-2]: uncertainty_columns[0],
                excel.columns[-1]: uncertainty_columns[1],
            },
            inplace=True,
        )
    else:
        for col in uncertainty_columns:
            excel.loc[:, col] = excel.loc[:, 2050]

    swap_patterns = [
        "technical life",
        "efficiency",
        "Hydrogen output, at LHV",
    ]  # cases where bigger is better
    swap = [any(term in idx.lower() for term in swap_patterns) for idx in excel.index]
    tmp = excel.loc[swap, "2050-pessimist"]
    excel.loc[swap, "2050-pessimist"] = excel.loc[swap, "2050-optimist"]
    excel.loc[swap, "2050-optimist"] = tmp

    if expectation:
        excel.loc[:, 2050] = excel.loc[:, f"2050-{expectation}"].combine_first(
            excel.loc[:, 2050]
        )
    excel.drop(columns=uncertainty_columns, inplace=True)

    # fix for battery with different excel sheet format
    if tech == "battery":
        excel.rename(columns={"Technology": 2040}, inplace=True)

    if expectation:
        excel = excel.loc[:, [2020, 2050]]

    parameters = [
        "efficiency",
        "investment",
        "Fixed O&M",
        "Variable O&M",
        "production capacity for one unit",
        "Output capacity expansion cost",
        "Hydrogen output",
        "Hydrogen (% total input_e (MWh / MWh))",
        " - hereof recoverable for district heating (%-points of heat loss)",
        "Cb coefficient",
        "Cv coefficient",
        "Distribution network costs",
        "Technical life",
        "Energy storage expansion cost",
        "Output capacity expansion cost (M€2015/MW)",
        "Heat input",
        "Heat  input",
        "Electricity input",
        "Eletricity input",
        "Heat out",
        "capture rate",
        "FT Liquids Output, MWh/MWh Total Input",
        " - hereof recoverable for district heating (%-points of heat loss)",
        "Bio SNG (% of fuel input)",
        "Total O&M",
    ]

    df = pd.DataFrame()
    for para in parameters:
        # attr = excel[excel.index.str.contains(para)]
        attr = excel[[para in index for index in excel.index]]
        if len(attr) != 0:
            df = pd.concat([df, attr])
    df.index = df.index.str.replace("€", "EUR")

    df = df.reindex(columns=df.columns[df.columns.isin(years)])
    df = df[~df.index.duplicated(keep="first")]

    # replace missing data
    df.replace("-", np.nan, inplace=True)
    # average data  in format "lower_value-upper_value"
    df = df.applymap(
        lambda x: (float((x).split("-")[0]) + float((x).split("-")[1])) / 2
        if (type(x) == str and "-" in x)
        else x
    )
    # remove symbols "~", ">", "<" and " "
    for sym in ["~", ">", "<", " "]:
        df = df.applymap(lambda x: x.replace(sym, "") if type(x) == str else x)

    df = df.astype(float)
    df = df.mask(
        df.apply(pd.to_numeric, errors="coerce").isnull(),
        df.astype(str).apply(lambda x: x.str.strip()),
    )
    # print(df)

    # ## Modify data loaded from DEA on a per-technology case
    # if (tech == "offwind") and snakemake.config['offwind_no_gridcosts']:
    #     df.loc['Nominal investment (*total) [MEUR/MW_e, 2020]'] -= excel.loc['Nominal investment (installation: grid connection) [M€/MW_e, 2020]']  # TODO generalise

    # Exlucde indirect costs for centralised system with additional piping.
    if tech.startswith("industrial heat pump"):
        df = df.drop("Indirect investments cost (MEUR per MW)")

    if tech == "biogas plus hydrogen":
        df.drop(df.loc[df.index.str.contains("GJ SNG")].index, inplace=True)

    if tech == "BtL":
        df.drop(df.loc[df.index.str.contains("1,000 t FT Liquids")].index, inplace=True)

    if tech == "methanolisation":
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech == "Fischer-Tropsch":
        df.drop(df.loc[df.index.str.contains("l FT Liquids")].index, inplace=True)

    if tech == "biomass boiler":
        df.drop(
            df.loc[df.index.str.contains("Possible additional")].index, inplace=True
        )
        df.drop(df.loc[df.index.str.contains("Total efficiency")].index, inplace=True)

    if tech == "Haber-Bosch":
        df.drop(
            df.loc[
                df.index.str.contains("Specific investment mark-up factor optional ASU")
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains(
                    "Specific investment (MEUR /TPD Ammonia output", regex=False
                )
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Fixed O&M (MEUR /TPD Ammonia", regex=False)
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Variable O&M (EUR /t Ammonia)", regex=False)
            ].index,
            inplace=True,
        )

    if tech == "air separation unit":
        # Bugfix: DEA renewable fuels 04/2022 has wrong unit (MEUR instead of kEUR)
        df.index = df.index.str.replace(
            "Fixed O&M (MEUR /TPD Ammonia)",
            "Fixed O&M (kEUR /TPD Ammonia)",
            regex=False,
        )

        # Calculate ASU cost separate to HB facility in terms of t N2 output
        df.loc[
            [
                "Specific investment (MEUR /TPD Ammonia output)",
                "Fixed O&M (kEUR /TPD Ammonia)",
                "Variable O&M (EUR /t Ammonia)",
            ]
        ] *= (
            df.loc["Specific investment mark-up factor optional ASU"] - 1.0
        ) / excel.loc[
            "N2 Consumption, t/t Ammonia"
        ]
        # Convert output to hourly generation
        df.loc[
            [
                "Specific investment (MEUR /TPD Ammonia output)",
                "Fixed O&M (kEUR /TPD Ammonia)",
            ]
        ] *= 24

        # Rename costs for correct units
        df.index = df.index.str.replace("MEUR /TPD Ammonia output", "MEUR/t_N2/h")
        df.index = df.index.str.replace("kEUR /TPD Ammonia", "kEUR/t_N2/h/year")
        df.index = df.index.str.replace("EUR /t Ammonia", "EUR/t_N2")

        df.drop(
            df.loc[
                df.index.str.contains("Specific investment mark-up factor optional ASU")
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains(
                    "Specific investment (MEUR /MW Ammonia output)", regex=False
                )
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Fixed O&M (kEUR/MW Ammonia/year)", regex=False)
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Variable O&M (EUR/MWh Ammonia)", regex=False)
            ].index,
            inplace=True,
        )

    if "solid biomass power" in tech:
        df.index = df.index.str.replace("EUR/MWeh", "EUR/MWh")

    df_final = pd.DataFrame(index=df.index, columns=years)

    # [RTD-interpolation-example]
    for index in df_final.index:
        values = np.interp(
            x=years,
            xp=df.columns.values.astype(float),
            fp=df.loc[index, :].values.astype(float),
        )
        df_final.loc[index, :] = values

    # if year-specific data is missing and not fixed by interpolation fill forward with same values
    df_final = df_final.fillna(method="ffill", axis=1)

    df_final["source"] = source_dict["DEA"] + ", " + excel_file.replace("inputs/", "")
    if tech in new_format and (tech != "electrolysis"):
        for attr in ["investment", "Fixed O&M"]:
            to_drop = df[
                df.index.str.contains(attr) & ~df.index.str.contains("\(\*total\)")
            ].index
            df_final.drop(to_drop, inplace=True)

        df_final["unit"] = df_final.rename(
            index=lambda x: x[x.rfind("[") + 1 : x.rfind("]")]
        ).index.values
    else:
        df_final["unit"] = df_final.rename(
            index=lambda x: x[x.rfind("(") + 1 : x.rfind(")")]
        ).index.values
    df_final.index = df_final.index.str.replace(r" \(.*\)", "", regex=True)

    return df_final


def order_data(tech_data):
    """
    check if the units of different variables are conform
    -> print warning if not
    return a pd.Dataframe 'data' in pypsa tech data syntax (investment, FOM,
    VOM, efficiency)
    """

    clean_df = {}
    for tech in tech_data.index.levels[0]:
        clean_df[tech] = pd.DataFrame()
        switch = False
        df = tech_data.loc[tech]

        # --- investment ----
        investment = df[
            (
                df.index.str.contains("investment")
                | df.index.str.contains("Distribution network costs")
            )
            & (
                (df.unit == "EUR/MW")
                | (df.unit == "EUR/MW_e")
                | (df.unit == "EUR/MW_th - heat output")
                | (df.unit == "EUR/MW_th excluding drive energy")
                | (df.unit == "EUR/MW_th")
                | (df.unit == "EUR/MW_MeOH")
                | (df.unit == "EUR/MW_FT/year")
                | (df.unit == "EUR/MW_NH3")
                | (df.unit == "EUR/MWhCapacity")
                | (df.unit == "EUR/MWh")
                | (df.unit == "EUR/MW_CH4")
                | (df.unit == "EUR/MWh/year")
                | (df.unit == "EUR/MW_e, 2020")
                | (df.unit == "EUR/MW input")
                | (df.unit == "EUR/t_N2/h")
            )  # air separation unit
        ].copy()
        if len(investment) != 1:
            switch = True
            print(
                "check investment: ",
                tech,
                " ",
                df[df.index.str.contains("investment")].unit,
            )
        else:
            investment["parameter"] = "investment"
            clean_df[tech] = investment

        # ---- FOM ----------------
        if len(investment):
            fixed = df[
                (
                    df.index.str.contains("Fixed O&M")
                    | df.index.str.contains("Total O&M")
                )
                & (
                    (df.unit == investment.unit[0] + "/year")
                    | (df.unit == "EUR/MW/km/year")
                    | (df.unit == "EUR/MW/year")
                    | (df.unit == "EUR/MW_e/y, 2020")
                    | (df.unit == "EUR/MW_e/y")
                    | (df.unit == "EUR/MW_FT/year")
                    | (df.unit == "EUR/MW_CH4/year")
                    | (df.unit == "% of specific investment/year")
                    | (df.unit == investment.unit.str.split(" ")[0][0] + "/year")
                )
            ].copy()
            if (len(fixed) != 1) and (len(df[df.index.str.contains("Fixed O&M")]) != 0):
                switch = True
                print(
                    "check FOM: ",
                    tech,
                    " ",
                    df[df.index.str.contains("Fixed O&M")].unit,
                )
            if len(fixed) == 1:
                fixed["parameter"] = "fixed"
                clean_df[tech] = pd.concat([clean_df[tech], fixed])
                fom = pd.DataFrame(columns=fixed.columns)
                if not any(fixed.unit.str.contains("% of specific investment/year")):
                    fom[years] = fixed[years] / investment[years].values * 100
                else:
                    fom[years] = fixed[years]
                fom["parameter"] = "FOM"
                fom["unit"] = "%/year"
                fom["source"] = fixed["source"]
                clean_df[tech] = pd.concat([clean_df[tech], fom])

        # ---- VOM -----
        vom = df[
            df.index.str.contains("Variable O&M")
            & (
                (df.unit == "EUR/MWh")
                | (df.unit == "EUR/MWh_e")
                | (df.unit == "EUR/MWh_th")
                | (df.unit == "EUR/MWh_FT")
                | (df.unit == "EUR/MWh_NH3")
                | (df.unit == "EUR/MWh_MeOH")
                | (df.unit == "EUR/MWh/year")
                | (df.unit == "EUR/MWh/km")
                | (df.unit == "EUR/MWh")
                | (df.unit == "EUR/MWhoutput")
                | (tech == "biogas upgrading")
            )
        ].copy()
        if len(vom) == 1:
            vom.loc[:, "parameter"] = "VOM"
            clean_df[tech] = pd.concat([clean_df[tech], vom])

        elif len(vom) != 1 and len(df[df.index.str.contains("Variable O&M")]) != 0:
            switch = True
            print(
                "check VOM: ", tech, " ", df[df.index.str.contains("Variable O&M")].unit
            )

        # ----- lifetime --------
        lifetime = df[
            df.index.str.contains("Technical life") & (df.unit == "years")
        ].copy()
        if len(lifetime) != 1:
            switch = True
            print(
                "check lifetime: ",
                tech,
                " ",
                df[df.index.str.contains("Technical life")].unit,
            )
        else:
            lifetime["parameter"] = "lifetime"
            clean_df[tech] = pd.concat([clean_df[tech], lifetime])

        # ----- efficiencies ------
        efficiency = df[
            (
                df.index.str.contains("efficiency")
                | (df.index.str.contains("Hydrogen output, at LHV"))
                | (df.index.str.contains("FT Liquids Output, MWh/MWh Total Input"))
                | (df.index.str.contains("hereof recoverable for district heating"))
                | (df.index.str.contains("Bio SNG"))
                | (df.index == ("Hydrogen"))
            )
            & (
                (df.unit == "%")
                | (df.unit == "% total size")
                | (df.unit == "% of fuel input")
                | (df.unit == "MWh_H2/MWh_e")
                | (df.unit == "%-points of heat loss")
                | df.unit.str.contains("MWh_FT/MWh_H2")
            )
            & (~df.index.str.contains("name plate"))
        ].copy()

        if tech == "Fischer-Tropsch":
            efficiency[years] *= 100

        # take annual average instead of name plate efficiency
        if any(efficiency.index.str.contains("annual average")):
            efficiency = efficiency[efficiency.index.str.contains("annual average")]

        # hydrogen electrolysiswith recoverable heat
        heat_recovery_label = "hereof recoverable for district heating"
        with_heat_recovery = efficiency.index.str.contains(heat_recovery_label)
        if with_heat_recovery.any():
            efficiency_heat = efficiency[with_heat_recovery].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
            efficiency_h2 = efficiency[efficiency.index.str.contains("Hydrogen")].copy()
            efficiency_h2["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_h2])

        # check if electric and heat efficiencies are given
        if any(["Electric" in ind for ind in efficiency.index]) and any(
            ["Heat" in ind for ind in efficiency.index]
        ):
            efficiency_heat = efficiency[efficiency.index.str.contains("Heat")].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
            efficiency = efficiency[efficiency.index.str.contains("Electric")].copy()
            efficiency["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency])

        elif len(efficiency) != 1:
            switch = True
            if not any(efficiency.index.str.contains("Round trip")):
                print(
                    "check efficiency: ",
                    tech,
                    " ",
                    df[df.index.str.contains("efficiency")].unit,
                )
        else:
            efficiency["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency])

        # add c_v and c_b coefficient
        if "Cb coefficient" in df.index:
            c_b = df.loc[df.index.str.contains("Cb coefficient")].dropna().copy()
            if len(c_b):
                c_b["parameter"] = "c_b"
                clean_df[tech] = pd.concat([clean_df[tech], c_b])
        if "Cv coefficient" in df.index:
            c_v = df.loc[df.index.str.contains("Cv coefficient")].dropna().copy()
            if len(c_v):
                c_v["parameter"] = "c_v"
                clean_df[tech] = pd.concat([clean_df[tech], c_v])

        if switch:
            print("---------------------------------------")

    # concat data
    data = (
        pd.concat(clean_df)
        .reset_index()
        .rename(columns={"level_0": "technology", "level_1": "further description"})
        .set_index(["technology", "parameter"])
    )

    # # add water tank charger/ discharger
    # charger = tech_data.loc[("central water tank storage", "Round trip efficiency")].copy()
    # charger["further description"] = "efficiency from sqr(Round trip efficiency)"
    # charger[years] = charger[years]**0.5*10
    # charger.rename(index={"Round trip efficiency": "efficiency"},
    #                level=1, inplace=True)
    # charger.rename(index={'central water tank storage':"water tank charger"},
    #                level=0, inplace=True)
    # data = pd.concat([data, charger], sort=True)
    # charger.rename(index={"water tank charger": "water tank discharger"},
    #                level=0, inplace=True)
    # data = pd.concat([data, charger], sort=True)

    return data


def add_mean_solar_rooftop(data):
    # take mean of rooftop commercial and residential
    rooftop = (
        data.loc[data.index.get_level_values(0).str.contains("solar-rooftop")][years]
        .astype(float)
        .groupby(level=1)
        .mean()
    )
    for col in data.columns[~data.columns.isin(years)]:
        rooftop[col] = data.loc["solar-rooftop residential"][col]
    # set multi index
    rooftop = pd.concat([rooftop], keys=["solar-rooftop"])
    rooftop["source"] = "Calculated. See 'further description'."
    rooftop[
        "further description"
    ] = "Mixed investment costs based on average of 50% 'solar-rooftop commercial' and 50% 'solar-rooftop residential'"
    # add to data
    data = pd.concat([data, rooftop])
    # add solar assuming 50% utility and 50% rooftop
    solar = (
        (data.loc[["solar-rooftop", "solar-utility"]][years])
        .astype(float)
        .groupby(level=1)
        .mean()
    )
    for col in data.columns[~data.columns.isin(years)]:
        solar[col] = data.loc["solar-rooftop residential"][col]
    solar["source"] = "Calculated. See 'further description'."
    solar[
        "further description"
    ] = "Mixed investment costs based on average of 50% 'solar-rooftop' and 50% 'solar-utility'"
    # set multi index
    solar = pd.concat([solar], keys=["solar"])
    return pd.concat([data, solar])


def add_description(data):
    """
    add as a column to the tech data the excel sheet name,
    add comment for offwind connection costs
    """
    # add excel sheet names to data frame
    wished_order = list(years) + ["unit", "source", "further description"]
    data = data.reindex(columns=wished_order)
    sheets = data.reset_index()["technology"].map(sheet_names).fillna("")
    sheets.index = data.index
    data["further description"] = sheets + ":  " + data["further description"]

    # add comment for offwind investment
    if snakemake.config["offwind_no_gridcosts"]:
        data.loc[
            ("offwind", "investment"), "further description"
        ] += " grid connection costs substracted from investment costs"

    return data


def convert_units(data):
    """
    convert investment and efficiency units to be align with old pypsa
    assumptions
    """
    # convert efficiency from % -> per unit
    data.loc[
        data.index.get_level_values(1).isin(["efficiency", "efficiency-heat"]), years
    ] /= 100
    data.loc[
        data.index.get_level_values(1).isin(["efficiency", "efficiency-heat"]), "unit"
    ] = "per unit"

    # convert MW -> kW
    to_convert = data.index.get_level_values(1).isin(
        ["fixed", "investment"]
    ) & data.unit.str.contains("/MW")
    data.loc[to_convert, years] /= 1e3
    data.loc[to_convert, "unit"] = data.loc[to_convert, "unit"].str.replace(
        "/MW", "/kW"
    )

    return data


def add_country_column(costs: pd.DataFrame, country: str) -> pd.DataFrame:
    """Adds country column to cost dataframe.

    Parameters
    ----------
    costs : pd.DataFrame
        Cost dataframe.
    country : str
        country that needs to be added.
    """
    costs["country"] = country
    return costs


def get_tech_data(data_in, expectation, sheet_names):
    d_by_tech = get_data_from_DEA(
        data_in=data_in,
        expectation=expectation,
        sheet_names=sheet_names,
    )
    tech_data = pd.concat(d_by_tech).sort_index()
    return tech_data


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os

        from _helpers import mock_snakemake

        os.chdir(os.path.join(os.getcwd(), "scripts"))
        snakemake = mock_snakemake("compile_global_costs")

    years = snakemake.config["years"]

    excel_files = [v for k, v in snakemake.input.items() if "dea" in k]

    data_in = get_excel_sheets(excel_files)

    tech_data_idn = (
        get_tech_data(
            data_in=data_in,
            expectation=snakemake.config["expectation"],
            sheet_names=dea_indonesia,
        )
        .pipe(clean_up_units, years)
        .pipe(lambda x: x[x[years].ne(0).sum(axis=1) != 0])
        .pipe(add_mean_solar_rooftop)
        .pipe(order_data)
        .pipe(convert_units)
        .pipe(add_country_column, "ID")
        .set_index("country", append=True)
    )

    tech_data_vn = (
        get_tech_data(
            data_in=data_in,
            expectation=snakemake.config["expectation"],
            sheet_names=dea_vietnam,
        )
        .pipe(clean_up_units, years)
        .pipe(lambda x: x[x[years].ne(0).sum(axis=1) != 0])
        .pipe(add_mean_solar_rooftop)
        .pipe(order_data)
        .pipe(convert_units)
        .pipe(add_country_column, "VN")
        .set_index("country", append=True)
    )

    data = pd.concat([tech_data_idn, tech_data_vn])
    # add excel sheet names and further description
    # data = add_description(data)
    # convert efficiency from %-> per unit and investment from MW->kW to compare

    for year in years:
        costs = data[[year, "unit", "source", "further description"]].rename(
            columns={year: "value"}
        )
        costs["value"] = costs["value"].astype(float)

        costs.sort_index(inplace=True)
        costs.loc[:, "value"] = round(
            costs.value.astype(float), snakemake.config.get("ndigits", 2)
        )
        costs.to_csv([v for v in snakemake.output if str(year) in v][0])
