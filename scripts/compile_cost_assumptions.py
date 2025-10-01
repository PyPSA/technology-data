# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8
"""
Script creates cost csv for chosen years from different source (source_dict).
The data is standardized for uniform:
    - cost years (depending on the rate of inflation )
    - technology names
    - units

Technology data from the Danish Energy Agency Technology Database are preferred.
If data are missing from all sources, these are taken from the old PyPSA cost
assumptions (with a printed warning)

The script is structured as follows:

    (1) DEA data:
        (a) read + convert units to same base
        (b) specify assumptions for  certain technologies
        (c) convert to pypsa cost syntax (investment, FOM, VOM, efficiency)
    (2) read data from other sources which need additional formatting:
        (a) old pypsa cost assumptions
        (b) Fraunhofer ISE cost assumptions
    (3) merge data from all sources for every year and save it as a csv

@author: Marta, Lisa
"""

import logging
from datetime import date

import numpy as np
import pandas as pd
from currency_converter import ECB_URL, CurrencyConverter
from scipy import interpolate

from scripts._helpers import (
    adjust_for_inflation,
    configure_logging,
    get_relative_fn,
    mock_snakemake,
    prepare_inflation_rate,
)

logger = logging.getLogger(__name__)

try:
    pd.set_option("future.no_silent_downcasting", True)
except Exception:
    pass
# ---------- sources -------------------------------------------------------
source_dict = {
    "DEA": "Danish Energy Agency",
    # solar utility
    "Vartiaien": "Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity",
    # solar rooftop
    "ETIP": "European PV Technology and Innovation Platform",
    # fuel cost
    "zappa": "Is a 100% renewable European power system feasible by 2050?",
    # co2 intensity
    "co2": "Entwicklung der spezifischen Kohlendioxid-Emissionen des deutschen Strommix in den Jahren 1990 - 2018",
    # gas pipeline costs
    "ISE": "WEGE ZU EINEM KLIMANEUTRALEN ENERGIESYSEM, Anhang zur Studie, Fraunhofer-Institut für Solare Energiesysteme ISE, Freiburg",
    # Water desalination costs
    "Caldera2016": "Caldera et al 2016: Local cost of seawater RO desalination based on solar PV and windenergy: A global estimate. (https://doi.org/10.1016/j.desal.2016.02.004)",
    "Caldera2017": "Caldera et al 2017: Learning Curve for Seawater Reverse Osmosis Desalination Plants: Capital Cost Trend of the Past, Present, and Future (https://doi.org/10.1002/2017WR021402)",
    # home battery storage and inverter investment costs
    "EWG": "Global Energy System based on 100% Renewable Energy, Energywatchgroup/LTU University, 2019",
    "HyNOW": "Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014",
    # efficiencies + lifetime SMR / SMR + CC
    "IEA": "IEA Global average levelised cost of hydrogen production by energy source and technology, 2019 and 2050 (2020), https://www.iea.org/data-and-statistics/charts/global-average-levelised-cost-of-hydrogen-production-by-energy-source-and-technology-2019-and-2050",
    # SMR capture rate
    "Timmerberg": "Hydrogen and hydrogen-derived fuels through methane decomposition of natural gas – GHG emissions and costs Timmerberg et al. (2020), https://doi.org/10.1016/j.ecmx.2020.100043",
    # geothermal (enhanced geothermal systems)
    "Aghahosseini2020": "Aghahosseini, Breyer 2020: From hot rock to useful energy: A global estimate of enhanced geothermal systems potential, https://www.sciencedirect.com/science/article/pii/S0306261920312551",
    # review of existing deep geothermal projects
    "Breede2015": "Breede et al. 2015: Overcoming challenges in the classification of deep geothermal potential, https://eprints.gla.ac.uk/169585/",
    # Study of deep geothermal systems in the Northern Upper Rhine Graben
    "Frey2022": "Frey et al. 2022: Techno-Economic Assessment of Geothermal Resources in the Variscan Basement of the Northern Upper Rhine Graben",
    # vehicles
    "vehicles": "PATHS TO A CLIMATE-NEUTRAL ENERGY SYSTEM The German energy transformation in its social context. https://www.ise.fraunhofer.de/en/publications/studies/paths-to-a-climate-neutral-energy-system.html",
}

# [DEA-sheet-names]
dea_sheet_names = {
    "onwind": "20 Onshore turbines",
    "offwind": "21 Offshore turbines",
    "solar-utility": "22 Utility-scale PV",
    "solar-utility single-axis tracking": "22 Utility-scale PV tracker",
    "solar-rooftop residential": "22 Rooftop PV residential",
    "solar-rooftop commercial": "22 Rooftop PV commercial",
    "OCGT": "52 OCGT - Natural gas",
    "CCGT": "05 Gas turb. CC, steam extract.",
    "oil": "50 Diesel engine farm",
    "biomass CHP": "09c Straw, Large, 40 degree",
    "biomass EOP": "09c Straw, Large, 40 degree",
    "biomass HOP": "09c Straw HOP",
    "central coal CHP": "01 Coal CHP",
    "central gas CHP": "04 Gas turb. simple cycle, L",
    "central gas CHP CC": "04 Gas turb. simple cycle, L",
    "central solid biomass CHP": "09a Wood Chips, Large 50 degree",
    "central solid biomass CHP CC": "09a Wood Chips, Large 50 degree",
    "central solid biomass CHP powerboost CC": "09a Wood Chips, Large 50 degree",
    "central air-sourced heat pump": "40 Comp. hp, airsource 3 MW",
    "central geothermal heat source": "45.1.b Geothermal DH, 2000m, E",
    "central excess-heat-sourced heat pump": "40 Comp. hp, excess heat 10 MW",
    "central ground-sourced heat pump": "40 Absorption heat pump, DH",
    "central resistive heater": "41 Electric Boilers",
    "central gas boiler": "44 Natural Gas DH Only",
    "decentral gas boiler": "202 Natural gas boiler",
    "direct firing gas": "312.a Direct firing Natural Gas",
    "direct firing gas CC": "312.a Direct firing Natural Gas",
    "direct firing solid fuels": "312.b Direct firing Sold Fuels",
    "direct firing solid fuels CC": "312.b Direct firing Sold Fuels",
    "decentral ground-sourced heat pump": "207.7 Ground source existing",
    "decentral air-sourced heat pump": "207.3 Air to water existing",
    "central water pit storage": "140 PTES seasonal",
    "central water tank storage": "141 Large hot water tank",
    "decentral water tank storage": "142 Small scale hot water tank",
    "fuel cell": "12 LT-PEMFC CHP",
    "hydrogen storage underground": "151c Hydrogen Storage - Caverns",
    "hydrogen storage tank type 1 including compressor": "151a Hydrogen Storage - Tanks",
    "micro CHP": "219 LT-PEMFC mCHP - natural gas",
    "biogas": "81 Biogas, Basic plant, small",
    "biogas CC": "81 Biogas, Basic plant, small",
    "biogas upgrading": "82 Upgrading 3,000 Nm3 per h",
    "battery": "180 Lithium Ion Battery",
    "industrial heat pump medium temperature": "302.a High temp. hp Up to 125 C",
    "industrial heat pump high temperature": "302.b High temp. hp Up to 150",
    "electric boiler steam": "310.1 Electric boiler steam  ",
    "gas boiler steam": "311.1c Steam boiler Gas",
    "solid biomass boiler steam": "311.1e Steam boiler Wood",
    "solid biomass boiler steam CC": "311.1e Steam boiler Wood",
    "biomass boiler": "204 Biomass boiler, automatic",
    "electrolysis": "86 AEC 100 MW",
    "direct air capture": "403.a Direct air capture",
    "biomass CHP capture": "401.a Post comb - small CHP",
    "cement capture": "401.c Post comb - Cement kiln",
    "BioSNG": "84 Gasif. CFB, Bio-SNG",
    "BtL": "85 Gasif. Ent. Flow FT, liq fu ",
    "biomass-to-methanol": "97 Methanol from biomass gasif.",
    "biogas plus hydrogen": "99 SNG from methan. of biogas",
    "methanolisation": "98 Methanol from hydrogen",
    "Fischer-Tropsch": "102 Hydrogen to Jet",
    "central hydrogen CHP": "12 LT-PEMFC CHP",
    "Haber-Bosch": "103 Hydrogen to Ammonia",
    "air separation unit": "103 Hydrogen to Ammonia",
    "waste CHP": "08 WtE CHP, Large, 50 degree",
    "waste CHP CC": "08 WtE CHP, Large, 50 degree",
    "biochar pyrolysis": "105 Slow pyrolysis, Straw",
    "electrolysis small": "86 AEC 10 MW",
    "gas storage": "150 Underground Storage of Gas",
}
# [DEA-sheet-names]

uncrtnty_lookup = {
    "onwind": "J:K",
    "offwind": "J:K",
    "solar-utility": "J:K",
    "solar-utility single-axis tracking": "J:K",
    "solar-rooftop residential": "J:K",
    "solar-rooftop commercial": "J:K",
    "OCGT": "I:J",
    "CCGT": "I:J",
    "oil": "I:J",
    "biomass CHP": "I:J",
    "biomass EOP": "I:J",
    "biomass HOP": "I:J",
    "central coal CHP": "",
    "central gas CHP": "I:J",
    "central gas CHP CC": "I:J",
    "central hydrogen CHP": "I:J",
    "central solid biomass CHP": "I:J",
    "central solid biomass CHP CC": "I:J",
    "central solid biomass CHP powerboost CC": "I:J",
    "solar": "",
    "central air-sourced heat pump": "J:K",
    "central geothermal heat source": "H:K",
    "central excess-heat-sourced heat pump": "H:K",
    "central ground-sourced heat pump": "I:J",
    "central resistive heater": "I:J",
    "central gas boiler": "I:J",
    "decentral gas boiler": "I:J",
    "direct firing gas": "H:I",
    "direct firing gas CC": "H:I",
    "direct firing solid fuels": "H:I",
    "direct firing solid fuels CC": "H:I",
    "decentral ground-sourced heat pump": "I:J",
    "decentral air-sourced heat pump": "I:J",
    "central water pit storage": "I:L",
    "central water tank storage": "J:K",
    "decentral water tank storage": "J:K",
    "fuel cell": "I:J",
    "hydrogen storage underground": "J:K",
    "hydrogen storage tank type 1 including compressor": "J:K",
    "micro CHP": "I:J",
    "biogas": "I:J",
    "biogas CC": "I:J",
    "biogas upgrading": "I:J",
    "electrolysis": "I:J",
    "battery": "H,K",
    "direct air capture": "I:J",
    "cement capture": "I:J",
    "biomass CHP capture": "I:J",
    "BioSNG": "I:J",
    "BtL": "J:K",
    "biomass-to-methanol": "J:K",
    "biogas plus hydrogen": "J:K",
    "industrial heat pump medium temperature": "H:I",
    "industrial heat pump high temperature": "H:I",
    "electric boiler steam": "H:I",
    "gas boiler steam": "H:I",
    "solid biomass boiler steam": "H:I",
    "solid biomass boiler steam CC": "H:I",
    "biomass boiler": "I:J",
    "Fischer-Tropsch": "I:J",
    "Haber-Bosch": "I:J",
    "air separation unit": "I:J",
    "methanolisation": "J:K",
    "waste CHP": "I:J",
    "waste CHP CC": "I:J",
    "biochar pyrolysis": "J:K",
    "biomethanation": "J:K",
    "electrolysis small": "I:J",
    "gas storage": "",
}

# since February 2022 DEA uses a new format for the technology data
# all Excel sheets of updated technologies have a different layout and are
# given in EUR_2020 money (instead of EUR_2015)
cost_year_2020 = [
    "solar-utility",
    "solar-utility single-axis tracking",
    "solar-rooftop residential",
    "solar-rooftop commercial",
    "offwind",
    "electrolysis",
    "biogas",
    "biogas CC",
    "biogas upgrading",
    "direct air capture",
    "biomass CHP capture",
    "cement capture",
    "BioSNG",
    "BtL",
    "biomass-to-methanol",
    "biogas plus hydrogen",
    "methanolisation",
    "Fischer-Tropsch",
    "biochar pyrolysis",
    "biomethanation",
    "electrolysis small",
    "central water pit storage",
    "central water tank storage",
    "decentral water tank storage",
    "hydrogen storage underground",
    "hydrogen storage tank type 1 including compressor",
    "battery",
    "gas storage",
]

manual_cost_year_assignments_2020 = [
    "central water pit charger",
    "central water pit discharger",
    "central water tank charger",
    "central water tank discharger",
    "decentral water tank charger",
    "decentral water tank discharger",
    "battery storage",
    "battery inverter",
    "gas storage charger",
    "gas storage discharger",
    "gas storage discharger",
]

cost_year_2019 = [
    "direct firing gas",
    "direct firing gas CC",
    "direct firing solid fuels",
    "direct firing solid fuels CC",
    "industrial heat pump medium temperature",
    "industrial heat pump high temperature",
    "electric boiler steam",
    "gas boiler steam",
    "solid biomass boiler steam",
    "solid biomass boiler steam CC",
]


# -------- FUNCTIONS ---------------------------------------------------


def get_excel_sheets(list_of_excel_files: list) -> dict:
    """
    The function reads Excel files and returns them in a dictionary.
    The dictionary has the files names as keys and the lists of sheet names as values.

    Parameters
    ----------
    list_of_excel_files : list
        Excel files to process

    Returns
    -------
    Dictionary
        data from DEA
    """

    excel_sheets_dictionary = {}
    for entry in list_of_excel_files:
        if entry[-5:] == ".xlsx":
            excel_sheets_dictionary[entry] = pd.ExcelFile(entry).sheet_names
    logger.info(f"found {len(excel_sheets_dictionary)} excel sheets: ")
    for key in excel_sheets_dictionary.keys():
        logger.info(f"* {key}")
    return excel_sheets_dictionary


def get_sheet_location(
    tech_name: str, sheet_names_dict: dict, input_data_dict: dict
) -> str:
    """
    The function returns a dictionary. The dictionary has the technology names as keys and
    the Excel file names where the technology is saved as values

    Parameters
    ----------
    tech_name : str
        technology name
    sheet_names_dict : dict
        dictionary having the technology name as keys and Excel sheet names as values
    input_data_dict : dict
        dictionary having the files names as keys and the lists of sheet names as values

    Returns
    -------
    str
        Excel file name where the technology is present
    """

    key_list = [
        key
        for key, value in input_data_dict.items()
        if any(sheet_names_dict[tech_name] in s for s in value)
    ]

    if len(key_list) == 1:
        return key_list[0]
    elif len(key_list) > 1:
        logger.info(f"{tech_name} appears in more than one sheet name")
        return "Multiple sheets found"
    else:
        logger.info(
            f"tech {tech_name} with sheet name {sheet_names_dict[tech_name]} not found in excel sheets. "
        )
        return "Sheet not found"


def get_dea_maritime_data(
    fn: str, years: list, input_data_df: pd.DataFrame
) -> pd.DataFrame:
    """
    The function returns a dataframe containing the technology data for shipping from the DEA database.

    Parameters
    ----------
    fn : str
        path to DEA input data file for shipping
    years : list
        years for which a cost assumption is provided
    input_data_df : pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        technology data cost assumptions enriched with shipping data from DEA
    """

    dea_maritime_data_sheet_names = [
        "Container feeder, diesel",
        "Container feeder, methanol",
        "Container feeder, ammonia",
        "Container, diesel",
        "Container, methanol",
        "Container, ammonia",
        "Tank&bulk, diesel",
        "Tank&bulk, methanol",
        "Tankbulk, ammonia",
    ]

    excel = pd.read_excel(
        fn,
        sheet_name=dea_maritime_data_sheet_names,
        index_col=[0, 1],
        usecols="A:F",
        na_values="N/A",
        engine="calamine",
    )

    wished_index = [
        "Typical ship lifetime (years)",
        "Upfront ship cost (mill. €)",
        "Fixed O&M (€/year)",
        "Variable O&M (€/nm)",
    ]

    for sheet in excel.keys():
        df = excel[sheet]
        df = df.iloc[1:, :].set_axis(df.iloc[0], axis=1)

        assert "Typical operational speed" in df.index.get_level_values(1)[22]
        # in unit GJ/nm
        efficiency = df.iloc[22]

        df = df[df.index.get_level_values(1).isin(wished_index)]
        df = df.droplevel(level=0)
        df.loc["efficiency (GJ/nm)"] = efficiency
        df = df.reindex(columns=pd.Index(years).union(df.columns))
        df = df.astype(float)
        df = df.interpolate(axis=1, limit_direction="both")
        df = df[years]

        # dropna
        df = df.dropna(how="all", axis=0)
        # add column for units
        df["unit"] = df.rename(
            index=lambda x: x[x.rfind("(") + 1 : x.rfind(")")]
        ).index.values
        df["unit"] = df.unit.str.replace("€", "EUR")
        # remove units from index
        df.index = df.index.str.replace(r" \(.*\)", "", regex=True)

        # convert million Euro -> Euro
        df_i = df[df.unit == "mill. EUR"].index
        df.loc[df_i, years] *= 1e6
        df.loc[df_i, "unit"] = "EUR"

        # convert FOM in % of investment/year
        if "Fixed O&M" in df.index:
            df.loc["Fixed O&M", years] /= df.loc["Upfront ship cost", years] * 100
            df.loc["Fixed O&M", "unit"] = "%/year"

        # convert nm in km
        # 1 Nautical Mile (nm) = 1.852 Kilometers (km)
        df_i = df[df.unit.str.contains("/nm")].index
        df.loc[df_i, years] /= 1.852
        df.loc[df_i, "unit"] = df.loc[df_i, "unit"].str.replace("/nm", "/km")

        # 1 GJ = 1/3600 * 1e9 Wh = 1/3600 * 1e3 MWh
        df_i = df[df.unit.str.contains("GJ")].index
        df.loc[df_i, years] *= 1e3 / 3600
        df.loc[df_i, "unit"] = df.loc[df_i, "unit"].str.replace("GJ", "MWh")

        # add source + cost year
        df["source"] = f"Danish Energy Agency, {get_relative_fn(fn)}"
        # cost year is 2023 p.10
        df["currency_year"] = 2023
        # add sheet name
        df["further description"] = sheet

        # FOM, VOM,efficiency, lifetime, investment
        rename = {
            "Typical ship lifetime": "lifetime",
            "Upfront ship cost": "investment",
            "Fixed O&M": "FOM",
            "Variable O&M": "VOM",
        }

        df = df.rename(index=rename)

        df = pd.concat([df], keys=[sheet], names=["technology", "parameter"])

        input_data_df = pd.concat([input_data_df, df])

    return input_data_df


def get_dea_vehicle_data(
    fn: str, years: list, technology_dataframe: pd.DataFrame
) -> pd.DataFrame:
    """
    The function gets heavy-duty vehicle data from DEA.

    Parameters
    ----------
    fn : str
        path to DEA input data file for shipping
    years : list
        years for which a cost assumption is provided
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        technology data cost assumptions enriched with shipping data from DEA
    """

    dea_vehicle_data_sheet_names = [
        "Diesel L1",
        "Diesel L2",
        "Diesel L3",
        "Diesel B1",
        "Diesel B2",
        "BEV L1",
        "BEV L2",
        "BEV L3",
        "BEV B1",
        "BEV B2",
        "FCV L1",
        "FCV L2",
        "FCV L3",
        "FCV B1",
        "FCV B2",
    ]
    excel = pd.read_excel(
        fn,
        sheet_name=dea_vehicle_data_sheet_names,
        index_col=0,
        usecols="A:F",
        na_values="no data",
        engine="calamine",
    )

    wished_index = [
        "Typical vehicle lifetime (years)",
        "Upfront vehicle cost (€)",
        "Fixed maintenance cost (€/year)",
        "Variable maintenance cost (€/km)",
        "Motor size (kW)",
    ]

    # clarify DEA names
    types = {
        "L1": "Truck Solo max 26 tons",
        "L2": "Truck Trailer max 56 tons",
        "L3": "Truck Semi-Trailer max 50 tons",
        "B1": "Bus city",
        "B2": "Coach",
    }

    for sheet in excel.keys():
        df = excel[sheet]
        tech = sheet.split()[0] + " " + types.get(sheet.split()[1], "")
        df = df.iloc[1:, :].set_axis(df.iloc[0], axis=1)
        # "Fuel energy - typical load (MJ/km)"
        # represents efficiency for average weight vehicle carries during normal
        # operation, currently assuming mean between urban, regional and long haul
        assert df.index[27] == "Fuel energy - typical load (MJ/km)"
        efficiency = df.iloc[28:31].mean()
        df = df[df.index.isin(wished_index)]
        df.loc["efficiency (MJ/km)"] = efficiency
        df = df.reindex(columns=pd.Index(years).union(df.columns))
        df = df.interpolate(axis=1, limit_direction="both")
        df = df[years]

        # add column for units
        df["unit"] = df.rename(
            index=lambda x: x[x.rfind("(") + 1 : x.rfind(")")]
        ).index.values
        df["unit"] = df.unit.str.replace("€", "EUR")
        # remove units from index
        df.index = df.index.str.replace(r" \(.*\)", "", regex=True)

        # convert MJ in kWh -> 1 kWh = 3.6 MJ
        df_i = df.index[df.unit == "MJ/km"]
        df.loc[df_i, years] /= 3.6
        df.loc[df_i, "unit"] = "kWh/km"

        # convert FOM in % of investment/year
        df.loc["Fixed maintenance cost", years] /= (
            df.loc["Upfront vehicle cost", years] * 100
        )
        df.loc["Fixed maintenance cost", "unit"] = "%/year"

        # clarify costs are per vehicle
        df.loc["Upfront vehicle cost", "unit"] += "/vehicle"

        # add source + cost year
        df["source"] = f"Danish Energy Agency, {get_relative_fn(fn)}"
        # cost year is 2022 p.12
        df["currency_year"] = 2022
        # add sheet name
        df["further description"] = sheet

        # FOM, VOM,efficiency, lifetime, investment
        rename = {
            "Typical vehicle lifetime": "lifetime",
            "Upfront vehicle cost": "investment",
            "Fixed maintenance cost": "FOM",
            "Variable maintenance cost": "VOM",
        }

        df = df.rename(index=rename)

        to_keep = ["Motor size", "lifetime", "FOM", "VOM", "efficiency", "investment"]
        df = df[df.index.isin(to_keep)]

        df = pd.concat([df], keys=[tech], names=["technology", "parameter"])

        technology_dataframe = pd.concat([technology_dataframe, df])

    return technology_dataframe


def get_data_DEA(
    years: list,
    tech_name: str,
    sheet_names_dict: dict,
    input_data_dict: dict,
    offwind_no_grid_costs_flag: bool = True,
    expectation: str = None,
) -> pd.DataFrame:
    """
    The function interpolates costs for a given technology from DEA database sheet and
    stores technology data from DEA in a dictionary.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    tech_name : str
        technology name
    sheet_names_dict : dict
        dictionary having the technology name as keys and Excel sheet names as values
    input_data_dict : dict
        dictionary where the keys are the path to the DEA inputs and the values are the sheet names
    offwind_no_grid_costs_flag : bool
        flag to remove grid connection costs from DEA for offwind. Such costs are calculated separately in pypsa-eur
    expectation : str
        tech data uncertainty. The possible options are [None, "optimist", "pessimist"]

    Returns
    -------
    pandas.DataFrame
        technology data from DEA
    """

    excel_file = get_sheet_location(tech_name, sheet_names_dict, input_data_dict)
    if excel_file == "Sheet not found" or excel_file == "Multiple sheets found":
        logger.info(f"excel file not found for technology: {tech_name}")
        return pd.DataFrame()

    if tech_name in [
        "direct air capture",
        "cement capture",
        "biomass CHP capture",
    ]:
        usecols = "A:F"
    elif tech_name in [
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
    elif tech_name in [
        "Fischer-Tropsch",
        "Haber-Bosch",
        "air separation unit",
        "gas storage",
    ]:
        usecols = "B:F"
    elif tech_name in [
        "central water pit storage",
    ]:
        usecols = "B:H"
    else:
        usecols = "B:G"

    usecols += f",{uncrtnty_lookup[tech_name]}"

    if (
        (tech_name in cost_year_2019)
        or (tech_name in cost_year_2020)
        or ("renewable_fuels" in excel_file)
    ):
        skiprows = [0]
    else:
        skiprows = [0, 1]

    excel = pd.read_excel(
        excel_file,
        sheet_name=sheet_names_dict[tech_name],
        index_col=0,
        usecols=usecols,
        skiprows=skiprows,
        na_values="N.A",
        engine="calamine",
    )

    excel.dropna(axis=1, how="all", inplace=True)

    excel.index = excel.index.fillna(" ")
    excel.index = excel.index.astype(str)
    excel.dropna(axis=0, how="all", inplace=True)

    if 2020 not in excel.columns:
        selection = excel[excel.isin([2020])].dropna(how="all").index
        excel.columns = excel.loc[selection].iloc[0, :].fillna("Technology", limit=1)
        excel.drop(selection, inplace=True)

    uncertainty_columns = ["2050-optimist", "2050-pessimist"]
    if uncrtnty_lookup[tech_name]:
        # hydrogen storage sheets have reverse order of lower/upper estimates
        if tech_name in [
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
        # drop duplicates
        excel = excel[~excel.index.duplicated()]
        excel.loc[:, 2050] = excel.loc[:, f"2050-{expectation}"].combine_first(
            excel.loc[:, 2050]
        )
    excel.drop(columns=uncertainty_columns, inplace=True)

    if expectation:
        excel = excel.loc[:, [2020, 2050]]

    parameters = [
        "efficiency",
        "investment",
        "Fixed O&M",
        "Variable O&M",
        "production capacity for one unit",
        "Output capacity expansion cost",
        "Hydrogen Output",
        "Hydrogen (% total input_e (MWh / MWh))",
        "Hydrogen [% total input_e",
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
        "FT Liquids Output, [MWh/MWh Total Input]",
        " - hereof recoverable for district heating [%-points of heat loss]",
        " - hereof recoverable for district heating (%-points of heat loss)",
        "Bio SNG Output [% of fuel input]",
        "Methanol Output",
        "District heat  Output",
        "Electricity Output",
        "Total O&M",
        "Biochar Output",  # biochar pyrolysis
        "Pyrolysis oil Output",  # biochar pyrolysis
        "Pyrolysis gas Output",  # biochar pyrolysis
        "Heat Output",  # biochar pyrolysis
        "Specific energy content [GJ/ton] biochar",  # biochar pyrolysis
        "Electricity Consumption",
        "Feedstock Consumption",  # biochar pyrolysis
        "Methane Output",
        "CO2 Consumption",
        "Hydrogen Consumption",
        " - of which is equipment excluding heat pump",
        " - of which is heat pump including its installation",
        "Input capacity",
        "Output capacity",
        "Energy storage capacity",
        "Typical temperature difference in storage [hot/cold, K]",
        "Max. storage temperature, hot",
        "Storage temperature, discharged",
        "Energy losses during storage",
    ]

    # this is not good at all but requires significant changes to `test_compile_cost_assumptions` otherwise
    if tech_name == "central geothermal heat source":
        parameters += [
            " - of which is installation",
            "Heat generation capacity for one unit (MW)",
            "Heat generation from geothermal heat (MJ/s)",
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
    df = df.apply(
        lambda row: row.apply(
            lambda x: (
                (float(x.split("-")[0]) + float(x.split("-")[1])) / 2
                if isinstance(x, str) and "-" in x
                else x
            )
        ),
        axis=1,
    )

    # remove symbols "~", ">", "<" and " "
    for sym in ["~", ">", "<", " "]:
        df = df.apply(
            lambda col: col.apply(
                lambda x: x.replace(sym, "") if isinstance(x, str) else x
            )
        )

    df = df.astype(float)
    df = df.mask(
        df.apply(pd.to_numeric, errors="coerce").isnull(),
        df.astype(str).apply(lambda x: x.str.strip()),
    )

    # Modify data loaded from DEA on a per-technology case
    if (tech_name == "offwind") and offwind_no_grid_costs_flag:
        df.loc["Nominal investment (*total) [MEUR/MW_e, 2020]"] -= excel.loc[
            "Nominal investment (installation: grid connection) [M€/MW_e, 2020]"
        ]

    # Exclude indirect costs for centralised system with additional piping.
    if tech_name.startswith("industrial heat pump"):
        df = df.drop("Indirect investments cost (MEUR per MW)")

    if tech_name == "biogas plus hydrogen":
        df.drop(df.loc[df.index.str.contains("GJ SNG")].index, inplace=True)

    if tech_name == "BtL":
        df.drop(df.loc[df.index.str.contains("1,000 t FT Liquids")].index, inplace=True)

    if tech_name == "biomass-to-methanol":
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech_name == "methanolisation":
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech_name == "Fischer-Tropsch":
        df.drop(df.loc[df.index.str.contains("l FT Liquids")].index, inplace=True)

    if tech_name == "biomass boiler":
        df.drop(
            df.loc[df.index.str.contains("Possible additional")].index, inplace=True
        )
        df.drop(df.loc[df.index.str.contains("Total efficiency")].index, inplace=True)

    if tech_name == "Haber-Bosch":
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

    if tech_name == "air separation unit":
        divisor = (
            (df.loc["Specific investment mark-up factor optional ASU"] - 1.0)
            / excel.loc["N2 Consumption, [t/t] Ammonia"]
        ).astype(float)

        # Calculate ASU cost separate to HB facility in terms of t N2 output
        df.loc[
            [
                "Specific investment [MEUR /TPD Ammonia output]",
                "Fixed O&M [kEUR /TPD Ammonia]",
                "Variable O&M [EUR /t Ammonia]",
            ]
        ] *= divisor
        # Convert output to hourly generation
        df.loc[
            [
                "Specific investment [MEUR /TPD Ammonia output]",
                "Fixed O&M [kEUR /TPD Ammonia]",
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
                    "Specific investment [MEUR /MW Ammonia output]", regex=False
                )
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Fixed O&M [kEUR/MW Ammonia/year]", regex=False)
            ].index,
            inplace=True,
        )
        df.drop(
            df.loc[
                df.index.str.contains("Variable O&M [EUR/MWh Ammonia]", regex=False)
            ].index,
            inplace=True,
        )

    if "solid biomass power" in tech_name:
        df.index = df.index.str.replace("EUR/MWeh", "EUR/MWh")

    if "biochar pyrolysis" in tech_name:
        df = biochar_pyrolysis_harmonise_dea(df)

    elif tech_name == "central geothermal heat source":
        # we need to convert from costs per MW of the entire system (including heat pump)
        # to costs per MW of the geothermal heat source only
        # heat_source_costs [MEUR/MW_heat_source] = heat_source_costs [MEUR/MW_entire_system] * MW_entire_system / MW_heat_source
        df.loc["Nominal investment (MEUR per MW)"] = (
            (
                df.loc[" - of which is equipment excluding heat pump"]
                + df.loc[" - of which is installation"]
            )
            * df.loc["Heat generation capacity for one unit (MW)"]
            / df.loc["Heat generation from geothermal heat (MJ/s)"]
        )

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
    df_final = df_final.ffill(axis=1)

    df_final["source"] = f"{source_dict['DEA']}, {get_relative_fn(excel_file)}"
    if (
        tech_name in cost_year_2020
        and ("for_carbon_capture_transport_storage" not in excel_file)
        and ("renewable_fuels" not in excel_file)
        and ("for_energy_storage" not in excel_file)
    ):
        for attr in ["investment", "Fixed O&M"]:
            to_drop = df[
                df.index.str.contains(attr) & ~df.index.str.contains(r"\(\*total\)")
            ].index
            df_final.drop(to_drop, inplace=True)

        df_final["unit"] = df_final.rename(
            index=lambda x: x[x.rfind("[") + 1 : x.rfind("]")]
        ).index.values
    else:
        df_final.index = df_final.index.str.replace(r"\[", "(", regex=True).str.replace(
            r"\]", ")", regex=True
        )
        df_final["unit"] = df_final.rename(
            index=lambda x: x[x.rfind("(") + 1 : x.rfind(")")]
        ).index.values
    df_final.index = df_final.index.str.replace(r" \(.*\)", "", regex=True)

    return df_final


def add_desalination_data(cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds technology data for seawater desalination (SWRO) and water storage.

    Parameters
    ----------
    cost_dataframe : pandas.DataFrame
        cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    # Interpolate cost based on historic costs/cost projection to fitting year
    cs = [2070, 1917, 1603, 1282, 1025]  # in USD/(m^3/d)
    ys = [2015, 2022, 2030, 2040, 2050]
    c = np.interp(year, ys, cs)
    c *= 24  # in USD/(m^3/h)
    c /= 1.17  # in EUR/(m^3/h)

    tech_name = "seawater desalination"

    cost_dataframe.loc[(tech_name, "investment"), "value"] = c
    cost_dataframe.loc[(tech_name, "investment"), "unit"] = "EUR/(m^3-H2O/h)"
    cost_dataframe.loc[(tech_name, "investment"), "source"] = (
        source_dict["Caldera2017"] + ", Table 4."
    )
    cost_dataframe.loc[(tech_name, "investment"), "currency_year"] = 2015

    cost_dataframe.loc[(tech_name, "FOM"), "value"] = 4.0
    cost_dataframe.loc[(tech_name, "FOM"), "unit"] = "%/year"
    cost_dataframe.loc[(tech_name, "FOM"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )
    cost_dataframe.loc[(tech_name, "FOM"), "currency_year"] = 2015

    cost_dataframe.loc[(tech_name, "FOM"), "value"] = 4.0
    cost_dataframe.loc[(tech_name, "FOM"), "unit"] = "%/year"
    cost_dataframe.loc[(tech_name, "FOM"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )

    cost_dataframe.loc[(tech_name, "lifetime"), "value"] = 30
    cost_dataframe.loc[(tech_name, "lifetime"), "unit"] = "years"
    cost_dataframe.loc[(tech_name, "lifetime"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )

    salinity = snakemake.config["desalination"]["salinity"]
    cost_dataframe.loc[(tech_name, "electricity-input"), "value"] = (
        0.0003 * salinity**2 + 0.0018 * salinity + 2.6043
    )
    cost_dataframe.loc[(tech_name, "electricity-input"), "unit"] = "kWh/m^3-H2O"
    cost_dataframe.loc[(tech_name, "electricity-input"), "source"] = (
        source_dict["Caldera2016"] + ", Fig. 4."
    )

    tech_name = "clean water tank storage"
    cost_dataframe.loc[(tech_name, "investment"), "value"] = 65
    cost_dataframe.loc[(tech_name, "investment"), "unit"] = "EUR/m^3-H2O"
    cost_dataframe.loc[(tech_name, "investment"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )
    cost_dataframe.loc[(tech_name, "investment"), "currency_year"] = 2013

    cost_dataframe.loc[(tech_name, "FOM"), "value"] = 2
    cost_dataframe.loc[(tech_name, "FOM"), "unit"] = "%/year"
    cost_dataframe.loc[(tech_name, "FOM"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )
    cost_dataframe.loc[(tech_name, "FOM"), "currency_year"] = 2013

    cost_dataframe.loc[(tech_name, "lifetime"), "value"] = 30
    cost_dataframe.loc[(tech_name, "lifetime"), "unit"] = "years"
    cost_dataframe.loc[(tech_name, "lifetime"), "source"] = (
        source_dict["Caldera2016"] + ", Table 1."
    )

    return cost_dataframe


def add_co2_intensity(cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds CO2 intensity for the carriers.

    Parameters
    ----------
    cost_dataframe : pandas.DataFrame
        cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    TJ_to_MWh = 277.78
    cost_dataframe.loc[("gas", "CO2 intensity"), "value"] = (
        55827 / 1e3 / TJ_to_MWh
    )  # Erdgas
    cost_dataframe.loc[("coal", "CO2 intensity"), "value"] = (
        93369 / 1e3 / TJ_to_MWh
    )  # Steinkohle
    cost_dataframe.loc[("lignite", "CO2 intensity"), "value"] = (
        113031 / 1e3 / TJ_to_MWh
    )  # Rohbraunkohle Rheinland
    cost_dataframe.loc[("oil", "CO2 intensity"), "value"] = (
        74020 / 1e3 / TJ_to_MWh
    )  # Heizöl, leicht
    cost_dataframe.loc[("methanol", "CO2 intensity"), "value"] = (
        0.2482  # t_CO2/MWh_th, based on stochiometric composition.
    )
    cost_dataframe.loc[("solid biomass", "CO2 intensity"), "value"] = 0.3

    oil_specific_energy = 44  # GJ/t
    CO2_CH2_mass_ratio = 44 / 14  # kg/kg (1 mol per mol)
    CO2_C_mass_ratio = 44 / 12  # kg/kg
    methane_specific_energy = 50  # GJ/t
    CO2_CH4_mass_ratio = 44 / 16  # kg/kg (1 mol per mol)
    biomass_specific_energy = 18  # GJ/t LHV
    biomass_carbon_content = 0.5
    cost_dataframe.loc[("oil", "CO2 intensity"), "value"] = (
        (1 / oil_specific_energy) * 3.6 * CO2_CH2_mass_ratio
    )  # tCO2/MWh
    cost_dataframe.loc[("gas", "CO2 intensity"), "value"] = (
        (1 / methane_specific_energy) * 3.6 * CO2_CH4_mass_ratio
    )  # tCO2/MWh
    cost_dataframe.loc[("solid biomass", "CO2 intensity"), "value"] = (
        biomass_carbon_content * (1 / biomass_specific_energy) * 3.6 * CO2_C_mass_ratio
    )  # tCO2/MWh

    cost_dataframe.loc[("oil", "CO2 intensity"), "source"] = (
        "Stoichiometric calculation with 44 GJ/t diesel and -CH2- approximation of diesel"
    )
    cost_dataframe.loc[("gas", "CO2 intensity"), "source"] = (
        "Stoichiometric calculation with 50 GJ/t CH4"
    )
    cost_dataframe.loc[("solid biomass", "CO2 intensity"), "source"] = (
        "Stoichiometric calculation with 18 GJ/t_DM LHV and 50% C-content for solid biomass"
    )
    cost_dataframe.loc[("coal", "CO2 intensity"), "source"] = source_dict["co2"]
    cost_dataframe.loc[("lignite", "CO2 intensity"), "source"] = source_dict["co2"]

    cost_dataframe.loc[pd.IndexSlice[:, "CO2 intensity"], "unit"] = "tCO2/MWh_th"

    return cost_dataframe


def add_solar_from_other(years: list, cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds solar from other sources than DEA (since the lifetime assumed in
    DEA is very optimistic).

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    cost_dataframe : pandas.DataFrame
        costs

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    # solar utility from Vartiaian 2019
    interpolated_data = np.interp(
        x=years, xp=[2020, 2030, 2040, 2050], fp=[431, 275, 204, 164]
    )
    # the paper says 'In this report, all results are given in real 2019
    # money.'
    interpolated_data = interpolated_data / (
        1 + snakemake.config["rate_inflation"]
    ) ** (2019 - snakemake.config["eur_year"])
    solar_uti = pd.Series(data=interpolated_data, index=years)

    # solar rooftop from ETIP 2019
    interpolated_data = np.interp(x=years, xp=[2020, 2030, 2050], fp=[1150, 800, 550])
    # using 2016 money in page 10
    interpolated_data = interpolated_data / (
        1 + snakemake.config["rate_inflation"]
    ) ** (2016 - snakemake.config["eur_year"])
    solar_roof = pd.Series(data=interpolated_data, index=years)

    # solar utility from Vartiaian 2019
    if snakemake.config["solar_utility_from_vartiaien"]:
        cost_dataframe.loc[("solar-utility", "investment"), "value"] = solar_uti[year]
        cost_dataframe.loc[("solar-utility", "investment"), "source"] = source_dict[
            "Vartiaien"
        ]
        cost_dataframe.loc[("solar-utility", "investment"), "currency_year"] = 2019

        cost_dataframe.loc[("solar-utility", "lifetime"), "value"] = 30
        cost_dataframe.loc[("solar-utility", "lifetime"), "source"] = source_dict[
            "Vartiaien"
        ]
        cost_dataframe.loc[("solar-utility", "lifetime"), "currency_year"] = 2019

    if snakemake.config["solar_rooftop_from_etip"]:
        # solar rooftop from ETIP 2019
        cost_dataframe.loc[("solar-rooftop", "investment"), "value"] = solar_roof[year]
        cost_dataframe.loc[("solar-rooftop", "investment"), "source"] = source_dict[
            "ETIP"
        ]
        cost_dataframe.loc[("solar-rooftop", "investment"), "currency_year"] = 2019

        cost_dataframe.loc[("solar-rooftop", "lifetime"), "value"] = 30
        cost_dataframe.loc[("solar-rooftop", "lifetime"), "source"] = source_dict[
            "ETIP"
        ]
        cost_dataframe.loc[("solar-rooftop", "lifetime"), "currency_year"] = 2019

    # lifetime & efficiency for solar
    cost_dataframe.loc[("solar", "lifetime"), "value"] = cost_dataframe.loc[
        (["solar-rooftop", "solar-utility"], "lifetime"), "value"
    ].mean()
    cost_dataframe.loc[("solar", "lifetime"), "unit"] = "years"
    cost_dataframe.loc[("solar", "lifetime"), "currency_year"] = 2019
    cost_dataframe.loc[("solar", "lifetime"), "source"] = (
        "Assuming 50% rooftop, 50% utility"
    )

    return cost_dataframe


# [add-h2-from-other]
def add_h2_from_other(cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function assumes higher efficiency for electrolysis (0.8) and fuel cell (0.58).

    Parameters
    ----------
    cost_dataframe : pandas.DataFrame
        costs

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    cost_dataframe.loc[("electrolysis", "efficiency"), "value"] = 0.8
    cost_dataframe.loc[("fuel cell", "efficiency"), "value"] = 0.58
    cost_dataframe.loc[("electrolysis", "efficiency"), "source"] = "budischak2013"
    cost_dataframe.loc[("electrolysis", "efficiency"), "currency_year"] = 2013
    cost_dataframe.loc[("fuel cell", "efficiency"), "source"] = "budischak2013"
    cost_dataframe.loc[("fuel cell", "efficiency"), "currency_year"] = 2013

    return cost_dataframe


# [unify-diw-inflation]
def unify_diw(cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds currency year for the DIW costs from 2010.

    Parameters
    ----------
    cost_dataframe : pandas.DataFrame
        costs

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    cost_dataframe.loc[("PHS", "investment"), "currency_year"] = 2010
    cost_dataframe.loc[("ror", "investment"), "currency_year"] = 2010
    cost_dataframe.loc[("hydro", "investment"), "currency_year"] = 2010

    return cost_dataframe


def biochar_pyrolysis_harmonise_dea(df: pd.DataFrame) -> pd.DataFrame:
    """
    The function harmonises biochar and pyrolysis costs.

    Parameters
    ----------
    df : pandas.DataFrame
        costs

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    # data for 2020 not available
    if 2020 in df.columns:
        df.drop(columns=2020, inplace=True)
    # normalize biochar and total heat output to feedstock input
    idx = df.index.str.contains("Total Input")
    idx2 = df.index.str.contains("Feedstock Consumption")
    df.loc[idx] = df.loc[idx].astype(float) / df.loc[idx2].values.astype(float)
    df.index = df.index.str.replace("Total Input", "feedstock")

    # all pyrolysis product except char are combusted for heat
    df_sum = pd.concat(
        (
            df.iloc[df.index.str.contains("Pyrolysis oil Output")],
            df.iloc[df.index.str.contains("Pyrolysis gas Output")],
            df.iloc[df.index.str.contains("Heat Output")],
        ),
        axis=0,
    ).sum(axis=0, skipna=False)
    df.iloc[df.index.str.contains("Heat Output")] = df_sum * 100

    to_drop = df[
        df.index.str.contains("Pyrolysis oil Output")
        | df.index.str.contains("Pyrolysis gas Output")
        | df.index.str.contains("Electricity Consumption")
        | df.index.str.contains("Feedstock Consumption")
    ].index
    df.drop(to_drop, inplace=True)

    # normalizing costs to biochar output
    df_divid = pd.concat(
        (
            df.iloc[df.index.str.contains("Biochar Output")],
            df.iloc[df.index.str.contains("Heat Output")],
        ),
        axis=0,
    ).sum(axis=0, skipna=False)
    biochar_totoutput = df.iloc[df.index.str.contains("Biochar Output")] / df_divid
    idx3 = df.index.str.contains("EUR")
    df.loc[idx3] = df.loc[idx3].values.astype(float) / biochar_totoutput.values.astype(
        float
    )
    df.index = df.index.str.replace(" output from pyrolysis process", "", regex=True)

    # rename units
    df.rename(
        index={
            df.loc[df.index.str.contains("Specific investment")].index[0]: df.loc[
                df.index.str.contains("Specific investment")
            ].index.str.replace("MW", "MW_biochar")[0],
            df.loc[df.index.str.contains("Fixed O&M")].index[0]: df.loc[
                df.index.str.contains("Fixed O&M")
            ].index.str.replace("MW", "MW_biochar")[0],
            df.loc[df.index.str.contains("Variable O&M")].index[0]: df.loc[
                df.index.str.contains("Variable O&M")
            ].index.str.replace("MWh", "MWh_biochar")[0],
        },
        inplace=True,
    )

    df_div = (
        df.iloc[df.index.str.contains("Specific energy content")].astype(float) / 3.6
    )
    df.iloc[df.index.str.contains("Specific energy content")] = df.iloc[
        df.index.str.contains("Biochar Output")
    ].astype(float) / df_div.values.astype(float)

    df.rename(
        index={
            df.loc[df.index.str.contains("Specific energy content")].index.values[
                0
            ]: "yield biochar [ton biochar/MWh_feedstock]",
            df.loc[df.index.str.contains("Biochar Output")].index.values[
                0
            ]: "efficiency biochar [MWh_biochar/MWh_feedstock]",
            df.loc[df.index.str.contains("Heat Output")].index.values[
                0
            ]: "efficiency heat [% MWh_feedstock]",
        },
        inplace=True,
    )

    return df


def get_data_from_DEA(
    years: list,
    sheet_names_dict: dict,
    input_data_dictionary: dict,
    offwind_no_grid_costs: bool = True,
    expectation: str = None,
) -> dict:
    """
    The function stores technology data from DEA in a dictionary.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    sheet_names_dict : dict
        dictionary having the technology name as keys and Excel sheet names as values
    input_data_dictionary : dict
        dictionary where the keys are the path to the DEA inputs and the values are the sheet names
    offwind_no_grid_costs : bool
        flag to remove grid connection costs from DEA for offwind. Such costs are calculated separately in pypsa-eur
    expectation : str
        tech data uncertainty. The possible options are [None, "optimist", "pessimist"]

    Returns
    -------
    Dictionary
        technology data from DEA
    """

    data_by_tech_dict = {}

    for tech_name, dea_tech in sheet_names_dict.items():
        logger.info(f"{tech_name} in PyPSA corresponds to {dea_tech} in DEA database.")
        df = get_data_DEA(
            years,
            tech_name,
            sheet_names_dict,
            input_data_dictionary,
            offwind_no_grid_costs,
            expectation,
        ).fillna(0)
        data_by_tech_dict[tech_name] = df

    return data_by_tech_dict


def clean_up_units(
    technology_dataframe: pd.DataFrame, value_column: str = "", source: str = ""
) -> pd.DataFrame:
    """
    The function converts units of an input dataframe. Namely, it converts:
        - power: Mega Watt (MW)
        - energy: Mega-Watt-hour (MWh)
        - currency: Euro (EUR)

    Parameters
    ----------
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions
    value_column : str
        column to modify
    source : str
        either empty string or 'dea'

    Returns
    -------
    pandas.DataFrame
        technology data with converted units
    """

    # Currency conversion
    REPLACEMENTS = [
        ("€", "EUR"),
        ("$", "USD"),
        ("₤", "GBP"),
    ]
    # Download the full history, this will be up-to-date. Current value is:
    # https://www.ecb.europa.eu/stats/eurofxref/eurofxref-hist.zip
    c = CurrencyConverter(ECB_URL, fallback_on_missing_rate=True)

    for old, new in REPLACEMENTS:
        technology_dataframe.unit = technology_dataframe.unit.str.replace(
            old, new, regex=False
        )
        technology_dataframe.loc[
            technology_dataframe.unit.str.contains(new), value_column
        ] *= c.convert(1, new, "EUR", date=date(2020, 1, 1))
        technology_dataframe.unit = technology_dataframe.unit.str.replace(new, "EUR")

    technology_dataframe.unit = technology_dataframe.unit.str.replace(" per ", "/")
    technology_dataframe.unit = technology_dataframe.unit.str.replace(" / ", "/")
    technology_dataframe.unit = technology_dataframe.unit.str.replace(" /", "/")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("J/s", "W")

    # units
    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("MEUR"), value_column
    ] *= 1e6
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MEUR", "EUR")

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("mio EUR"), value_column
    ] *= 1e6
    technology_dataframe.unit = technology_dataframe.unit.str.replace("mio EUR", "EUR")

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("mill. EUR"), value_column
    ] *= 1e6
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "mill. EUR", "EUR"
    )

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("1000EUR"), value_column
    ] *= 1e3
    technology_dataframe.unit = technology_dataframe.unit.str.replace("1000EUR", "EUR")

    technology_dataframe.unit = technology_dataframe.unit.str.replace("k EUR", "kEUR")
    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("kEUR"), value_column
    ] *= 1e3
    technology_dataframe.unit = technology_dataframe.unit.str.replace("kEUR", "EUR")

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("/kW"), value_column
    ] *= 1e3

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("kW")
        & ~technology_dataframe.unit.str.contains("/kW"),
        value_column,
    ] /= 1e3
    technology_dataframe.unit = technology_dataframe.unit.str.replace("kW", "MW")

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("/GWh"), value_column
    ] /= 1e3
    technology_dataframe.unit = technology_dataframe.unit.str.replace("/GWh", "/MWh")

    technology_dataframe.loc[
        technology_dataframe.unit.str.contains("/GJ"), value_column
    ] *= 3.6
    technology_dataframe.unit = technology_dataframe.unit.str.replace("/GJ", "/MWh")

    # Harmonise individual units so that they can be handled later
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        " a year", "/year"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace("2015EUR", "EUR")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("2015-EUR", "EUR")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("2020-EUR", "EUR")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("EUR2015", "EUR")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("EUR-2015", "EUR")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MWe", "MW_e")
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "EUR/MW of total input_e", "EUR/MW_e"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        r"MWh/MWh\)", "MWh_H2/MWh_e", regex=True
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MWth", "MW_th")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MWheat", "MW_th")
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MWhth", "MWh_th")
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWhheat", "MWh_th"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWH Liquids", "MWh_FT"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW Liquids", "MW_FT"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW Methanol", "MW_MeOH"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace("MW output", "MW")
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW/year FT Liquids/year", "MW_FT/year"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW/year Methanol", "MW_MeOH/year"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWh FT Liquids/year", "MWh_FT"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWh methanol", "MWh_MeOH"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW/year SNG", "MW_CH4/year"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWh SNG", "MWh_CH4"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW-methanol", "MW_MeOH"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW SNG", "MW_CH4"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "EUR/MWh of total input", "EUR/MWh_e"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "EUR/MWeh", "EUR/MWh_e"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "% -points of heat loss", "MWh_th/MWh_el"
    )
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "FT Liquids Output, MWh/MWh Total Input", "MWh_FT/MWh_H2"
    )
    # biomass-to-methanol-specific
    if isinstance(technology_dataframe.index, pd.MultiIndex):
        technology_dataframe.loc[
            technology_dataframe.index.get_level_values(1) == "Methanol Output,", "unit"
        ] = "MWh_MeOH/MWh_th"
        technology_dataframe.loc[
            technology_dataframe.index.get_level_values(1) == "District heat  Output,",
            "unit",
        ] = "MWh_th/MWh_th"
        technology_dataframe.loc[
            technology_dataframe.index.get_level_values(1) == "Electricity Output,",
            "unit",
        ] = "MWh_e/MWh_th"

    # Ammonia-specific
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW Ammonia output", "MW_NH3"
    )  # specific investment
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MW Ammonia", "MW_NH3"
    )  # fom
    technology_dataframe.unit = technology_dataframe.unit.str.replace(
        "MWh Ammonia", "MWh_NH3"
    )  # vom
    technology_dataframe.loc[technology_dataframe.unit == "EUR/MW/y", "unit"] = (
        "EUR/MW/year"
    )

    # convert per unit costs to MW
    cost_per_unit = technology_dataframe.unit.str.contains("/unit")
    technology_dataframe.loc[cost_per_unit, value_column] = technology_dataframe.loc[
        cost_per_unit, value_column
    ].apply(
        lambda val_x: (
            val_x
            / technology_dataframe.loc[
                (val_x.name[0], "Heat production capacity for one unit")
            ][value_column]
        ).iloc[0, :],
        axis=1,
    )
    technology_dataframe.loc[cost_per_unit, "unit"] = technology_dataframe.loc[
        cost_per_unit, "unit"
    ].str.replace("/unit", "/MW_th")

    if source == "dea":
        # clarify MW -> MW_th
        # see on p.278 of docu: "However, the primary purpose of the heat pumps in the
        # technology catalogue is heating. In this chapter the unit MW is referring to
        # the heat output (also MJ/s) unless otherwise noted"
        techs_mwth = [
            "central air-sourced heat pump",
            "central gas boiler",
            "central resistive heater",
            "decentral air-sourced heat pump",
            "decentral gas boiler",
            "decentral ground-sourced heat pump",
        ]
        technology_dataframe.loc[techs_mwth, "unit"] = technology_dataframe.loc[
            techs_mwth, "unit"
        ].replace(
            {
                "EUR/MW": "EUR/MW_th",
                "EUR/MW/year": "EUR/MW_th/year",
                "EUR/MWh": "EUR/MWh_th",
                "MW": "MW_th",
            }
        )

        # clarify MW -> MW_e
        techs_e = ["fuel cell"]
        technology_dataframe.loc[techs_e, "unit"] = technology_dataframe.loc[
            techs_e, "unit"
        ].replace(
            {
                "EUR/MW": "EUR/MW_e",
                "EUR/MW/year": "EUR/MW_e/year",
                "EUR/MWh": "EUR/MWh_e",
                "MW": "MW_e",
            }
        )

    technology_dataframe.unit = technology_dataframe.unit.str.replace(r"\)", "")
    return technology_dataframe


def set_specify_assumptions(
    years: list, technology_dataframe: pd.DataFrame
) -> pd.DataFrame:
    """
    The function implements more specific investment and efficiency assumptions for the following technologies:
        - central resistive heater (investment costs for large > 10 MW generators are assumed)
        - decentral gas boiler (grid connection costs)
        - biogas upgrading (include grid connection costs)
        - heat pumps (efficiencies for radiators assumed)

    Furthermore, to avoid duplicates some investment + efficiency, rows are dropped for:
        - decentral gas boilers (drop duplicated efficiency)
        - PV module (drop efficiency)

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology dataframe
    """

    # for central resistive heater there are investment costs for small (1-5MW)
    # and large (>10 MW) generators, assume the costs for large generators
    to_drop = [
        ("central resistive heater", "Nominal investment, 400/690 V; 1-5 MW"),
        ("decentral gas boiler", "Heat efficiency, annual average, net"),
    ]

    # for decentral gas boilers total and heat efficiency given, the values are
    # the same, drop one of the rows to avoid duplicates

    # for decentral gas boilers there are investment costs and possible
    # additional investments which apply for grid connection if the house is
    # not connected yet those costs are added as an extra row since the
    # lifetime of the branch pipe is assumed to be  50 years (see comment K in
    # Excel sheet)
    boiler_connect = technology_dataframe.loc[
        [
            ("decentral gas boiler", "Possible additional specific investment"),
            ("decentral gas boiler", "Technical lifetime"),
        ]
    ]
    boiler_connect.loc[("decentral gas boiler", "Technical lifetime"), years] = 50.0
    boiler_connect.rename(
        index={"decentral gas boiler": "decentral gas boiler connection"}, inplace=True
    )
    technology_dataframe = pd.concat([technology_dataframe, boiler_connect])
    to_drop.append(("decentral gas boiler", "Possible additional specific investment"))

    # biogas upgrading investment costs should include grid injection costs
    index = technology_dataframe.loc["biogas upgrading"].index.str.contains(
        "investment"
    )
    name = "investment (upgrading, methane redution and grid injection)"
    inv = (
        technology_dataframe.loc["biogas upgrading"]
        .loc[index]
        .groupby(["unit", "source"])
        .sum()
        .reset_index()
    )
    new = pd.concat(
        [technology_dataframe.loc["biogas upgrading"].loc[~index], inv]
    ).rename({0: name})
    new.index = pd.MultiIndex.from_product([["biogas upgrading"], new.index.to_list()])
    technology_dataframe.drop("biogas upgrading", level=0, inplace=True)
    technology_dataframe = pd.concat([technology_dataframe, new])

    # drop PV module conversion efficiency
    technology_dataframe = technology_dataframe.drop(
        "PV module conversion efficiency [p.u.]", level=1
    )

    # heat pump efficiencies are assumed the one's for existing building,
    # in the DEA they do differ between heating the floor area or heating with
    # radiators, since most households heat with radiators and there
    # efficiencies are lower (conservative approach) those are assumed
    # furthermore the total efficiency is assumed which includes auxiliary electricity
    # consumption
    name = "Heat efficiency, annual average, net, radiators"
    techs_radiator = technology_dataframe.xs(name, level=1).index
    for tech_name in techs_radiator:
        df = technology_dataframe.loc[tech_name]
        df = df[(~df.index.str.contains("efficiency")) | (df.index == name)]
        df.rename(index={name: name + ", existing one family house"}, inplace=True)
        df.index = pd.MultiIndex.from_product([[tech_name], df.index.to_list()])
        technology_dataframe.drop(tech_name, level=0, inplace=True)
        technology_dataframe = pd.concat([technology_dataframe, df])

    technology_dataframe = technology_dataframe.drop(to_drop)

    return technology_dataframe.sort_index()


def set_round_trip_efficiency(
    years: list, technology_dataframe: pd.DataFrame
) -> pd.DataFrame:
    """
    The function get round trip efficiency for hydrogen and battery storage.
    It assumes for battery sqrt(DC efficiency) and it splits it into inverter + storage.
    Finally, it renames investment rows for easier sorting.

    Parameters
    ----------
    years: list
        years for which a cost assumption is provided
    technology_dataframe: pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology dataframe
    """

    technology_dataframe.loc[
        ("hydrogen storage underground", "Round trip efficiency"), years
    ] *= 100.0
    technology_dataframe.loc[
        ("hydrogen storage tank type 1 including compressor", "Round trip efficiency"),
        years,
    ] *= 100.0

    # battery split into inverter and storage, assume for efficiency sqr(round trip DC)
    df = technology_dataframe.loc["battery"]
    inverter = df.loc[
        [
            "Round trip efficiency DC",
            "Output capacity expansion cost",
            "Technical lifetime",
            "Fixed O&M",
        ]
    ]

    inverter.rename(
        index={
            "Output capacity expansion cost": "Output capacity expansion cost investment"
        },
        inplace=True,
    )

    # Manual correction based on footnote.
    inverter.loc["Technical lifetime", years] = 10.0
    inverter.loc["Technical lifetime", "source"] += ", Note K."

    inverter.index = pd.MultiIndex.from_product(
        [["battery inverter"], inverter.index.to_list()]
    )

    storage = df.reindex(index=["Technical lifetime", "Energy storage expansion cost"])
    storage.rename(
        index={
            "Energy storage expansion cost": "Energy storage expansion cost investment"
        },
        inplace=True,
    )
    storage.index = pd.MultiIndex.from_product(
        [["battery storage"], storage.index.to_list()]
    )
    technology_dataframe.drop("battery", level=0, inplace=True)
    technology_dataframe = pd.concat([technology_dataframe, inverter, storage])

    return technology_dataframe.sort_index()


def order_data(years: list, technology_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function check if the units of different variables are conform and logs warnings if not.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        technology data in pypsa tech data syntax (investment, FOM,VOM, efficiency)
    """

    clean_df = {}
    for tech_name in technology_dataframe.index.get_level_values(0).unique():
        clean_df[tech_name] = pd.DataFrame()
        switch = False
        df = technology_dataframe.loc[tech_name]
        if tech_name == "methanolisation":
            print(df)
            input('')

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
                | (df.unit == "EUR/MWh Capacity")
                | (df.unit == "EUR/MWh")
                | (df.unit == "EUR/MW_CH4")
                | (df.unit == "EUR/MWh/year")
                | (df.unit == "EUR/MW_e, 2020")
                | (df.unit == "EUR/MW input")
                | (df.unit == "EUR/MW-methanol")
                | (df.unit == "EUR/t_N2/h")  # air separation unit
                | (df.unit == "EUR/MW_biochar")
            )
        ].copy()

        if len(investment) != 1:
            switch = True
            if df[df.index.str.contains("investment")].unit.empty:
                logger.info(f"check investment: {str(tech_name)} is not available")
            else:
                logger.info(
                    f"check investment: {str(tech_name)} {str(df[df.index.str.contains('investment')].unit)}"
                )
        else:
            investment["parameter"] = "investment"
            clean_df[tech_name] = investment

        # ---- FOM ----------------
        if len(investment):
            fixed = df[
                (
                    df.index.str.contains("Fixed O&M")
                    | df.index.str.contains("Total O&M")
                )
                & (
                    (df.unit == investment.unit.iloc[0] + "/year")
                    | (df.unit == "EUR/MW/km/year")
                    | (df.unit == "EUR/MW/year")
                    | (df.unit == "EUR/MW_e/y, 2020")
                    | (df.unit == "EUR/MW_e/y")
                    | (df.unit == "EUR/MW_FT/year")
                    | (df.unit == "EUR/MWh_FT")
                    | (df.unit == "EUR/MW_MeOH/year")
                    | (df.unit == "EUR/MW_CH4/year")
                    | (df.unit == "EUR/MW_biochar/year")
                    | (df.unit == "% of specific investment/year")
                    | (df.unit == investment.unit.str.split(" ").iloc[0][0] + "/year")
                )
            ].copy()

            if (len(fixed) != 1) and (len(df[df.index.str.contains("Fixed O&M")]) != 0):
                switch = True
                if df[df.index.str.contains("Fixed O&M")].unit.empty:
                    logger.info("check FOM: ", str(tech_name), " is not available")
                else:
                    logger.info(
                        f"check FOM: {str(tech_name)} {str(df[df.index.str.contains('Fixed O&M')].unit)}",
                    )
            if tech_name == "central water pit storage":
                # For current data, the FOM values for central water pit storage are too high by a factor of 1000.
                # See issue: https://github.com/PyPSA/technology-data/issues/203
                fixed[years] /= 1000  # in €/MWhCapacity/year
            if tech_name == "Fischer-Tropsch":
                fixed[years] *= 8000 # conversion from €/MWh to €/MW/year, assuming 8000 full load hours
            if len(fixed) == 1:
                fixed["parameter"] = "fixed"
                clean_df[tech_name] = pd.concat([clean_df[tech_name], fixed])
                fom = pd.DataFrame(columns=fixed.columns)
                if not any(fixed.unit.str.contains("% of specific investment/year")):
                    investment[investment == 0] = float("nan")
                    investment = investment.ffill(axis=1).fillna(0)
                    fom[years] = fixed[years] / investment[years].values * 100
                else:
                    fom[years] = fixed[years]
                fom["parameter"] = "FOM"
                fom["unit"] = "%/year"
                fom["source"] = fixed["source"]
                clean_df[tech_name] = pd.concat([clean_df[tech_name], fom])

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
                | (df.unit == "EUR/MWh_CH4")
                | (df.unit == "EUR/MWh_biochar")
                | (tech_name == "biogas upgrading")
            )
        ].copy()
        if len(vom) == 1:
            vom.loc[:, "parameter"] = "VOM"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], vom])

        elif len(vom) != 1 and len(df[df.index.str.contains("Variable O&M")]) != 0:
            switch = True
            if df[df.index.str.contains("Variable O&M")].unit.empty:
                logger.info(f"check VOM: {str(tech_name)} is not available")
            else:
                logger.info(
                    f"check VOM: {str(tech_name)} {str(df[df.index.str.contains('Variable O&M')].unit)}"
                )

        # ----- lifetime --------
        lifetime = df[
            df.index.str.contains("Technical life") & (df.unit == "years")
        ].copy()
        if len(lifetime) != 1:
            switch = True
            if df[df.index.str.contains("Technical life")].unit.empty:
                logger.info(f"check lifetime: {tech_name} is not available")
            else:
                logger.info(
                    f"check lifetime: {tech_name} {str(df[df.index.str.contains('Technical life')].unit)}"
                )
        else:
            lifetime["parameter"] = "lifetime"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], lifetime])

        if tech_name == "gas storage":
            lifetime_value = 100
            gas_storage_lifetime_index = "estimation: most underground storage are already built, they do have a long lifetime"
            gas_storage_lifetime_df = pd.DataFrame(index=[gas_storage_lifetime_index])

            for year in years:
                gas_storage_lifetime_df[year] = [lifetime_value]

            gas_storage_lifetime_df["parameter"] = "lifetime"
            gas_storage_lifetime_df["source"] = "TODO no source"
            gas_storage_lifetime_df["unit"] = "years"

            logger.info(f"Lifetime for {tech_name} manually set to {lifetime_value}")
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], gas_storage_lifetime_df]
            )

        # ----- efficiencies ------
        efficiency = df[
            (
                (df.index.str.contains("efficiency"))
                | (df.index.str.contains("Hydrogen output, at LHV"))
                | (df.index.str.contains("Hydrogen Output"))
                | (df.index.str.contains("FT Liquids Output"))
                | (df.index.str.contains("Methanol Output"))
                | (df.index.str.contains("District heat  Output"))
                | (df.index.str.contains("Electricity Output"))
                | (df.index.str.contains("hereof recoverable for district heating"))
                | (df.index.str.contains("Bio SNG"))
                | (df.index.str.contains("biochar"))
                | (df.index == ("Hydrogen"))
            )
            & (
                (df.unit == "%")
                | (df.unit == "% total size")
                | (df.unit == "% of fuel input")
                | (df.unit == "MWh_H2/MWh_e")
                | (df.unit == "%-points of heat loss")
                | (df.unit == "MWh_MeOH/MWh_th")
                | (df.unit == "MWh_e/MWh_th")
                | (df.unit == "MWh_th/MWh_th")
                | (df.unit == "MWh/MWh Total Input")
                | df.unit.str.contains("MWh_FT/MWh_H2")
                | df.unit.str.contains("MWh/MWh Total Input")
                | df.unit.str.contains("MWh_biochar/MWh_feedstock")
                | df.unit.str.contains("ton biochar/MWh_feedstock")
                | df.unit.str.contains("MWh_CH4/MWh_H2")
                | df.unit.str.contains("% MWh_feedstock")
            )
        ].copy()

        if tech_name == "Fischer-Tropsch":
            efficiency[years] *= 100

        # take annual average instead of name plate efficiency, unless central air-sourced heat pump
        if (
            any(efficiency.index.str.contains("annual average"))
            and tech_name != "central air-sourced heat pump"
        ):
            efficiency = efficiency[efficiency.index.str.contains("annual average")]
        elif any(efficiency.index.str.contains("name plate")):
            efficiency = efficiency[efficiency.index.str.contains("name plate")]

        # hydrogen electrolysiswith recoverable heat
        heat_recovery_label = "hereof recoverable for district heating"
        with_heat_recovery = efficiency.index.str.contains(heat_recovery_label)
        if with_heat_recovery.any():
            efficiency_heat = efficiency[with_heat_recovery].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_heat])
            efficiency_h2 = efficiency[
                efficiency.index.str.contains("Hydrogen Output")
            ].copy()
            efficiency_h2["parameter"] = "efficiency"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_h2])

        # check if electric and heat efficiencies are given
        if any(["Electric" in ind for ind in efficiency.index]) and any(
            ["Heat" in ind for ind in efficiency.index]
        ):
            efficiency_heat = efficiency[efficiency.index.str.contains("Heat")].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_heat])
            efficiency = efficiency[efficiency.index.str.contains("Electric")].copy()
            efficiency["parameter"] = "efficiency"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency])

        elif tech_name == "biomass-to-methanol":
            efficiency_heat = efficiency[
                efficiency.index.str.contains("District heat")
            ].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            efficiency_heat.loc[:, years] *= 100  # in %
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_heat])
            efficiency_elec = efficiency[
                efficiency.index.str.contains("Electric")
            ].copy()
            efficiency_elec["parameter"] = "efficiency-electricity"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_elec])
            efficiency_meoh = efficiency[
                efficiency.index.str.contains("Methanol")
            ].copy()
            efficiency_meoh["parameter"] = "efficiency"
            efficiency_meoh.loc[:, years] *= 100  # in %
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_meoh])

        elif tech_name == "biochar pyrolysis":
            efficiency_biochar = efficiency[
                efficiency.index.str.contains("efficiency biochar")
            ].copy()
            efficiency_biochar["parameter"] = "efficiency-biochar"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_biochar])
            efficiency_biochar_mass = efficiency[
                efficiency.index.str.contains("yield biochar")
            ].copy()
            efficiency_biochar_mass["parameter"] = "yield-biochar"
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], efficiency_biochar_mass]
            )
            efficiency_heat = efficiency[
                efficiency.index.str.contains("efficiency heat")
            ].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency_heat])

        elif len(efficiency) != 1:
            switch = True
            if not any(efficiency.index.str.contains("Round trip")):
                if df[df.index.str.contains("efficiency")].unit.empty:
                    logger.info(f"check efficiency: {str(tech_name)} is not available")
                else:
                    logger.info(
                        f"check efficiency: {str(tech_name)} {df[df.index.str.contains('efficiency')].unit}"
                    )
        else:
            efficiency["parameter"] = "efficiency"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], efficiency])

        # add storage temperature for TES
        if tech_name == "central water pit storage":
            top_storage_temp_ptes = df.loc[
                df.index.str.contains("Max. storage temperature, hot")
            ].copy()
            top_storage_temp_ptes["parameter"] = "Top storage temperature"
            top_storage_temp_ptes.rename(
                index={
                    "Max. storage temperature, hot": "Typical max. storage temperature"
                },
                inplace=True,
            )
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], top_storage_temp_ptes]
            )

            bottom_storage_temp_ptes = df.loc[
                df.index.str.contains("Storage temperature, discharged")
            ].copy()
            bottom_storage_temp_ptes["parameter"] = "Bottom storage temperature"
            bottom_storage_temp_ptes.rename(
                index={
                    "Storage temperature, discharged": "Typical bottom storage temperature"
                },
                inplace=True,
            )
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], bottom_storage_temp_ptes]
            )
            energy_loss = df.loc[
                df.index.str.contains("Energy losses during storage")
            ].copy()
            energy_loss["parameter"] = "standing losses"
            energy_loss.loc[("Energy losses during storage", years)] = (
                energy_loss.loc[("Energy losses during storage", years)]
                / (
                    78
                    - bottom_storage_temp_ptes.loc[
                        ("Typical bottom storage temperature", years)
                    ]
                )
                * 100
                / 24
            )  # 78°C is the average temperature for ptes
            energy_loss["unit"] = "%/hour"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], energy_loss])

        if tech_name == "central water tank storage":
            temp_difference_central_ttes = df.loc[
                df.index.str.contains("Typical temperature difference in storage")
            ].copy()
            temp_difference_central_ttes["parameter"] = "temperature difference"
            temp_difference_central_ttes.rename(
                index={
                    "Typical temperature difference in storage": "Typical temperature difference"
                },
                inplace=True,
            )
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], temp_difference_central_ttes]
            )
            energy_loss = df.loc[
                df.index.str.contains("Energy losses during storage")
            ].copy()
            energy_loss["parameter"] = "standing losses"
            energy_loss[years] = energy_loss[years] / 24
            energy_loss["unit"] = "%/hour"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], energy_loss])

        if tech_name == "decentral water tank storage":
            temp_difference_decentral_ttes = df.loc[
                df.index.str.contains("Typical temperature difference in storage")
            ].copy()
            temp_difference_decentral_ttes["parameter"] = "temperature difference"
            temp_difference_decentral_ttes.rename(
                index={
                    "Typical temperature difference in storage": "Typical temperature difference"
                },
                inplace=True,
            )
            clean_df[tech_name] = pd.concat(
                [clean_df[tech_name], temp_difference_decentral_ttes]
            )
            energy_loss = df.loc[
                df.index.str.contains("Energy losses during storage")
            ].copy()
            energy_loss["parameter"] = "standing losses"
            energy_loss[years] = energy_loss[years]
            energy_loss["unit"] = "%/hour"
            clean_df[tech_name] = pd.concat([clean_df[tech_name], energy_loss])

        # add c_v and c_b coefficient
        if "Cb coefficient" in df.index:
            c_b = df.loc[df.index.str.contains("Cb coefficient")].dropna().copy()
            if len(c_b):
                c_b["parameter"] = "c_b"
                clean_df[tech_name] = pd.concat([clean_df[tech_name], c_b])
        if "Cv coefficient" in df.index:
            c_v = df.loc[df.index.str.contains("Cv coefficient")].dropna().copy()
            if len(c_v):
                c_v["parameter"] = "c_v"
                clean_df[tech_name] = pd.concat([clean_df[tech_name], c_v])

        if switch:
            logger.info("---------------------------------------")

    # concat data
    output_data_dataframe = (
        pd.concat(clean_df)
        .reset_index()
        .rename(columns={"level_0": "technology", "level_1": "further description"})
        .set_index(["technology", "parameter"])
    )

    # add central water tank charger/ discharger
    charger_tank = technology_dataframe.loc[
        ("central water tank storage", " - Charge efficiency")
    ].copy()
    charger_tank["further description"] = "Charger efficiency"
    charger_tank.rename(
        index={" - Charge efficiency": "efficiency"}, level=1, inplace=True
    )
    charger_tank.rename(
        index={"central water tank storage": "central water tank charger"},
        level=0,
        inplace=True,
    )
    output_data_dataframe = pd.concat([output_data_dataframe, charger_tank], sort=True)
    charger_tank.rename(
        index={"central water tank charger": "central water tank discharger"},
        level=0,
        inplace=True,
    )
    charger_tank["further description"] = "Discharger efficiency"

    output_data_dataframe = pd.concat([output_data_dataframe, charger_tank], sort=True)

    # add decentral water tank charger/ discharger
    charger_tank = technology_dataframe.loc[
        ("decentral water tank storage", " - Charge efficiency")
    ].copy()
    charger_tank["further description"] = "Charger efficiency"
    charger_tank.rename(
        index={" - Charge efficiency": "efficiency"}, level=1, inplace=True
    )
    charger_tank.rename(
        index={"decentral water tank storage": "decentral water tank charger"},
        level=0,
        inplace=True,
    )
    output_data_dataframe = pd.concat([output_data_dataframe, charger_tank], sort=True)
    charger_tank.rename(
        index={"decentral water tank charger": "decentral water tank discharger"},
        level=0,
        inplace=True,
    )
    charger_tank["further description"] = "Discharger efficiency"

    output_data_dataframe = pd.concat([output_data_dataframe, charger_tank], sort=True)

    # add water pit charger/ discharger
    charger_pit = technology_dataframe.loc[
        ("central water pit storage", " - Charge efficiency")
    ].copy()
    charger_pit[years] *= 100
    charger_pit["further description"] = "Charger efficiency"

    charger_pit.rename(
        index={" - Charge efficiency": "efficiency"}, level=1, inplace=True
    )
    charger_pit.rename(
        index={"central water pit storage": "central water pit charger"},
        level=0,
        inplace=True,
    )
    output_data_dataframe = pd.concat([output_data_dataframe, charger_pit], sort=True)
    charger_pit.rename(
        index={"central water pit charger": "central water pit discharger"},
        level=0,
        inplace=True,
    )
    charger_pit["further description"] = "Discharger efficiency"
    output_data_dataframe = pd.concat([output_data_dataframe, charger_pit], sort=True)

    # add energy to power ratio for central water tank storage
    power_ratio_tank = (
        technology_dataframe.loc[
            ("central water tank storage", "Input capacity for one unit")
        ]
        .copy()
        .squeeze()
    )
    storage_capacity_tank = (
        technology_dataframe.loc[
            ("central water tank storage", "Energy storage capacity for one unit")
        ]
        .copy()
        .squeeze()
    )

    power_ratio_tank[years] = storage_capacity_tank[years].div(power_ratio_tank[years])
    power_ratio_tank["further description"] = (
        "Ratio between energy storage and input capacity"
    )
    power_ratio_tank["unit"] = "h"
    power_ratio_tank = power_ratio_tank.to_frame().T
    power_ratio_tank.rename(
        index={"Input capacity for one unit": "energy to power ratio"},
        level=1,
        inplace=True,
    )
    output_data_dataframe = pd.concat(
        [output_data_dataframe, power_ratio_tank], sort=True
    )

    # add energy to power ratio for decentral water tank storage
    power_ratio_tank = (
        technology_dataframe.loc[
            ("decentral water tank storage", "Input capacity for one unit")
        ]
        .copy()
        .squeeze()
    )
    storage_capacity_tank = (
        technology_dataframe.loc[
            ("decentral water tank storage", "Energy storage capacity for one unit")
        ]
        .copy()
        .squeeze()
    )

    power_ratio_tank[years] = storage_capacity_tank[years].div(power_ratio_tank[years])
    power_ratio_tank["further description"] = (
        "Ratio between energy storage and input capacity"
    )
    power_ratio_tank["unit"] = "h"
    power_ratio_tank = power_ratio_tank.to_frame().T
    power_ratio_tank.rename(
        index={"Input capacity for one unit": "energy to power ratio"},
        level=1,
        inplace=True,
    )
    output_data_dataframe = pd.concat(
        [output_data_dataframe, power_ratio_tank], sort=True
    )

    # add energy to power ratio for water pit storage
    power_ratio_pit = (
        technology_dataframe.loc[
            ("central water pit storage", "Input capacity for one unit")
        ]
        .copy()
        .squeeze()
    )
    storage_capacity_pit = (
        technology_dataframe.loc[
            ("central water pit storage", "Energy storage capacity for one unit")
        ]
        .copy()
        .squeeze()
    )

    power_ratio_pit[years] = storage_capacity_pit[years].div(power_ratio_pit[years])
    power_ratio_pit["further description"] = (
        "Ratio between energy storage and input capacity"
    )
    power_ratio_pit["unit"] = "h"
    power_ratio_pit = power_ratio_pit.to_frame().T
    power_ratio_pit.rename(
        index={"Input capacity for one unit": "energy to power ratio"},
        level=1,
        inplace=True,
    )
    output_data_dataframe = pd.concat(
        [output_data_dataframe, power_ratio_pit], sort=True
    )

    # add gas storage charger/ discharger
    # process equipment, injection (2200MW) withdrawal (6600MW)
    # assuming half of investment costs for injection, half for withdrawal
    investment_gas_storage_charger = technology_dataframe.loc[
        ("gas storage", "Total investment cost")
    ].copy()
    investment_gas_storage_charger[years] = (
        investment_gas_storage_charger[years] / 2 / 2200 / 1e3
    )
    investment_gas_storage_charger.loc[
        ("gas storage", "Total investment cost"), "unit"
    ] = "EUR/kW"

    investment_gas_storage_charger.rename(
        index={"Total investment cost": "investment"}, level=1, inplace=True
    )
    investment_gas_storage_charger.rename(
        index={"gas storage": "gas storage charger"},
        level=0,
        inplace=True,
    )
    output_data_dataframe = pd.concat(
        [output_data_dataframe, investment_gas_storage_charger], sort=True
    )

    investment_gas_storage_discharger = technology_dataframe.loc[
        ("gas storage", "Total investment cost")
    ].copy()
    investment_gas_storage_discharger[years] = (
        investment_gas_storage_discharger[years] / 2 / 6600 / 1e3
    )
    investment_gas_storage_discharger.loc[
        ("gas storage", "Total investment cost"), "unit"
    ] = "EUR/kW"

    investment_gas_storage_discharger.rename(
        index={"Total investment cost": "investment"}, level=1, inplace=True
    )
    investment_gas_storage_discharger.rename(
        index={"gas storage": "gas storage discharger"},
        level=0,
        inplace=True,
    )
    output_data_dataframe = pd.concat(
        [output_data_dataframe, investment_gas_storage_discharger], sort=True
    )

    return output_data_dataframe


def add_description(
    years: list,
    technology_dataframe: pd.DataFrame,
    offwind_no_grid_costs_flag: bool = True,
) -> pd.DataFrame:
    """
    The function adds the Excel sheet name as a column to the tech data and adds comments for offwind connection costs.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions
    offwind_no_grid_costs_flag : bool
        flag to remove grid connection costs from DEA for offwind. Such costs are calculated separately in pypsa-eur

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # add Excel sheet names to technology_dataframe frame
    wished_order = years + ["unit", "source", "further description"]
    technology_dataframe = technology_dataframe.reindex(columns=wished_order)
    technology_dataframe.index.set_names(["technology", "parameter"], inplace=True)
    sheets = (
        technology_dataframe.reset_index()["technology"].map(dea_sheet_names).fillna("")
    )
    sheets.index = technology_dataframe.index
    technology_dataframe["further description"] = (
        sheets + ":  " + technology_dataframe["further description"]
    )

    # add comment for offwind investment
    if offwind_no_grid_costs_flag:
        technology_dataframe.loc[("offwind", "investment"), "further description"] += (
            " grid connection costs subtracted from investment costs"
        )

    return technology_dataframe


def convert_units(years: list, technology_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function converts investment and efficiency units to be aligned with old pypsa assumptions.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe : pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # convert efficiency from % -> per unit
    technology_dataframe.loc[
        technology_dataframe.index.get_level_values(1).isin(
            ["efficiency", "efficiency-heat"]
        ),
        years,
    ] /= 100
    technology_dataframe.loc[
        technology_dataframe.index.get_level_values(1).isin(
            ["efficiency", "efficiency-heat"]
        ),
        "unit",
    ] = "per unit"

    # convert MW -> kW
    to_convert = technology_dataframe.index.get_level_values(1).isin(
        ["fixed", "investment"]
    ) & technology_dataframe.unit.str.contains("/MW")
    technology_dataframe.loc[to_convert, years] /= 1e3
    technology_dataframe.loc[to_convert, "unit"] = technology_dataframe.loc[
        to_convert, "unit"
    ].str.replace("/MW", "/kW")

    return technology_dataframe


def add_carbon_capture(
    years: list,
    sheet_names_dict: dict,
    new_technology_dataframe: pd.DataFrame,
    technology_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """
    The function adds carbon capture rates.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    sheet_names_dict : dict
        dictionary having the technology name as keys and Excel sheet names as values
    new_technology_dataframe:
        updated technology data cost assumptions
    technology_dataframe : pandas.DataFrame
        existing technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    for tech_name in ["cement capture", "biomass CHP capture"]:
        new_technology_dataframe.loc[(tech_name, "capture_rate"), years] = (
            technology_dataframe.loc[
                (tech_name, "Ax) CO2 capture rate, net"), years
            ].values[0]
            / 100
        )
        new_technology_dataframe.loc[(tech_name, "capture_rate"), "unit"] = "per unit"

    for tech_name in ["direct air capture", "cement capture", "biomass CHP capture"]:
        new_technology_dataframe.loc[(tech_name, "investment"), years] = (
            technology_dataframe.loc[(tech_name, "Specific investment"), years].values[
                0
            ]
            * 1e6
        )
        new_technology_dataframe.loc[(tech_name, "investment"), "unit"] = "EUR/(tCO2/h)"

        new_technology_dataframe.loc[(tech_name, "FOM"), years] = (
            technology_dataframe.loc[(tech_name, "Fixed O&M"), years].values[0]
            / technology_dataframe.loc[
                (tech_name, "Specific investment"), years
            ].values[0]
            * 100
        )
        new_technology_dataframe.loc[(tech_name, "FOM"), "unit"] = "%/year"

        name_list = [
            ("C2) Eletricity input ", "electricity-input"),
            ("C1) Heat  input ", "heat-input"),
            ("C1) Heat out ", "heat-output"),
            (
                "CO₂ compression and dehydration - Electricity input",
                "compression-electricity-input",
            ),
            ("CO₂ compression and dehydration - Heat out", "compression-heat-output"),
        ]

        for dea_name, our_name in name_list:
            new_technology_dataframe.loc[(tech_name, our_name), years] = (
                technology_dataframe.loc[(tech_name, dea_name), years].values[0]
            )
            new_technology_dataframe.loc[(tech_name, our_name), "unit"] = "MWh/tCO2"

        new_technology_dataframe.loc[tech_name, "source"] = (
            new_technology_dataframe.loc[(tech_name, "lifetime"), "source"]
        )
        new_technology_dataframe.loc[tech_name, "further description"] = (
            sheet_names_dict[tech_name]
        )

    return new_technology_dataframe


def rename_pypsa_old(cost_dataframe_pypsa: pd.DataFrame) -> pd.DataFrame:
    """
    The function renames old technology names to new ones to compare converts units from water tanks to compare.

    Parameters
    ----------
    cost_dataframe_pypsa: pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    to_drop = ["retrofitting I", "retrofitting II"]
    cost_dataframe_pypsa.drop(to_drop, level=0, inplace=True)

    # rename to new names
    cost_dataframe_pypsa.rename({"central CHP": "central gas CHP"}, inplace=True)
    cost_dataframe_pypsa.rename(
        {"hydrogen underground storage": "hydrogen storage underground"}, inplace=True
    )

    # convert EUR/m^3 to EUR/kWh for 40 K diff and 1.17 kWh/m^3/K
    cost_dataframe_pypsa.loc[
        ("decentral water tank storage", "investment"), "value"
    ] /= 1.17 * 40
    cost_dataframe_pypsa.loc[("decentral water tank storage", "investment"), "unit"] = (
        "EUR/kWh"
    )

    return cost_dataframe_pypsa


def add_manual_input(technology_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds input from manual_input.csv.

    Parameters
    ----------
    technology_dataframe: pandas.DataFrame
        technology data cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    df = pd.read_csv(
        snakemake.input["manual_input"], quotechar='"', sep=",", keep_default_na=False
    )
    df = df.rename(columns={"further_description": "further description"})

    content_list = []
    for tech_name in df["technology"].unique():
        c0 = df[df["technology"] == tech_name]
        for param in c0["parameter"].unique():
            queried_df = df.query("technology == @tech_name and parameter == @param")

            row_series = pd.Series(
                index=snakemake.config["years"],
                data=np.interp(
                    snakemake.config["years"], queried_df["year"], queried_df["value"]
                ),
                name=param,
            )
            row_series["parameter"] = param
            row_series["technology"] = tech_name
            try:
                row_series["currency_year"] = int(queried_df["currency_year"].values[0])
            except ValueError:
                row_series["currency_year"] = np.nan
            for col in ["unit", "source", "further description"]:
                row_series[col] = "; and ".join(queried_df[col].unique().astype(str))
            row_series = row_series.rename(
                {"further_description": "further description"}
            )  # match column name between manual_input and original TD workflow
            content_list.append(row_series)

    new_df = pd.DataFrame(content_list).set_index(["technology", "parameter"])
    technology_dataframe.index.set_names(["technology", "parameter"], inplace=True)
    # overwrite DEA data with manual input
    technology_dataframe = new_df.combine_first(technology_dataframe)

    return technology_dataframe


def rename_ISE(cost_dataframe_ise: pd.DataFrame) -> pd.DataFrame:
    """
    The function renames ISE costs to fit to tech data.

    Parameters
    ----------
    cost_dataframe_ise: pandas.DataFrame
        ISE cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    cost_dataframe_ise.rename(
        index={
            "Investition": "investment",
            "Lebensdauer": "lifetime",
            "M/O-Kosten": "FOM",
        },
        columns={
            "Einheit": "unit",
            "2020": 2020,
            "2025": 2025,
            "2030": 2030,
            "2035": 2035,
            "2040": 2040,
            "2045": 2045,
            "2050": 2050,
        },
        inplace=True,
    )
    cost_dataframe_ise.index.names = ["technology", "parameter"]
    cost_dataframe_ise["unit"] = cost_dataframe_ise.unit.replace(
        {"a": "years", "% Invest": "%"}
    )
    cost_dataframe_ise["source"] = source_dict["ISE"]
    cost_dataframe_ise["further description"] = cost_dataframe_ise.reset_index()[
        "technology"
    ].values
    # could not find specific currency year in report, assume year of publication
    cost_dataframe_ise["currency_year"] = 2020

    return cost_dataframe_ise


def rename_ISE_vehicles(costs_vehicles_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function renames ISE vehicles costs to fit to tech data.
    energy

    Parameters
    ----------
    costs_vehicles_dataframe: pandas.DataFrame
        vehicles ISE cost assumptions

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    costs_vehicles_dataframe.rename(
        index={
            "Investition": "investment",
            "Lebensdauer": "lifetime",
            "M/O-Kosten": "FOM",
            "Wirkungsgrad*": "efficiency",
            "PKW Batterie-Elektromotor": "Battery electric (passenger cars)",
            "LKW Batterie-Elektromotor": "Battery electric (trucks)",
            "LKW H2- Brennstoffzelle": "Hydrogen fuel cell (trucks)",
            "PKW H2- Brennstoffzelle": "Hydrogen fuel cell (passenger cars)",
            "LKW ICE- Flï¿½ssigtreibstoff": "Liquid fuels ICE (trucks)",
            "PKW ICE- Flï¿½ssigtreibstoff": "Liquid fuels ICE (passenger cars)",
            "LKW Ladeinfrastruktur Brennstoffzellen Fahrzeuge * LKW": "Charging infrastructure fuel cell vehicles trucks",
            "PKW Ladeinfrastruktur Brennstoffzellen Fahrzeuge * PKW": "Charging infrastructure fuel cell vehicles passenger cars",
            "PKW Ladeinfrastruktur schnell  (reine) Batteriefahrzeuge*": "Charging infrastructure fast (purely) battery electric vehicles passenger cars",
            "Ladeinfrastruktur langsam (reine) Batteriefahrzeuge*": "Charging infrastructure slow (purely) battery electric vehicles passenger cars",
        },
        columns={
            "Einheit": "unit",
            "2020": 2020,
            "2025": 2025,
            "2030": 2030,
            "2035": 2035,
            "2040": 2040,
            "2045": 2045,
            "2050": 2050,
        },
        inplace=True,
    )
    costs_vehicles_dataframe.index.names = ["technology", "parameter"]
    costs_vehicles_dataframe["unit"] = costs_vehicles_dataframe.unit.replace(
        {"a": "years", "% Invest": "%"}
    )
    costs_vehicles_dataframe["source"] = source_dict["vehicles"]
    # could not find specific currency year in report, assume year of publication
    costs_vehicles_dataframe["currency_year"] = 2020
    costs_vehicles_dataframe["further description"] = (
        costs_vehicles_dataframe.reset_index()["technology"].values
    )
    return costs_vehicles_dataframe


def carbon_flow(
    years: list, cost_dataframe: pd.DataFrame, year_to_use: int
) -> pd.DataFrame:
    """
    The function renames ISE vehicles costs to fit to tech data.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    cost_dataframe: pandas.DataFrame
        cost dataframe
    year_to_use: int
        year to use

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # NB: This requires some digits of accuracy; rounding to two digits creates carbon inbalances when scaling up
    c_in_char = 0  # Carbon ending up in char: zero avoids inbalace -> assumed to be circulated back and eventually end up in one of the other output streams
    medium_out = ""
    CH4_specific_energy = 50  # GJ/t methane

    btlcost_data = np.interp(x=years, xp=[2020, 2050], fp=[3500, 2000])
    btl_cost = pd.Series(data=btlcost_data, index=years)

    bmH2cost_data = np.interp(x=years, xp=[2020, 2050], fp=[4000, 2500])
    bmH2_cost = pd.Series(data=bmH2cost_data, index=years)

    btleta_data = np.interp(x=years, xp=[2020, 2050], fp=[0.35, 0.45])
    btl_eta = pd.Series(data=btleta_data, index=years)

    # Adding pelletizing cost to biomass boiler
    cost_dataframe.loc[("biomass boiler", "pelletizing cost"), "value"] = 9
    cost_dataframe.loc[("biomass boiler", "pelletizing cost"), "unit"] = (
        "EUR/MWh_pellets"
    )
    cost_dataframe.loc[("biomass boiler", "pelletizing cost"), "currency_year"] = 2019
    cost_dataframe.loc[("biomass boiler", "pelletizing cost"), "source"] = (
        "Assumption based on doi:10.1016/j.rser.2019.109506"
    )

    for tech_name in [
        "Fischer-Tropsch",
        "methanolisation",
        "BtL",
        "biomass-to-methanol",
        "BioSNG",
        "biogas",
        "biogas CC",
        "digestible biomass to hydrogen",
        "solid biomass to hydrogen",
        "electrobiofuels",
    ]:
        inv_cost = 0
        eta = 0
        lifetime = 0
        FOM = 0
        VOM = 0
        currency_year = np.nan
        source = "TODO"
        co2_capture_rate = 0.90

        if (tech_name, "capture rate") not in cost_dataframe.index:
            cost_dataframe.loc[(tech_name, "capture rate"), "value"] = co2_capture_rate
            cost_dataframe.loc[(tech_name, "capture rate"), "unit"] = "per unit"
            cost_dataframe.loc[(tech_name, "capture rate"), "source"] = (
                "Assumption based on doi:10.1016/j.biombioe.2015.01.006"
            )

        if tech_name == "BtL":
            inv_cost = btl_cost[year_to_use]
            medium_out = "oil"
            eta = btl_eta[year_to_use]
            source = "doi:10.1016/j.enpol.2017.05.013"
            currency_year = 2017

        if tech_name == "biomass-to-methanol":
            medium_out = "methanol"

        elif tech_name == "BioSNG":
            medium_out = "gas"
            lifetime = 25

        elif tech_name in ["biogas", "biogas CC"]:
            eta = 1
            source = "Assuming input biomass is already given in biogas output"
            AD_CO2_share = 0.4  # volumetric share in biogas (rest is CH4)

        elif tech_name == "biogas plus hydrogen":
            # NB: this falls between power to gas and biogas and should be used with care, due to possible minor
            # differences in resource use etc. which may tweak results in favour of one tech or another
            eta = 1.6
            H2_in = 0.46

            heat_out = 0.19
            source = "Calculated from data in Danish Energy Agency, data_sheets_for_renewable_fuels.xlsx"
            cost_dataframe.loc[(tech_name, "hydrogen input"), "value"] = H2_in
            cost_dataframe.loc[(tech_name, "hydrogen input"), "unit"] = "MWh_H2/MWh_CH4"
            cost_dataframe.loc[(tech_name, "hydrogen input"), "source"] = source

            cost_dataframe.loc[(tech_name, "heat output"), "value"] = heat_out
            cost_dataframe.loc[(tech_name, "heat output"), "unit"] = "MWh_th/MWh_CH4"
            cost_dataframe.loc[(tech_name, "heat output"), "source"] = source
            currency_year = cost_dataframe.loc[
                ("biogas plus hydrogen", "VOM"), "currency_year"
            ]

            # TODO: this needs to be refined based on e.g. stoichiometry:
            AD_CO2_share = 0.1  # volumetric share in biogas (rest is CH4).

        elif tech_name == "digestible biomass to hydrogen":
            inv_cost = bmH2_cost[year_to_use]
            eta = 0.39
            FOM = 4.25
            currency_year = 2014
            costs.loc[(tech_name, "FOM"), "currency_year"] = 2014
            source = "Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014"  # source_dict('HyNOW')

        elif tech_name == "solid biomass to hydrogen":
            inv_cost = bmH2_cost[year_to_use]
            eta = 0.56
            FOM = 4.25
            currency_year = 2014
            cost_dataframe.loc[(tech_name, "FOM"), "currency_year"] = 2014
            source = "Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014"  # source_dict('HyNOW')

        if eta > 0:
            cost_dataframe.loc[(tech_name, "efficiency"), "value"] = eta
            cost_dataframe.loc[(tech_name, "efficiency"), "unit"] = "per unit"
            cost_dataframe.loc[(tech_name, "efficiency"), "source"] = source

        if tech_name in ["BioSNG", "BtL", "biomass-to-methanol"]:
            input_CO2_intensity = cost_dataframe.loc[
                ("solid biomass", "CO2 intensity"), "value"
            ]

            cost_dataframe.loc[(tech_name, "C in fuel"), "value"] = (
                cost_dataframe.loc[(tech_name, "efficiency"), "value"]
                * cost_dataframe.loc[(medium_out, "CO2 intensity"), "value"]
                / input_CO2_intensity
            )
            cost_dataframe.loc[(tech_name, "C stored"), "value"] = (
                1 - cost_dataframe.loc[(tech_name, "C in fuel"), "value"] - c_in_char
            )
            cost_dataframe.loc[(tech_name, "CO2 stored"), "value"] = (
                input_CO2_intensity
                * cost_dataframe.loc[(tech_name, "C stored"), "value"]
            )

            cost_dataframe.loc[(tech_name, "C in fuel"), "unit"] = "per unit"
            cost_dataframe.loc[(tech_name, "C stored"), "unit"] = "per unit"
            cost_dataframe.loc[(tech_name, "CO2 stored"), "unit"] = "tCO2/MWh_th"

            cost_dataframe.loc[(tech_name, "C in fuel"), "source"] = (
                "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            )
            cost_dataframe.loc[(tech_name, "C stored"), "source"] = (
                "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            )
            cost_dataframe.loc[(tech_name, "CO2 stored"), "source"] = (
                "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            )

        elif tech_name in ["electrobiofuels"]:
            input_CO2_intensity = cost_dataframe.loc[
                ("solid biomass", "CO2 intensity"), "value"
            ]
            oil_CO2_intensity = cost_dataframe.loc[("oil", "CO2 intensity"), "value"]

            cost_dataframe.loc[("electrobiofuels", "C in fuel"), "value"] = (
                cost_dataframe.loc[("BtL", "C in fuel"), "value"]
                + cost_dataframe.loc[("BtL", "C stored"), "value"]
                * cost_dataframe.loc[("Fischer-Tropsch", "capture rate"), "value"]
            )
            cost_dataframe.loc[("electrobiofuels", "C in fuel"), "unit"] = "per unit"
            cost_dataframe.loc[("electrobiofuels", "C in fuel"), "source"] = (
                "Stoichiometric calculation"
            )

            cost_dataframe.loc[("electrobiofuels", "efficiency-biomass"), "value"] = (
                cost_dataframe.loc[("electrobiofuels", "C in fuel"), "value"]
                * input_CO2_intensity
                / oil_CO2_intensity
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-biomass"), "unit"] = (
                "per unit"
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-biomass"), "source"] = (
                "Stoichiometric calculation"
            )

            efuel_scale_factor = (
                cost_dataframe.loc[("BtL", "C stored"), "value"]
                * cost_dataframe.loc[("Fischer-Tropsch", "capture rate"), "value"]
            )

            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "value"] = (
                cost_dataframe.loc[("Fischer-Tropsch", "efficiency"), "value"]
                / efuel_scale_factor
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "unit"] = (
                "per unit"
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "source"] = (
                "Stoichiometric calculation"
            )

            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "value"] = 1 / (
                1
                / cost_dataframe.loc[
                    ("electrobiofuels", "efficiency-hydrogen"), "value"
                ]
                + 1
                / cost_dataframe.loc[("electrobiofuels", "efficiency-biomass"), "value"]
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "unit"] = (
                "per unit"
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "source"] = (
                "Stoichiometric calculation"
            )

            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "value"] = (
                cost_dataframe.loc[("Fischer-Tropsch", "efficiency"), "value"]
                / efuel_scale_factor
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "unit"] = (
                "per unit"
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-hydrogen"), "source"] = (
                "Stoichiometric calculation"
            )

            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "value"] = 1 / (
                1
                / cost_dataframe.loc[
                    ("electrobiofuels", "efficiency-hydrogen"), "value"
                ]
                + 1
                / cost_dataframe.loc[("electrobiofuels", "efficiency-biomass"), "value"]
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "unit"] = (
                "per unit"
            )
            cost_dataframe.loc[("electrobiofuels", "efficiency-tot"), "source"] = (
                "Stoichiometric calculation"
            )

            inv_cost = (
                btl_cost[year_to_use]
                + cost_dataframe.loc[("Fischer-Tropsch", "investment"), "value"]
                * efuel_scale_factor
            )
            VOM = (
                cost_dataframe.loc[("BtL", "VOM"), "value"]
                + cost_dataframe.loc[("Fischer-Tropsch", "VOM"), "value"]
                * efuel_scale_factor
            )
            FOM = cost_dataframe.loc[("BtL", "FOM"), "value"]
            medium_out = "oil"
            currency_year = cost_dataframe.loc[
                ("Fischer-Tropsch", "investment"), "currency_year"
            ]
            cost_dataframe.loc[(tech_name, "FOM"), "currency_year"] = 2015
            source = "combination of BtL and electrofuels"

        elif tech_name in ["biogas", "biogas CC", "biogas plus hydrogen"]:
            CH4_density = 0.657  # kg/Nm3
            CO2_density = 1.98  # kg/Nm3
            CH4_vol_energy_density = (
                (1 - AD_CO2_share) * CH4_specific_energy * CH4_density / (1000 * 3.6)
            )  # MJ/Nm3 -> MWh/Nm3
            CO2_weight_share = (
                AD_CO2_share * CO2_density
            )  # TODO: what value is used for AD_CO2_share in this if branch?

            cost_dataframe.loc[(tech_name, "CO2 stored"), "value"] = (
                CO2_weight_share / CH4_vol_energy_density / 1000
            )  # tCO2/MWh,in (NB: assuming the input is already given in the biogas potential and cost
            cost_dataframe.loc[(tech_name, "CO2 stored"), "unit"] = "tCO2/MWh_th"
            cost_dataframe.loc[(tech_name, "CO2 stored"), "source"] = (
                "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            )

        if inv_cost > 0:
            cost_dataframe.loc[(tech_name, "investment"), "value"] = inv_cost
            cost_dataframe.loc[(tech_name, "investment"), "unit"] = "EUR/kW_th"
            cost_dataframe.loc[(tech_name, "investment"), "source"] = source
            cost_dataframe.loc[(tech_name, "investment"), "currency_year"] = (
                currency_year
            )

        if lifetime > 0:
            cost_dataframe.loc[(tech_name, "lifetime"), "value"] = lifetime
            cost_dataframe.loc[(tech_name, "lifetime"), "unit"] = "years"
            cost_dataframe.loc[(tech_name, "lifetime"), "source"] = source

        if FOM > 0:
            cost_dataframe.loc[(tech_name, "FOM"), "value"] = FOM
            cost_dataframe.loc[(tech_name, "FOM"), "unit"] = "%/year"
            cost_dataframe.loc[(tech_name, "FOM"), "source"] = source

        if VOM > 0:
            cost_dataframe.loc[(tech_name, "VOM"), "value"] = VOM
            cost_dataframe.loc[(tech_name, "VOM"), "unit"] = "EUR/MWh_th"
            cost_dataframe.loc[(tech_name, "VOM"), "source"] = source
            cost_dataframe.loc[(tech_name, "VOM"), "currency_year"] = currency_year

    return cost_dataframe


def energy_penalty(cost_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds energy penalty for biomass carbon capture.

    Parameters
    ----------
    cost_dataframe: pandas.DataFrame
        cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # Need to take steam production for CC into account, assumed with the main feedstock,
    # e.g. the input biomass is used also for steam, and the efficiency for el and heat is scaled down accordingly

    for tech_name in [
        "central solid biomass CHP CC",
        "waste CHP CC",
        "solid biomass boiler steam CC",
        "direct firing solid fuels CC",
        "direct firing gas CC",
        "biogas CC",
    ]:
        if "powerboost" in tech_name:
            boiler = "electric boiler steam"
            feedstock = "solid biomass"
            co2_capture = cost_dataframe.loc[(feedstock, "CO2 intensity"), "value"]
        elif "gas" in tech_name:
            boiler = "gas boiler steam"
            feedstock = "gas"
            co2_capture = cost_dataframe.loc[(feedstock, "CO2 intensity"), "value"]
        elif "biogas" in tech_name:
            boiler = "gas boiler steam"
            co2_capture = cost_dataframe.loc[(tech_name, "CO2 stored"), "value"]
        else:
            boiler = "solid biomass boiler steam"
            feedstock = "solid biomass"
            co2_capture = cost_dataframe.loc[(feedstock, "CO2 intensity"), "value"]

        # Scaling biomass input to account for heat demand of carbon capture
        scalingFactor = 1 / (
            1
            + co2_capture
            * cost_dataframe.loc[("biomass CHP capture", "heat-input"), "value"]
            / cost_dataframe.loc[(boiler, "efficiency"), "value"]
        )

        eta_steam = (1 - scalingFactor) * cost_dataframe.loc[
            (boiler, "efficiency"), "value"
        ]
        eta_old = cost_dataframe.loc[(tech_name, "efficiency"), "value"]

        eta_main = (
            cost_dataframe.loc[(tech_name, "efficiency"), "value"] * scalingFactor
        )

        # Adapting investment share of tech due to steam boiler addition. Investment per MW_el.
        cost_dataframe.loc[(tech_name, "investment"), "value"] = (
            cost_dataframe.loc[(tech_name, "investment"), "value"] * eta_old / eta_main
            + cost_dataframe.loc[(boiler, "investment"), "value"] * eta_steam / eta_main
        )
        cost_dataframe.loc[(tech_name, "investment"), "source"] = (
            "Combination of " + tech_name + " and " + boiler
        )
        cost_dataframe.loc[(tech_name, "investment"), "further description"] = ""

        if cost_dataframe.loc[(tech_name, "VOM"), "value"]:
            break
        else:
            cost_dataframe.loc[(tech_name, "VOM"), "value"] = 0.0

        cost_dataframe.loc[(tech_name, "VOM"), "value"] = (
            cost_dataframe.loc[(tech_name, "VOM"), "value"] * eta_old / eta_main
            + cost_dataframe.loc[(boiler, "VOM"), "value"] * eta_steam / eta_main
        )
        cost_dataframe.loc[(tech_name, "VOM"), "source"] = (
            "Combination of " + tech_name + " and " + boiler
        )
        cost_dataframe.loc[(tech_name, "VOM"), "further description"] = ""

        cost_dataframe.loc[(tech_name, "efficiency"), "value"] = eta_main
        cost_dataframe.loc[(tech_name, "efficiency"), "source"] = (
            "Combination of " + tech_name + " and " + boiler
        )
        cost_dataframe.loc[(tech_name, "efficiency"), "further description"] = ""

        if "CHP" in tech_name:
            cost_dataframe.loc[(tech_name, "efficiency-heat"), "value"] = (
                cost_dataframe.loc[(tech_name, "efficiency-heat"), "value"]
                * scalingFactor
                + cost_dataframe.loc[("solid biomass", "CO2 intensity"), "value"]
                * (
                    cost_dataframe.loc[("biomass CHP capture", "heat-output"), "value"]
                    + cost_dataframe.loc[
                        ("biomass CHP capture", "compression-heat-output"), "value"
                    ]
                )
            )
            cost_dataframe.loc[(tech_name, "efficiency-heat"), "source"] = (
                "Combination of " + tech_name + " and " + boiler
            )
            cost_dataframe.loc[
                (tech_name, "efficiency-heat"), "further description"
            ] = ""

        if "biogas CC" in tech_name:
            cost_dataframe.loc[(tech_name, "VOM"), "value"] = 0
            cost_dataframe.loc[(tech_name, "VOM"), "unit"] = "EUR/MWh"

        cost_dataframe.loc[(tech_name, "VOM"), "value"] = (
            cost_dataframe.loc[(tech_name, "VOM"), "value"] * eta_old / eta_main
            + cost_dataframe.loc[(boiler, "VOM"), "value"] * eta_steam / eta_main
        )
        cost_dataframe.loc[(tech_name, "VOM"), "source"] = (
            "Combination of " + tech_name + " and " + boiler
        )
        cost_dataframe.loc[(tech_name, "VOM"), "further description"] = ""

    return cost_dataframe


def add_egs_data(technology_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds enhanced geothermal systems cost assumptions.
    Data taken from Aghahosseini, Breyer 2020: From hot rock to useful energy...

    Parameters
    ----------
    technology_dataframe: pandas.DataFrame
        technology data

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    parameters = [
        "CO2 intensity",
        "lifetime",
        "efficiency residential heat",
        "efficiency electricity",
        "FOM",
    ]
    tech_name_list = ["geothermal"]
    multi_i = pd.MultiIndex.from_product(
        [tech_name_list, parameters], names=["technology", "parameter"]
    )
    geoth_df = pd.DataFrame(index=multi_i, columns=data.columns)
    years = [col for col in data.columns if isinstance(col, int)]

    # lifetime
    geoth_df.loc[("geothermal", "lifetime"), years] = 30  # years
    geoth_df.loc[("geothermal", "lifetime"), "unit"] = "years"
    geoth_df.loc[("geothermal", "lifetime"), "source"] = source_dict["Aghahosseini2020"]

    # co2 emissions
    geoth_df.loc[("geothermal", "CO2 intensity"), years] = 0.12  # tCO2/MWh_el
    geoth_df.loc[("geothermal", "CO2 intensity"), "unit"] = "tCO2/MWh_el"
    geoth_df.loc[("geothermal", "CO2 intensity"), "source"] = source_dict[
        "Aghahosseini2020"
    ]
    geoth_df.loc[("geothermal", "CO2 intensity"), "further description"] = (
        "Likely to be improved; Average of 85 percent of global egs power plant capacity"
    )

    # efficiency for heat generation using organic rankine cycle
    geoth_df.loc[("geothermal", "efficiency residential heat"), years] = 0.8
    geoth_df.loc[("geothermal", "efficiency residential heat"), "unit"] = "per unit"
    geoth_df.loc[("geothermal", "efficiency residential heat"), "source"] = (
        "{}; {}".format(source_dict["Aghahosseini2020"], source_dict["Breede2015"])
    )
    geoth_df.loc[
        ("geothermal", "efficiency residential heat"), "further description"
    ] = "This is a rough estimate, depends on local conditions"

    # efficiency for electricity generation using organic rankine cycle
    geoth_df.loc[("geothermal", "efficiency electricity"), years] = 0.1
    geoth_df.loc[("geothermal", "efficiency electricity"), "unit"] = "per unit"
    geoth_df.loc[("geothermal", "efficiency electricity"), "source"] = "{}; {}".format(
        source_dict["Aghahosseini2020"], source_dict["Breede2015"]
    )
    geoth_df.loc[("geothermal", "efficiency electricity"), "further description"] = (
        "This is a rough estimate, depends on local conditions"
    )

    # relative additional capital cost of using residual heat for district heating (25 percent)
    geoth_df.loc[("geothermal", "district heating cost"), years] = 0.25
    geoth_df.loc[("geothermal", "district heating cost"), "unit"] = "%"
    geoth_df.loc[("geothermal", "district heating cost"), "source"] = "{}".format(
        source_dict["Frey2022"]
    )
    geoth_df.loc[("geothermal", "district heating cost"), "further description"] = (
        "If capital cost of electric generation from EGS is 100%, district heating adds additional 25%"
    )

    # fixed operational costs
    geoth_df.loc[("geothermal", "FOM"), years] = 2.0
    geoth_df.loc[("geothermal", "FOM"), "unit"] = "%/year"
    geoth_df.loc[("geothermal", "FOM"), "source"] = source_dict["Aghahosseini2020"]
    geoth_df.loc[("geothermal", "FOM"), "further description"] = (
        "Both for flash, binary and ORC plants. See Supplemental Material for details"
    )

    geoth_df = geoth_df.dropna(axis=1, how="all")

    return pd.concat([technology_dataframe, geoth_df])


def annuity(n: float, r: float = 0.07) -> float:
    """
    The function calculates the annuity factor for an asset with lifetime n years and discount rate of r.

    Parameters
    ----------
    n: float:
        lifetime
    r: float
        discount rate

    Returns
    -------
    float
        annuity
    """

    if r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def add_home_battery_costs(
    ewg_cost_file_name: str, years: list, cost_dataframe: pd.DataFrame
) -> pd.DataFrame:
    """
    The function adds investment costs for home battery storage and inverter.
    Since home battery costs are not part of the DEA catalogue, utility-scale
    costs are multiplied by a factor determined by data from the EWG study.

    Parameters
    ----------
    ewg_cost_file_name: str
        file name for the cost assumptions from the EWG study
    years : list
        years for which a cost assumption is provided
    cost_dataframe: pandas.DataFrame
        existing cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # get DEA assumptions for utility scale
    home_battery = data.loc[["battery storage", "battery inverter"]].rename(
        index=lambda x: "home " + x, level=0
    )

    # get EWG cost assumptions
    costs_ewg = pd.read_csv(ewg_cost_file_name, index_col=list(range(2))).sort_index()
    v = costs_ewg.unstack()[[str(year) for year in years]].swaplevel(axis=1)

    # annualise EWG cost assumptions
    fixed = (annuity(v["lifetime"]) + v["FOM"] / 100.0) * v["investment"]

    # battery storage index in EWG --------------
    battery_store_i = [
        "Battery PV prosumer - commercial storage",
        "Battery PV prosumer - industrial storage",
        "Battery PV prosumer - residential storage",
        "Battery storage",
    ]

    battery_store_ewg = fixed.loc[battery_store_i].T

    def get_factor(df, cols, utility_col):
        """Get factor by which costs are increasing for home installations"""
        return (
            df[cols]
            .div(df[utility_col], axis=0)
            .mean(axis=1)
            .rename(index=lambda x: float(x))
        )

    # take mean of cost increase for commercial and residential storage compared to utility-scale
    home_cols = [
        "Battery PV prosumer - commercial storage",
        "Battery PV prosumer - residential storage",
    ]
    factor = get_factor(battery_store_ewg, home_cols, "Battery storage")

    home_cost = (
        home_battery.loc[("home battery storage", "investment"), years] * factor
    ).values
    home_battery.loc[("home battery storage", "investment"), years] = home_cost

    # battery inverter index in EWG -----------------------
    battery_inverter_i = [
        "Battery PV prosumer - commercial interface",
        "Battery PV prosumer - industrial interface PHES",
        "Battery PV prosumer - residential interface",
        "Battery interface",
    ]

    battery_inverter_ewg = fixed.loc[battery_inverter_i].T

    home_cols = [
        "Battery PV prosumer - commercial interface",
        "Battery PV prosumer - residential interface",
    ]
    factor = get_factor(battery_inverter_ewg, home_cols, "Battery interface")
    home_cost = (
        home_battery.loc[("home battery inverter", "investment"), years] * factor
    ).values
    home_battery.loc[("home battery inverter", "investment"), years] = home_cost

    # adjust source
    home_battery["source"] = home_battery["source"].apply(
        lambda x: source_dict["EWG"] + ", " + x
    )

    return pd.concat([cost_dataframe, home_battery])


def add_SMR_data(years: list, technology_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    The function adds steam methane reforming (SMR) technology data. The data set used are:
    - investment cost : currently no cost reduction for investment costs of SMR CC assumed.
        - IEA (2020) [1]: assumes cost reduction -19.2% from 2019-2050
        - Agora [2]: no cost reduction
    - carbon capture rate:
        - IEA (2020) [1]: 0.9
        - Agora [2]: 0.9
        - [3]: 0.9
        - Timmerberg et al.: 0.56-0.9
    - efficiency:
        - Agora SMR + CC (LHV/LHV) 0.58

    [1] IEA (2020) https://www.iea.org/data-and-statistics/charts/global-average-levelised-cost-of-hydrogen-production-by-energy-source-and-technology-2019-and-2050
    [2] Agora (2021) p.52 https://static.agora-energiewende.de/fileadmin/Projekte/2021/2021_02_EU_H2Grid/A-EW_203_No-regret-hydrogen_WEB.pdf
    [3] p.12 https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1011506/Hydrogen_Production_Costs_2021.pdf

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe: pandas.DataFrame
        technology cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    parameters = ["FOM", "investment", "lifetime", "efficiency"]
    techs = ["SMR", "SMR CC"]
    multi_i = pd.MultiIndex.from_product(
        [techs, parameters], names=["technology", "parameter"]
    )
    SMR_df = pd.DataFrame(index=multi_i, columns=technology_dataframe.columns)

    # efficiencies per unit in LHV (stays constant 2019 to 2050)
    SMR_df.loc[("SMR", "efficiency"), years] = 0.76
    SMR_df.loc[("SMR CC", "efficiency"), years] = 0.69
    SMR_df.loc[(techs, "efficiency"), "source"] = source_dict["IEA"]
    SMR_df.loc[(techs, "efficiency"), "unit"] = "per unit (in LHV)"

    # lifetime
    SMR_df.loc[(techs, "lifetime"), years] = 30
    SMR_df.loc[(techs, "lifetime"), "source"] = source_dict["IEA"]
    SMR_df.loc[(techs, "lifetime"), "unit"] = "years"

    # FOM
    SMR_df.loc[(techs, "FOM"), years] = 5
    SMR_df.loc[(techs, "FOM"), "source"] = source_dict["DEA"]
    SMR_df.loc[(techs, "FOM"), "unit"] = "%/year"
    SMR_df.loc[(techs, "FOM"), "currency_year"] = 2015
    SMR_df.loc[(techs, "FOM"), "further description"] = (
        "Technology data for renewable fuels, in pdf on table 3 p.311"
    )

    # investment
    # investment given in unit EUR/kg H_2/h -> convert to EUR/MW_CH4
    # lower heating value (LHV) of H2
    LHV_H2 = 33.33  # unit kWh/kg
    SMR = 12500 / LHV_H2 * 1e3 * 1 / SMR_df.loc[("SMR", "efficiency"), years]
    SMR_CCS = 14500 / LHV_H2 * 1e3 * 1 / SMR_df.loc[("SMR", "efficiency"), years]

    SMR_df.loc[("SMR", "investment"), years] = SMR
    SMR_df.loc[("SMR CC", "investment"), years] = SMR_CCS
    SMR_df.loc[(techs, "investment"), "source"] = source_dict["DEA"]
    SMR_df.loc[(techs, "investment"), "unit"] = "EUR/MW_CH4"
    SMR_df.loc[(techs, "investment"), "currency_year"] = 2015
    SMR_df.loc[(techs, "investment"), "further description"] = (
        "Technology data for renewable fuels, in pdf on table 3 p.311"
    )

    # carbon capture rate
    SMR_df.loc[("SMR CC", "capture_rate"), years] = 0.9
    SMR_df.loc[("SMR CC", "capture_rate"), "source"] = source_dict["IEA"]
    SMR_df.loc[("SMR CC", "capture_rate"), "unit"] = "per unit"
    SMR_df.loc[("SMR CC", "capture_rate"), "further description"] = (
        "wide range: capture rates between 54%-90%"
    )

    SMR_df = SMR_df.dropna(axis=1, how="all")

    return pd.concat([technology_dataframe, SMR_df])


def add_mean_solar_rooftop(
    years: list, technology_dataframe: pd.DataFrame
) -> pd.DataFrame:
    """
    The function adds costs for solar rooftop.

    Parameters
    ----------
    years : list
        years for which a cost assumption is provided
    technology_dataframe: pandas.DataFrame
        technology cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated technology data
    """

    # take mean of rooftop commercial and residential
    rooftop = (
        technology_dataframe.loc[
            technology_dataframe.index.get_level_values(0).str.contains("solar-rooftop")
        ][years]
        .astype(float)
        .groupby(level=1)
        .mean()
    )
    for col in technology_dataframe.columns[~technology_dataframe.columns.isin(years)]:
        rooftop[col] = technology_dataframe.loc["solar-rooftop residential"][col]
    # set multi index
    rooftop = pd.concat([rooftop], keys=["solar-rooftop"])
    rooftop["source"] = "Calculated. See 'further description'."
    rooftop["further description"] = (
        "Mixed investment costs based on average of 50% 'solar-rooftop commercial' and 50% 'solar-rooftop residential'"
    )
    # add to technology_dataframe
    rooftop.index.names = technology_dataframe.index.names
    technology_dataframe = pd.concat([technology_dataframe, rooftop])
    # add solar assuming 50% utility and 50% rooftop
    solar = (
        (technology_dataframe.loc[["solar-rooftop", "solar-utility"]][years])
        .astype(float)
        .groupby(level=1)
        .mean()
    )
    for col in technology_dataframe.columns[~technology_dataframe.columns.isin(years)]:
        solar[col] = technology_dataframe.loc["solar-rooftop residential"][col]
    solar["source"] = "Calculated. See 'further description'."
    solar["further description"] = (
        "Mixed investment costs based on average of 50% 'solar-rooftop' and 50% 'solar-utility'"
    )
    # set multi index
    solar = pd.concat([solar], keys=["solar"])
    solar.index.names = technology_dataframe.index.names
    return pd.concat([technology_dataframe, solar])


def geometric_series(
    nominator: float, denominator: float = 1.0, number_of_terms: int = 1, start: int = 1
) -> float:
    """
    The function computes a geometric series. The geometric series is given with a constant ratio between successive terms.
    When moving to infinity the geometric series converges to a limit. https://en.wikipedia.org/wiki/Series_(mathematics)
    For example, for nominator = 1.0, denominator = 2.0, number_of_terms = 3 and start = 0 results in
    1/2**0 + 1/2**1 + 1/2**2 = 1 + 1/2 + 1/4 = 1.75. If number_of_terms grows, the sum converges to 2

    Parameters
    ----------
    nominator : float
        nominator of the ratio
    denominator : float
        denominator of the ratio
    number_of_terms : int
        number of terms in the sum
    start : int
        if the value is 0, it means it starts at the first term

    Returns
    -------
    float
        sum of the terms
    """
    return sum(
        [nominator / denominator**i for i in range(start, start + number_of_terms)]
    )


def add_energy_storage_database(
    pnnl_storage_file_name: str,
    pnnl_energy_storage_dict: dict,
    cost_dataframe: pd.DataFrame,
    data_year: int,
) -> (pd.DataFrame, pd.Series):
    """
    The function adds energy storage database compiled.

    Learning rate drop. For example, the nominal DC SB learning rate for RFBs is set at
    4.5%, 1.5% for lead-acid batteries, compared to 10% for Li-ion batteries, corresponding to cost drops of
    17%, 6%, and 35%, respectively. For the rest of the categories for battery-based systems, the learning
    rates were kept the same for all batteries as described in the ESGC 2020 report.

    Fix cost drop. Due to the uncertainties in both anticipated deployments and the correct learning rate to use during the
    initial phase, this work assumes a fixed-cost drop for zinc batteries, gravity, and thermal storage
    systems. For example, a 20% cost drop in DC SB and 10% drop in DCBOS was assumed for zinc batteries,
    while keeping the cost drops for power equipment in line with Li-ion BESS, while system integration,
    EPC, and project development costs are maintained at 90% of Li-ion BESS 2030 values.

    Parameters
    ----------
    pnnl_storage_file_name: str
        PNNL storage file name
    pnnl_energy_storage_dict: dict
        PNNL storage configuration dictionary
    cost_dataframe: pandas.DataFrame
        existing cost dataframe
    data_year: int
        year to consider

    Returns
    -------
    tuple with DataFrame and Series
        updated cost dataframe and technologies
    """

    logger.info(f"Add energy storage database compiled for year {data_year}")
    # a) Import csv file
    df = pd.read_excel(
        pnnl_storage_file_name,
        sheet_name="energy-storage-database",
        dtype={
            "technology": str,
            "type": str,
            "carrier": str,
            "parameter": str,
            "year": int,
            "value": float,
            "unit": str,
            "source": str,
            "note": str,
            "reference": str,
            "ref_size_MW": float,
            "EP_ratio_h": float,
        },
        engine="calamine",
    )
    df = df.drop(columns=["ref_size_MW", "EP_ratio_h"])
    df = df.fillna(df.dtypes.replace({"float64": 0.0, "O": "NULL"}))
    df.loc[:, "unit"] = df.unit.str.replace("NULL", "per unit")

    # b) Change data to PyPSA format (aggregation of components, units, currency, etc.)
    df = clean_up_units(df, "value")  # base clean up

    # rewrite technology to be charger, store, discharger, bidirectional-charger
    df.loc[:, "carrier"] = df.carrier.str.replace("NULL", "")
    df.loc[:, "carrier"] = df["carrier"].apply(lambda x: x.split("-"))
    carrier_list_len = df["carrier"].apply(lambda x: len(x))
    carrier_str_len = df["carrier"].apply(lambda x: len(x[0]))
    carrier_first_item = df["carrier"].apply(lambda x: x[0])
    carrier_last_item = df["carrier"].apply(lambda x: x[-1])
    bicharger_filter = carrier_list_len == 3
    charger_filter = (carrier_list_len == 2) & (carrier_first_item == "elec")
    discharger_filter = (carrier_list_len == 2) & (carrier_last_item == "elec")
    store_filter = (carrier_list_len == 1) & (carrier_str_len > 0)
    reference_filter = (carrier_list_len == 1) & (
        carrier_first_item == "reference_value"
    )
    df = df[~reference_filter]  # remove reference values
    df.loc[bicharger_filter, "technology_type"] = "bicharger"
    df.loc[charger_filter, "technology_type"] = "charger"
    df.loc[discharger_filter, "technology_type"] = "discharger"
    df.loc[store_filter, "technology_type"] = "store"
    df.loc[df.unit == "EUR/MWh-year", "technology_type"] = "store"
    # Some investment inputs need to be distributed between charger and discharger
    for tech_name in df.technology.unique():
        nan_filter = (
            (df.technology == tech_name)
            & (carrier_str_len == 0)
            & (df.parameter == "investment")
        )
        store_filter = nan_filter & (df.unit == "EUR/MWh")
        if not df.loc[store_filter].empty:
            df.loc[store_filter, "technology_type"] = (
                "store"  # value will be aggregated later in the groupby
            )
        # charger and discharger with 50% distribution e.g. in case of Hydrogen
        power_filter = nan_filter & (df.unit == "EUR/MW")
        if not df.loc[power_filter].empty:
            agg = (
                df.loc[power_filter]
                .groupby(["technology", "year"])
                .sum(numeric_only=True)
            )
            charger_investment_filter = (
                charger_filter
                & (df.technology == tech_name)
                & (df.parameter == "investment")
            )
            discharger_investment_filter = (
                discharger_filter
                & (df.technology == tech_name)
                & (df.parameter == "investment")
            )
            df.loc[charger_investment_filter & df.year == 2021, "value"] += (
                agg.loc[(tech_name, 2021)] / 2
            )
            df.loc[charger_investment_filter & df.year == 2030, "value"] += (
                agg.loc[(tech_name, 2030)] / 2
            )
            df.loc[discharger_investment_filter & df.year == 2021, "value"] += (
                agg.loc[(tech_name, 2021)] / 2
            )
            df.loc[discharger_investment_filter & df.year == 2030, "value"] += (
                agg.loc[(tech_name, 2030)] / 2
            )
    df.loc[:, "technology"] = df["technology"] + "-" + df["technology_type"]

    # aggregate technology_type and unit
    df = (
        df.groupby(["technology", "unit", "year"])
        .agg(
            {
                "technology": "first",
                "year": "first",
                "parameter": "first",
                "value": "sum",
                "unit": "first",
                "type": "first",
                "carrier": "first",
                "technology_type": "first",
                "source": "first",
                "note": "first",
                "reference": "first",
            }
        )
        .reset_index(drop=True)
    )

    # calculate %/year FOM on aggregated values
    for tech_name in df.technology.unique():
        for year in df.year.unique():
            df_tech = df.loc[(df.technology == tech_name) & (df.year == year)].copy()
            a = df_tech.loc[df_tech.unit == "EUR/MW-year", "value"].values
            b = df_tech.loc[df_tech.unit == "EUR/MW", "value"].values
            df.loc[df_tech.loc[df_tech.unit == "EUR/MW-year"].index, "value"] = (
                a / b * 100
            )  # EUR/MW-year / EUR/MW = %/year
            c = df_tech.loc[df_tech.unit == "EUR/MWh-year", "value"].values
            d = df_tech.loc[df_tech.unit == "EUR/MWh", "value"].values
            df.loc[df_tech.loc[df_tech.unit == "EUR/MWh-year"].index, "value"] = (
                c / d * 100
            )  # EUR/MWh-year / EUR/MWh = %/year

    df.loc[:, "unit"] = df.unit.str.replace("EUR/MW-year", "%/year")
    df.loc[:, "unit"] = df.unit.str.replace("EUR/MWh-year", "%/year")

    # c) Linear Inter/Extrapolation
    # data available for 2021 and 2030, but value for given "year" passed by function needs to be calculated
    for tech_name in df.technology.unique():
        for param in df.parameter.unique():
            filter = (df.technology == tech_name) & (df.parameter == param)
            y = df.loc[filter, "value"]
            if y.empty:
                continue  # nothing to interpolate
            elif y.iloc[0] == y.iloc[1] or param == "efficiency" or param == "lifetime":
                ynew = y.iloc[1]  # assume new value is the same as 2030
            elif y.iloc[0] != y.iloc[1]:
                x = df.loc[filter, "year"]  # both values 2021+2030
                first_segment_diff = y.iloc[0] - y.iloc[1]
                endp_first_segment = y.iloc[1]

                # Below we create linear segments between 2021-2030
                # While the first segment is known, the others are defined by the initial segments with a accumulating quadratic decreasing gradient
                other_segments_points = [2034, 2039, 2044, 2049, 2054, 2059]

                if (
                    tech_name == "Hydrogen-discharger"
                    or tech_name == "Pumped-Heat-store"
                ):
                    x1 = pd.concat(
                        [x, pd.DataFrame(other_segments_points)], ignore_index=True
                    )
                    y1 = y
                    factor = 5
                    for i in range(
                        len(other_segments_points)
                    ):  # -1 because of segments
                        cost_at_year = endp_first_segment - geometric_series(
                            nominator=first_segment_diff,
                            denominator=factor,
                            number_of_terms=i + 1,
                        )
                        y1 = pd.concat(
                            [y1, pd.DataFrame([cost_at_year])], ignore_index=True
                        )
                    f = interpolate.interp1d(
                        x1.squeeze(),
                        y1.squeeze(),
                        kind="linear",
                        fill_value="extrapolate",
                    )
                elif tech_name == "Hydrogen-charger":
                    x2 = pd.concat(
                        [x, pd.DataFrame(other_segments_points)], ignore_index=True
                    )
                    y2 = y
                    factor = 6.5
                    for i in range(len(other_segments_points)):
                        cost_at_year = endp_first_segment - geometric_series(
                            nominator=first_segment_diff,
                            denominator=factor,
                            number_of_terms=i + 1,
                        )
                        y2 = pd.concat(
                            [y2, pd.DataFrame([cost_at_year])], ignore_index=True
                        )
                    f = interpolate.interp1d(
                        x2.squeeze(),
                        y2.squeeze(),
                        kind="linear",
                        fill_value="extrapolate",
                    )
                else:
                    x3 = pd.concat(
                        [x, pd.DataFrame(other_segments_points)], ignore_index=True
                    )
                    y3 = y
                    factor = 2
                    for i in range(len(other_segments_points)):
                        cost_at_year = endp_first_segment - geometric_series(
                            nominator=first_segment_diff,
                            denominator=factor,
                            number_of_terms=i + 1,
                        )
                        y3 = pd.concat(
                            [y3, pd.DataFrame([cost_at_year])], ignore_index=True
                        )
                    f = interpolate.interp1d(
                        x3.squeeze(),
                        y3.squeeze(),
                        kind="linear",
                        fill_value="extrapolate",
                    )

                option = pnnl_energy_storage_dict
                if option.get("approx_beyond_2030") == ["geometric_series"]:
                    ynew = f(data_year)
                if option.get("approx_beyond_2030") == ["same_as_2030"]:
                    if data_year <= 2030:
                        # apply linear interpolation
                        ynew = f(data_year)
                    if data_year > 2030:
                        # apply same value as 2030
                        ynew = y.iloc[1]  # assume new value is the same as 2030

            df_new = pd.DataFrame(
                [
                    {
                        "technology": tech_name,
                        "year": data_year,
                        "parameter": param,
                        "value": ynew,
                        "unit": df.loc[filter, "unit"].unique().item(),
                        "source": df.loc[filter, "source"].unique().item(),
                        "carrier": df.loc[filter, "carrier"].iloc[1],
                        "technology_type": df.loc[filter, "technology_type"]
                        .unique()
                        .item(),
                        "type": df.loc[filter, "type"].unique().item(),
                        "note": df.loc[filter, "note"].iloc[1],
                        "reference": df.loc[filter, "reference"].iloc[1],
                    }
                ]
            )
            # not concat if df year is 2021 or 2030 (otherwise duplicate)
            if data_year == 2021 or data_year == 2030:
                continue
            else:
                df = pd.concat([df, df_new], ignore_index=True)

    # d) Combine metadata and add to cost database
    df.loc[:, "source"] = df["source"] + ", " + df["reference"]
    for i in df.index:
        df.loc[i, "further description"] = str(
            {
                "carrier": df.loc[i, "carrier"],
                "technology_type": [df.loc[i, "technology_type"]],
                "type": [df.loc[i, "type"]],
                "note": [df.loc[i, "note"]],
            }
        )
    # keep only relevant columns
    df = df.loc[
        df.year == data_year,
        ["technology", "parameter", "value", "unit", "source", "further description"],
    ]
    tech_names = df.technology.unique()
    df = df.set_index(["technology", "parameter"])

    return pd.concat([cost_dataframe, df]), tech_names


# %% *************************************************************************
#  ---------- MAIN ------------------------------------------------------------
if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions")

    configure_logging(snakemake)

    years_list = list(snakemake.config["years"])
    inflation_rate = prepare_inflation_rate(snakemake.input.inflation_rate)

    # p.77 Figure 51 share of vehicle-km driven by truck

    # (1) DEA data
    # (a)-------- get data from DEA excel sheets ----------------------------------
    # read excel sheet names of all excel files
    excel_files = [v for k, v in snakemake.input.items() if "dea" in k.casefold()]
    data_in = get_excel_sheets(excel_files)
    # create dictionary with raw data from DEA sheets
    d_by_tech = get_data_from_DEA(
        years_list,
        dea_sheet_names,
        data_in,
        snakemake.config["offwind_no_gridcosts"],
        expectation=snakemake.config["expectation"],
    )
    # concat into pd.Dataframe
    tech_data = pd.concat(d_by_tech).sort_index()
    # clean up units
    tech_data = clean_up_units(tech_data, years_list, source="dea")

    # (b) ------ specific assumptions for some technologies -----------------------

    # specify investment and efficiency assumptions for:
    # resistive heater, decentral gas boiler, biogas upgrading and heat pumps
    tech_data = set_specify_assumptions(years_list, tech_data)

    # round trip efficiency for hydrogen + battery storage
    tech_data = set_round_trip_efficiency(years_list, tech_data)

    # drop all rows which only contains zeros
    tech_data = tech_data.loc[(tech_data[years_list] != 0).sum(axis=1) != 0]

    # (c) -----  get tech data in pypsa syntax -----------------------------------
    # make categories: investment, FOM, VOM, efficiency, c_b, c_v
    data = order_data(years_list, tech_data)
    # add Excel sheet names and further description
    data = add_description(years_list, data, snakemake.config["offwind_no_gridcosts"])
    # convert efficiency from %-> per unit and investment from MW->kW to compare
    data = convert_units(years_list, data)
    # add carbon capture
    data = add_carbon_capture(years_list, dea_sheet_names, data, tech_data)

    # adjust for inflation
    for x in data.index.get_level_values("technology"):
        if x in cost_year_2020 or x in manual_cost_year_assignments_2020:
            data.at[x, "currency_year"] = 2020
        elif x in cost_year_2019:
            data.at[x, "currency_year"] = 2019
        else:
            data.at[x, "currency_year"] = 2015

    # add heavy-duty assumptions, cost year is 2022
    data = get_dea_vehicle_data(snakemake.input.dea_vehicles, years_list, data)

    # add shipping data
    data = get_dea_maritime_data(snakemake.input.dea_ship, years_list, data)

    # %% (2) -- get data from other sources which need formatting -----------------
    # (a)  ---------- get old pypsa costs ---------------------------------------
    costs_pypsa = pd.read_csv(
        snakemake.input.pypsa_costs, index_col=[0, 2]
    ).sort_index()
    # rename some techs and convert units
    costs_pypsa = rename_pypsa_old(costs_pypsa)

    # (b1) ------- add vehicle costs from Fraunhofer vehicle study ------------------------
    costs_vehicles = pd.read_csv(
        snakemake.input.fraunhofer_vehicles_costs,
        engine="python",
        index_col=[0, 1],
        encoding="ISO-8859-1",
    )
    # rename + reorder to fit to other data
    costs_vehicles = rename_ISE_vehicles(costs_vehicles)
    if "NT" in costs_vehicles.index:
        costs_vehicles.drop(["NT"], axis=0, inplace=True, level=0)
    costs_vehicles = convert_units(years_list, costs_vehicles)
    # add costs for vehicles
    data = pd.concat([data, costs_vehicles], sort=True)

    # (b) ------- add costs from Fraunhofer ISE study --------------------------
    costs_ISE = pd.read_csv(
        snakemake.input.fraunhofer_costs,
        engine="python",
        index_col=[0, 1],
        encoding="ISO-8859-1",
    )
    # rename + reorder to fit to other data
    costs_ISE = rename_ISE(costs_ISE)
    # add costs for gas pipelines
    data = pd.concat([data, costs_ISE.loc[["Gasnetz"]]], sort=True)

    data = add_manual_input(data)
    # add costs for home batteries

    if snakemake.config["energy_storage_database"].get("ewg_home_battery", True):
        data = add_home_battery_costs(snakemake.input.EWG_costs, years_list, data)
    # add SMR assumptions
    data = add_SMR_data(years_list, data)
    # add solar rooftop costs by taking the mean of commercial and residential
    data = add_mean_solar_rooftop(years_list, data)

    data.index.names = ["technology", "parameter"]
    # %% (3) ------ add additional sources and save cost as csv ------------------
    # [RTD-target-multiindex-df]
    for year in years_list:
        costs = data[
            [year, "unit", "source", "further description", "currency_year"]
        ].rename(columns={year: "value"})
        costs["value"] = costs["value"].astype(float)

        # biomass is differentiated by biomass CHP and HOP
        costs.loc[("solid biomass", "fuel"), "value"] = 12
        costs.loc[("solid biomass", "fuel"), "unit"] = "EUR/MWh_th"
        costs.loc[("solid biomass", "fuel"), "source"] = (
            "JRC ENSPRESO ca avg for MINBIOWOOW1 (secondary forest residue wood chips), ENS_Ref for 2040"
        )
        costs.loc[("solid biomass", "fuel"), "currency_year"] = 2010

        costs.loc[("digestible biomass", "fuel"), "value"] = 15
        costs.loc[("digestible biomass", "fuel"), "unit"] = "EUR/MWh_th"
        costs.loc[("digestible biomass", "fuel"), "source"] = (
            "JRC ENSPRESO ca avg for MINBIOAGRW1, ENS_Ref for 2040"
        )
        costs.loc[("digestible biomass", "fuel"), "currency_year"] = 2010

        # add solar data from other source than DEA
        if any(
            [
                snakemake.config["solar_utility_from_vartiaien"],
                snakemake.config["solar_rooftop_from_etip"],
            ]
        ):
            costs = add_solar_from_other(years_list, costs)

        # add desalination and clean water tank storage
        costs = add_desalination_data(costs)
        # add energy storage database
        if snakemake.config["energy_storage_database"]["pnnl_energy_storage"].get(
            "add_data", True
        ):
            costs, tech = add_energy_storage_database(
                snakemake.input["pnnl_energy_storage"],
                snakemake.config["energy_storage_database"]["pnnl_energy_storage"],
                costs,
                year,
            )
            costs.loc[tech, "currency_year"] = 2020

        # add electrolyzer and fuel cell efficiency from other source than DEA
        if snakemake.config["energy_storage_database"].get("h2_from_budischak", True):
            costs = add_h2_from_other(costs)

        # CO2 intensity
        costs = add_co2_intensity(costs)

        # carbon balances
        costs = carbon_flow(years_list, costs, year)

        # energy penalty of carbon capture
        costs = energy_penalty(costs)

        # include old pypsa costs
        check = pd.concat([costs_pypsa, costs], sort=True)

        # missing technologies
        missing = costs_pypsa.index.levels[0].difference(costs.index.levels[0])
        if len(missing) & (year == years_list[0]):
            logger.info("************************************************************")
            logger.info("warning, in new cost assumptions the following components: ")
            for i in range(len(missing)):
                logger.info(f"{i + 1} {missing[i]}")
            logger.info(" are missing and the old cost assumptions are assumed.")
            logger.info("************************************************************")

        to_add = costs_pypsa.loc[missing].drop("year", axis=1)
        to_add.loc[:, "further description"] = " from old pypsa cost assumptions"
        # TODO check currency year from old pypsa cost assumptions
        to_add["currency_year"] = 2015
        costs_tot = pd.concat([costs, to_add], sort=False)

        # single components missing
        comp_missing = costs_pypsa.index.difference(costs_tot.index)
        if year == years_list[0]:
            logger.info(
                "single parameters of technologies are missing, using old PyPSA assumptions: "
            )
            logger.info(comp_missing)
            logger.info("old c_v and c_b values are assumed where given")
        to_add = costs_pypsa.loc[comp_missing].drop("year", axis=1)
        to_add.loc[:, "further description"] = " from old pypsa cost assumptions"
        # more data on geothermal is added downstream, so old assumptions are redundant
        to_add = to_add.drop("geothermal")
        # TODO check currency year from old pypsa cost assumptions
        to_add["currency_year"] = 2015
        costs_tot = pd.concat([costs_tot, to_add], sort=False)

        # unify the cost from DIW2010
        costs_tot = unify_diw(costs_tot)
        costs_tot.drop("fixed", level=1, inplace=True)

        # adjust for inflation
        techs = costs_tot.index.get_level_values(0).unique()
        costs_tot["currency_year"] = costs_tot.currency_year.astype(float)
        costs_tot = adjust_for_inflation(
            inflation_rate, costs_tot, techs, snakemake.config["eur_year"], "value"
        )

        # format and sort
        costs_tot.sort_index(inplace=True)
        costs_tot.loc[:, "value"] = round(
            costs_tot.value.astype(float), snakemake.config.get("ndigits", 2)
        )
        costs_tot.to_csv([v for v in snakemake.output if str(year) in v][0])
