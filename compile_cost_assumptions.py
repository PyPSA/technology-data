#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:21:20 2020

@author: bw0928
"""

import pandas as pd
import numpy as np
import os

# %% -------- PARAMETER ------------------------------------------------------
path_in = "/home/ws/bw0928/Dokumente/compile_costs_new/technology_data/inputs/"

years = np.arange(2020, 2055, 5)
rate_inflation = 0.02
solar_from_DEA = False  # add solar data from DEA if false from Vartiaien/ETIP
h2_from_budischak = False  # add fuel cell/electrolysis efficiencies from budischak
# remove grid connection costs from DEA for offwind because they are calculated
# seperately in pypsa-eur
offwind_no_gridcosts = True
# ---------- sources -------------------------------------------------------
source_DEA = 'DEA'
# solar
source_Vartiainen = 'Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity'
source_ETIP = 'European PV Technology and Innovation Platform'
# nuclear, coal, lignite from Lazards
source_Lazards = 'Lazard s Levelized Cost of Energy Analysis - Version 13.0'
# and the fuel cost is from Zappa's paper
zappa_paper = 'Is a 100% renewable European power system feasible by 2050?'
# co2 intensity
source_co2 = 'Entwicklung der spezifischen Kohlendioxid-Emissionen des deutschen Strommix in den Jahren 1990 - 2018'
# ---------------------------------------------------------------------------

sheet_names = {'onwind': '20 Onshore turbines',
               'offwind': '21 Offshore turbines',
               'solar-utility': '22 Photovoltaics Large',
               'solar-rooftop': '22 Photovoltaics Small',
               'OCGT': '52 OCGT - Natural gas',
               'CCGT': '05 Gas turb. CC, steam extract.',
               'oil': '50 Diesel engine farm',
               'biomass CHP': '09c Straw, Large, 40 degree',
               'biomass EOP': '09c Straw, Large, 40 degree',
               'biomass HOP': '09c Straw HOP',
               'central coal CHP': '01 Coal CHP',
               'central gas CHP': '04 Gas turb. simple cycle, L',
               'central solid biomass CHP': '09b Wood Pellets, Medium',
               'solar': '22 Photovoltaics Medium',
               'central air-sourced heat pump': '40 Comp. hp, airsource 3 MW',
               'central ground-sourced heat pump': '40 Absorption heat pump, DH',
               'central resistive heater': '41 Electric Boilers',
               'central gas boiler': '44 Natural Gas DH Only',
               'decentral gas boiler': '202 Natural gas boiler',
               'decentral ground-sourced heat pump': '207.7 Ground source existing',
               'decentral air-sourced heat pump': '207.3 Air to water existing',
               # 'decentral resistive heater': '216 Electric heating',
               'central water tank storage': '140 PTES seasonal',
               # 'decentral water tank storage': '142 Small scale hot water tank',
               'fuel cell': '12 LT-PEMFC CHP',
               'hydrogen storage underground': '151c Hydrogen Storage - Caverns',
               'hydrogen storage tank': '151a Hydrogen Storage - Tanks',
               'micro CHP': '219 LT-PEMFC mCHP - natural gas',
               'biogas upgrading': '82 Biogas, upgrading',
               'battery': '180 Lithium Ion Battery',
               'electrolysis': '88 Alkaline Electrolyser',
               # 'electricity distribution rural': '101 2 el distri Rural',
               # 'electricity distribution urban': '101 4 el distri  city',
               # 'gas distribution rural': '102 7 gas  Rural',
               # 'gas distribution urban': '102 9 gas City',
               # 'DH distribution rural': '103_12 DH_Distribu Rural',
               # 'DH distribution urban': '103_14 DH_Distribu City',
               # 'DH distribution low T': '103_16 DH_Distr New area LTDH',
               # 'gas pipeline': '102 6 gas Main distri line',
               # "DH main transmission": "103_11 DH transmission",
               }

# %% -------- FUNCTIONS ---------------------------------------------------

def get_excel_sheets(path_in):
    """"
    read all excel sheets of a given input path (path_in) and return
    them as a dictionary (data_in)
    """
    data_in = {}
    for entry in os.listdir(path_in):
        if entry[-5:] == ".xlsx":
            data_in[entry] = pd.ExcelFile(path_in + entry).sheet_names
    print("found ", len(data_in), " excel sheets: ")
    for key in data_in.keys():
        print("* ", key)
    return data_in


def get_sheet_location(tech, sheet_names, data_in):
    """
    looks up in which excel file technology is saved
    """
    for key in data_in:
        if sheet_names[tech] in data_in[key]:
            return key
    print("******* warning *************")
    print("tech ", tech, " with sheet name ", sheet_names[tech],
          "  not found in excel sheets.")
    print("****************************")
    return None

#
def get_data_DEA(tech):
    """
    interpolate cost for a given technology from DEA database sheet
    """
    excel_file = get_sheet_location(tech, sheet_names, data_in)
    if excel_file is None:
        print("excel file not found for tech ", tech)
        return None

    excel = pd.read_excel(path_in + excel_file,
                          sheet_name=sheet_names[tech],
                          index_col=0,
                          usecols='B:G', skiprows=[0, 1])

    excel.index = excel.index.fillna(" ")
    excel.dropna(axis=0, how="all", inplace=True)

    if 2020 not in excel.columns:
        excel.reset_index(inplace=True)
        excel.columns = (excel.loc[excel[excel == 2020].dropna(
            how="all").index] .iloc[0, :].fillna("Technology", limit=1))
        excel.drop(excel[excel == 2020].dropna(how="all").index, inplace=True)
        excel.set_index(excel.columns[0], inplace=True)

    parameters = ["efficiency", "investment", "Fixed O&M",
                  "Variable O&M", "production capacity for one unit",
                  "Output capacity expansion cost",
                  "Hydrogen output", "Cb coefficient",
                  "Cv coefficient",
                  "Distribution network costs", "Technical life",
                  "Energy storage expansion cost",
                  'Output capacity expansion cost (M€2015/MW)']

    df = pd.DataFrame()
    for para in parameters:
        attr = excel[excel.index.str.contains(para)]
        if len(attr) != 0:
            df = df.append(attr)
    df.index = df.index.str.replace('€', 'EUR')

    df = df.loc[:, excel.columns.isin(years)]

    # replace missing data
    df.replace("-", np.nan, inplace=True)
    # average data  in format "lower_value-upper_value"
    df = df.applymap(lambda x: (float((x).split("-")[0])+float((x).split("-")[1]))/2 if (type(x)==str and "-" in x) else x)
    # remove approx. symbol "~"
    df = df.applymap(lambda x: float(x.replace("~","")) if type(x)==str else x)

    df = df.astype(float)

    if (tech == "offwind") & offwind_no_gridcosts:
        df.loc['Nominal investment (MEUR/MW)'] -= excel.loc[' - of which grid connection']


    df_final = pd.DataFrame(index=df.index, columns=years)

    for index in df_final.index:
        values = np.interp(x=years, xp=df.columns.values.astype(float), fp=df.loc[index, :].values.astype(float))
        df_final.loc[index, :] = values

    df_final["source"] = source_DEA + ", " + excel_file
    df_final["unit"] = (df_final.rename(index=lambda x:
                                        x[x.rfind("(")+1: x.rfind(")")]).index.values)
    df_final.index = df_final.index.str.replace(r" \(.*\)","")

    return df_final

#
def add_conventional_data(costs):
    """"
    add technology data for convetional carriers from Lazards, DIW and BP
    """
    # nuclear from Lazards
    costs.loc[('nuclear', 'investment'), 'value'] = 8595 / \
        (1 + rate_inflation)**(2019 - 2015)
    costs.loc[('nuclear', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('nuclear', 'investment'), 'source'] = source_Lazards

    costs.loc[('nuclear', 'FOM'), 'value'] = 1.4
    costs.loc[('nuclear', 'FOM'), 'unit'] = "%/year"
    costs.loc[('nuclear', 'FOM'), 'source'] = source_Lazards

    costs.loc[('nuclear', 'VOM'), 'value'] = 3.5
    costs.loc[('nuclear', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('nuclear', 'VOM'), 'source'] = source_Lazards

    costs.loc[('nuclear', 'efficiency'), 'value'] = 0.33
    costs.loc[('nuclear', 'efficiency'), 'unit'] = "per unit"
    costs.loc[('nuclear', 'efficiency'), 'source'] = source_Lazards

    costs.loc[('nuclear', 'fuel'), 'value'] = 2.6
    costs.loc[('nuclear', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('nuclear', 'fuel'), 'source'] = source_Lazards
    costs.loc[('uranium', 'fuel'), 'value'] = 2.6
    costs.loc[('uranium', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('uranium', 'fuel'), 'source'] = source_Lazards

    costs.loc[('nuclear', 'lifetime'), 'value'] = 40
    costs.loc[('nuclear', 'lifetime'), 'unit'] = "years"
    costs.loc[('nuclear', 'lifetime'), 'source'] = source_Lazards

    # coal from Lazards and BP 2019
    costs.loc[('coal', 'investment'), 'value'] = 4162.5 / \
        (1 + rate_inflation)**(2019 - 2015)
    costs.loc[('coal', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('coal', 'investment'), 'source'] = source_Lazards

    costs.loc[('coal', 'FOM'), 'value'] = 1.6
    costs.loc[('coal', 'FOM'), 'unit'] = "%/year"
    costs.loc[('coal', 'FOM'), 'source'] = source_Lazards

    costs.loc[('coal', 'VOM'), 'value'] = 3.5
    costs.loc[('coal', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('coal', 'VOM'), 'source'] = source_Lazards

    costs.loc[('coal', 'efficiency'), 'value'] = 0.33
    costs.loc[('coal', 'efficiency'), 'unit'] = "per unit"
    costs.loc[('coal', 'efficiency'), 'source'] = source_Lazards

    costs.loc[('coal', 'fuel'), 'value'] = 8.15
    costs.loc[('coal', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('coal', 'fuel'), 'source'] = 'BP 2019'
    costs.loc[('gas', 'fuel'), 'value'] = 20.1
    costs.loc[('gas', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('gas', 'fuel'), 'source'] = 'BP 2019'

    costs.loc[('coal', 'lifetime'), 'value'] = 40
    costs.loc[('coal', 'lifetime'), 'unit'] = "years"
    costs.loc[('coal', 'lifetime'), 'source'] = source_Lazards

    # lignite from Lazards and DIW
    costs.loc[('lignite', 'investment'), 'value'] = 4162.5 / \
        (1 + rate_inflation)**(2019 - 2015)
    costs.loc[('lignite', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('lignite', 'investment'), 'source'] = source_Lazards

    costs.loc[('lignite', 'FOM'), 'value'] = 1.6
    costs.loc[('lignite', 'FOM'), 'unit'] = "%/year"
    costs.loc[('lignite', 'FOM'), 'source'] = source_Lazards

    costs.loc[('lignite', 'VOM'), 'value'] = 3.5
    costs.loc[('lignite', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('lignite', 'VOM'), 'source'] = source_Lazards

    costs.loc[('lignite', 'efficiency'), 'value'] = 0.33
    costs.loc[('lignite', 'efficiency'), 'unit'] = 'per unit'
    costs.loc[('lignite', 'efficiency'), 'source'] = source_Lazards

    costs.loc[('lignite', 'fuel'), 'value'] = 2.9
    costs.loc[('lignite', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('lignite', 'fuel'), 'source'] = 'DIW'

    costs.loc[('lignite', 'lifetime'), 'value'] = 40
    costs.loc[('lignite', 'lifetime'), 'unit']  = "years"
    costs.loc[('lignite', 'lifetime'), 'source'] = source_Lazards

    return costs


def add_co2_intensity(costs):
    """"
    add CO2 intensity for the carriers
    """
    TJ_to_MWh = 277.78
    costs.loc[('gas', 'CO2 intensity'), 'value'] = 55827 / 1e3 / TJ_to_MWh  # Erdgas
    costs.loc[('coal', 'CO2 intensity'), 'value'] = 93369 / 1e3 / TJ_to_MWh  # Steinkohle
    costs.loc[('lignite', 'CO2 intensity'), 'value'] = 113031 / 1e3 / TJ_to_MWh  # Rohbraunkohle Rheinland
    costs.loc[('oil', 'CO2 intensity'), 'value'] = 74020 / 1e3 / TJ_to_MWh  # Heizöl, leicht
    costs.at[('solid biomass', 'CO2 intensity'), 'value'] = 0.3

    costs.loc[('gas', 'CO2 intensity'), 'source'] = source_co2
    costs.loc[('coal', 'CO2 intensity'), 'source'] = source_co2
    costs.loc[('lignite', 'CO2 intensity'), 'source'] = source_co2
    costs.loc[('oil', 'CO2 intensity'), 'source'] = source_co2
    costs.at[('solid biomass', 'CO2 intensity'), 'source'] = "TODO"

    costs.loc[pd.IndexSlice[:, "CO2 intensity"], "unit"] = "tCO2/MWh_th"

    return costs


def add_solar_from_other(costs):
    """"
    add solar from other sources than DEA (since they are very optimistic)
    """
    # solar utility from Vartiaian 2019
    data = np.interp(x=years, xp=[2020, 2030, 2040, 2050],
                     fp=[431, 275, 204, 164])
    # the paper says 'In this report, all results are given in real 2019
    # money.'
    data = data / (1 + rate_inflation)**(2019 - 2015)
    solar_uti = pd.Series(data=data, index=years)

    # solar rooftop from ETIP 2019
    data = np.interp(x=years, xp=[2020, 2030, 2050], fp=[1150, 800, 550])
    # using 2016 money in page 10
    data = data / (1 + rate_inflation)**(2016 - 2015)
    solar_roof = pd.Series(data=data, index=years)

    # solar utility from Vartiaian 2019
    costs.loc[('solar-utility', 'investment'), 'value'] = solar_uti[year]
    costs.loc[('solar-utility', 'investment'), 'source'] = source_Vartiainen

    costs.loc[('solar-utility', 'lifetime'), 'value'] = 30
    costs.loc[('solar-utility', 'lifetime'), 'source'] = source_Vartiainen

    # solar rooftop from ETIP 2019 - convert EUR/kW -> EUR/MW
    costs.loc[('solar-rooftop', 'investment'), 'value'] = solar_roof[year]
    costs.loc[('solar-rooftop', 'investment'), 'source'] = source_ETIP

    costs.loc[('solar-rooftop', 'lifetime'), 'value'] = 30
    costs.loc[('solar-rooftop', 'lifetime'), 'source'] = source_ETIP

    # lifetime&efficiency for solar
    costs.loc[('solar', 'lifetime'), 'value'] = costs.loc[(
        ['solar-rooftop', 'solar-utility'], 'lifetime'), 'value'].mean()
    costs.loc[('solar', 'lifetime'), 'unit'] = 'years'
    costs.loc[('solar', 'lifetime'),
              'source'] = 'Assuming 50% rooftop, 50% utility'
    # costs.loc[('solar', 'efficiency'), 'value'] = 1
    # costs.loc[('solar', 'efficiency'), 'unit'] = 'per unit'

    return costs


def add_h2_from_other(costs):
    """
    assume higher efficiency for electrolysis(0.8) and fuel cell(0.58)
    """
    costs.loc[('electrolysis', 'efficiency'), 'value'] = 0.8
    costs.loc[('fuel cell', 'efficiency'), 'value'] = 0.58
    costs.loc[('electrolysis', 'efficiency'), 'source'] = 'budischak2013'
    costs.loc[('fuel cell', 'efficiency'), 'source'] = 'budischak2013'

    return costs


def add_costs_ccs(costs, techs_ccs=["central solid biomass CHP",
                                    "central gas CHP"]  # SMR
                  ):
    """"
    add costs and efficiencies from DIW for CCS for technologies 'techs_css'
    """
    for tech_ccs in techs_ccs:
        name = tech_ccs + " CCS"
        costs = costs.append(costs.loc[tech_ccs].set_index(
            pd.MultiIndex.from_product([[name], costs.loc[tech_ccs].index])))
        costs.loc[(name, 'efficiency'), 'value'] *= 0.9
        # costs extra for CCS from DIW
        costs.loc[(name, 'investment'), 'value'] += 600
        costs.loc[(name, 'investment'), 'source'] += " , DIW (CCS)"

    return costs


def unify_diw(costs):
    """"
    include inflation for the DIW costs from 2010
    """
    inflation = (1 + rate_inflation)**(2010 - 2015)
    costs.loc[('PHS', 'investment'), 'value'] /= inflation
    costs.loc[('ror', 'investment'), 'value'] /= inflation
    costs.loc[('hydro', 'investment'), 'value'] /= inflation

    return costs
# %% ---------- MAIN -------------------------------------------------------
# --------- get data from DEA excel sheets -----------------------------------
data_in = get_excel_sheets(path_in)
d_by_tech = {}

for tech in sheet_names.keys():
    print(tech + ' in PyPSA corresponds to ' + sheet_names[tech] +
          ' in DEA database.')
    df = get_data_DEA(tech=tech).fillna(0)
    # df.fillna(value=0, inplace=True)
    d_by_tech[tech] = df
# %% --------- clean up units --------------------------------
tech_data = pd.concat(d_by_tech)

tech_data.unit = tech_data.unit.str.replace(" per ", "/")
tech_data.unit = tech_data.unit.str.replace(" / ", "/")
tech_data.unit = tech_data.unit.str.replace("J/s", "W")

# units
tech_data.loc[tech_data.unit.str.contains("MEUR"), years] *= 1e6
tech_data.unit = tech_data.unit.str.replace("MEUR", "EUR")

tech_data.loc[tech_data.unit.str.contains("1000EUR"), years] *= 1e3
tech_data.unit = tech_data.unit.str.replace("1000EUR", "EUR")

tech_data.loc[tech_data.unit.str.contains("kEUR"), years] *= 1e3
tech_data.unit = tech_data.unit.str.replace("kEUR", "EUR")

tech_data.loc[tech_data.unit.str.contains("kW"), years] /= 1e3
tech_data.unit = tech_data.unit.str.replace("kW", "MW")

tech_data.loc[tech_data.unit.str.contains("/GWh"), years] /= 1e3
tech_data.unit = tech_data.unit.str.replace("/GWh", "/MWh")

tech_data.loc[tech_data.unit.str.contains("/GJ"), years] *= 3.6
tech_data.unit = tech_data.unit.str.replace("/GJ", "/MWh")

tech_data.unit = tech_data.unit.str.replace(" a year", "/year")
tech_data.unit = tech_data.unit.str.replace("2015EUR", "EUR")
tech_data.unit = tech_data.unit.str.replace("2015-EUR", "EUR")
tech_data.unit = tech_data.unit.str.replace("EUR2015", "EUR")
tech_data.unit = tech_data.unit.str.replace("EUR-2015", "EUR")
tech_data.unit = tech_data.unit.str.replace("MWe", "MW_e")
tech_data.unit = tech_data.unit.str.replace("MWth", "MW_th")
tech_data.unit = tech_data.unit.str.replace("MWheat", "MW_th")
tech_data.unit = tech_data.unit.str.replace("MWhheat", "MWh_th")
tech_data.loc[tech_data.unit=='EUR/MW/y', "unit"] = 'EUR/MW/year'

# convert per unit costs to MW
techs_per_unit = tech_data.xs("Heat production capacity for one unit",
                             level=1).index
for tech in techs_per_unit:
    df = tech_data.loc[tech]
    cap = df.loc["Heat production capacity for one unit"]
    df.loc[df.unit.str.contains("/unit"), years] /= cap.loc[years]
    df.loc[df.unit.str.contains("/unit"), "unit"] = df.loc[df.unit.str.contains("/unit"), "unit"].str.replace("/unit", "/"+cap.unit+"_th")


# clarify MW -> MW_th
# see on p.278 of docu: "However, the primary purpose of the heat pumps in the
# technology catalogue is heating. In this chapter the unit MW is referring to
# the heat output (also MJ/s) unless otherwise noted"
techs_mwth = ['central air-sourced heat pump', 'central gas boiler',
              'central resistive heater', 'decentral air-sourced heat pump',
              'decentral gas boiler', 'decentral ground-sourced heat pump' ]
tech_data.loc[techs_mwth, "unit"] = (tech_data.loc[techs_mwth, "unit"]
                                     .replace({"EUR/MW": "EUR/MW_th",
                                               "EUR/MW/year": "EUR/MW_th/year",
                                               'EUR/MWh':'EUR/MWh_th',
                                               "MW": "MW_th"}))

# clarify MW -> MW_e
techs_e = ['fuel cell']
tech_data.loc[techs_e, "unit"] = (tech_data.loc[techs_e, "unit"]
                                  .replace({"EUR/MW": "EUR/MW_e",
                                            "EUR/MW/year": "EUR/MW_e/year",
                                            'EUR/MWh':'EUR/MWh_e',
                                             "MW": "MW_e"}))



# %% ------ specific assumptions for some technologies ------------------------

# for central resistive heater there are investment costs for small (1-5MW) and
# large (>10 MW) generators, assume the costs for large generators
to_drop = [("central resistive heater", 'Nominal investment, 400/690 V; 1-5 MW')]

# for decentral gas boilers total and heat efficiency given, the values are
# the same, drop one of the rows to avoid duplicates
to_drop.append(("decentral gas boiler", "Heat efficiency, annual average, net"))

# for decentral gas boilers there are investment costs and possbile additional
# investments which apply for grid connection if the house is not connected yet
# those costs are currently excluded with the assumption that there are no new
# decentral gas boilers build in houses with no gas grid connection
to_drop.append(("decentral gas boiler", "Specific investment"))

# hydrogen storage assume round trip efficiency
to_drop.append(("hydrogen storage tank", ' - Charge efficiency'))
to_drop.append(("hydrogen storage tank", ' - Discharge efficiency'))
to_drop.append(("hydrogen storage underground", ' - Charge efficiency'))
to_drop.append(("hydrogen storage underground", ' - Discharge efficiency'))
tech_data.loc[("hydrogen storage underground", "Round trip efficiency"), years] *= 100
tech_data.loc[("hydrogen storage tank", "Round trip efficiency"), years] *= 100
# drop PV module conversion efficiency
tech_data = tech_data.drop("PV module conversion efficiency", level=1)

# biogas upgrading investment costs should include grid injection costs
index = tech_data.loc["biogas upgrading"].index.str.contains("investment")
name = 'investment (upgrading, methane redution and grid injection)'
inv = tech_data.loc["biogas upgrading"].loc[index].groupby(["unit", "source"]).sum().reset_index()
new = pd.concat([tech_data.loc["biogas upgrading"].loc[~index],
                                               inv]).rename({0:name})
new.index = pd.MultiIndex.from_product([["biogas upgrading"],
                                        new.index.to_list()])
tech_data.drop("biogas upgrading", level=0, inplace=True)
tech_data = pd.concat([tech_data, new])

# heat pump efficiencies are assumed the one's for existing building,
# in the DEA they do differ between heating the floor area or heating with radiators,
# since most households heat with radiators and there efficiencies are lower
# (conservative approach) those are assumed
# furthermore the total efficiency is assumed which includes auxilary electricity
# consumption
name = 'Heat efficiency, annual average, net, radiators'
techs_radiator = tech_data.xs(name, level=1).index
for tech in techs_radiator:
    df = tech_data.loc[tech]
    df = df[(~df.index.str.contains("efficiency")) | (df.index==name)]
    df.rename(index={name: name + ", existing one family house"}, inplace=True)
    df.index = pd.MultiIndex.from_product([[tech], df.index.to_list()])
    tech_data.drop(tech, level=0, inplace=True)
    tech_data = pd.concat([tech_data, df])

# battery split into inverter and storage, assume for efficiency sqr(round trip DC)
df = tech_data.loc["battery"]
inverter = df.loc[['Round trip efficiency DC',
                   'Output capacity expansion cost',
                   'Technical lifetime', 'Fixed O&M']]

inverter.rename(index ={'Output capacity expansion cost':
                        'Output capacity expansion cost investment'},
                inplace=True)
inverter.index = pd.MultiIndex.from_product([["battery inverter"],
                                             inverter.index.to_list()])
storage = df.reindex(index=['Technical lifetime',
                            'Energy storage expansion cost'])
storage.rename(index={'Energy storage expansion cost':
                      'Energy storage expansion cost investment'}, inplace=True)
storage.index = pd.MultiIndex.from_product([["battery storage"],
                                            storage.index.to_list()])
tech_data.drop("battery", level=0, inplace=True)
tech_data = pd.concat([tech_data, inverter, storage])

tech_data = tech_data.drop(to_drop)

# drop all rows which only contains zeros
tech_data = tech_data.loc[(tech_data[years]!=0).sum(axis=1)!=0]

# %%
clean_df = {}
for tech in tech_data.index.levels[0]:
    clean_df[tech] = pd.DataFrame()
    switch = False
    df = tech_data.loc[tech]

    # --- investment ----
    investment = df[(df.index.str.contains("investment") |
                     df.index.str.contains("Distribution network costs"))
                 & ((df.unit=="EUR/MW")|
                    (df.unit=="EUR/MW_e")|
                    (df.unit=="EUR/MW_th - heat output")|
                    (df.unit=="EUR/MW_th excluding drive energy")|
                    (df.unit=="EUR/MW_th") |
                    (df.unit=="EUR/MWhCapacity") |
                    (df.unit=="EUR/MWh") |
                    (df.unit=="EUR/MWh/year") |
                    (df.unit=="EUR/MW input"))]
    if len(investment)!=1:
        switch = True
        print("check investment: ", tech, " ",
              df[df.index.str.contains("investment")].unit)
    else:
        investment.loc[:, "parameter"] = "investment"
        clean_df[tech] = investment

    # ---- FOM ----------------
    if len(investment):
        fixed = df[df.index.str.contains("Fixed O&M") &
                   ((df.unit==investment.unit[0]+"/year")|
                    (df.unit=="EUR/MW/km/year")|
                    (df.unit=="EUR/MW/year")|
                    (df.unit==investment.unit.str.split(" ")[0][0]+"/year"))]
        if (len(fixed)!=1) and (len(df[df.index.str.contains("Fixed O&M")])!=0):
            switch = True
            print("check FOM: ", tech, " ",
                  df[df.index.str.contains("Fixed O&M")].unit)
        if len(fixed) == 1:
            fixed.loc[:,"parameter"] = "fixed"
            clean_df[tech] = pd.concat([clean_df[tech], fixed])
            fom = pd.DataFrame(columns=fixed.columns)
            fom[years] = fixed[years]/investment[years].values*100
            fom["parameter"] = "FOM"
            fom["unit"] = "%/year"
            fom["source"] = fixed["source"]
            clean_df[tech] = pd.concat([clean_df[tech], fom])

    # ---- VOM -----
    vom = df[df.index.str.contains("Variable O&M")& ((df.unit=="EUR/MWh") |
                                                     (df.unit=="EUR/MWh_e") |
                                                     (df.unit=="EUR/MWh_th") |
                                                     (df.unit=="EUR/MWh/year") |
                                                     (df.unit=="EUR/MWh/km") |
                                                     (df.unit=="EUR/MWh") |
                                                     (df.unit=="EUR/MWhoutput") |
                                                     (tech == "biogas upgrading"))]
    if len(vom)!=1 and len(df[df.index.str.contains("Variable O&M")])!=0:
        switch = True
        print("check VOM: ", tech, " ",
              df[df.index.str.contains("Variable O&M")].unit)
    if len(vom)==1:
        vom.loc[:,"parameter"] = "VOM"
        clean_df[tech] = pd.concat([clean_df[tech], vom])

    # ----- lifetime --------
    lifetime = df[df.index.str.contains("Technical life") & (df.unit=="years")]
    if len(lifetime)!=1:
        switch  = True
        print("check lifetime: ", tech, " ",
              df[df.index.str.contains("Technical life")].unit)
    else:
        lifetime.loc[:,"parameter"] = "lifetime"
        clean_df[tech] = pd.concat([clean_df[tech], lifetime])


    # ----- efficiencies ------
    efficiency = df[(df.index.str.contains("efficiency") |
                     df.index.str.contains("Hydrogen output, at LHV"))
                     & ((df.unit=="%") |  (df.unit =="% total size"))
                     & (~df.index.str.contains("name plate"))]

    # take annual average instead of name plate efficiency
    if any(efficiency.index.str.contains("annual average")):
        efficiency = efficiency[efficiency.index.str.contains("annual average")]

    # check if electric and heat efficiencies are given
    if (any(["Electric" in ind for ind in efficiency.index]) and
        any(["Heat" in ind for ind in efficiency.index])):
        print("eff heat and electric in ", tech)
        efficiency_heat = efficiency[efficiency.index.str.contains("Heat")]
        efficiency_heat.loc[:,"parameter"] = "efficiency-heat"
        clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
        efficiency = efficiency[efficiency.index.str.contains("Electric")]
        efficiency.loc[:,"parameter"] = "efficiency"
        clean_df[tech] = pd.concat([clean_df[tech], efficiency])

    elif len(efficiency)!=1:
        switch  = True
        print("check efficiency: ", tech, " ",
              df[df.index.str.contains("efficiency")].unit)
    else:
        efficiency.loc[:,"parameter"] = "efficiency"
        clean_df[tech] = pd.concat([clean_df[tech], efficiency])

    # add c_v and c_b coefficient
    if "Cb coefficient" in df.index:
        c_b = df.loc[df.index.str.contains("Cb coefficient")].dropna()
        if len(c_b):
            c_b.loc[:, "parameter"] = "c_b"
            clean_df[tech] = pd.concat([clean_df[tech], c_b])
    if "Cv coefficient" in df.index:
        c_v = df.loc[df.index.str.contains("Cv coefficient")].dropna()
        if len(c_v):
            c_v.loc[:, "parameter"] = "c_v"
            clean_df[tech] = pd.concat([clean_df[tech], c_v])

    if switch:
        print("---------------------------------------")

# %% -------- concat data and convert units to compare with old--------------
data = (pd.concat(clean_df).reset_index().rename(columns={"level_0":"technology",
                                                          "level_1": "further description"})
        .set_index(["technology", "parameter"]))

# add water tank charger/ discharger
charger = tech_data.loc[("central water tank storage", "Round trip efficiency")]
charger["further description"] = "efficiency from sqr(Round trip efficiency)"
charger[years] = charger[years]**0.5*10
charger.rename(index={"Round trip efficiency": "efficiency"},
               level=1, inplace=True)
charger.rename(index={'central water tank storage':"water tank charger"},
               level=0, inplace=True)
data = pd.concat([data, charger])
charger.rename(index={"water tank charger": "water tank discharger"},
               level=0, inplace=True)
data = pd.concat([data, charger])

# add excel sheet names to data frame
wished_order = list(years) + ["unit", "source", "further description"]
data = data.reindex(columns=wished_order)
sheets = data.reset_index()["technology"].map(sheet_names).fillna("")
sheets.index = data.index
data["further description"] = sheets + ":  " + data["further description"]

# add comment for offwind investment
if offwind_no_gridcosts:
    data.loc[("offwind", "investment"),
             "further description"] += " grid connection costs substracted from investment costs"

# convert efficiency from % -> per unit
data.loc[data.index.get_level_values(1).isin(["efficiency", "efficiency-heat"])
         , years] /= 100
data.loc[data.index.get_level_values(1).isin(["efficiency", "efficiency-heat"])
         , "unit"] = "per unit"

# convert MW -> kW
to_convert = (data.index.get_level_values(1).isin(["fixed", "investment"]) &
              data.unit.str.contains("/MW"))
data.loc[to_convert, years] /= 1e3
data.loc[to_convert, "unit"] = (data.loc[to_convert, "unit"].str
                               .replace("/MW","/kW"))

# %% ------------ get old pypsa costs ---------------------------------------
costs_pypsa = pd.read_csv('inputs/costs_PyPSA.csv',
                          index_col=[0,2]).sort_index()
to_drop = [
           #  'hydrogen storage',
           # 'hydrogen underground storage',
           'retrofitting I', 'retrofitting II']
costs_pypsa.drop(to_drop, level=0, inplace=True)

# central CHP is gas-fired
costs_pypsa.rename({'central CHP': 'central gas CHP'}, inplace=True)
costs_pypsa.rename({'hydrogen storage': 'hydrogen storage tank'}, inplace=True)
costs_pypsa.rename({'hydrogen underground storage': 'hydrogen storage underground'}, inplace=True)

#convert EUR/m^3 to EUR/kWh for 40 K diff and 1.17 kWh/m^3/K
costs_pypsa.loc[('decentral water tank storage','investment'),
                'value'] /= 1.17*40
costs_pypsa.loc[('decentral water tank storage','investment'),'unit'] = 'EUR/kWh'

# %% ------ add additional sources and save cost.csv ------------------
for year in years:
    costs = (data[[year, "unit", "source", "further description"]]
             .rename(columns={year: "value"}))
    costs["value"] = costs["value"].astype(float)

    # biomass is differentiated by biomass CHP and HOP
    costs.loc[('solid biomass', 'fuel'), 'value'] = 25.2
    costs.loc[('solid biomass', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('solid biomass', 'fuel'), 'source'] = zappa_paper

    # add solar data from other source than DEA
    if not solar_from_DEA:
        costs = add_solar_from_other(costs)

    if h2_from_budischak:
        costs = add_h2_from_other(costs)

    # add data from conventional carriers
    costs = add_conventional_data(costs)
    # CO2 intensity
    costs = add_co2_intensity(costs)
    # CCS
    costs = add_costs_ccs(costs)

    # include old pypsa costs
    check = pd.concat([costs_pypsa, costs], axis=1)

    # missing technologies
    missing = costs_pypsa.index.levels[0].difference(costs.index.levels[0])
    if len(missing):
        print("************************************************************")
        print("warning, in new cost assumptions the following components: ")
        for i in range(len(missing)):
            print("    ", i + 1, missing[i])
        print(" are missing and the old cost assumptions are assumed.")
        print("************************************************************")

    to_add = costs_pypsa.loc[missing].drop("year", axis=1)
    to_add.loc[:,"further description"] = " from old pypsa cost assumptions"
    costs_tot = pd.concat([costs, to_add], sort=False)

    # single components missing
    comp_missing = costs_pypsa.index.difference(costs_tot.index)
    print("single components missing: ")
    print(comp_missing)
    to_add = costs_pypsa.loc[comp_missing].drop("year", axis=1)
    to_add.loc[:, "further description"] = " from old pypsa cost assumptions"
    costs_tot = pd.concat([costs_tot, to_add], sort=False)

    # take c_v and c_b values from old cost assumptions TODO check again!
    print("old c_v and c_b values are assumed where given")
    c_value_index =  costs_pypsa[costs_pypsa.index.get_level_values(1).isin(
                                ["c_b", "c_v"])].index
    costs_tot.loc[c_value_index,
                  ["value", "unit", "source"]] = costs_pypsa.loc[c_value_index,
                                                       ["value", "unit", "source"]]
    costs_tot.loc[c_value_index, "further description"] = " from old pypsa cost assumptions"

    # unify the cost from DIW2010
    costs_tot = unify_diw(costs_tot)
    costs_tot.drop("fixed", level=1, inplace=True)
    costs_tot.sort_index(inplace=True)
    costs_tot = round(costs_tot, ndigits=2)
    costs_tot.to_csv("outputs/costs_{}.csv".format(year))

# %% open questions:
# c_b and c_v values very different!
# battery inverter efficiency much higher in DEA and investment lower
# electrolyser and learning rates
# decentral resistive heater, decentral water tanks currently from old pypsa costs

# TODO
# unify units from old and new (convert EUR/m² to kW)
# c_b and c_v values
# check battery
# look up solar thermal, oil boiler, decentral resistive heater in DEA
# add offwind grid component connection costs (sea cable ...)

