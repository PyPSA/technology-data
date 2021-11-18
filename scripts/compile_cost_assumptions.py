#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script creates cost csv for choosen years from different source (source_dict).
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

import pandas as pd
import numpy as np

years = snakemake.config['years']

# ---------- sources -------------------------------------------------------
source_dict = {
                'DEA': 'Danish Energy Agency',
                # solar utility
                'Vartiaien': 'Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity',
                # solar rooftop
                'ETIP': 'European PV Technology and Innovation Platform',
                # nuclear, coal, lignite
                'Lazards': 'Lazard s Levelized Cost of Energy Analysis - Version 13.0',
                # fuel cost
                'zappa':  'Is a 100% renewable European power system feasible by 2050?',
                # co2 intensity
                "co2" :'Entwicklung der spezifischen Kohlendioxid-Emissionen des deutschen Strommix in den Jahren 1990 - 2018',
                # gas pipeline costs
                "ISE": "WEGE ZU EINEM KLIMANEUTRALEN ENERGIESYSEM, Anhang zur Studie, Fraunhofer-Institut für Solare Energiesysteme ISE, Freiburg",
                # Water desalination costs
                "Caldera2016": "Caldera et al 2016: Local cost of seawater RO desalination based on solar PV and windenergy: A global estimate. (https://doi.org/10.1016/j.desal.2016.02.004)",
                "Caldera2017": "Caldera et al 2017: Learning Curve for Seawater Reverse Osmosis Desalination Plants: Capital Cost Trend of the Past, Present, and Future (https://doi.org/10.1002/2017WR021402)",
                # home battery storage and inverter investment costs
                "EWG": "Global Energy System based on 100% Renewable Energy, Energywatchgroup/LTU University, 2019",
                "HyNOW" : "Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014"
                }

# [DEA-sheet-names]
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
               'central solid biomass CHP': '09a Wood Chips, Large 50 degree',#'09b Wood Pellets, Medium',
               # TODO: Dummy for now:
               # 'central solid biomass CHP CCS': '09b Wood Pellets, Medium',
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
               'hydrogen storage tank incl. compressor': '151a Hydrogen Storage - Tanks',
               'micro CHP': '219 LT-PEMFC mCHP - natural gas',
               'biogas' : '81 Biogas Plant, Basic conf.',
               'biogas upgrading': '82 Biogas, upgrading',
               'battery': '180 Lithium Ion Battery',
               'industrial heat pump medium temperature': '302.a High temp. hp Up to 125 C',
               'industrial heat pump high temperature': '302.b High temp. hp Up to 150',
               'electric boiler steam': '310.1 Electric boiler steam  ',
               'gas boiler steam': '311.1c Steam boiler Gas',
               'solid biomass boiler steam': '311.1e Steam boiler Wood',
               'electrolysis': '86 AEC 100MW', #'88 Alkaline Electrolyser',
               'direct air capture' : '403.a Direct air capture',
               'biomass CHP capture' : '401.a Post comb - small CHP',
               'cement capture' : '401.c Post comb - Cement kiln',
               'BioSNG' : '84 Gasif. CFB, Bio-SNG',
               'BtL' : '85 Gasif. Ent. Flow FT, liq fu ',
               'methanolisation': '98 Methanol from power',
               'Fischer-Tropsch': '102 Hydrogen to Jet',
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
# [DEA-sheet-names]

uncrtnty_lookup = {'onwind': 'J:K',
                    'offwind': 'J:K',
                    'solar-utility': 'J:K',
                    'solar-rooftop': '',
                    'OCGT': 'I:J',
                    'CCGT': 'I:J',
                    'oil': 'I:J',
                    'biomass CHP': 'I:J',
                    'biomass EOP': 'I:J',
                    'biomass HOP': 'I:J',
                    'central coal CHP': '',
                    'central gas CHP': 'I:J',
                    'central solid biomass CHP': 'I:J',
                    # 'central solid biomass CHP CCS': 'I:J',
                    'solar': '',
                    'central air-sourced heat pump': 'J:K',
                    'central ground-sourced heat pump': 'I:J',
                    'central resistive heater': 'I:J',
                    'central gas boiler': 'I:J',
                    'decentral gas boiler': 'I:J',
                    'decentral ground-sourced heat pump': 'I:J',
                    'decentral air-sourced heat pump': 'I:J',
                    'central water tank storage': 'J:K',
                    'fuel cell': 'I:J',
                    'hydrogen storage underground': 'J:K',
                    'hydrogen storage tank incl. compressor': 'J:K',
                    'micro CHP': 'I:J',
                    'biogas': 'I:J',
                    'biogas upgrading': 'I:J',
                    'electrolysis': 'I:J',
                    'battery': 'L,N',
                    'direct air capture': 'I:J',
                    'cement capture': 'I:J',
                    'biomass CHP capture': 'I:J',
                    'BioSNG' : 'I:J',
                    'BtL' : 'J:K',
                    'industrial heat pump medium temperature':'H:I',
                    'industrial heat pump high temperature':'H:I',
                    'electric boiler steam':'H:I',
                    'gas boiler steam':'H:I',
                    'solid biomass boiler steam':'H:I',
                    'Fischer-Tropsch': 'I:J',
                    'methanolisation': 'J:K',
}


# %% -------- FUNCTIONS ---------------------------------------------------

def get_excel_sheets(excel_files):
    """"
    read all excel sheets and return
    them as a dictionary (data_in)
    """
    data_in = {}
    for entry in excel_files:
        if entry[-5:] == ".xlsx":
            data_in[entry] = pd.ExcelFile(entry).sheet_names
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
def get_data_DEA(tech, data_in, expectation=None):
    """
    interpolate cost for a given technology from DEA database sheet

    uncertainty can be "optimist", "pessimist" or None|""
    """
    excel_file = get_sheet_location(tech, sheet_names, data_in)
    if excel_file is None:
        print("excel file not found for tech ", tech)
        return None

    if tech=="battery":
        usecols = "B:J"
    elif tech in ['direct air capture', 'cement capture', 'biomass CHP capture']:
        usecols = "A:F"
    elif tech in ['industrial heat pump medium temperature','industrial heat pump high temperature','electric boiler steam','gas boiler steam','solid biomass boiler steam']:
        usecols = "A:E"
    elif tech in ['Fischer-Tropsch']:
        usecols = "B:F"
    else:
        usecols = "B:G"

    usecols += f",{uncrtnty_lookup[tech]}"

    excel = pd.read_excel(excel_file,
                          sheet_name=sheet_names[tech],
                          index_col=0,
                          usecols=usecols,
                          skiprows=[0, 1],
                          na_values="N.A")

    excel.dropna(axis=1, how="all", inplace=True)


    excel.index = excel.index.fillna(" ")
    excel.index = excel.index.astype(str)
    excel.dropna(axis=0, how="all", inplace=True)

    print(excel)
    if 2020 not in excel.columns:
        selection = excel[excel.isin([2020])].dropna(how="all").index
        excel.columns = excel.loc[selection].iloc[0, :].fillna("Technology", limit=1)
        excel.drop(selection, inplace=True)

    uncertainty_columns = ["2050-optimist", "2050-pessimist"]
    if uncrtnty_lookup[tech]:
        # hydrogen storage sheets have reverse order of lower/upper estimates
        if tech in ["hydrogen storage tank incl. compressor", "hydrogen storage cavern"]:
            uncertainty_columns.reverse()
        excel.rename(columns={excel.columns[-2]: uncertainty_columns[0],
                                excel.columns[-1]: uncertainty_columns[1]
                                }, inplace=True)
    else:
        for col in uncertainty_columns:
            excel.loc[:,col] = excel.loc[:,2050]

    swap_patterns = ["technical life", "efficiency", "Hydrogen output, at LHV"] # cases where bigger is better
    swap = [any(term in idx.lower() for term in swap_patterns) for idx in excel.index]
    tmp = excel.loc[swap, "2050-pessimist"]
    excel.loc[swap, "2050-pessimist"] = excel.loc[swap, "2050-optimist"]
    excel.loc[swap, "2050-optimist"] = tmp

    if expectation:
        excel.loc[:,2050] = excel.loc[:,f"2050-{expectation}"].combine_first(excel.loc[:,2050])
    excel.drop(columns=uncertainty_columns, inplace=True)

    # fix for battery with different excel sheet format
    if tech == "battery":
        excel.rename(columns={"Technology":2040}, inplace=True)

    if expectation:
        excel = excel.loc[:,[2020,2050]]

    parameters = ["efficiency", "investment", "Fixed O&M",
                  "Variable O&M", "production capacity for one unit",
                  "Output capacity expansion cost",
                  "Hydrogen output",
                  "Hydrogen (% total input_e (MWh / MWh))",
                  "Cb coefficient",
                  "Cv coefficient",
                  "Distribution network costs", "Technical life",
                  "Energy storage expansion cost",
                  'Output capacity expansion cost (M€2015/MW)',
                  'Heat input', 'Heat  input', 'Electricity input', 'Eletricity input', 'Heat out',
                  'capture rate', 
                  "FT Liquids Output, MWh/MWh Total Input",
                  "Bio SNG (% of fuel input)",
                  "Total O&M"]


    df = pd.DataFrame()
    for para in parameters:
        # attr = excel[excel.index.str.contains(para)]
        attr = excel[[para in index for index in excel.index]]
        if len(attr) != 0:
            df = df.append(attr)
    df.index = df.index.str.replace('€', 'EUR')

    df = df.reindex(columns=df.columns[df.columns.isin(years)])
    df = df[~df.index.duplicated(keep='first')]

    # replace missing data
    df.replace("-", np.nan, inplace=True)
    # average data  in format "lower_value-upper_value"
    df = df.applymap(lambda x: (float((x).split("-")[0])
                                + float((x).split("-")[1]))/2 if (type(x)==str and "-" in x) else x)
    # remove symbols "~", ">", "<"
    for sym in ["~", ">", "<"]:
        df = df.applymap(lambda x: x.replace(sym,"") if type(x)==str else x)

    df = df.astype(float)

    if (tech == "offwind") and snakemake.config['offwind_no_gridcosts']:
        df.loc['Nominal investment (MEUR/MW)'] -= excel.loc[' - of which grid connection']

    # Exlucde indirect costs for centralised system with additional piping.
    if tech.startswith('industrial heat pump'):
        df = df.drop('Indirect investments cost (MEUR per MW)')

    if tech == 'methanolisation':
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech == 'Fischer-Tropsch':
        df.drop(df.loc[df.index.str.contains("l FT Liquids")].index, inplace=True)

    df_final = pd.DataFrame(index=df.index, columns=years)

    # [RTD-interpolation-example]
    for index in df_final.index:
        values = np.interp(x=years, xp=df.columns.values.astype(float), fp=df.loc[index, :].values.astype(float))
        df_final.loc[index, :] = values

    # if year-specific data is missing and not fixed by interpolation fill forward with same values
    df_final = df_final.fillna(method='ffill', axis=1)

    df_final["source"] = source_dict["DEA"] + ", " + excel_file.replace("inputs/","")
    df_final["unit"] = (df_final.rename(index=lambda x:
                                        x[x.rfind("(")+1: x.rfind(")")]).index.values)
    df_final.index = df_final.index.str.replace(r" \(.*\)","", regex=True)

    return df_final

def add_desalinsation_data(costs):
    """
    add technology data for sea water desalination (SWRO) and water storage.
    """

    # Interpolate cost based on historic costs/cost projection to fitting year
    cs = [2070,1917,1603,1282,1025] # in USD/(m^3/d)
    ys = [2015,2022,2030,2040,2050]
    c = np.interp(year, ys, cs)
    c *= 24                         # in USD/(m^3/h)
    c /= 1.17                       # in EUR/(m^3/h)

    tech = "seawater desalination"
    costs.loc[(tech, 'investment'), 'value'] = c
    costs.loc[(tech, 'investment'), 'unit'] = "EUR/(m^3-H2O/h)"
    costs.loc[(tech, 'investment'), 'source'] = source_dict['Caldera2017'] + ", Table 4."

    costs.loc[(tech, 'FOM'), 'value'] = 4.
    costs.loc[(tech, 'FOM'), 'unit'] = "%/year"
    costs.loc[(tech, 'FOM'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    costs.loc[(tech, 'lifetime'), 'value'] = 30
    costs.loc[(tech, 'lifetime'), 'unit'] = "years"
    costs.loc[(tech, 'lifetime'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    salinity = snakemake.config['desalination']['salinity']
    costs.loc[(tech, 'electricity-input'), 'value'] = (0.0003*salinity**2+0.0018*salinity+2.6043)
    costs.loc[(tech, 'electricity-input'), 'unit'] = "kWh/m^3-H2O"
    costs.loc[(tech, 'electricity-input'), 'source'] = source_dict['Caldera2016'] + ", Fig. 4."

    tech = "clean water tank storage"
    costs.loc[(tech, 'investment'), 'value'] = 65
    costs.loc[(tech, 'investment'), 'unit'] = "EUR/m^3-H2O"
    costs.loc[(tech, 'investment'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    costs.loc[(tech, 'FOM'), 'value'] = 2
    costs.loc[(tech, 'FOM'), 'unit'] = "%/year"
    costs.loc[(tech, 'FOM'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    costs.loc[(tech, 'lifetime'), 'value'] = 30
    costs.loc[(tech, 'lifetime'), 'unit'] = "years"
    costs.loc[(tech, 'lifetime'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    costs = adjust_for_inflation(costs, ['seawater desalination'], 2015)
    costs = adjust_for_inflation(costs, ['clean water tank storage'], 2013)

    return costs

def add_conventional_data(costs):
    """"
    add technology data for conventional carriers from Lazards, DIW and BP
    """
    # nuclear from Lazards
    costs.loc[('nuclear', 'investment'), 'value'] = 4000# / \
        #(1 + snakemake.config['rate_inflation'])**(2019 - snakemake.config['eur_year'])
    costs.loc[('nuclear', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('nuclear', 'investment'), 'source'] = "Assumption based on 10.1016/j.energy.2020.117015"

    costs.loc[('nuclear', 'FOM'), 'value'] = 1.4
    costs.loc[('nuclear', 'FOM'), 'unit'] = "%/year"
    costs.loc[('nuclear', 'FOM'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear', 'VOM'), 'value'] = 3.5
    costs.loc[('nuclear', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('nuclear', 'VOM'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear', 'efficiency'), 'value'] = 0.33
    costs.loc[('nuclear', 'efficiency'), 'unit'] = "per unit"
    costs.loc[('nuclear', 'efficiency'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear', 'fuel'), 'value'] = 2.6
    costs.loc[('nuclear', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('nuclear', 'fuel'), 'source'] = source_dict['Lazards']
    costs.loc[('uranium', 'fuel'), 'value'] = 2.6
    costs.loc[('uranium', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('uranium', 'fuel'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear', 'lifetime'), 'value'] = 60
    costs.loc[('nuclear', 'lifetime'), 'unit'] = "years"
    costs.loc[('nuclear', 'lifetime'), 'source'] = source_dict['Lazards']


    costs.loc[('nuclear_new', 'investment'), 'value'] = 6000
    costs.loc[('nuclear_new', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('nuclear_new', 'investment'), 'source'] = "Assumption based on 10.1016/j.energy.2020.117015"

    costs.loc[('nuclear_new', 'FOM'), 'value'] = 1.4
    costs.loc[('nuclear_new', 'FOM'), 'unit'] = "%/year"
    costs.loc[('nuclear_new', 'FOM'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear_new', 'VOM'), 'value'] = 3.5
    costs.loc[('nuclear_new', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('nuclear_new', 'VOM'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear_new', 'efficiency'), 'value'] = 0.33
    costs.loc[('nuclear_new', 'efficiency'), 'unit'] = "per unit"
    costs.loc[('nuclear_new', 'efficiency'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear_new', 'fuel'), 'value'] = 2.6
    costs.loc[('nuclear_new', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('nuclear_new', 'fuel'), 'source'] = source_dict['Lazards']

    costs.loc[('nuclear_new', 'lifetime'), 'value'] = 60
    costs.loc[('nuclear_new', 'lifetime'), 'unit'] = "years"
    costs.loc[('nuclear_new', 'lifetime'), 'source'] = "Assumption"


    # coal from Lazards and BP 2019
    costs.loc[('coal', 'investment'), 'value'] = 4162.5 / \
        (1 + snakemake.config['rate_inflation'])**(2019 - snakemake.config['eur_year'])
    costs.loc[('coal', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('coal', 'investment'), 'source'] = source_dict['Lazards']

    costs.loc[('coal', 'FOM'), 'value'] = 1.6
    costs.loc[('coal', 'FOM'), 'unit'] = "%/year"
    costs.loc[('coal', 'FOM'), 'source'] = source_dict['Lazards']

    costs.loc[('coal', 'VOM'), 'value'] = 3.5
    costs.loc[('coal', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('coal', 'VOM'), 'source'] = source_dict['Lazards']

    costs.loc[('coal', 'efficiency'), 'value'] = 0.33
    costs.loc[('coal', 'efficiency'), 'unit'] = "per unit"
    costs.loc[('coal', 'efficiency'), 'source'] = source_dict['Lazards']

    costs.loc[('coal', 'fuel'), 'value'] = 8.15
    costs.loc[('coal', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('coal', 'fuel'), 'source'] = 'BP 2019'
    costs.loc[('gas', 'fuel'), 'value'] = 20.1
    costs.loc[('gas', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('gas', 'fuel'), 'source'] = 'BP 2019'

    costs.loc[('coal', 'lifetime'), 'value'] = 40
    costs.loc[('coal', 'lifetime'), 'unit'] = "years"
    costs.loc[('coal', 'lifetime'), 'source'] = source_dict['Lazards']

    # lignite from Lazards and DIW
    costs.loc[('lignite', 'investment'), 'value'] = 4162.5 / \
        (1 + snakemake.config['rate_inflation'])**(2019 - snakemake.config['eur_year'])
    costs.loc[('lignite', 'investment'), 'unit'] = "EUR/kW_e"
    costs.loc[('lignite', 'investment'), 'source'] = source_dict['Lazards']

    costs.loc[('lignite', 'FOM'), 'value'] = 1.6
    costs.loc[('lignite', 'FOM'), 'unit'] = "%/year"
    costs.loc[('lignite', 'FOM'), 'source'] = source_dict['Lazards']

    costs.loc[('lignite', 'VOM'), 'value'] = 3.5
    costs.loc[('lignite', 'VOM'), 'unit'] = "EUR/MWh_e"
    costs.loc[('lignite', 'VOM'), 'source'] = source_dict['Lazards']

    costs.loc[('lignite', 'efficiency'), 'value'] = 0.33
    costs.loc[('lignite', 'efficiency'), 'unit'] = 'per unit'
    costs.loc[('lignite', 'efficiency'), 'source'] = source_dict['Lazards']

    costs.loc[('lignite', 'fuel'), 'value'] = 2.9
    costs.loc[('lignite', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('lignite', 'fuel'), 'source'] = 'DIW'

    costs.loc[('lignite', 'lifetime'), 'value'] = 40
    costs.loc[('lignite', 'lifetime'), 'unit']  = "years"
    costs.loc[('lignite', 'lifetime'), 'source'] = source_dict['Lazards']

    return costs


def add_co2_intensity(costs):
    """"
    add CO2 intensity for the carriers
    """
    TJ_to_MWh = 277.78
    costs.loc[('coal', 'CO2 intensity'), 'value'] = 93369 / 1e3 / TJ_to_MWh  # Steinkohle
    costs.loc[('lignite', 'CO2 intensity'), 'value'] = 113031 / 1e3 / TJ_to_MWh  # Rohbraunkohle Rheinland

    oil_specific_energy = 44 #GJ/t
    CO2_CH2_mass_ratio = 44/14 #kg/kg (1 mol per mol)
    CO2_C_mass_ratio = 44/12 #kg/kg
    methane_specific_energy = 50 #GJ/t
    CO2_CH4_mass_ratio = 44/16 #kg/kg (1 mol per mol)
    biomass_specific_energy = 18 #GJ/t LHV
    biomass_carbon_content = 0.5
    costs.loc[('oil', 'CO2 intensity'), 'value'] = (1/oil_specific_energy) * 3.6 * CO2_CH2_mass_ratio #tCO2/MWh
    costs.loc[('gas', 'CO2 intensity'), 'value'] = (1/methane_specific_energy) * 3.6 * CO2_CH4_mass_ratio #tCO2/MWh
    costs.loc[('solid biomass', 'CO2 intensity'), 'value'] = biomass_carbon_content * (1/biomass_specific_energy) * 3.6 * CO2_C_mass_ratio #tCO2/MWh

    costs.loc[('oil', 'CO2 intensity'), 'source'] = "Stoichiometric calculation with 44 GJ/t diesel and -CH2- approximation of diesel"
    costs.loc[('solid biomass', 'CO2 intensity'), 'source'] = "Stoichiometric calculation with 18 GJ/t_DM LHV and 50% C-content for solid biomass"

    costs.loc[('gas', 'CO2 intensity'), 'source'] = "Stoichiometric calculation with 50 GJ/t CH4"
    costs.loc[('coal', 'CO2 intensity'), 'source'] = source_dict["co2"]
    costs.loc[('lignite', 'CO2 intensity'), 'source'] = source_dict["co2"]
    costs.loc[pd.IndexSlice[:, "CO2 intensity"], "unit"] = "tCO2/MWh_th"

    return costs

# [add-solar-from-others]
def add_solar_from_other(costs):
    """"
    add solar from other sources than DEA (since the life time assumed in
    DEA is very optimistic)
    """

    # solar utility from Vartiaian 2019
    data = np.interp(x=years, xp=[2020, 2030, 2040, 2050],
                     fp=[431, 275, 204, 164])
    # the paper says 'In this report, all results are given in real 2019
    # money.'
    data = data / (1 + snakemake.config['rate_inflation'])**(2019 - snakemake.config['eur_year'])
    solar_uti = pd.Series(data=data, index=years)

    # solar rooftop from ETIP 2019
    data = np.interp(x=years, xp=[2020, 2030, 2050], fp=[1150, 800, 550])
    # using 2016 money in page 10
    data = data / (1 + snakemake.config['rate_inflation'])**(2016 - snakemake.config['eur_year'])
    solar_roof = pd.Series(data=data, index=years)

    # solar utility from Vartiaian 2019
    if snakemake.config['solar_utility_from_vartiaien']:
        costs.loc[('solar-utility', 'investment'), 'value'] = solar_uti[year]
        costs.loc[('solar-utility', 'investment'), 'source'] = source_dict['Vartiaien']

        costs.loc[('solar-utility', 'lifetime'), 'value'] = 30
        costs.loc[('solar-utility', 'lifetime'), 'source'] = source_dict['Vartiaien']

    if snakemake.config['solar_rooftop_from_etip']:
        # solar rooftop from ETIP 2019
        costs.loc[('solar-rooftop', 'investment'), 'value'] = solar_roof[year]
        costs.loc[('solar-rooftop', 'investment'), 'source'] = source_dict['ETIP']

        costs.loc[('solar-rooftop', 'lifetime'), 'value'] = 30
        costs.loc[('solar-rooftop', 'lifetime'), 'source'] = source_dict['ETIP']

    # lifetime&efficiency for solar
    costs.loc[('solar', 'lifetime'), 'value'] = costs.loc[(
        ['solar-rooftop', 'solar-utility'], 'lifetime'), 'value'].mean()
    costs.loc[('solar', 'lifetime'), 'unit'] = 'years'
    costs.loc[('solar', 'lifetime'),
              'source'] = 'Assuming 50% rooftop, 50% utility'
    # costs.loc[('solar', 'efficiency'), 'value'] = 1
    # costs.loc[('solar', 'efficiency'), 'unit'] = 'per unit'

    return costs

# [add-h2-from-other]
def add_h2_from_other(costs):
    """
    assume higher efficiency for electrolysis(0.8) and fuel cell(0.58)
    """
    costs.loc[('electrolysis', 'efficiency'), 'value'] = 0.8
    costs.loc[('fuel cell', 'efficiency'), 'value'] = 0.58
    costs.loc[('electrolysis', 'efficiency'), 'source'] = 'budischak2013'
    costs.loc[('fuel cell', 'efficiency'), 'source'] = 'budischak2013'

    return costs

# [unify-diw-inflation]
def unify_diw(costs):
    """"
    include inflation for the DIW costs from 2010
    """
    inflation = (1 + snakemake.config['rate_inflation'])**(2010 - snakemake.config['eur_year'])
    costs.loc[('PHS', 'investment'), 'value'] /= inflation
    costs.loc[('ror', 'investment'), 'value'] /= inflation
    costs.loc[('hydro', 'investment'), 'value'] /= inflation

    return costs


def get_data_from_DEA(data_in, expectation=None):
    """
    saves technology data from DEA in dictionary d_by_tech
    """
    d_by_tech = {}

    for tech, dea_tech in sheet_names.items():
        print(f'{tech} in PyPSA corresponds to {dea_tech} in DEA database.')
        df = get_data_DEA(tech, data_in, expectation).fillna(0)
        d_by_tech[tech] = df

    return d_by_tech

def adjust_for_inflation(costs, techs, ref_year):
    """
    adjust the investment costs for the specified techs for inflation.

    techs: str or list
        One or more techs in costs index for which the inflation adjustment is done.
    ref_year: int
        Reference year for which the costs are provided and based on which the inflation adjustment is done.
    costs: pd.Dataframe
        Dataframe containing the costs data with multiindex on technology and one index key 'investment'.
    """

    inflation = (1 + snakemake.config['rate_inflation'])**(ref_year - snakemake.config['eur_year'])
    costs.loc[(techs, 'investment'), 'value'] /= inflation

    return costs


def clean_up_units(tech_data):
    """
    converts units of a pd.Dataframe tech_data to match:
    power: Mega Watt (MW)
    energy: Mega-Watt-hour (MWh)
    currency: Euro (EUR)

    clarifies if MW_th or MW_e

    """

    tech_data.unit = tech_data.unit.str.replace(" per ", "/")
    tech_data.unit = tech_data.unit.str.replace(" / ", "/")
    tech_data.unit = tech_data.unit.str.replace(" /", "/")
    tech_data.unit = tech_data.unit.str.replace("J/s", "W")

    # units
    tech_data.loc[tech_data.unit.str.contains("MEUR"), years] *= 1e6
    tech_data.unit = tech_data.unit.str.replace("MEUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("mio EUR"), years] *= 1e6
    tech_data.unit = tech_data.unit.str.replace("mio EUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("1000EUR"), years] *= 1e3
    tech_data.unit = tech_data.unit.str.replace("1000EUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("kEUR"), years] *= 1e3
    tech_data.unit = tech_data.unit.str.replace("kEUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("/kW"), years] *= 1e3

    tech_data.loc[tech_data.unit.str.contains("kW")  & ~tech_data.unit.str.contains("/kW"), years] /= 1e3
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
    tech_data.unit = tech_data.unit.str.replace("EUR/MW of total input_e", "EUR/MW_e")
    tech_data.unit = tech_data.unit.str.replace("MWh/MWh\)", "MWh_H2/MWh_e", regex=True)
    tech_data.unit = tech_data.unit.str.replace("MWth", "MW_th")
    tech_data.unit = tech_data.unit.str.replace("MWheat", "MW_th")
    tech_data.unit = tech_data.unit.str.replace("MWhth", "MWh_th")
    tech_data.unit = tech_data.unit.str.replace("MWhheat", "MWh_th")
    tech_data.unit = tech_data.unit.str.replace("MWH Liquids", "MWh_FT")
    tech_data.unit = tech_data.unit.str.replace("MW Liquids", "MW_FT")
    tech_data.unit = tech_data.unit.str.replace("MW Methanol", "MW_MeOH")
    tech_data.unit = tech_data.unit.str.replace("MW output", "MW")
    tech_data.unit = tech_data.unit.str.replace("MW/year FT Liquids/year", "MW_FT/year")
    tech_data.unit = tech_data.unit.str.replace("MWh FT Liquids/year", "MWh_FT")
    tech_data.unit = tech_data.unit.str.replace("EUR/MWh of total input", "EUR/MWh_e")


    tech_data.unit = tech_data.unit.str.replace("FT Liquids Output, MWh/MWh Total Inpu", "MWh_FT/MWh_H2")
    tech_data.loc[tech_data.unit=='EUR/MW/y', "unit"] = 'EUR/MW/year'

    # convert per unit costs to MW
    cost_per_unit = tech_data.unit.str.contains("/unit")
    tech_data.loc[cost_per_unit, years] = tech_data.loc[cost_per_unit, years].apply(
                                                lambda x: (x / tech_data.loc[(x.name[0],
                                                                            "Heat production capacity for one unit")][years]).iloc[0,:],
                                              axis=1)
    tech_data.loc[cost_per_unit, "unit"] = tech_data.loc[cost_per_unit,
                                                         "unit"].str.replace("/unit", "/MW_th")

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

    tech_data.loc[('methanolisation', 'Variable O&M'), "unit"] = "EUR/MWh_MeOH"

    return tech_data


def set_specify_assumptions(tech_data):
    """
    for following technologies more specific investment and efficiency
    assumptions are taken:

        - central resistive heater (investment costs for large > 10 MW
                                    generators are assumed)
        - decentral gas boiler (grid connection costs)
        - biogas upgrading (include grid connection costs)
        - heat pumps (efficiencies for radiators assumed)

    to avoid duplicates some investment + efficiency data is dropped for:

        - decentral gas boilers (drop duplicated efficiency)
        - PV module (drop efficiency)

    """

    # for central resistive heater there are investment costs for small (1-5MW)
    # and large (>10 MW) generators, assume the costs for large generators
    to_drop = [("central resistive heater", 'Nominal investment, 400/690 V; 1-5 MW')]

    # for decentral gas boilers total and heat efficiency given, the values are
    # the same, drop one of the rows to avoid duplicates
    to_drop.append(("decentral gas boiler", "Heat efficiency, annual average, net"))

    # for decentral gas boilers there are investment costs and possible
    # additional investments which apply for grid connection if the house is
    # not connected yet those costs are added as an extra row since the
    # lifetime of the branchpipe is assumed to be  50 years (see comment K in
    # excel sheet)
    boiler_connect = tech_data.loc[[("decentral gas boiler",
                                     "Possible additional specific investment"),
                                    ("decentral gas boiler",
                                     "Technical lifetime")]]
    boiler_connect.loc[("decentral gas boiler", "Technical lifetime"), years] = 50
    boiler_connect.rename(index={"decentral gas boiler":
                                 "decentral gas boiler connection"}, inplace=True)
    tech_data = pd.concat([tech_data, boiler_connect])
    to_drop.append(("decentral gas boiler", "Possible additional specific investment"))

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

    # drop PV module conversion efficiency
    tech_data = tech_data.drop("PV module conversion efficiency", level=1)

    # heat pump efficiencies are assumed the one's for existing building,
    # in the DEA they do differ between heating the floor area or heating with
    # radiators, since most households heat with radiators and there
    # efficiencies are lower (conservative approach) those are assumed
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

    tech_data = tech_data.drop(to_drop)

    return tech_data.sort_index()


def set_round_trip_efficiency(tech_data):
    """
    get round trip efficiency for hydrogen and battery storage
    assume for battery sqrt(DC efficiency) and split into inverter + storage
    rename investment rows for easier sorting
    """

    # hydrogen storage
    to_drop = [("hydrogen storage tank incl. compressor", ' - Charge efficiency')]
    to_drop.append(("hydrogen storage tank incl. compressor", ' - Discharge efficiency'))
    to_drop.append(("hydrogen storage underground", ' - Charge efficiency'))
    to_drop.append(("hydrogen storage underground", ' - Discharge efficiency'))
    tech_data.loc[("hydrogen storage underground", "Round trip efficiency"), years] *= 100
    tech_data.loc[("hydrogen storage tank incl. compressor", "Round trip efficiency"), years] *= 100



    # battery split into inverter and storage, assume for efficiency sqr(round trip DC)
    df = tech_data.loc["battery"]
    inverter = df.loc[['Round trip efficiency DC',
                       'Output capacity expansion cost',
                       'Technical lifetime', 'Fixed O&M']]

    inverter.rename(index ={'Output capacity expansion cost':
                            'Output capacity expansion cost investment'},
                    inplace=True)

    # Manual correction based on footnote.
    inverter.loc['Technical lifetime', years] = 10.
    inverter.loc['Technical lifetime', 'source'] += ', Note K.'

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

    return tech_data.sort_index()


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
        investment = df[(df.index.str.contains("investment") |
                         df.index.str.contains("Distribution network costs"))
                     & ((df.unit=="EUR/MW")|
                        (df.unit=="EUR/MW_e")|
                        (df.unit=="EUR/MW_th - heat output")|
                        (df.unit=="EUR/MW_th excluding drive energy")|
                        (df.unit=="EUR/MW_th") |
                        (df.unit=="EUR/MW_MeOH") |
                        (df.unit=="EUR/MW_FT/year") |
                        (df.unit=="EUR/MWhCapacity") |
                        (df.unit=="EUR/MWh") |
                        (df.unit=="EUR/MWh/year") |
                        (df.unit=="EUR/MW input"))].copy()
        if len(investment)!=1:
            switch = True
            print("check investment: ", tech, " ",
                  df[df.index.str.contains("investment")].unit)
        else:
            investment["parameter"] = "investment"
            clean_df[tech] = investment

        # ---- FOM ----------------
        if len(investment):
            fixed = df[(df.index.str.contains("Fixed O&M") |
                        df.index.str.contains("Total O&M")) &
                       ((df.unit==investment.unit[0]+"/year")|
                        (df.unit=="EUR/MW/km/year")|
                        (df.unit=="EUR/MW/year")|
                        (df.unit=="EUR/MW_FT/year")|
                        (df.unit=='% of specific investment/year')|
                        (df.unit==investment.unit.str.split(" ")[0][0]+"/year"))].copy()
            if (len(fixed)!=1) and (len(df[df.index.str.contains("Fixed O&M")])!=0):
                switch = True
                print("check FOM: ", tech, " ",
                      df[df.index.str.contains("Fixed O&M")].unit)
            if len(fixed) == 1:
                fixed["parameter"] = "fixed"
                clean_df[tech] = pd.concat([clean_df[tech], fixed])
                fom = pd.DataFrame(columns=fixed.columns)
                if not any(fixed.unit.str.contains('% of specific investment/year')):
                    fom[years] = fixed[years]/investment[years].values*100
                else:
                    fom[years] = fixed[years]
                fom["parameter"] = "FOM"
                fom["unit"] = "%/year"
                fom["source"] = fixed["source"]
                clean_df[tech] = pd.concat([clean_df[tech], fom])

        # ---- VOM -----
        vom = df[df.index.str.contains("Variable O&M") & ((df.unit=="EUR/MWh") |
                                                         (df.unit=="EUR/MWh_e") |
                                                         (df.unit=="EUR/MWh_th") |
                                                         (df.unit=="EUR/MWh_FT") |
                                                         (df.unit=="EUR/MWh_MeOH") |
                                                         (df.unit=="EUR/MWh/year") |
                                                         (df.unit=="EUR/MWh/km") |
                                                         (df.unit=="EUR/MWh") |
                                                         (df.unit=="EUR/MWhoutput") |
                                                         (tech == "biogas upgrading"))].copy()
        if len(vom)==1:
            vom.loc[:,"parameter"] = "VOM"
            clean_df[tech] = pd.concat([clean_df[tech], vom])

        elif len(vom)!=1 and len(df[df.index.str.contains("Variable O&M")])!=0:
            switch = True
            print("check VOM: ", tech, " ",
                  df[df.index.str.contains("Variable O&M")].unit)

        # ----- lifetime --------
        lifetime = df[df.index.str.contains("Technical life") & (df.unit=="years")].copy()
        if len(lifetime)!=1:
            switch  = True
            print("check lifetime: ", tech, " ",
                  df[df.index.str.contains("Technical life")].unit)
        else:
            lifetime["parameter"] = "lifetime"
            clean_df[tech] = pd.concat([clean_df[tech], lifetime])


        # ----- efficiencies ------
        efficiency = df[(df.index.str.contains("efficiency") |
                         (df.index.str.contains("Hydrogen output, at LHV"))|
                         (df.index.str.contains("FT Liquids Output, MWh/MWh Total Input"))|
                         (df.index.str.contains("Bio SNG"))|
                         (df.index == ("Hydrogen")))
                         & ((df.unit=="%") |  (df.unit =="% total size") |  (df.unit =="% of fuel input") |  (df.unit =="MWh_H2/MWh_e") | df.unit.str.contains("MWh_FT/MWh_H2"))
                         & (~df.index.str.contains("name plate"))].copy()

        if tech == 'Fischer-Tropsch':
            efficiency[years] *= 100

        # take annual average instead of name plate efficiency
        if any(efficiency.index.str.contains("annual average")):
            efficiency = efficiency[efficiency.index.str.contains("annual average")]

        # check if electric and heat efficiencies are given
        if (any(["Electric" in ind for ind in efficiency.index]) and
            any(["Heat" in ind for ind in efficiency.index])):
            efficiency_heat = efficiency[efficiency.index.str.contains("Heat")].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
            efficiency = efficiency[efficiency.index.str.contains("Electric")].copy()
            efficiency["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency])

        elif len(efficiency)!=1:
            switch  = True
            if not any(efficiency.index.str.contains("Round trip")):
                print("check efficiency: ", tech, " ",
                       df[df.index.str.contains("efficiency")].unit)
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
    data = (pd.concat(clean_df).reset_index().rename(columns={"level_0":"technology",
                                                          "level_1": "further description"})
        .set_index(["technology", "parameter"]))

    # add water tank charger/ discharger
    charger = tech_data.loc[("central water tank storage", "Round trip efficiency")].copy()
    charger["further description"] = "efficiency from sqr(Round trip efficiency)"
    charger[years] = charger[years]**0.5*10
    charger.rename(index={"Round trip efficiency": "efficiency"},
                   level=1, inplace=True)
    charger.rename(index={'central water tank storage':"water tank charger"},
                   level=0, inplace=True)
    data = pd.concat([data, charger], sort=True)
    charger.rename(index={"water tank charger": "water tank discharger"},
                   level=0, inplace=True)
    data = pd.concat([data, charger], sort=True)

    return data


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
    if snakemake.config['offwind_no_gridcosts']:
        data.loc[("offwind", "investment"),
                 "further description"] += " grid connection costs substracted from investment costs"

    return data


def convert_units(data):
    """
    convert investment and efficiency units to be align with old pypsa
    assumptions
    """
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

    return data


def add_gas_storage(data):
    """
    add gas storage tech data, different methodolgy than other sheets and
    therefore added later
    """

    gas_storage = pd.read_excel(snakemake.input.dea_storage,
                                sheet_name="150 Underground Storage of Gas",
                                index_col=1)
    gas_storage.dropna(axis=1, how="all", inplace=True)

    # establishment of one cavern ~ 100*1e6 Nm3 = 1.1 TWh
    investment = gas_storage.loc['Total cost, 100 mio Nm3 active volume'][0]
    # convert million EUR/1.1 TWh -> EUR/kWh
    investment /= (1.1 * 1e3)
    data.loc[("gas storage", "investment"), years] = investment
    data.loc[("gas storage", "investment"), "source"] = source_dict["DEA"]
    data.loc[("gas storage", "investment"), "further description"] = "150 Underground Storage of Gas, Establishment of one cavern (units converted)"
    data.loc[("gas storage", "investment"), "unit"] = "EUR/kWh"
    data.loc[("gas storage", "lifetime"), years] = 100
    data.loc[("gas storage", "lifetime"), "source"] = "TODO no source"
    data.loc[("gas storage", "lifetime"), "further description"] = "estimation: most underground storage are already build, they do have a long lifetime"
    data.loc[("gas storage", "lifetime"), "unit"] = "years"

    # process equipment, injection (2200MW) withdrawl (6600MW)
    # assuming half of investment costs for injection, half for withdrawl
    investment_charge = gas_storage.loc["Total investment cost"].iloc[0,0]/2/2200*1e3
    investment_discharge = gas_storage.loc["Total investment cost"].iloc[0,0]/2/6600*1e3
    data.loc[("gas storage charger", "investment"), years] = investment_charge
    data.loc[("gas storage discharger", "investment"), years] = investment_discharge
    data.loc[("gas storage charger", "investment"), "source"] = source_dict["DEA"]
    data.loc[("gas storage charger", "investment"), "further description"] = "150 Underground Storage of Gas, Process equipment (units converted)"
    data.loc[("gas storage charger", "investment"), "unit"] = "EUR/kW"
    data.loc[("gas storage discharger", "investment"), "source"] = source_dict["DEA"]
    data.loc[("gas storage discharger", "investment"), "further description"] = "150 Underground Storage of Gas, Process equipment (units converted)"
    data.loc[("gas storage discharger", "investment"), "unit"] = "EUR/kW"

    # operation + maintenance 400-500 million m³ = 4.4-5.5 TWh
    FOM = gas_storage.loc["Total, incl. administration"].iloc[0] /(5.5*investment*1e3)*100
    data.loc[("gas storage", "FOM"), years] = FOM
    data.loc[("gas storage", "FOM"), "source"] = source_dict["DEA"]
    data.loc[("gas storage", "FOM"), "further description"] = "150 Underground Storage of Gas, Operation and Maintenace, salt cavern (units converted)"
    data.loc[("gas storage", "FOM"), "unit"] = "%"

    return data

def add_carbon_capture(data, tech_data):

    for tech in ['cement capture', 'biomass CHP capture']:
        data.loc[(tech,"capture_rate"), years] = tech_data.loc[(tech,'Ax) CO2 capture rate, net'), years].values[0]/100
        data.loc[(tech,"capture_rate"), 'unit'] = 'per unit'


    for tech in ['direct air capture', 'cement capture', 'biomass CHP capture']:
        print(tech, tech_data.loc[tech].index)

        data.loc[(tech,"investment"), years] = tech_data.loc[(tech,'Specific investment'), years].values[0]*1e6
        data.loc[(tech,"investment"), 'unit'] = 'EUR/(tCO2/h)'

        data.loc[(tech,"FOM"), years] = tech_data.loc[(tech,'Fixed O&M'), years].values[0]/tech_data.loc[(tech,'Specific investment'), years].values[0]*100
        data.loc[(tech,"FOM"), 'unit'] = '%/year'

        name_list = [('C2) Eletricity input ',"electricity-input"),
                     ('C1) Heat  input ',"heat-input"),
                     ('C1) Heat out ','heat-output'),
                     ('CO₂ compression and dehydration - Electricity input',"compression-electricity-input"),
                     ('CO₂ compression and dehydration - Heat out',"compression-heat-output")]

        for dea_name, our_name in name_list:
            data.loc[(tech,our_name), years] = tech_data.loc[(tech,dea_name), years].values[0]
            data.loc[(tech,our_name), 'unit'] = 'MWh/tCO2'

        data.loc[tech,'source'] = data.loc[(tech,'lifetime'),'source']
        data.loc[tech,'further description'] = sheet_names[tech]

    return data

def rename_pypsa_old(costs_pypsa):
    """
    renames old technology names to new ones to compare
    converts units from water tanks to compare
    """

    to_drop = ['retrofitting I', 'retrofitting II']
    costs_pypsa.drop(to_drop, level=0, inplace=True)

    # rename to new names
    costs_pypsa.rename({'central CHP': 'central gas CHP'}, inplace=True)
    costs_pypsa.rename({'hydrogen storage': 'hydrogen storage tank'}, inplace=True)
    costs_pypsa.rename({'hydrogen underground storage': 'hydrogen storage underground'},
                       inplace=True)

    #convert EUR/m^3 to EUR/kWh for 40 K diff and 1.17 kWh/m^3/K
    costs_pypsa.loc[('decentral water tank storage','investment'),
                    'value'] /= 1.17*40
    costs_pypsa.loc[('decentral water tank storage','investment'),'unit'] = 'EUR/kWh'

    return costs_pypsa

def add_manual_input(data):

    df = pd.read_csv(snakemake.input['manual_input'], quotechar='"',sep=',', keep_default_na=False)

    # Inflation adjustment for investment and VOM
    mask = df[df['parameter'].isin(['investment','VOM'])].index
    df.loc[mask, 'value'] /= (1+snakemake.config['rate_inflation'])**(df.loc[mask, 'currency_year']-snakemake.config['eur_year'])

    l = []
    for tech in df['technology'].unique():
        c0 = df[df['technology'] == tech]
        for param in c0['parameter'].unique():

            c = df.query('technology == @tech and parameter == @param')

            s = pd.Series(index=snakemake.config['years'],
                      data=np.interp(snakemake.config['years'], c['year'], c['value']),
                      name=param)
            s['parameter'] = param
            s['technology'] = tech
            for col in ['unit','source','further_description']:
                s[col] = "; and\n".join(c[col].unique().astype(str))
            
            l.append(s)

    new_df = pd.DataFrame(l).set_index(['technology','parameter'])
    data = data.combine_first(new_df)
    
    return data


def rename_ISE(costs_ISE):
    """
    rename ISE costs to fit to tech data
    """
    costs_ISE.rename(index = {"Investition": "investment",
                          "Lebensdauer": "lifetime",
                          "M/O-Kosten": "FOM"},
                 columns = {"Einheit": "unit",
                            "2020": 2020,
                            "2025": 2025,
                            "2030": 2030,
                            "2035": 2035,
                            "2040": 2040,
                            "2045": 2045,
                            "2050": 2050}, inplace=True)
    costs_ISE.index.names = ["technology", "parameter"]
    costs_ISE.unit.replace({"a": "years", "% Invest": "%"}, inplace=True)
    costs_ISE["source"] = source_dict["ISE"]
    costs_ISE['further description'] = costs_ISE.reset_index()["technology"].values

    return costs_ISE

def carbon_flow(costs):
    # NB: This requires some digits of accuracy; rounding to two digits creates carbon inbalances when scaling up
    c_in_char = 0 # Carbon ending up in char: zero avoids inbalace -> assumed to be circulated back and eventually end up in one of the other output streams
    medium_out = ''
    CH4_specific_energy = 50 #GJ/t methane

    for tech in ['BtL', 'BioSNG', 'methanation', 'Fischer-Tropsch', 'biogas', 'digestible biomass to hydrogen', 'solid biomass to hydrogen']:
        inv_cost = 0
        eta = 0
        lifetime = 0
        FOM = 0
        VOM = 0
        source = 'TODO'
        co2_capture_rate = 0.98

        if tech == 'BtL':
            inv_cost = 2000
            medium_out = 'oil'
            eta = 0.45
            source = "doi:10.1016/j.enpol.2017.05.013"

        elif tech == 'BioSNG':
            medium_out = 'gas'
            lifetime = 25

        elif tech == 'biogas':
            eta = 1
            source = "Assuming input biomass is already given in biogas output"

        elif tech == 'digestible biomass to hydrogen':
            inv_cost = 2500
            eta = 0.39
            FOM = 4.25
            source = 'Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014' #source_dict('HyNOW')

        elif tech == 'solid biomass to hydrogen':
            inv_cost = 2500
            eta = 0.56
            FOM = 4.25
            source = 'Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014' #source_dict('HyNOW')

        # elif tech == 'Fischer-Tropsch':

        # elif tech == 'methanation':

        if eta > 0:
            costs.loc[(tech, 'efficiency'), 'value'] = eta
            costs.loc[(tech, 'efficiency'), 'unit'] = "per unit"
            costs.loc[(tech, 'efficiency'), 'source'] = source

        costs.loc[(tech, 'capture rate'), 'value'] = co2_capture_rate
        costs.loc[(tech, 'capture rate'), 'unit'] = "per unit"
        costs.loc[(tech, 'capture rate'), 'source'] = "Assumption based on doi:10.1016/j.biombioe.2015.01.006"

        if inv_cost > 0:
            costs.loc[(tech, 'investment'), 'value'] = inv_cost
            costs.loc[(tech, 'investment'), 'unit'] = "EUR/kW_th"
            costs.loc[(tech, 'investment'), 'source'] = source

        if lifetime > 0:
            costs.loc[(tech, 'lifetime'), 'value'] = lifetime
            costs.loc[(tech, 'lifetime'), 'unit'] = "years"
            costs.loc[(tech, 'lifetime'), 'source'] = source

        if FOM > 0:
            costs.loc[(tech, 'FOM'), 'value'] = FOM
            costs.loc[(tech, 'FOM'), 'unit'] = "%/year"
            costs.loc[(tech, 'FOM'), 'source'] = source

        if VOM > 0:
            costs.loc[(tech, 'VOM'), 'value'] = VOM
            costs.loc[(tech, 'VOM'), 'unit'] = "EUR/MWh_th"
            costs.loc[(tech, 'VOM'), 'source'] = source

        if tech in ['BioSNG', 'BtL']:
            input_CO2_intensity = costs.loc[('solid biomass', 'CO2 intensity'), 'value']

            costs.loc[(tech, 'C in fuel'), 'value'] = costs.loc[(tech, 'efficiency'), 'value'] \
                                                  * costs.loc[(medium_out, 'CO2 intensity'), 'value'] \
                                                  / input_CO2_intensity
            costs.loc[(tech, 'C stored'), 'value'] = 1 - costs.loc[(tech, 'C in fuel'), 'value'] - c_in_char
            costs.loc[(tech, 'CO2 stored'), 'value'] = input_CO2_intensity * costs.loc[(tech, 'C stored'), 'value']

            costs.loc[(tech, 'C in fuel'), 'unit'] = "per unit"
            costs.loc[(tech, 'C stored'), 'unit'] = "per unit"
            costs.loc[(tech, 'CO2 stored'), 'unit'] = "tCO2/MWh_th"

            costs.loc[(tech, 'C in fuel'), 'source'] = "Stoichiometric calculation"
            costs.loc[(tech, 'C stored'), 'source'] = "Stoichiometric calculation"
            costs.loc[(tech, 'CO2 stored'), 'source'] = "Stoichiometric calculation"

        elif tech == 'biogas':
            AD_CO2_share = 0.4 #volumetric share in biogas (rest is CH4)

            CH4_density = 0.657 #kg/Nm3
            CO2_density = 1.98 #kg/Nm3
            CH4_vol_energy_density = CH4_specific_energy * CH4_density / (1000 * 3.6) #MJ/Nm3 -> MWh/Nm3
            CO2_weight_share = AD_CO2_share * CO2_density

            costs.loc[('biogas', 'CO2 stored'), 'value'] = CO2_weight_share / CH4_vol_energy_density / 1000 #tCO2/MWh,in (NB: assuming the input is already given in the biogas potential and cost
            costs.loc[('biogas', 'CO2 stored'), 'unit'] = "tCO2/MWh_th"
            costs.loc[('biogas', 'CO2 stored'), 'source'] = "Stoichiometric calculation"

    return costs


def add_home_battery_costs(costs):
    """
    adds investment costs for home battery storage and inverter.
    Since home battery costs are not part of the DEA cataloque, utility-scale
    costs are multiplied by a factor determined by data from the EWG study
    """
    # get DEA assumptions for utility scale
    home_battery = (data.loc[["battery storage", "battery inverter"]]
                    .rename(index=lambda x: "home " + x, level=0))

    # get EWG cost assumptions
    costs_ewg = pd.read_csv(snakemake.input.EWG_costs,
                            index_col=list(range(2))).sort_index()
    v = costs_ewg.unstack()[[str(year) for year in years]].swaplevel(axis=1)

    def annuity(n,r=0.07):
        """
        Calculate the annuity factor for an asset with lifetime n years and
        discount rate of r
        """
        if isinstance(r, pd.Series):
            return pd.Series(1/n, index=r.index).where(r == 0, r/(1. - 1./(1.+r)**n))
        elif r > 0:
            return r/(1. - 1./(1.+r)**n)
        else:
            return 1/n

    # annualise EWG cost assumptions
    fixed = (annuity(v["lifetime"])+v["FOM"]/100.) * v["investment"]

    # battery storage index in EWG --------------
    battery_store_i = [
                        'Battery PV prosumer - commercial storage',
                        'Battery PV prosumer - industrial storage',
                       'Battery PV prosumer - residential storage',
                       'Battery storage']

    battery_store_ewg = fixed.loc[battery_store_i].T

    def get_factor(df, cols, utility_col):
        """get factor by which costs are increasing for home installations"""
        return (df[cols].div(df[utility_col], axis=0).mean(axis=1)
                .rename(index=lambda x: float(x)))

    # take mean of cost increase for commercial and residential storage compared to utility-scale
    home_cols = ['Battery PV prosumer - commercial storage',
                 'Battery PV prosumer - residential storage']
    factor = get_factor(battery_store_ewg, home_cols, "Battery storage")

    home_cost =  (home_battery.loc[("home battery storage", "investment"), years] * factor).values
    home_battery.loc[("home battery storage", "investment"), years] = home_cost

    # battery inverter index in EWG -----------------------
    battery_inverter_i = [
                        'Battery PV prosumer - commercial interface',
                        'Battery PV prosumer - industrial interface PHES',
                        'Battery PV prosumer - residential interface',
                        'Battery interface']

    battery_inverter_ewg = fixed.loc[battery_inverter_i].T

    home_cols = ['Battery PV prosumer - commercial interface',
                 'Battery PV prosumer - residential interface']
    factor = get_factor(battery_inverter_ewg, home_cols, "Battery interface")
    home_cost = (home_battery.loc[("home battery inverter", "investment"), years] * factor).values
    home_battery.loc[("home battery inverter", "investment"), years] = home_cost

    # adjust source
    home_battery["source"] = home_battery["source"].apply(lambda x: source_dict["EWG"] + ", " + x)

    return pd.concat([costs, home_battery])
# %% *************************************************************************
#  ---------- MAIN ------------------------------------------------------------
# for testing
if 'snakemake' not in globals():
    from vresutils.snakemake import MockSnakemake
    snakemake = MockSnakemake(
        input=dict(
                    pypsa_costs = "inputs/costs_PyPSA.csv",
                    fraunhofer_costs = "inputs/Fraunhofer_ISE_costs.csv",
                    fraunhofer_energy_prices = "inputs/Fraunhofer_ISE_energy_prices.csv",
                    EWG_costs = "inputs/EWG_costs.csv",
                    dea_transport = "inputs/energy_transport_data_sheet_dec_2017.xlsx",
                    dea_renewable_fuels = "inputs/data_sheets_for_renewable_fuels.xlsx",
                    dea_storage = "inputs/technology_data_catalogue_for_energy_storage.xlsx",
                    dea_generation = "inputs/technology_data_for_el_and_dh_-_0009.xlsx",
                    dea_heating = "inputs/technologydatafor_heating_installations_marts_2018.xlsx",
                    dea_industrial = "inputs/technology_data_for_industrial_process_heat_0002.xlsx"
        ),
        output=["outputs/costs_{}.csv".format(year) for year in [2020, 2025, 2030, 2035, 2040, 2045, 2050]]

    )
    import yaml
    with open('config.yaml', encoding='utf8') as f:
        snakemake.config = yaml.safe_load(f)

# (1) DEA data
# (a)-------- get data from DEA excel sheets ----------------------------------
# read excel sheet names of all excel files
excel_files = [v for k,v in snakemake.input.items() if "dea" in k]
data_in = get_excel_sheets(excel_files)
# create dictionary with raw data from DEA sheets
d_by_tech = get_data_from_DEA(data_in, expectation=snakemake.config["expectation"])
# concat into pd.Dataframe
tech_data = pd.concat(d_by_tech).sort_index()
# clean up units
tech_data = clean_up_units(tech_data)
# (b) ------ specific assumptions for some technologies -----------------------

# specify investment and efficiency assumptions for:
# resistive heater, decentral gas boiler, biogas upgrading and heat pumps
tech_data = set_specify_assumptions(tech_data)

# round trip efficiency for hydrogen + battery storage
tech_data = set_round_trip_efficiency(tech_data)

# drop all rows which only contains zeros
tech_data = tech_data.loc[(tech_data[years]!=0).sum(axis=1)!=0]

# (c) -----  get tech data in pypsa syntax -----------------------------------
# make categories: investment, FOM, VOM, efficiency, c_b, c_v
data = order_data(tech_data)
# add excel sheet names and further description
data = add_description(data)
# convert efficiency from %-> per unit and investment from MW->kW to compare
data = convert_units(data)
# add gas storage (different methodology than other sheets)
data = add_gas_storage(data)
# add carbon capture
data = add_carbon_capture(data, tech_data)


# %% (2) -- get data from other sources which need formatting -----------------
# (a)  ---------- get old pypsa costs ---------------------------------------
costs_pypsa = pd.read_csv(snakemake.input.pypsa_costs,
                          index_col=[0,2]).sort_index()
# rename some techs and convert units
costs_pypsa = rename_pypsa_old(costs_pypsa)

# (b) ------- add costs from Fraunhofer ISE study --------------------------
costs_ISE = pd.read_csv(snakemake.input.fraunhofer_costs,
                        engine="python",
                        index_col=[0,1],
                        encoding = "ISO-8859-1")
# rename + reorder to fit to other data
costs_ISE = rename_ISE(costs_ISE)
# add costs for gas pipelines
data = pd.concat([data, costs_ISE.loc[["Gasnetz"]]], sort=True)

data = add_manual_input(data)
# add costs for home batteries
data = add_home_battery_costs(data)
# %% (3) ------ add additional sources and save cost as csv ------------------
# [RTD-target-multiindex-df]
for year in years:
    costs = (data[[year, "unit", "source", "further description"]]
             .rename(columns={year: "value"}))
    costs["value"] = costs["value"].astype(float)

    # biomass is differentiated by biomass CHP and HOP
    costs.loc[('solid biomass', 'fuel'), 'value'] = 12
    costs.loc[('solid biomass', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('solid biomass', 'fuel'), 'source'] = "JRC ENSPRESO ca avg for MINBIOWOOW1 (secondary forest residue wood chips), ENS_Ref for 2040"

    costs.loc[('digestible biomass', 'fuel'), 'value'] = 15
    costs.loc[('digestible biomass', 'fuel'), 'unit'] = 'EUR/MWh_th'
    costs.loc[('digestible biomass', 'fuel'), 'source'] = "JRC ENSPRESO ca avg for MINBIOAGRW1, ENS_Ref for 2040"

    # add solar data from other source than DEA
    if any([snakemake.config['solar_utility_from_vartiaien'], snakemake.config['solar_rooftop_from_etip']]):
        costs = add_solar_from_other(costs)
    
    # add desalination and clean water tank storage
    costs = add_desalinsation_data(costs)

    # add electrolyzer and fuel cell efficiency from other source than DEA
    if snakemake.config['h2_from_budischak']:
        costs = add_h2_from_other(costs)

    # add data from conventional carriers
    costs = add_conventional_data(costs)
    # CO2 intensity
    costs = add_co2_intensity(costs)

    #carbon balances
    costs = carbon_flow(costs)

    # include old pypsa costs
    check = pd.concat([costs_pypsa, costs], sort=True, axis=1)

    # missing technologies
    missing = costs_pypsa.index.levels[0].difference(costs.index.levels[0])
    if (len(missing) & (year==years[0])):
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
    if (year==years[0]):
        print("single parameters of technologies are missing, using old PyPSA assumptions: ")
        print(comp_missing)
        print("old c_v and c_b values are assumed where given")
    to_add = costs_pypsa.loc[comp_missing].drop("year", axis=1)
    to_add.loc[:, "further description"] = " from old pypsa cost assumptions"
    costs_tot = pd.concat([costs_tot, to_add], sort=False)

    # unify the cost from DIW2010
    costs_tot = unify_diw(costs_tot)
    costs_tot.drop("fixed", level=1, inplace=True)
    costs_tot.sort_index(inplace=True)
    costs_tot = round(costs_tot, ndigits=snakemake.config.get("ndigits", 5))
    costs_tot.to_csv([v for v in snakemake.output if str(year) in v][0])