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
try:
    pd.set_option('future.no_silent_downcasting', True)
except Exception:
    pass
# ---------- sources -------------------------------------------------------
source_dict = {
                'DEA': 'Danish Energy Agency',
                # solar utility
                'Vartiaien': 'Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity',
                # solar rooftop
                'ETIP': 'European PV Technology and Innovation Platform',
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
                "HyNOW" : "Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014",
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
		        "vehicles" : "PATHS TO A CLIMATE-NEUTRAL ENERGY SYSTEM The German energy transformation in its social context. https://www.ise.fraunhofer.de/en/publications/studies/paths-to-a-climate-neutral-energy-system.html"
                }

# [DEA-sheet-names]
sheet_names = {'onwind': '20 Onshore turbines',
               'offwind': '21 Offshore turbines',
               'solar-utility': '22 Utility-scale PV',
               'solar-utility single-axis tracking': '22 Utility-scale PV tracker',
               'solar-rooftop residential': '22 Rooftop PV residential',
               'solar-rooftop commercial': '22 Rooftop PV commercial',
               'OCGT': '52 OCGT - Natural gas',
               'CCGT': '05 Gas turb. CC, steam extract.',
               'oil': '50 Diesel engine farm',
               'biomass CHP': '09c Straw, Large, 40 degree',
               'biomass EOP': '09c Straw, Large, 40 degree',
               'biomass HOP': '09c Straw HOP',
               'central coal CHP': '01 Coal CHP',
               'central gas CHP': '04 Gas turb. simple cycle, L',
               'central gas CHP CC': '04 Gas turb. simple cycle, L',
               'central solid biomass CHP': '09a Wood Chips, Large 50 degree',
               'central solid biomass CHP CC': '09a Wood Chips, Large 50 degree',
               'central solid biomass CHP powerboost CC': '09a Wood Chips, Large 50 degree',
               # 'solid biomass power': '09a Wood Chips extract. plant',
               # 'solid biomass power CC': '09a Wood Chips extract. plant',
               'central air-sourced heat pump': '40 Comp. hp, airsource 3 MW',
               'central ground-sourced heat pump': '40 Absorption heat pump, DH',
               'central resistive heater': '41 Electric Boilers',
               'central gas boiler': '44 Natural Gas DH Only',
               'decentral gas boiler': '202 Natural gas boiler',
               'direct firing gas': '312.a Direct firing Natural Gas',
               'direct firing gas CC': '312.a Direct firing Natural Gas',
               'direct firing solid fuels': '312.b Direct firing Sold Fuels',
               'direct firing solid fuels CC': '312.b Direct firing Sold Fuels',
               'decentral ground-sourced heat pump': '207.7 Ground source existing',
               'decentral air-sourced heat pump': '207.3 Air to water existing',
               # 'decentral resistive heater': '216 Electric heating',
               'central water tank storage': '140 PTES seasonal',
               # 'decentral water tank storage': '142 Small scale hot water tank',
               'fuel cell': '12 LT-PEMFC CHP',
               'hydrogen storage underground': '151c Hydrogen Storage - Caverns',
               'hydrogen storage tank type 1 including compressor': '151a Hydrogen Storage - Tanks',
               'micro CHP': '219 LT-PEMFC mCHP - natural gas',
               'biogas' : '81 Biogas, Basic plant, small',
               'biogas CC' : '81 Biogas, Basic plant, small',
               'biogas upgrading': '82 Upgrading 3,000 Nm3 per h',
               'battery': '180 Lithium Ion Battery',
               'industrial heat pump medium temperature': '302.a High temp. hp Up to 125 C',
               'industrial heat pump high temperature': '302.b High temp. hp Up to 150',
               'electric boiler steam': '310.1 Electric boiler steam  ',
               'gas boiler steam': '311.1c Steam boiler Gas',
               'solid biomass boiler steam': '311.1e Steam boiler Wood',
               'solid biomass boiler steam CC': '311.1e Steam boiler Wood',
               'biomass boiler': '204 Biomass boiler, automatic',
               'electrolysis': '86 AEC 100 MW', 
               'direct air capture': '403.a Direct air capture',
               'biomass CHP capture': '401.a Post comb - small CHP',
               'cement capture': '401.c Post comb - Cement kiln',
               'BioSNG': '84 Gasif. CFB, Bio-SNG',
               'BtL': '85 Gasif. Ent. Flow FT, liq fu ',
               'biomass-to-methanol': '97 Methanol from biomass gasif.',
               'biogas plus hydrogen': '99 SNG from methan. of biogas',
               'methanolisation': '98 Methanol from hydrogen',
               'Fischer-Tropsch': '102 Hydrogen to Jet',
               'central hydrogen CHP': '12 LT-PEMFC CHP',
               'Haber-Bosch': '103 Hydrogen to Ammonia',
               'air separation unit': '103 Hydrogen to Ammonia',
               'waste CHP': '08 WtE CHP, Large, 50 degree',
               'waste CHP CC': '08 WtE CHP, Large, 50 degree',
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
                    'solar-utility single-axis tracking': 'J:K',
                    'solar-rooftop residential':  'J:K',
                    'solar-rooftop commercial':  'J:K',
                    'OCGT': 'I:J',
                    'CCGT': 'I:J',
                    'oil': 'I:J',
                    'biomass CHP': 'I:J',
                    'biomass EOP': 'I:J',
                    'biomass HOP': 'I:J',
                    'central coal CHP': '',
                    'central gas CHP': 'I:J',
                    'central gas CHP CC': 'I:J',
                    'central hydrogen CHP': 'I:J',
                    'central solid biomass CHP': 'I:J',
                    'central solid biomass CHP CC': 'I:J',
                    'central solid biomass CHP powerboost CC': 'I:J',
                    # 'solid biomass power': 'J:K',
                    # 'solid biomass power CC': 'J:K',
                    'solar': '',
                    'central air-sourced heat pump': 'J:K',
                    'central ground-sourced heat pump': 'I:J',
                    'central resistive heater': 'I:J',
                    'central gas boiler': 'I:J',
                    'decentral gas boiler': 'I:J',
                    'direct firing gas': 'H:I',
                    'direct firing gas CC': 'H:I',
                    'direct firing solid fuels': 'H:I',
                    'direct firing solid fuels CC': 'H:I',
                    'decentral ground-sourced heat pump': 'I:J',
                    'decentral air-sourced heat pump': 'I:J',
                    'central water tank storage': 'J:K',
                    'fuel cell': 'I:J',
                    'hydrogen storage underground': 'J:K',
                    'hydrogen storage tank type 1 including compressor': 'J:K',
                    'micro CHP': 'I:J',
                    'biogas': 'I:J',
                    'biogas CC': 'I:J',
                    'biogas upgrading': 'I:J',
                    'electrolysis': 'I:J',
                    'battery': 'L,N',
                    'direct air capture': 'I:J',
                    'cement capture': 'I:J',
                    'biomass CHP capture': 'I:J',
                    'BioSNG' : 'I:J',
                    'BtL' : 'J:K',
                    'biomass-to-methanol' : 'J:K',
                    'biogas plus hydrogen' : 'J:K',
                    'industrial heat pump medium temperature':'H:I',
                    'industrial heat pump high temperature':'H:I',
                    'electric boiler steam':'H:I',
                    'gas boiler steam':'H:I',
                    'solid biomass boiler steam':'H:I',
                    'solid biomass boiler steam CC':'H:I',
                    'biomass boiler': 'I:J',
                    'Fischer-Tropsch': 'I:J',
                    'Haber-Bosch': 'I:J',
                    'air separation unit': 'I:J',
                    'methanolisation': 'J:K',
                    'waste CHP': 'I:J',
                    'waste CHP CC': 'I:J',
}

# since February 2022 DEA uses a new format for the technology data
# all excel sheets of updated technologies have a different layout and are
# given in EUR_2020 money (instead of EUR_2015)
cost_year_2020 = ['solar-utility',
              'solar-utility single-axis tracking',
              'solar-rooftop residential',
              'solar-rooftop commercial',
              'offwind',
              'electrolysis',
              'biogas',
              'biogas CC',
              'biogas upgrading',
              'direct air capture',
              'biomass CHP capture',
              'cement capture',
              'BioSNG',
              'BtL',
              'biomass-to-methanol',
              'biogas plus hydrogen',
              'methanolisation',
              'Fischer-Tropsch'
              ]

cost_year_2019 = ['direct firing gas',
                'direct firing gas CC',
                'direct firing solid fuels',
                'direct firing solid fuels CC',
                'industrial heat pump medium temperature',
                'industrial heat pump high temperature',
                'electric boiler steam',
                'gas boiler steam',
                'solid biomass boiler steam',
                'solid biomass boiler steam CC',
                ]


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
    elif tech in ['industrial heat pump medium temperature', 'industrial heat pump high temperature',
                  'electric boiler steam', "gas boiler steam", "solid biomass boiler steam", "solid biomass boiler steam CC", "direct firing gas", "direct firing gas CC", "direct firing solid fuels", "direct firing solid fuels CC"]:
        usecols = "A:E"
    elif tech in ['Fischer-Tropsch', 'Haber-Bosch', 'air separation unit']:
        usecols = "B:F"
    else:
        usecols = "B:G"

    usecols += f",{uncrtnty_lookup[tech]}"


    if ((tech in cost_year_2019) or (tech in cost_year_2020) or ("renewable_fuels" in excel_file)):
        skiprows = [0]
    else:
        skiprows = [0,1]

    excel = pd.read_excel(excel_file,
                          sheet_name=sheet_names[tech],
                          index_col=0,
                          usecols=usecols,
                          skiprows=skiprows,
                          na_values="N.A")
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
        if tech in ["hydrogen storage tank type 1 including compressor", "hydrogen storage cavern"]:
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
                  "Hydrogen Output",
                  "Hydrogen (% total input_e (MWh / MWh))",
                  "Hydrogen [% total input_e",
                  " - hereof recoverable for district heating (%-points of heat loss)",
                  "Cb coefficient",
                  "Cv coefficient",
                  "Distribution network costs", "Technical life",
                  "Energy storage expansion cost",
                  'Output capacity expansion cost (M€2015/MW)',
                  'Heat input', 'Heat  input', 'Electricity input', 'Eletricity input', 'Heat out',
                  'capture rate',
                  "FT Liquids Output, MWh/MWh Total Input",
                  " - hereof recoverable for district heating [%-points of heat loss]",
                  " - hereof recoverable for district heating (%-points of heat loss)",
                  "Bio SNG Output [% of fuel input]", 
                  "Methanol Output", 
                  "District heat  Output",
                  "Electricity Output",
                  "Total O&M"]


    df = pd.DataFrame()
    for para in parameters:
        # attr = excel[excel.index.str.contains(para)]
        attr = excel[[para in index for index in excel.index]]
        if len(attr) != 0:
            df = pd.concat([df, attr])
    df.index = df.index.str.replace('€', 'EUR')

    df = df.reindex(columns=df.columns[df.columns.isin(years)])
    df = df[~df.index.duplicated(keep='first')]

    # replace missing data
    df.replace("-", np.nan, inplace=True)
    # average data  in format "lower_value-upper_value"
    df = df.apply(lambda row: row.apply(lambda x: (float(x.split("-")[0])
                                                   + float(x.split("-")[1]))
                                        / 2 if isinstance(x, str) and "-" in x else x),
                  axis=1)

    # remove symbols "~", ">", "<" and " "
    for sym in ["~", ">", "<", " "]:
        df = df.apply(lambda col: col.apply(lambda x: x.replace(sym, "")
                                            if isinstance(x, str) else x))


    df = df.astype(float)
    df = df.mask(df.apply(pd.to_numeric, errors='coerce').isnull(), df.astype(str).apply(lambda x: x.str.strip()))
    # print(df)

    ## Modify data loaded from DEA on a per-technology case
    if (tech == "offwind") and snakemake.config['offwind_no_gridcosts']:
        df.loc['Nominal investment (*total) [MEUR/MW_e, 2020]'] -= excel.loc['Nominal investment (installation: grid connection) [M€/MW_e, 2020]']

    # Exlucde indirect costs for centralised system with additional piping.
    if tech.startswith('industrial heat pump'):
        df = df.drop('Indirect investments cost (MEUR per MW)')

    if tech == 'biogas plus hydrogen':
        df.drop(df.loc[df.index.str.contains("GJ SNG")].index, inplace=True)

    if tech == 'BtL':
        df.drop(df.loc[df.index.str.contains("1,000 t FT Liquids")].index, inplace=True)

    if tech == "biomass-to-methanol":
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech == 'methanolisation':
        df.drop(df.loc[df.index.str.contains("1,000 t Methanol")].index, inplace=True)

    if tech == 'Fischer-Tropsch':
        df.drop(df.loc[df.index.str.contains("l FT Liquids")].index, inplace=True)

    if tech == 'biomass boiler':
        df.drop(df.loc[df.index.str.contains("Possible additional")].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Total efficiency")].index, inplace=True)

    if tech == "Haber-Bosch":
        df.drop(df.loc[df.index.str.contains("Specific investment mark-up factor optional ASU")].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Specific investment (MEUR /TPD Ammonia output", regex=False)].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Fixed O&M (MEUR /TPD Ammonia", regex=False)].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Variable O&M (EUR /t Ammonia)", regex=False)].index, inplace=True)

    if tech == "air separation unit":
        divisor = ((df.loc["Specific investment mark-up factor optional ASU"] - 1.0)
                    / excel.loc["N2 Consumption, [t/t] Ammonia"]).astype(float)
        
        # Calculate ASU cost separate to HB facility in terms of t N2 output
        df.loc[[
            "Specific investment [MEUR /TPD Ammonia output]",
            "Fixed O&M [kEUR /TPD Ammonia]",
            "Variable O&M [EUR /t Ammonia]"
            ]] *= divisor
        # Convert output to hourly generation
        df.loc[[
            "Specific investment [MEUR /TPD Ammonia output]",
            "Fixed O&M [kEUR /TPD Ammonia]",
            ]] *= 24

        # Rename costs for correct units
        df.index = df.index.str.replace("MEUR /TPD Ammonia output", "MEUR/t_N2/h")
        df.index = df.index.str.replace("kEUR /TPD Ammonia", "kEUR/t_N2/h/year")
        df.index = df.index.str.replace("EUR /t Ammonia", "EUR/t_N2")

        df.drop(df.loc[df.index.str.contains("Specific investment mark-up factor optional ASU")].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Specific investment [MEUR /MW Ammonia output]", regex=False)].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Fixed O&M [kEUR/MW Ammonia/year]", regex=False)].index, inplace=True)
        df.drop(df.loc[df.index.str.contains("Variable O&M [EUR/MWh Ammonia]", regex=False)].index, inplace=True)
        
    if "solid biomass power" in tech:
        df.index = df.index.str.replace("EUR/MWeh", "EUR/MWh")

    df_final = pd.DataFrame(index=df.index, columns=years)

    # [RTD-interpolation-example]
    for index in df_final.index:
        values = np.interp(x=years, xp=df.columns.values.astype(float), fp=df.loc[index, :].values.astype(float))
        df_final.loc[index, :] = values

    # if year-specific data is missing and not fixed by interpolation fill forward with same values
    df_final = df_final.ffill(axis=1)

    df_final["source"] = source_dict["DEA"] + ", " + excel_file.replace("inputs/","")
    if tech in cost_year_2020 and (not ("for_carbon_capture_transport_storage" in excel_file)) and (not ("renewable_fuels" in excel_file)):
        for attr in ["investment", "Fixed O&M"]:
            to_drop = df[df.index.str.contains(attr) &
                         ~df.index.str.contains("\(\*total\)")].index
            df_final.drop(to_drop, inplace=True)

        df_final["unit"] = (df_final.rename(index=lambda x:
                                            x[x.rfind("[")+1: x.rfind("]")]).index.values)
    else:
        df_final.index = df_final.index.str.replace("\[", "(", regex=True).str.replace("\]", ")", regex=True)
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
    costs.loc[(tech, 'investment'), 'currency_year'] = 2015
    
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
    costs.loc[(tech, 'investment'), 'currency_year'] = 2013

    costs.loc[(tech, 'FOM'), 'value'] = 2
    costs.loc[(tech, 'FOM'), 'unit'] = "%/year"
    costs.loc[(tech, 'FOM'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

    costs.loc[(tech, 'lifetime'), 'value'] = 30
    costs.loc[(tech, 'lifetime'), 'unit'] = "years"
    costs.loc[(tech, 'lifetime'), 'source'] = source_dict['Caldera2016'] + ", Table 1."

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
    costs.loc[('methanol', 'CO2 intensity'), 'value'] = 0.2482 # t_CO2/MWh_th, based on stochiometric composition.
    costs.loc[('solid biomass', 'CO2 intensity'), 'value'] = 0.3

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
    costs.loc[('gas', 'CO2 intensity'), 'source'] = "Stoichiometric calculation with 50 GJ/t CH4"
    costs.loc[('solid biomass', 'CO2 intensity'), 'source'] = "Stoichiometric calculation with 18 GJ/t_DM LHV and 50% C-content for solid biomass"
    costs.loc[('coal', 'CO2 intensity'), 'source'] = source_dict["co2"]
    costs.loc[('lignite', 'CO2 intensity'), 'source'] = source_dict["co2"]


    costs.loc[pd.IndexSlice[:, "CO2 intensity"], "unit"] = "tCO2/MWh_th"

    return costs


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
        costs.loc[('solar-utility', 'investment'), 'currency_year'] = 2019

        costs.loc[('solar-utility', 'lifetime'), 'value'] = 30
        costs.loc[('solar-utility', 'lifetime'), 'source'] = source_dict['Vartiaien']
        costs.loc[('solar-utility', 'lifetime'), 'currency_year'] = 2019

    if snakemake.config['solar_rooftop_from_etip']:
        # solar rooftop from ETIP 2019
        costs.loc[('solar-rooftop', 'investment'), 'value'] = solar_roof[year]
        costs.loc[('solar-rooftop', 'investment'), 'source'] = source_dict['ETIP']
        costs.loc[('solar-rooftop', 'investment'), 'currency_year'] = 2019

        costs.loc[('solar-rooftop', 'lifetime'), 'value'] = 30
        costs.loc[('solar-rooftop', 'lifetime'), 'source'] = source_dict['ETIP']
        costs.loc[('solar-rooftop', 'lifetime'), 'currency_year'] = 2019

    # lifetime&efficiency for solar
    costs.loc[('solar', 'lifetime'), 'value'] = costs.loc[(
        ['solar-rooftop', 'solar-utility'], 'lifetime'), 'value'].mean()
    costs.loc[('solar', 'lifetime'), 'unit'] = 'years'
    costs.loc[('solar', 'lifetime'), 'currency_year'] = 2019
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
    costs.loc[('electrolysis', 'efficiency'), 'currency_year'] =  2013
    costs.loc[('fuel cell', 'efficiency'), 'source'] = 'budischak2013'
    costs.loc[('fuel cell', 'efficiency'), 'currency_year'] =  2013

    return costs

# [unify-diw-inflation]
def unify_diw(costs):
    """"
    add currency year for the DIW costs from 2010
    """

    costs.loc[('PHS', 'investment'), 'currency_year'] = 2010
    costs.loc[('ror', 'investment'), 'currency_year'] = 2010
    costs.loc[('hydro', 'investment'), 'currency_year'] = 2010

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

def adjust_for_inflation(inflation_rate, costs, techs, ref_year, col):
    """
    adjust the investment costs for the specified techs for inflation.

    techs: str or list
        One or more techs in costs index for which the inflation adjustment is done.
    ref_year: int
        Reference year for which the costs are provided and based on which the inflation adjustment is done.
    costs: pd.Dataframe
        Dataframe containing the costs data with multiindex on technology and one index key 'investment'.
    """
    
    def get_factor(inflation_rate, ref_year, eur_year):
        if (pd.isna(ref_year)) or (ref_year<1900): return np.nan
        if ref_year == eur_year: return 1
        mean = inflation_rate.mean()
        if ref_year< eur_year:
            new_index = np.arange(ref_year+1, eur_year+1)
            df = 1 + inflation_rate.reindex(new_index).fillna(mean)    
            return df.cumprod().loc[eur_year]
        else:
            new_index = np.arange(eur_year+1, ref_year+1)
            df = 1 + inflation_rate.reindex(new_index).fillna(mean)
            return 1/df.cumprod().loc[ref_year]
    
    inflation = costs.currency_year.apply(lambda x: get_factor(inflation_rate, x, snakemake.config['eur_year']))

    paras = ["investment", "VOM", "fuel"]
    filter_i = costs.index.get_level_values(0).isin(techs) & costs.index.get_level_values(1).isin(paras) 
    costs.loc[filter_i, col] = costs.loc[filter_i, col].mul(inflation.loc[filter_i], axis=0)


    return costs


def clean_up_units(tech_data, value_column="", source=""):
    """
    converts units of a pd.Dataframe tech_data to match:
    power: Mega Watt (MW)
    energy: Mega-Watt-hour (MWh)
    currency: Euro (EUR)

    clarifies if MW_th or MW_e
    """
    from currency_converter import CurrencyConverter
    from datetime import date
    from currency_converter import ECB_URL

    # Currency conversion
    REPLACEMENTS = [
        ('€', 'EUR'),
        ('$', 'USD'),
        ('₤', 'GBP'),
    ]
    # Download the full history, this will be up to date. Current value is:
    # https://www.ecb.europa.eu/stats/eurofxref/eurofxref-hist.zip
    c = CurrencyConverter(ECB_URL)
    c = CurrencyConverter(fallback_on_missing_rate=True)

    for old, new in REPLACEMENTS:
        tech_data.unit = tech_data.unit.str.replace(old, new, regex=False)
        tech_data.loc[tech_data.unit.str.contains(new), value_column] *= c.convert(1, new, "EUR", date=date(2020, 1, 1))
        tech_data.unit = tech_data.unit.str.replace(new, "EUR")

    tech_data.unit = tech_data.unit.str.replace(" per ", "/")
    tech_data.unit = tech_data.unit.str.replace(" / ", "/")
    tech_data.unit = tech_data.unit.str.replace(" /", "/")
    tech_data.unit = tech_data.unit.str.replace("J/s", "W")

    # units
    tech_data.loc[tech_data.unit.str.contains("MEUR"), value_column] *= 1e6
    tech_data.unit = tech_data.unit.str.replace("MEUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("mio EUR"), value_column] *= 1e6
    tech_data.unit = tech_data.unit.str.replace("mio EUR", "EUR")
    
    tech_data.loc[tech_data.unit.str.contains("mill. EUR"), value_column] *= 1e6
    tech_data.unit = tech_data.unit.str.replace("mill. EUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("1000EUR"), value_column] *= 1e3
    tech_data.unit = tech_data.unit.str.replace("1000EUR", "EUR")

    tech_data.unit = tech_data.unit.str.replace("k EUR", "kEUR")
    tech_data.loc[tech_data.unit.str.contains("kEUR"), value_column] *= 1e3
    tech_data.unit = tech_data.unit.str.replace("kEUR", "EUR")

    tech_data.loc[tech_data.unit.str.contains("/kW"), value_column] *= 1e3

    tech_data.loc[tech_data.unit.str.contains("kW")  & ~tech_data.unit.str.contains("/kW"), value_column] /= 1e3
    tech_data.unit = tech_data.unit.str.replace("kW", "MW")

    tech_data.loc[tech_data.unit.str.contains("/GWh"), value_column] /= 1e3
    tech_data.unit = tech_data.unit.str.replace("/GWh", "/MWh")

    tech_data.loc[tech_data.unit.str.contains("/GJ"), value_column] *= 3.6
    tech_data.unit = tech_data.unit.str.replace("/GJ", "/MWh")

    # Harmonise individual units so that they can be handled later
    tech_data.unit = tech_data.unit.str.replace(" a year", "/year")
    tech_data.unit = tech_data.unit.str.replace("2015EUR", "EUR")
    tech_data.unit = tech_data.unit.str.replace("2015-EUR", "EUR")
    tech_data.unit = tech_data.unit.str.replace("2020-EUR", "EUR")
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
    tech_data.unit = tech_data.unit.str.replace("MW/year Methanol", "MW_MeOH/year")
    tech_data.unit = tech_data.unit.str.replace("MWh FT Liquids/year", "MWh_FT")
    tech_data.unit = tech_data.unit.str.replace("MWh methanol", "MWh_MeOH")
    tech_data.unit = tech_data.unit.str.replace("MW/year SNG", "MW_CH4/year")
    tech_data.unit = tech_data.unit.str.replace("MWh SNG", "MWh_CH4")
    tech_data.unit = tech_data.unit.str.replace("MW SNG", "MW_CH4")
    tech_data.unit = tech_data.unit.str.replace("EUR/MWh of total input", "EUR/MWh_e")
    tech_data.unit = tech_data.unit.str.replace("EUR/MWeh", "EUR/MWh_e")
    tech_data.unit = tech_data.unit.str.replace("% -points of heat loss", "MWh_th/MWh_el")


    tech_data.unit = tech_data.unit.str.replace("FT Liquids Output, MWh/MWh Total Inpu", "MWh_FT/MWh_H2")
    # biomass-to-methanol-specific
    if isinstance(tech_data.index, pd.MultiIndex):
        tech_data.loc[tech_data.index.get_level_values(1)=="Methanol Output,", "unit"] = "MWh_MeOH/MWh_th"
        tech_data.loc[tech_data.index.get_level_values(1)=='District heat  Output,', "unit"] =  "MWh_th/MWh_th"
        tech_data.loc[tech_data.index.get_level_values(1)=='Electricity Output,', "unit"] =  "MWh_e/MWh_th"
       
    # Ammonia-specific
    tech_data.unit = tech_data.unit.str.replace("MW Ammonia output", "MW_NH3") #specific investment
    tech_data.unit = tech_data.unit.str.replace("MW Ammonia", "MW_NH3") #fom
    tech_data.unit = tech_data.unit.str.replace("MWh Ammonia", "MWh_NH3") #vom
    tech_data.loc[tech_data.unit=='EUR/MW/y', "unit"] = 'EUR/MW/year'

    # convert per unit costs to MW
    cost_per_unit = tech_data.unit.str.contains("/unit")
    tech_data.loc[cost_per_unit, value_column] = tech_data.loc[cost_per_unit, value_column].apply(
                                                lambda x: (x / tech_data.loc[(x.name[0],
                                                                            "Heat production capacity for one unit")][value_column]).iloc[0,:],
                                              axis=1)
    tech_data.loc[cost_per_unit, "unit"] = tech_data.loc[cost_per_unit,
                                                         "unit"].str.replace("/unit", "/MW_th")

    if source == "dea":
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

        if "methanolisation" in tech_data.index:
            tech_data = tech_data.sort_index()
            tech_data.loc[('methanolisation', 'Variable O&M'), "unit"] = "EUR/MWh_MeOH"
    
    tech_data.unit = tech_data.unit.str.replace("\)", "")
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
    tech_data = tech_data.drop("PV module conversion efficiency [p.u.]", level=1)

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
    to_drop = [("hydrogen storage tank type 1 including compressor", ' - Charge efficiency')]
    to_drop.append(("hydrogen storage tank type 1 including compressor", ' - Discharge efficiency'))
    to_drop.append(("hydrogen storage underground", ' - Charge efficiency'))
    to_drop.append(("hydrogen storage underground", ' - Discharge efficiency'))
    tech_data.loc[("hydrogen storage underground", "Round trip efficiency"), years] *= 100
    tech_data.loc[("hydrogen storage tank type 1 including compressor", "Round trip efficiency"), years] *= 100



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
    for tech in tech_data.index.get_level_values(0).unique():
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
                        (df.unit=="EUR/MW_NH3") |
                        (df.unit=="EUR/MWhCapacity") |
                        (df.unit=="EUR/MWh") |
                        (df.unit=="EUR/MW_CH4") |
                        (df.unit=="EUR/MWh/year") |
                        (df.unit=="EUR/MW_e, 2020") |
                        (df.unit=="EUR/MW input") |
                        (df.unit=='EUR/MW-methanol') |
                        (df.unit=="EUR/t_N2/h")) # air separation unit
                    ].copy()
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
                       ((df.unit==investment.unit.iloc[0]+"/year")|
                        (df.unit=="EUR/MW/km/year")|
                        (df.unit=="EUR/MW/year")|
                        (df.unit=="EUR/MW_e/y, 2020")|
                        (df.unit=="EUR/MW_e/y")|
                        (df.unit=="EUR/MW_FT/year")|
                        (df.unit=="EUR/MWh_FT")|
                        (df.unit=="EUR/MW_MeOH/year")|
                        (df.unit=="EUR/MW_CH4/year")|
                        (df.unit=='% of specific investment/year')|
                        (df.unit==investment.unit.str.split(" ").iloc[0][0]+"/year"))].copy()
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
                                                         (df.unit=="EUR/MWh_NH3") |
                                                         (df.unit=="EUR/MWh_MeOH") |
                                                         (df.unit=="EUR/MWh/year") |
                                                         (df.unit=="EUR/MWh/km") |
                                                         (df.unit=="EUR/MWh") |
                                                         (df.unit=="EUR/MWhoutput") |
                                                         (df.unit=="EUR/MWh_CH4") |
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
                         (df.index.str.contains("Hydrogen Output"))|
                         (df.index.str.contains("FT Liquids Output, MWh/MWh Total Input"))|
                         (df.index.str.contains("Methanol Output"))|
                         (df.index.str.contains("District heat  Output"))|
                         (df.index.str.contains("Electricity Output"))|
                         (df.index.str.contains("hereof recoverable for district heating"))|
                         (df.index.str.contains("Bio SNG"))|
                         (df.index == ("Hydrogen")))
                         & ((df.unit=="%") |  (df.unit =="% total size") |
                           (df.unit =="% of fuel input") |
                           (df.unit =="MWh_H2/MWh_e") |
                           (df.unit =="%-points of heat loss") |
                           (df.unit =="MWh_MeOH/MWh_th") |
                           (df.unit =="MWh_e/MWh_th") |
                           (df.unit =="MWh_th/MWh_th") |
                           (df.unit =='MWh/MWh Total Input') |
                           df.unit.str.contains("MWh_FT/MWh_H2"))
                         & (~df.index.str.contains("name plate"))].copy()

        if tech == 'Fischer-Tropsch':
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
            efficiency_h2 = efficiency[efficiency.index.str.contains("Hydrogen Output")].copy()
            efficiency_h2["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_h2])

        # check if electric and heat efficiencies are given
        if (any(["Electric" in ind for ind in efficiency.index]) and
            any(["Heat" in ind for ind in efficiency.index])):
            efficiency_heat = efficiency[efficiency.index.str.contains("Heat")].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
            efficiency = efficiency[efficiency.index.str.contains("Electric")].copy()
            efficiency["parameter"] = "efficiency"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency])

        elif tech == "biomass-to-methanol":
            efficiency_heat = efficiency[efficiency.index.str.contains("District heat")].copy()
            efficiency_heat["parameter"] = "efficiency-heat"
            efficiency_heat.loc[:,years] *= 100 # in %
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_heat])
            efficiency_elec = efficiency[efficiency.index.str.contains("Electric")].copy()
            efficiency_elec["parameter"] = "efficiency-electricity"
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_elec])
            efficiency_meoh = efficiency[efficiency.index.str.contains("Methanol")].copy()
            efficiency_meoh["parameter"] = "efficiency"
            efficiency_meoh.loc[:,years] *= 100 # in %
            clean_df[tech] = pd.concat([clean_df[tech], efficiency_meoh])

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
    data.index.set_names(["technology", "parameter"], inplace=True)
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
    investment = gas_storage.loc['Total cost, 100 mio Nm3 active volume'].iloc[0]
    # convert million EUR/1.1 TWh -> EUR/kWh
    investment /= (1.1 * 1e3)
    data.loc[("gas storage", "investment"), years] = investment
    data.loc[("gas storage", "investment"), "source"] = source_dict["DEA"]
    data.loc[("gas storage", "investment"), "further description"] = "150 Underground Storage of Gas, Establishment of one cavern (units converted)"
    data.loc[("gas storage", "investment"), "unit"] = "EUR/kWh"
    data.loc[("gas storage", "investment"), "currency_year"] = 2015
    
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
    data.loc[("gas storage charger", "investment"), "currency_year"] = 2015

    
    data.loc[("gas storage discharger", "investment"), "source"] = source_dict["DEA"]
    data.loc[("gas storage discharger", "investment"), "further description"] = "150 Underground Storage of Gas, Process equipment (units converted)"
    data.loc[("gas storage discharger", "investment"), "unit"] = "EUR/kW"
    data.loc[("gas storage charger", "investment"), "currency_year"] = 2015

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
    costs_pypsa.rename({'hydrogen underground storage': 'hydrogen storage underground'},
                       inplace=True)

    #convert EUR/m^3 to EUR/kWh for 40 K diff and 1.17 kWh/m^3/K
    costs_pypsa.loc[('decentral water tank storage','investment'),
                    'value'] /= 1.17*40
    costs_pypsa.loc[('decentral water tank storage','investment'),'unit'] = 'EUR/kWh'

    return costs_pypsa

def add_manual_input(data):

    df = pd.read_csv(snakemake.input['manual_input'], quotechar='"',sep=',', keep_default_na=False)
    df = df.rename(columns={"further_description": "further description"})

   
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
            try:
                s["currency_year"] = int(c["currency_year"].values[0]) 
            except ValueError:
                s["currency_year"] = np.nan
            for col in ['unit','source','further description']:
                s[col] = "; and\n".join(c[col].unique().astype(str))
            s = s.rename({"further_description":"further description"}) # match column name between manual_input and original TD workflow
            l.append(s)

    new_df = pd.DataFrame(l).set_index(['technology','parameter'])
    data.index.set_names(["technology", "parameter"], inplace=True)
    # overwrite DEA data with manual input
    data = new_df.combine_first(data)

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
    costs_ISE["unit"] = costs_ISE.unit.replace({"a": "years", "% Invest": "%"})
    costs_ISE["source"] = source_dict["ISE"]
    costs_ISE['further description'] = costs_ISE.reset_index()["technology"].values
    # could not find specific currency year in report, assume year of publication
    costs_ISE['currency_year'] = 2020

    return costs_ISE


def rename_ISE_vehicles(costs_vehicles):
    """
    rename ISE_vehicles costs to fit to tech data
    """

    costs_vehicles.rename(index = {"Investition": "investment",
                          "Lebensdauer": "lifetime",
                          "M/O-Kosten": "FOM",
			"Wirkungsgrad*" : "efficiency",
			"PKW Batterie-Elektromotor" : "Battery electric (passenger cars)",
			"LKW Batterie-Elektromotor" : "Battery electric (trucks)",
			"LKW H2- Brennstoffzelle": "Hydrogen fuel cell (trucks)",
			"PKW H2- Brennstoffzelle": "Hydrogen fuel cell (passenger cars)",
			"LKW ICE- Flï¿½ssigtreibstoff": "Liquid fuels ICE (trucks)",
			"PKW ICE- Flï¿½ssigtreibstoff": "Liquid fuels ICE (passenger cars)",
			"LKW Ladeinfrastruktur Brennstoffzellen Fahrzeuge * LKW": "Charging infrastructure fuel cell vehicles trucks",
			"PKW Ladeinfrastruktur Brennstoffzellen Fahrzeuge * PKW": "Charging infrastructure fuel cell vehicles passenger cars",
			"PKW Ladeinfrastruktur schnell  (reine) Batteriefahrzeuge*" : "Charging infrastructure fast (purely) battery electric vehicles passenger cars",
			"Ladeinfrastruktur langsam (reine) Batteriefahrzeuge*" : "Charging infrastructure slow (purely) battery electric vehicles passenger cars"},
                 columns = {"Einheit": "unit",
                            "2020": 2020,
                            "2025": 2025,
                            "2030": 2030,
                            "2035": 2035,
                            "2040": 2040,
                            "2045": 2045,
                            "2050": 2050}, inplace=True)
    costs_vehicles.index.names = ["technology", "parameter"]
    costs_vehicles["unit"] = costs_vehicles.unit.replace({"a": "years", "% Invest": "%"})
    costs_vehicles["source"] = source_dict["vehicles"]
    # could not find specific currency year in report, assume year of publication
    costs_vehicles["currency_year"] = 2020
    costs_vehicles['further description'] =  costs_vehicles.reset_index()["technology"].values
    return costs_vehicles

def carbon_flow(costs,year):
    # NB: This requires some digits of accuracy; rounding to two digits creates carbon inbalances when scaling up
    c_in_char = 0 # Carbon ending up in char: zero avoids inbalace -> assumed to be circulated back and eventually end up in one of the other output streams
    medium_out = ''
    CH4_specific_energy = 50 #GJ/t methane

    btlcost_data = np.interp(x=years, xp=[2020, 2050], fp=[3500, 2000])
    btl_cost = pd.Series(data=btlcost_data, index=years)

    bmH2cost_data = np.interp(x=years, xp=[2020, 2050], fp=[4000, 2500])
    bmH2_cost = pd.Series(data=bmH2cost_data, index=years)

    btleta_data = np.interp(x=years, xp=[2020, 2050], fp=[0.35, 0.45])
    btl_eta = pd.Series(data=btleta_data, index=years)

    #Adding pelletizing cost to biomass boiler
    costs.loc[('biomass boiler', 'pelletizing cost'), 'value'] = 9
    costs.loc[('biomass boiler', 'pelletizing cost'), 'unit'] = "EUR/MWh_pellets"
    costs.loc[('biomass boiler', 'pelletizing cost'), 'currency_year'] = 2019
    costs.loc[('biomass boiler', 'pelletizing cost'), 'source'] = "Assumption based on doi:10.1016/j.rser.2019.109506"


    for tech in ['Fischer-Tropsch', 'methanolisation', 'BtL', 'biomass-to-methanol', 'BioSNG', 'biogas',
                 'biogas CC', 'digestible biomass to hydrogen',
                 'solid biomass to hydrogen', 'electrobiofuels']:
        inv_cost = 0
        eta = 0
        lifetime = 0
        FOM = 0
        VOM = 0
        currency_year = np.nan
        source = 'TODO'
        co2_capture_rate = 0.90

        if not (tech, "capture rate") in costs.index:
            costs.loc[(tech, 'capture rate'), 'value'] = co2_capture_rate
            costs.loc[(tech, 'capture rate'), 'unit'] = "per unit"
            costs.loc[(tech, 'capture rate'), 'source'] = "Assumption based on doi:10.1016/j.biombioe.2015.01.006"


        if tech == 'BtL':
            inv_cost = btl_cost[year]
            medium_out = 'oil'
            eta = btl_eta[year]
            source = "doi:10.1016/j.enpol.2017.05.013"
            currency_year = 2017

        if tech == 'biomass-to-methanol':
            medium_out = 'methanol'

        elif tech == 'BioSNG':
            medium_out = 'gas'
            lifetime = 25

        elif tech in ['biogas', 'biogas CC']:
            eta = 1
            source = "Assuming input biomass is already given in biogas output"
            AD_CO2_share = 0.4 #volumetric share in biogas (rest is CH4)


        elif tech == 'biogas plus hydrogen':
            #NB: this falls between power to gas and biogas and should be used with care, due to possible minor
            # differences in resource use etc. which may tweak results in favour of one tech or another
            eta = 1.6
            H2_in = 0.46

            heat_out = 0.19
            source = "Calculated from data in Danish Energy Agency, data_sheets_for_renewable_fuels.xlsx"
            costs.loc[(tech, 'hydrogen input'), 'value'] = H2_in
            costs.loc[(tech, 'hydrogen input'), 'unit'] = "MWh_H2/MWh_CH4"
            costs.loc[(tech, 'hydrogen input'), 'source'] = source

            costs.loc[(tech, 'heat output'), 'value'] = heat_out
            costs.loc[(tech, 'heat output'), 'unit'] = "MWh_th/MWh_CH4"
            costs.loc[(tech, 'heat output'), 'source'] = source
            currency_year = costs.loc[('biogas plus hydrogen', 'VOM'), "currency_year"]

            #TODO: this needs to be refined based on e.g. stoichiometry:
            AD_CO2_share = 0.1 #volumetric share in biogas (rest is CH4).

        elif tech == 'digestible biomass to hydrogen':
            inv_cost = bmH2_cost[year]
            eta = 0.39
            FOM = 4.25
            currency_year = 2014
            source = 'Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014' #source_dict('HyNOW')

        elif tech == 'solid biomass to hydrogen':
            inv_cost = bmH2_cost[year]
            eta = 0.56
            FOM = 4.25
            currency_year = 2014
            source = 'Zech et.al. DBFZ Report Nr. 19. Hy-NOW - Evaluierung der Verfahren und Technologien für die Bereitstellung von Wasserstoff auf Basis von Biomasse, DBFZ, 2014' #source_dict('HyNOW')

        if eta > 0:
            costs.loc[(tech, 'efficiency'), 'value'] = eta
            costs.loc[(tech, 'efficiency'), 'unit'] = "per unit"
            costs.loc[(tech, 'efficiency'), 'source'] = source

        if tech in ['BioSNG', 'BtL', 'biomass-to-methanol']:
            input_CO2_intensity = costs.loc[('solid biomass', 'CO2 intensity'), 'value']

            costs.loc[(tech, 'C in fuel'), 'value'] = costs.loc[(tech, 'efficiency'), 'value'] \
                                                  * costs.loc[(medium_out, 'CO2 intensity'), 'value'] \
                                                  / input_CO2_intensity
            costs.loc[(tech, 'C stored'), 'value'] = 1 - costs.loc[(tech, 'C in fuel'), 'value'] - c_in_char
            costs.loc[(tech, 'CO2 stored'), 'value'] = input_CO2_intensity * costs.loc[(tech, 'C stored'), 'value']

            costs.loc[(tech, 'C in fuel'), 'unit'] = "per unit"
            costs.loc[(tech, 'C stored'), 'unit'] = "per unit"
            costs.loc[(tech, 'CO2 stored'), 'unit'] = "tCO2/MWh_th"

            costs.loc[(tech, 'C in fuel'), 'source'] = "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            costs.loc[(tech, 'C stored'), 'source'] = "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"
            costs.loc[(tech, 'CO2 stored'), 'source'] = "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"

        elif tech in ['electrobiofuels']:

            input_CO2_intensity = costs.loc[('solid biomass', 'CO2 intensity'), 'value']
            oil_CO2_intensity = costs.loc[('oil', 'CO2 intensity'), 'value']

            costs.loc[('electrobiofuels', 'C in fuel'), 'value'] = (costs.loc[('BtL', 'C in fuel'), 'value']
                                                                    + costs.loc[('BtL', 'C stored'), 'value']
                                                                    * costs.loc[('Fischer-Tropsch', 'capture rate'), 'value'])
            costs.loc[('electrobiofuels', 'C in fuel'), 'unit'] = 'per unit'
            costs.loc[('electrobiofuels', 'C in fuel'), 'source'] = 'Stoichiometric calculation'

            costs.loc[('electrobiofuels', 'efficiency-biomass'), 'value'] = costs.loc[('electrobiofuels', 'C in fuel'), 'value'] \
                                                                            * input_CO2_intensity / oil_CO2_intensity
            costs.loc[('electrobiofuels', 'efficiency-biomass'), 'unit'] = 'per unit'
            costs.loc[('electrobiofuels', 'efficiency-biomass'), 'source'] = 'Stoichiometric calculation'


            efuel_scale_factor = costs.loc[('BtL', 'C stored'), 'value']* costs.loc[('Fischer-Tropsch', 'capture rate'), 'value']

            costs.loc[('electrobiofuels', 'efficiency-hydrogen'), 'value'] = costs.loc[('Fischer-Tropsch', 'efficiency'), 'value']\
                                                                             / efuel_scale_factor
            costs.loc[('electrobiofuels', 'efficiency-hydrogen'), 'unit'] = 'per unit'
            costs.loc[('electrobiofuels', 'efficiency-hydrogen'), 'source'] = 'Stoichiometric calculation'

            costs.loc[('electrobiofuels', 'efficiency-tot'), 'value'] = (1 /
                                                                         (1 / costs.loc[('electrobiofuels', 'efficiency-hydrogen'), 'value'] +
                                                                          1 / costs.loc[('electrobiofuels', 'efficiency-biomass'), 'value']))
            costs.loc[('electrobiofuels', 'efficiency-tot'), 'unit'] = 'per unit'
            costs.loc[('electrobiofuels', 'efficiency-tot'), 'source'] = 'Stoichiometric calculation'

            inv_cost = btl_cost[year] + costs.loc[('Fischer-Tropsch', 'investment'), 'value'] * efuel_scale_factor
            VOM = costs.loc[('BtL', 'VOM'), 'value'] + costs.loc[('Fischer-Tropsch', 'VOM'), 'value'] * efuel_scale_factor
            FOM = costs.loc[('BtL', 'FOM'), 'value']
            medium_out = 'oil'
            currency_year = costs.loc[('Fischer-Tropsch', 'investment'), "currency_year"]
            source = "combination of BtL and electrofuels"

        elif tech in ['biogas', 'biogas CC', 'biogas plus hydrogen']:
            CH4_density = 0.657 #kg/Nm3
            CO2_density = 1.98 #kg/Nm3
            CH4_vol_energy_density = CH4_specific_energy * CH4_density / (1000 * 3.6) #MJ/Nm3 -> MWh/Nm3
            CO2_weight_share = AD_CO2_share * CO2_density

            costs.loc[(tech, 'CO2 stored'), 'value'] = CO2_weight_share / CH4_vol_energy_density / 1000 #tCO2/MWh,in (NB: assuming the input is already given in the biogas potential and cost
            costs.loc[(tech, 'CO2 stored'), 'unit'] = "tCO2/MWh_th"
            costs.loc[(tech, 'CO2 stored'), 'source'] = "Stoichiometric calculation, doi:10.1016/j.apenergy.2022.120016"

        if inv_cost > 0:
            costs.loc[(tech, 'investment'), 'value'] = inv_cost
            costs.loc[(tech, 'investment'), 'unit'] = "EUR/kW_th"
            costs.loc[(tech, 'investment'), 'source'] = source
            costs.loc[(tech, 'investment'), 'currency_year'] = currency_year

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
            costs.loc[(tech, 'VOM'), 'currency_year'] = currency_year

    return costs

def energy_penalty(costs):

    # Energy penalty for biomass carbon capture
    # Need to take steam production for CC into account, assumed with the main feedstock,
    # e.g. the input biomass is used also for steam, and the efficiency for el and heat is scaled down accordingly

    for tech in ['central solid biomass CHP CC', 'waste CHP CC', 'solid biomass boiler steam CC', 'direct firing solid fuels CC', 'direct firing gas CC', 'biogas CC']:

        if 'powerboost' in tech:
            boiler = 'electric boiler steam'
            feedstock = 'solid biomass'
            co2_capture = costs.loc[(feedstock, 'CO2 intensity'), 'value']
        elif 'gas' in tech:
            boiler = 'gas boiler steam'
            feedstock = 'gas'
            co2_capture = costs.loc[(feedstock, 'CO2 intensity'), 'value']
        elif 'biogas' in tech:
            boiler = 'gas boiler steam'
            co2_capture = costs.loc[(tech, 'CO2 stored'), 'value']
        else:
            boiler = 'solid biomass boiler steam'
            feedstock = 'solid biomass'
            co2_capture = costs.loc[(feedstock, 'CO2 intensity'), 'value']

        #Scaling biomass input to account for heat demand of carbon capture
        scalingFactor = 1 / (1 + co2_capture * costs.loc[
            ('biomass CHP capture', 'heat-input'), 'value']
                             / costs.loc[(boiler, 'efficiency'), 'value'])

        eta_steam = (1 - scalingFactor) * costs.loc[(boiler, 'efficiency'), 'value']
        eta_old = costs.loc[(tech, 'efficiency'), 'value']

        temp = costs.loc[(tech, 'efficiency'), 'value']
        eta_main = costs.loc[(tech, 'efficiency'), 'value'] * scalingFactor

        # Adapting investment share of tech due to steam boiler addition. Investment per MW_el.
        costs.loc[(tech, 'investment'), 'value'] = costs.loc[(tech, 'investment'), 'value'] * eta_old / eta_main \
            + costs.loc[(boiler, 'investment'), 'value'] * eta_steam / eta_main
        costs.loc[(tech, 'investment'), 'source'] = 'Combination of ' + tech + ' and ' + boiler
        costs.loc[(tech, 'investment'), 'further description'] = ''

        if costs.loc[(tech, 'VOM'), 'value']:
            break
        else:
            costs.loc[(tech, 'VOM'), 'value'] = 0.

        costs.loc[(tech, 'VOM'), 'value'] = costs.loc[(tech, 'VOM'), 'value'] * eta_old / eta_main \
            + costs.loc[(boiler, 'VOM'), 'value'] * eta_steam / eta_main
        costs.loc[(tech, 'VOM'), 'source'] = 'Combination of ' + tech + ' and ' + boiler
        costs.loc[(tech, 'VOM'), 'further description'] = ''

        costs.loc[(tech, 'efficiency'), 'value'] = eta_main
        costs.loc[(tech, 'efficiency'), 'source'] = 'Combination of ' + tech + ' and ' + boiler
        costs.loc[(tech, 'efficiency'), 'further description'] = ''

        if 'CHP' in tech:
            costs.loc[(tech, 'efficiency-heat'), 'value'] = \
                costs.loc[(tech, 'efficiency-heat'), 'value'] * scalingFactor \
                    + costs.loc[('solid biomass', 'CO2 intensity'), 'value'] * \
                (costs.loc[('biomass CHP capture', 'heat-output'), 'value'] +
                 costs.loc[('biomass CHP capture', 'compression-heat-output'), 'value'])
            costs.loc[(tech, 'efficiency-heat'), 'source'] = 'Combination of ' + tech + ' and ' + boiler
            costs.loc[(tech, 'efficiency-heat'), 'further description'] = ''

        if 'biogas CC' in tech:
            costs.loc[(tech, 'VOM'), 'value'] = 0
            costs.loc[(tech, 'VOM'), 'unit'] = 'EUR/MWh'

        costs.loc[(tech, 'VOM'), 'value'] = costs.loc[(tech, 'VOM'), 'value'] * eta_old / eta_main \
            + costs.loc[(boiler, 'VOM'), 'value'] * eta_steam / eta_main
        costs.loc[(tech, 'VOM'), 'source'] = 'Combination of ' + tech + ' and ' + boiler
        costs.loc[(tech, 'VOM'), 'further description'] = ''

    return costs

def add_egs_data(data):
    """
    Adds data of enhanced geothermal systems.

    Data taken from Aghahosseini, Breyer 2020: From hot rock to useful energy...
    
    """ 
    parameters = ["CO2 intensity", "lifetime", "efficiency residential heat", "efficiency electricity", "FOM"]
    techs = ["geothermal"]
    multi_i = pd.MultiIndex.from_product([techs, parameters], names=["technology", "parameter"])
    geoth_df = pd.DataFrame(index=multi_i, columns=data.columns)
    years = [col for col in data.columns if isinstance(col, int)]

    # lifetime
    geoth_df.loc[("geothermal", "lifetime"), years] = 30 #years
    geoth_df.loc[("geothermal", "lifetime"), "unit"] = "years"
    geoth_df.loc[("geothermal", "lifetime"), "source"] = source_dict["Aghahosseini2020"]

    # co2 emissions
    geoth_df.loc[("geothermal", "CO2 intensity"), years] = 0.12 # tCO2/MWh_el
    geoth_df.loc[("geothermal", "CO2 intensity"), "unit"] = "tCO2/MWh_el"
    geoth_df.loc[("geothermal", "CO2 intensity"), "source"] = source_dict["Aghahosseini2020"]
    geoth_df.loc[("geothermal", "CO2 intensity"), "further description"] = "Likely to be improved; Average of 85 percent of global egs power plant capacity"

    # efficiency for heat generation using organic rankine cycle
    geoth_df.loc[("geothermal", "efficiency residential heat"), years] = 0.8
    geoth_df.loc[("geothermal", "efficiency residential heat"), "unit"] = "per unit"
    geoth_df.loc[("geothermal", "efficiency residential heat"), "source"] = "{}; {}".format(source_dict["Aghahosseini2020"], source_dict["Breede2015"]) 
    geoth_df.loc[("geothermal", "efficiency residential heat"), "further description"] = "This is a rough estimate, depends on local conditions"

    # efficiency for electricity generation using organic rankine cycle
    geoth_df.loc[("geothermal", "efficiency electricity"), years] = 0.1
    geoth_df.loc[("geothermal", "efficiency electricity"), "unit"] = "per unit"
    geoth_df.loc[("geothermal", "efficiency electricity"), "source"] = "{}; {}".format(source_dict["Aghahosseini2020"], source_dict["Breede2015"]) 
    geoth_df.loc[("geothermal", "efficiency electricity"), "further description"] = "This is a rough estimate, depends on local conditions"

    # relative additional capital cost of using residual heat for district heating (25 percent)
    geoth_df.loc[("geothermal", "district heating cost"), years] = 0.25
    geoth_df.loc[("geothermal", "district heating cost"), "unit"] = "%"
    geoth_df.loc[("geothermal", "district heating cost"), "source"] = "{}".format(source_dict["Frey2022"]) 
    geoth_df.loc[("geothermal", "district heating cost"), "further description"] = "If capital cost of electric generation from EGS is 100%, district heating adds additional 25%"

    # fixed operational costs
    geoth_df.loc[("geothermal", "FOM"), years] = 2.
    geoth_df.loc[("geothermal", "FOM"), "unit"] = "%/year"
    geoth_df.loc[("geothermal", "FOM"), "source"] = source_dict["Aghahosseini2020"] 
    geoth_df.loc[("geothermal", "FOM"), "further description"] = "Both for flash, binary and ORC plants. See Supplemental Material for details"
    
    geoth_df = geoth_df.dropna(axis=1, how='all')
    
    return pd.concat([data, geoth_df])


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


def add_SMR_data(data):
    """Add steam  methane reforming (SMR)  technology data.

    investment cost :
        Currently no cost reduction for investment costs of SMR CC assumed.

        - IEA (2020) [1]: assumes cost reduction -19.2% from 2019-2050
        - Agora [2]: no cost reduction

    carbon capture rate:
        - IEA (2020) [1]: 0.9
        - Agora [2]: 0.9
        - [3]: 0.9
        - Timmerberg et al.: 0.56-0.9

    efficiency:
        - Agora SMR + CC (LHV/LHV) 0.58

    [1] IEA (2020) https://www.iea.org/data-and-statistics/charts/global-average-levelised-cost-of-hydrogen-production-by-energy-source-and-technology-2019-and-2050
    [2] Agora (2021) p.52 https://static.agora-energiewende.de/fileadmin/Projekte/2021/2021_02_EU_H2Grid/A-EW_203_No-regret-hydrogen_WEB.pdf
    [3] p.12 https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1011506/Hydrogen_Production_Costs_2021.pdf
    """
    parameters = ["FOM", "investment", "lifetime", "efficiency"]
    techs = ["SMR", "SMR CC"]
    multi_i = pd.MultiIndex.from_product([techs, parameters], names=["technology", "parameter"])
    SMR_df = pd.DataFrame(index=multi_i, columns=data.columns)

    # efficiencies per unit in LHV (stays constant 2019 to 2050)
    SMR_df.loc[("SMR", "efficiency"), years] =  0.76
    SMR_df.loc[("SMR CC", "efficiency"), years] =  0.69
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
    SMR_df.loc[(techs, "FOM"), "further description"] = "Technology data for renewable fuels, in pdf on table 3 p.311"

    # investment
    # investment given in unit EUR/kg H_2/h -> convert to EUR/MW_CH4
    # lower heating value (LHV) of H2
    LHV_H2  = 33.33 # unit kWh/kg
    SMR = 12500 / LHV_H2 * 1e3 * 1/SMR_df.loc[("SMR", "efficiency"), years]
    SMR_CCS = 14500 / LHV_H2 * 1e3 * 1/SMR_df.loc[("SMR", "efficiency"), years]

    SMR_df.loc[("SMR", "investment"), years] = SMR
    SMR_df.loc[("SMR CC", "investment"), years] = SMR_CCS
    SMR_df.loc[(techs, "investment"), "source"] = source_dict["DEA"]
    SMR_df.loc[(techs, "investment"), "unit"] = "EUR/MW_CH4"
    SMR_df.loc[(techs, "investment"), "currency_year"] = 2015
    SMR_df.loc[(techs, "investment"), "further description"] = "Technology data for renewable fuels, in pdf on table 3 p.311"

    # carbon capture rate
    SMR_df.loc[("SMR CC", "capture_rate"), years] = 0.9
    SMR_df.loc[("SMR CC", "capture_rate"), "source"] = source_dict["IEA"]
    SMR_df.loc[("SMR CC", "capture_rate"), "unit"] = "EUR/MW_CH4"
    SMR_df.loc[("SMR CC", "capture_rate"), "further description"] = "wide range: capture rates betwen 54%-90%"
    
    SMR_df = SMR_df.dropna(axis=1, how='all')
    
    return pd.concat([data, SMR_df])


def add_mean_solar_rooftop(data):
    # take mean of rooftop commercial and residential
    rooftop = (data.loc[data.index.get_level_values(0)
                       .str.contains("solar-rooftop")][years]
               .astype(float).groupby(level=1).mean())
    for col in data.columns[~data.columns.isin(years)]:
        rooftop[col] = data.loc["solar-rooftop residential"][col]
    # set multi index
    rooftop = pd.concat([rooftop], keys=["solar-rooftop"])
    rooftop["source"] = "Calculated. See 'further description'."
    rooftop["further description"] = "Mixed investment costs based on average of 50% 'solar-rooftop commercial' and 50% 'solar-rooftop residential'"
    # add to data
    rooftop.index.names = data.index.names
    data = pd.concat([data, rooftop])
    # add solar assuming 50% utility and 50% rooftop
    solar = (data.loc[["solar-rooftop", "solar-utility"]][years]).astype(float).groupby(level=1).mean()
    for col in data.columns[~data.columns.isin(years)]:
        solar[col] = data.loc["solar-rooftop residential"][col] 
    solar["source"] = "Calculated. See 'further description'."
    solar["further description"] = "Mixed investment costs based on average of 50% 'solar-rooftop' and 50% 'solar-utility'"
    # set multi index
    solar = pd.concat([solar], keys=["solar"])
    solar.index.names = data.index.names
    return pd.concat([data, solar])


def add_energy_storage_database(costs, data_year):
    """Add energy storage database compiled 
    
    Learning rate drop. For example, the nominal DC SB learning rate for RFBs is set at
    4.5%, 1.5% for lead-acid batteries, compared to 10% for Li-ion batteries, corresponding to cost drops of
    17%, 6%, and 35%, respectively. For the rest of the categories for battery-based systems, the learning
    rates were kept the same for all batteries as described in the ESGC 2020 report.

    Fix cost drop. Due to the uncertainties in both anticipated deployments and the correct learning rate to use during the
    initial phase, this work assumes a fixed-cost drop for zinc batteries, gravity, and thermal storage
    systems. For example, a 20% cost drop in DC SB and 10% drop in DCBOS was assumed for zinc batteries,
    while keeping the cost drops for power equipment in line with Li-ion BESS, while system integration,
    EPC, and project development costs are maintained at 90% of Li-ion BESS 2030 values.
    """
    from scipy import interpolate

    print(f"Add energy storage database compiled for year {data_year}")
    # a) Import csv file
    df = pd.read_excel(
        snakemake.input["pnnl_energy_storage"],
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
    )
    df = df.drop(columns=["ref_size_MW", "EP_ratio_h"])
    df = df.fillna(df.dtypes.replace({"float64": 0.0, "O": "NULL"}))
    df.loc[:,"unit"] = df.unit.str.replace("NULL", "per unit")

    # b) Change data to PyPSA format (aggregation of components, units, currency, etc.) 
    df = clean_up_units(df, "value")  # base clean up

    # rewrite technology to be charger, store, discharger, bidirectional-charger
    df.loc[:,"carrier"] = df.carrier.str.replace("NULL", "")
    df.loc[:,"carrier"] = df["carrier"].apply(lambda x: x.split('-'))
    carrier_list_len = df["carrier"].apply(lambda x: len(x))
    carrier_str_len = df["carrier"].apply(lambda x: len(x[0]))
    carrier_first_item =  df["carrier"].apply(lambda x: x[0])
    carrier_last_item =  df["carrier"].apply(lambda x: x[-1])
    bicharger_filter = (carrier_list_len == 3)
    charger_filter = (carrier_list_len == 2) & (carrier_first_item == "elec")
    discharger_filter = (carrier_list_len == 2) & (carrier_last_item == "elec")
    store_filter = (carrier_list_len == 1) & (carrier_str_len > 0)
    reference_filter = (carrier_list_len == 1) & (carrier_first_item == "reference_value")
    df = df[~reference_filter]  # remove reference values
    df.loc[bicharger_filter,"technology_type"] = "bicharger"
    df.loc[charger_filter,"technology_type"] = "charger"
    df.loc[discharger_filter,"technology_type"] = "discharger"
    df.loc[store_filter,"technology_type"] = "store"
    df.loc[df.unit=="EUR/MWh-year", "technology_type"] = "store"
    # Some investment inputs need to be distributed between charger and discharger 
    for tech in df.technology.unique():
        nan_filter = (df.technology==tech) & (carrier_str_len==0) & (df.parameter=="investment")
        store_filter = nan_filter & (df.unit=="EUR/MWh")
        if not df.loc[store_filter].empty:
            df.loc[store_filter, "technology_type"] = "store"  # value will be aggregated later in the groupby
        # charger and discharger with 50% distribution e.g. in case of Hydrogen
        power_filter = nan_filter & (df.unit=="EUR/MW")
        if not df.loc[power_filter].empty:
            agg = df.loc[power_filter].groupby(["technology", "year"]).sum(numeric_only=True)
            charger_investment_filter = charger_filter & (df.technology==tech) & (df.parameter=="investment")
            discharger_investment_filter = discharger_filter & (df.technology==tech) & (df.parameter=="investment")
            df.loc[charger_investment_filter & df.year==2021, "value"] += agg.loc[(tech, 2021)]/2
            df.loc[charger_investment_filter & df.year==2030, "value"] += agg.loc[(tech, 2030)]/2
            df.loc[discharger_investment_filter & df.year==2021, "value"] += agg.loc[(tech, 2021)]/2
            df.loc[discharger_investment_filter & df.year==2030, "value"] += agg.loc[(tech, 2030)]/2
    df.loc[:,"technology"] = df["technology"] + "-" + df["technology_type"]

    # aggregate technology_type and unit
    df = df.groupby(["technology", "unit", "year"]).agg({
        'technology': 'first',
        'year': 'first',
        'parameter': 'first',
        'value': 'sum',
        'unit': 'first',
        'type': 'first',
        'carrier': 'first',
        'technology_type': 'first',
        'source': 'first',
        'note': 'first',
        'reference': 'first',
    }).reset_index(drop=True)

    # calculate %/year FOM on aggregated values
    for tech in df.technology.unique():
        for year in df.year.unique():
            df_tech = df.loc[(df.technology == tech) & (df.year == year)].copy()
            a = df_tech.loc[df_tech.unit=="EUR/MW-year", "value"].values
            b = df_tech.loc[df_tech.unit=="EUR/MW", "value"].values
            df.loc[df_tech.loc[df_tech.unit=="EUR/MW-year"].index, "value"] = a / b * 100 # EUR/MW-year / EUR/MW = %/year
            c = df_tech.loc[df_tech.unit=="EUR/MWh-year", "value"].values
            d = df_tech.loc[df_tech.unit=="EUR/MWh", "value"].values
            df.loc[df_tech.loc[df_tech.unit=="EUR/MWh-year"].index, "value"] = c / d * 100 # EUR/MWh-year / EUR/MWh = %/year

    df.loc[:,"unit"] = df.unit.str.replace("EUR/MW-year", "%/year")
    df.loc[:,"unit"] = df.unit.str.replace("EUR/MWh-year", "%/year")

    # c) Linear Inter/Extrapolation
    # data available for 2021 and 2030, but value for given "year" passed by function needs to be calculated
    for tech in df.technology.unique():
        for param in df.parameter.unique():
            filter = (df.technology == tech) & (df.parameter == param)
            y = df.loc[filter, "value"]
            if y.empty:
                continue  # nothing to interpolate
            elif y.iloc[0]==y.iloc[1] or param=="efficiency" or param=="lifetime":
                ynew = y.iloc[1]  # assume new value is the same as 2030
            elif y.iloc[0]!=y.iloc[1]:
                x = df.loc[filter, "year"] # both values 2021+2030
                first_segment_diff = y.iloc[0]-y.iloc[1]
                endp_first_segment = y.iloc[1]
                
                # Below we create linear segments between 2021-2030
                # While the first segment is known, the others are defined by the initial segments with a accumulating quadratic descreasing gradient
                other_segments_points = [2034, 2039, 2044, 2049, 2054, 2059]
                
                def geometric_series(nominator, denominator=1, number_of_terms=1, start=1):
                    """
                    A geometric series is a series with a constant ratio between successive terms.
                    When moving to infinity the geometric series converges to a limit.
                    https://en.wikipedia.org/wiki/Series_(mathematics)

                    Example:
                    --------
                    nominator = 1
                    denominator = 2
                    number_of_terms = 3
                    start = 0  # 0 means it starts at the first term
                    result = 1/1**0 + 1/2**1 + 1/2**2 = 1 + 1/2 + 1/4 = 1.75

                    If moving to infinity the result converges to 2
                    """
                    return sum([nominator/denominator**i for i in range(start, start+number_of_terms)])

                if  tech=="Hydrogen-discharger" or tech=="Pumped-Heat-store":
                    x1 = pd.concat([x,pd.DataFrame(other_segments_points)], ignore_index=True)
                    y1 = y
                    factor = 5
                    for i in range(len(other_segments_points)): # -1 because of segments
                        cost_at_year = endp_first_segment - geometric_series(nominator=first_segment_diff, denominator=factor, number_of_terms=i+1)
                        y1 = pd.concat([y1, pd.DataFrame([cost_at_year])], ignore_index=True)
                    f = interpolate.interp1d(x1.squeeze(), y1.squeeze(), kind='linear', fill_value="extrapolate")
                elif tech=="Hydrogen-charger":
                    x2 = pd.concat([x,pd.DataFrame(other_segments_points)], ignore_index=True)
                    y2 = y
                    factor = 6.5
                    for i in range(len(other_segments_points)):
                        cost_at_year = endp_first_segment - geometric_series(nominator=first_segment_diff, denominator=factor, number_of_terms=i+1)
                        y2 = pd.concat([y2, pd.DataFrame([cost_at_year])], ignore_index=True)
                    f = interpolate.interp1d(x2.squeeze(), y2.squeeze(), kind='linear', fill_value="extrapolate")  
                else:
                    x3 = pd.concat([x,pd.DataFrame(other_segments_points)], ignore_index=True)
                    y3 = y
                    factor = 2
                    for i in range(len(other_segments_points)):
                        cost_at_year = endp_first_segment - geometric_series(nominator=first_segment_diff, denominator=factor, number_of_terms=i+1)
                        y3 = pd.concat([y3, pd.DataFrame([cost_at_year])], ignore_index=True)
                    f = interpolate.interp1d(x3.squeeze(), y3.squeeze(), kind='linear', fill_value="extrapolate")
                
                option = snakemake.config['energy_storage_database']['pnnl_energy_storage']
                if option.get('approx_beyond_2030') == ["geometric_series"]:
                    ynew = f(data_year)
                if option.get('approx_beyond_2030') == ["same_as_2030"]:
                    if data_year <= 2030:
                        # apply linear interpolation
                        ynew = f(data_year)
                    if data_year > 2030:
                        # apply same value as 2030
                        ynew = y.iloc[1]  # assume new value is the same as 2030

            df_new = pd.DataFrame([{
                "technology": tech,
                "year": data_year,
                "parameter": param,
                "value": ynew, 
                "unit": df.loc[filter, "unit"].unique().item(),
                "source": df.loc[filter, "source"].unique().item(),
                'carrier': df.loc[filter, "carrier"].iloc[1],
                'technology_type': df.loc[filter, "technology_type"].unique().item(),
                'type': df.loc[filter, "type"].unique().item(),
                'note': df.loc[filter, "note"].iloc[1],
                'reference': df.loc[filter, "reference"].iloc[1],
            }])
            # not concat if df year is 2021 or 2030 (otherwhise duplicate)
            if data_year == 2021 or data_year == 2030:
                continue
            else:
                df = pd.concat([df, df_new], ignore_index=True)

    # d) Combine metadata and add to cost database
    df.loc[:,"source"] = df["source"] + ", " + df["reference"]
    for i in df.index:
        df.loc[i,"further description"] = str(
            {
                "carrier": df.loc[i,"carrier"],
                "technology_type": [df.loc[i,"technology_type"]],
                "type": [df.loc[i,"type"]],
                "note": [df.loc[i,"note"]],
            }
        )
    # keep only relevant columns
    df = df.loc[df.year == data_year,["technology", "parameter", "value", "unit", "source", "further description"]]
    tech = df.technology.unique()
    df = df.set_index(['technology', 'parameter'])

    return pd.concat([costs, df]), tech


def prepare_inflation_rate(fn):
    """read in annual inflation rate from Eurostat
    https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/dataflow/ESTAT/prc_hicp_aind/1.0?references=descendants&detail=referencepartial&format=sdmx_2.1_generic&compressed=true
    """
    inflation_rate = pd.read_excel(fn,
                                   sheet_name="Sheet 1", index_col=0,
                                   header=[8])
    inflation_rate = (inflation_rate.loc["European Union - 27 countries (from 2020)"]
                      .dropna()).loc["2001"::]
    inflation_rate.rename(index=lambda x: int(x), inplace=True)
    inflation_rate = inflation_rate.astype(float)
    
    inflation_rate /= 100
        
    return inflation_rate
    
# %% *************************************************************************
#  ---------- MAIN ------------------------------------------------------------
if __name__ == "__main__":
    if 'snakemake' not in globals():
        import os
        from _helpers import mock_snakemake
        #os.chdir(os.path.join(os.getcwd(), "scripts"))
        snakemake = mock_snakemake("compile_cost_assumptions")

    years = snakemake.config['years']
    inflation_rate = prepare_inflation_rate(snakemake.input.inflation_rate)
    
  

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
    tech_data = clean_up_units(tech_data, years, source="dea")

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
    
    # adjust for inflation
    for x in data.index.get_level_values("technology"):
        if x in cost_year_2020:
            data.at[x, "currency_year"] = 2020
        elif x in cost_year_2019:
            data.at[x, "currency_year"] = 2019
        else:
            data.at[x, "currency_year"] = 2015

    # %% (2) -- get data from other sources which need formatting -----------------
    # (a)  ---------- get old pypsa costs ---------------------------------------
    costs_pypsa = pd.read_csv(snakemake.input.pypsa_costs,
                            index_col=[0,2]).sort_index()
    # rename some techs and convert units
    costs_pypsa = rename_pypsa_old(costs_pypsa)

    # (b1) ------- add vehicle costs from Fraunhofer vehicle study ------------------------
    costs_vehicles = pd.read_csv(snakemake.input.fraunhofer_vehicles_costs,
                            engine="python",
                            index_col=[0,1],
                            encoding="ISO-8859-1")
    # rename + reorder to fit to other data
    costs_vehicles = rename_ISE_vehicles(costs_vehicles)
    if 'NT' in costs_vehicles.index:
    	costs_vehicles.drop(['NT'], axis=0, inplace=True, level=0)
    costs_vehicles = convert_units(costs_vehicles)
    # add costs for vehicles
    data = pd.concat([data, costs_vehicles], sort=True)


    # (b) ------- add costs from Fraunhofer ISE study --------------------------
    costs_ISE = pd.read_csv(snakemake.input.fraunhofer_costs,
                            engine="python",
                            index_col=[0,1],
                            encoding = "ISO-8859-1")
    # rename + reorder to fit to other data
    costs_ISE = rename_ISE(costs_ISE)   
    # add costs for gas pipelines
    data = pd.concat([data, costs_ISE.loc[["Gasnetz"]]], sort=True)

    # add (enhanced) geothermal systems data
    data = add_egs_data(data) 

    data = add_manual_input(data)
    # add costs for home batteries

    if snakemake.config["energy_storage_database"].get("ewg_home_battery", True):
        data = add_home_battery_costs(data)
    # add SMR assumptions
    data = add_SMR_data(data)
    # add solar rooftop costs by taking the mean of commercial and residential
    data = add_mean_solar_rooftop(data)
    # %% (3) ------ add additional sources and save cost as csv ------------------
    # [RTD-target-multiindex-df]
    for year in years:
        costs = (data[[year, "unit", "source", "further description",
                       "currency_year"]]
                .rename(columns={year: "value"}))
        costs["value"] = costs["value"].astype(float)

        # biomass is differentiated by biomass CHP and HOP
        costs.loc[('solid biomass', 'fuel'), 'value'] = 12
        costs.loc[('solid biomass', 'fuel'), 'unit'] = 'EUR/MWh_th'
        costs.loc[('solid biomass', 'fuel'), 'source'] = "JRC ENSPRESO ca avg for MINBIOWOOW1 (secondary forest residue wood chips), ENS_Ref for 2040"
        costs.loc[('solid biomass', 'fuel'), 'currency_year'] = 2010 
        
        costs.loc[('digestible biomass', 'fuel'), 'value'] = 15
        costs.loc[('digestible biomass', 'fuel'), 'unit'] = 'EUR/MWh_th'
        costs.loc[('digestible biomass', 'fuel'), 'source'] = "JRC ENSPRESO ca avg for MINBIOAGRW1, ENS_Ref for 2040"
        costs.loc[('digestible biomass', 'fuel'), 'currency_year'] = 2010 
        
        # add solar data from other source than DEA
        if any([snakemake.config['solar_utility_from_vartiaien'], snakemake.config['solar_rooftop_from_etip']]):
            costs = add_solar_from_other(costs)

        # add desalination and clean water tank storage
        costs = add_desalinsation_data(costs)
        # add energy storage database
        if snakemake.config['energy_storage_database']['pnnl_energy_storage'].get("add_data", True):
            costs, tech = add_energy_storage_database(costs, year)
            costs.loc[tech, "currency_year"] = 2020

        # add electrolyzer and fuel cell efficiency from other source than DEA
        if snakemake.config["energy_storage_database"].get("h2_from_budischak", True):
            costs = add_h2_from_other(costs)

        # CO2 intensity
        costs = add_co2_intensity(costs)

        # carbon balances
        costs = carbon_flow(costs,year)

        # energy penalty of carbon capture
        costs = energy_penalty(costs)

        # include old pypsa costs
        check = pd.concat([costs_pypsa, costs], sort=True)

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
        # TODO check currency year from old pypsa cost assumptions
        to_add["currency_year"] = 2015
        costs_tot = pd.concat([costs, to_add], sort=False)

        # single components missing
        comp_missing = costs_pypsa.index.difference(costs_tot.index)
        if (year==years[0]):
            print("single parameters of technologies are missing, using old PyPSA assumptions: ")
            print(comp_missing)
            print("old c_v and c_b values are assumed where given")
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
        costs_tot = adjust_for_inflation(inflation_rate, costs_tot, techs,
                                         costs_tot.currency_year, ["value"])
        
        # format and sort
        costs_tot.sort_index(inplace=True)
        costs_tot.loc[:,'value'] = round(costs_tot.value.astype(float),
                                         snakemake.config.get("ndigits", 2))
        costs_tot.to_csv([v for v in snakemake.output if str(year) in v][0])