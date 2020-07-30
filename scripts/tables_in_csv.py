# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:10:15 2019

@author: Marta
"""

import pandas as pd
import numpy as np


"""
Latex table including FOM, efficiencies and lifetimes
"""

#write latex table
idx = pd.IndexSlice
costs = pd.read_csv('outputs/costs_2030.csv',
                    index_col=list(range(2))).sort_index()

# filename='tables/costs.tex'

# file = open(filename, 'w')
# %%
technologies=['onwind', 'offwind', 'solar-utility', 'solar-rooftop', 'OCGT',
               'hydro', 'ror', 'PHS',
              'central gas CHP', 'central solid biomass CHP',
              'HVDC overhead', 'HVDC inverter pair',
              'battery storage',
              'battery inverter', 'electrolysis', 'fuel cell',
              'hydrogen storage underground', 'hydrogen storage tank',
              'DAC', 'methanation', 'helmeth',
              'central gas boiler', 'decentral gas boiler',
              'central resistive heater', 'decentral resistive heater',
              'central water tank storage', 'decentral water tank storage',
              'water tank charger',
              'decentral air-sourced heat pump',
              'central air-sourced heat pump',
              'decentral ground-sourced heat pump',
              'Gasnetz', 'H2 pipeline',
              'SMR', 'biogas upgrading',
              'micro CHP',
              'decentral solar thermal', 'central solar thermal',
              'electricity distribution grid', 'electricity grid connection',
              'gas storage', 'gas storage charger', 'gas storage discharger',
              ]

name={'onwind' : 'Onshore Wind',
      'offwind' : 'Offshore Wind',
      'solar-utility' : 'Solar PV (utility-scale)',
      'solar-rooftop' : 'Solar PV (rooftop)',
      'OCGT': 'OCGT',
      'CCGT': 'CCGT',
      'coal':  'Coal power plant',
      'lignite': 'Lignite',
      'nuclear': 'Nuclear',
      'hydro':'Reservoir hydro',
      'ror':'Run of river',
      'PHS':'PHS',
      'battery inverter': 'Battery inverter',
      'battery storage': 'Battery storage',
      'hydrogen storage underground': 'H$_2$ storage underground',
      'hydrogen storage tank': 'H$_2$ storage tank',
      'electrolysis': 'Electrolysis',
      'fuel cell': 'Fuel cell',
      'methanation': 'Methanation',
      'DAC': 'DAC (direct-air capture)',
      'central gas boiler': 'Gas boiler central',
      'decentral gas boiler': 'Gas boiler decentral',
      'central resistive heater':'Resistive heater central',
      'decentral resistive heater':'Resistive heater decentral',
      'central gas CHP': 'Gas CHP',
      'central coal CHP': 'Coal CHP',
      'biomass CHP':'Biomass CHP',
      'biomass EOP':'Biomass power plant',
      'biomass HOP':'Biomass central heat plant',
      'central water tank storage': 'Water tank storage central',
      'decentral water tank storage': 'Water tank storage decentral',
      'water tank charger': 'Water tank charger/discharger',
      'HVDC overhead':'HVDC overhead',
      'HVDC inverter pair':'HVDC inverter pair',
      'decentral air-sourced heat pump': 'Air-sourced heat pump decentral',
      'central air-sourced heat pump': 'Air-sourced heat pump central',
      'central ground-sourced heat pump': 'Ground-sourced heat pump central',
      'decentral air-sourced heat pump': 'Air-sourced heat pump decentral',
      'decentral ground-sourced heat pump':  'Ground-sourced heat pump decentral',
      'Gasnetz': 'Gas pipeline',
      'micro CHP': 'Micro CHP',
      'central solid biomass CHP': 'Solid biomass CHP central',
      'helmeth': 'Helmeth (Power to SNG, KIT project)',
      'H2 pipeline': 'H2 pipeline',
      'SMR': 'Steam Methane Reforming (SMR)',
      'biogas upgrading': 'Biogas upgrading',
      'decentral solar thermal': 'Solar thermal central',
      'central solar thermal': 'Solar thermal decentral',
      'electricity distribution grid': 'Electricity distribution grid',
       'electricity grid connection': 'Electricity grid connection',
       'gas storage': 'Gas storage (underground cavern)',
       'gas storage charger': 'Gas storage injection',
       'gas storage discharger': 'Gas storage withdrawl',
      }

dic_ref = {'Technology Data for Energy Plants for Electricity and District heating generation':'DEA_2019',
           'Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity': 'Vartiainen_2019',
           'European PV Technology and Innovation Platform' : 'Vartiainen_2017',
           'Lazard’s Levelized Cost of Energy Analysis - Version 13.0': 'Lazard_2019',
           #'budischak2013':'Budischak_2013',
           #'NREL http://www.nrel.gov/docs/fy09osti/45873.pdf; budischak2013': 'Steward_2009b, Budischak_2013',
           'Schaber thesis':'Schaber_2013',
           'Hagspiel':'Hagspiel_2014',
           'Fasihi':'Fasihi_2017',
           'HP' : ' ',
           'DIW DataDoc http://hdl.handle.net/10419/80348' : 'Schroeder_2013',
            888 : 'water tank charger',
           'BP 2019':'BP_2019',
           'https://www.eia.gov/environment/emissions/co2_vol_mass.php' : 'EIA_emission_coefficients',
           'DIW': 'Schroeder_2013',
           'IEA2011b' : 'BP_2019',
           'Is a 100% renewable European power system feasible by 2050?': 'Zappa_2019, JRC_biomass',
           'Entwicklung der spezifischen Kohlendioxid-Emissionen des deutschen Strommix in den Jahren 1990 - 2018': 'German_Environment_Agency',
}

# Solar thermal collector decentral & 270 & m$^{2}$ & 1.3 & 20 & variable & \cite{Henning20141003} \\
# Solar thermal collector central & 140 & m$^{2}$ & 1.4 & 20 & variable & \cite{Henning20141003} \\
# Building retrofitting\tnote{f} & see text &  & 1 & 50 & 1 & \cite{Henning20141003,PalzerThesis} \\
# High-density district heating network\tnote{f} & 220 & kW\th  & 1 & 40  & 1 & \cite{IEESWV} \\
# Gas distribution network\tnote{f} & 387 & kW\th & 2 & 40 & 1 & based on \cite{bnetza2017} \\
#
# %%
relevant = costs.unstack().loc[technologies]
# costs from 2050 for direct air capture assumed
relevant.loc["DAC", ("value", "investment")] = 210.50
quantities = ["FOM", "VOM", "efficiency", "investment", "lifetime"]
table = relevant["value"][quantities]
table["investment"] = (table["investment"].astype(str) + " "
                       + relevant[("unit", "investment")])

table["VOM"] = table["VOM"].astype(str) + " " + relevant[("unit", "VOM")]

# table.rename(columns={"FOM": "FOM [% per year]",
#                       "lifetime": "lifetime [years]",
#                       "efficiency": "efficiency [per unit]"}, inplace=True)
table.fillna(" ", inplace=True)
table["source"] = relevant["source"][quantities].apply(lambda row:
                                                       row.dropna().unique(),
                                                       axis=1)
table["source"] = table.source.apply(lambda x: ', '.join(x))
table.source.str.contains("DEA")
# shorten source names
table.loc[table.source.str.contains("DEA"), "source"] = "DEA"
table.loc[table.source.str.contains("DIW"), "source"] = "DIW"
table.loc[table.source.str.contains("Fasihi"), "source"] = "Fasihi"
table.loc[table.source.str.contains("Welder"), "source"] = "Welder"
table.loc[table.source.str.contains("gov.uk"), "source"] = "GOV UK"
table.loc[table.source.str.contains("Fraunhofer"), "source"] = "FA ISE"
table.rename(index=name, inplace=True)

table.replace({'_th': '$_{th}$',
               "_el": "$_{el}$",
               "_e": "$_{el}$",
               "kWel": "kW$_{el}$",
               "kWGas": "kW$_{Gas}$"}, regex=True, inplace=True)
# add costs for gas distribution
table.loc['Gas distribution grid'] = table.loc['Electricity distribution grid']


table.sort_index().to_csv("/home/ws/bw0928/Dokumente/own_projects/retrofitting_paper/data/costs2030_copy.csv")
# %% biomass transport
transport_costs = pd.read_csv("/home/ws/bw0928/Dokumente/pypsa-eur-sec/data/biomass/biomass_transport_costs.csv",
                              index_col=0)
countries = ['AL', 'AT', 'BA', 'BE', 'BG', 'CH', 'CZ', 'DE', 'DK', 'EE',
       'ES', 'FI', 'FR', 'GB', 'GR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU',
       'LV', 'ME', 'MK', 'NL', 'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI',
       'SK']

transport_costs = round(transport_costs.reindex(index=countries), ndigits=2)
transport_costs.index.name = "Country"
transport_costs.columns = ["cost"]
transport_costs.to_csv("/home/ws/bw0928/Dokumente/own_projects/retrofitting_paper/data/biomass_transport_costs.csv")
# %%
# add retrofitting cost assumptions
retro = pd.read_csv("/home/ws/bw0928/Dokumente/pypsa-eur-sec/data/retro/retro_cost_germany.csv",
                    index_col=0, usecols=[0,1,2,3])
components = ["wall", "floor", "roof", "window"]
retro = retro.reindex(index=components).rename(index = lambda x: "Building retrofitting " +x)
retro.columns = ["fixed costs [EUR/m$^2$]",
                 "variable costs for additional insulation material [EUR/cm]",
                 "lifetime [years]"]
source_retro = "IWU: Kosten energierelevanter Bau- und Anlagenteile bei der energetischen Modernisierung von Altbauten, Darmstadt, 2015"
source_short = "IWU2015"
retro["source"] = source_short
retro.to_csv("/home/ws/bw0928/Dokumente/own_projects/retrofitting_paper/data/retro_components_costs.csv")

# %% add costs dE
retro_de = pd.read_csv("/home/ws/bw0928/Dokumente/pypsa-eur-sec/resources/retro_cost_elec_s_38.csv",
                       skiprows=[1,2])
retro_de.rename(columns={'Unnamed: 0': "Country",
                         'Unnamed: 1': "Sector"}, inplace=True)
retro_de = round(retro_de, ndigits=2)
retro_de.set_index(["Country", "Sector"], inplace=True)
retro_de.unstack(level=-1).to_csv("/home/ws/bw0928/Dokumente/own_projects/retrofitting_paper/data/retro_costs.csv")
# retro_de.rename(columns={"dE": "dE moderate",
#                          "dE.1": "dE ambitious",
#                          "cost": "cost moderate",
#                          "cost.1": "cost ambitious"})
retro_de = retro_de[["dE", "cost", "dE.1", "cost.1"]]
columns = [["moderate", "ambitious"], ["dE", "costs [EUR/m2]"]]
index = pd.MultiIndex.from_product(columns, names=['strength', 'type'])
retro_de.columns = index
retro_de.rename(index={"tot":"total"}, level=1, inplace=True)
a=retro_de.unstack(level=-1).to_latex(bold_rows=True)
# %%