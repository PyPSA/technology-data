# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:10:15 2019

@author: Marta
"""
#%%
import pandas as pd
import numpy as np
import os

"""
Latex table including FOM, efficiencies and lifetimes
"""

#write latex table
# read 2020 costs
idx = pd.IndexSlice
root_path = os.getcwd()
costs = pd.read_csv(os.path.join(root_path, 'outputs', 'costs_2060.csv'),index_col=list(range(2))).sort_index()

filename='table_inputs.tex'

file = open(filename, 'w')

technologies=['onwind', 'offwind', 'solar-utility', 'solar-rooftop', 'OCGT',
              'CCGT', 'coal', 'lignite', 'nuclear', 'hydro', 'ror', 'PHS',
              'central gas CHP',
              'biomass CHP',
              #'central coal CHP',
              #'biomass HOP',
              #'biomass EOP',
              'HVDC overhead', 'HVDC inverter pair',
              'battery storage',
              'battery inverter',
              'home battery storage',
              'home battery inverter',
              'electrolysis',
              'fuel cell',
              'hydrogen storage underground',
              'hydrogen storage tank type 1',
              'direct air capture',
              'methanation',
              'central gas boiler',
              'decentral gas boiler',
              'central resistive heater',
              'decentral resistive heater',
              'central water tank storage',
              'decentral water tank storage',
              'water tank charger',
              'decentral air-sourced heat pump',
              'central air-sourced heat pump',
              'decentral ground-sourced heat pump',
              'biomass CHP capture',
              'Fischer-Tropsch',
              'SMR',
              'SMR CC',
              'BioSNG',
              'BtL',
              'biogas plus hydrogen',
              'industrial heat pump medium temperature',
              'industrial heat pump high temperature',
              'electric boiler steam',
              'gas boiler steam',
              'solid biomass boiler steam',
              'methanolisation',
              'Compressed-Air-Adiabatic-bicharger',
              'Compressed-Air-Adiabatic-store', 'Concrete-charger',
              'Concrete-discharger', 'Concrete-store', 'Gravity-Brick-bicharger',
              'Gravity-Brick-store', 'Gravity-Water-Aboveground-bicharger',
              'Gravity-Water-Aboveground-store',
              'Gravity-Water-Underground-bicharger',
              'Gravity-Water-Underground-store', 'HighT-Molten-Salt-charger',
              'HighT-Molten-Salt-discharger', 'HighT-Molten-Salt-store',
              'Hydrogen-charger', 'Hydrogen-discharger', 'Hydrogen-store',
              'Lead-Acid-bicharger', 'Lead-Acid-store', 'Liquid-Air-charger',
              'Liquid-Air-discharger', 'Liquid-Air-store',
              'Lithium-Ion-LFP-bicharger', 'Lithium-Ion-LFP-store',
              'Lithium-Ion-NMC-bicharger', 'Lithium-Ion-NMC-store',
              'LowT-Molten-Salt-charger', 'LowT-Molten-Salt-discharger',
              'LowT-Molten-Salt-store', 'Ni-Zn-bicharger', 'Ni-Zn-store',
              'Pumped-Heat-charger', 'Pumped-Heat-discharger',
              'Pumped-Heat-store', 'Pumped-Storage-Hydro-bicharger',
              'Pumped-Storage-Hydro-store', 'Sand-charger', 'Sand-discharger',
              'Sand-store', 'Vanadium-Redox-Flow-bicharger',
              'Vanadium-Redox-Flow-store', 'Zn-Air-bicharger', 'Zn-Air-store',
              'Zn-Br-Flow-bicharger', 'Zn-Br-Flow-store',
              'Zn-Br-Nonflow-bicharger', 'Zn-Br-Nonflow-store'
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
      'home battery inverter': 'Home battery inverter',
      'home battery storage': 'Home battery storage',
      'hydrogen storage underground': 'H$_2$ storage underground',
      'hydrogen storage tank type 1': 'H$_2$ storage tank',
      'electrolysis': 'Electrolysis',
      'fuel cell': 'Fuel cell',
      'methanation': 'Methanation',
      'direct air capture': 'direct air capture',
      'central gas boiler': 'Central gas boiler',
      'decentral gas boiler': 'Domestic gas boiler',
      'central resistive heater':'Central resistive heater',
      'decentral resistive heater':'Domestic resistive heater',
      'central gas CHP':' Gas CHP',
      'central coal CHP':' Coal CHP',
      'biomass CHP':'Biomass CHP',
      'biomass EOP':'Biomass power plant',
      'biomass HOP':'Biomass central heat plant',
      'central water tank storage': 'Central water tank storage',
      'decentral water tank storage': 'Domestic water tank storage',
      'water tank charger': 'Water tank charger/discharger',
      'HVDC overhead':'HVDC overhead',
      'HVDC inverter pair':'HVDC inverter pair',
      #'central heat pump': 'Central heat pump',
      #'decentral heat pump': 'Decentral heat pump',
      #'central ground-sourced heat pump': 'Central ground-sourced heat pump',
      'central air-sourced heat pump': 'Central air-sourced heat pump',
      'decentral air-sourced heat pump': 'Domestic air-sourced heat pump',
      'decentral ground-sourced heat pump':  'Domestic ground-sourced heat pump',
      'biomass CHP capture':'CO$_2$ capture in CHP',
      'Fischer-Tropsch':'Fischer-Tropsch',
      'SMR': 'Steam Methane Reforming',
      'SMR CC': 'Steam Methane Reforming with CC',
      'BioSNG': 'BioSNG',
      'BtL': 'BtL',
      'biogas plus hydrogen': 'biogas plus hydrogen',
      'industrial heat pump medium temperature': 'industrial heat pump medium temperature',
      'industrial heat pump high temperature': 'industrial heat pump high temperature',
      'electric boiler steam': 'electric boiler steam',
      'gas boiler steam': 'gas boiler steam',
      'solid biomass boiler steam': 'solid biomass boiler steam',
      'methanolisation': 'methanolisation'
      'Compressed-Air-Adiabatic-bicharger': 'Compressed-Air-Adiabatic-bicharger',
      'Compressed-Air-Adiabatic-store': 'Compressed-Air-Adiabatic-store',
      'Concrete-charger': 'Concrete-charger',
      'Concrete-discharger': 'Concrete-discharger',
      'Concrete-store': 'Concrete-store',
      'Gravity-Brick-bicharger': 'Gravity-Brick-bicharger',
      'Gravity-Brick-store': 'Gravity-Brick-store',
      'Gravity-Water-Aboveground-bicharger': 'Gravity-Water-Aboveground-bicharger',
      'Gravity-Water-Aboveground-store': 'Gravity-Water-Aboveground-store',
      'Gravity-Water-Underground-bicharger': 'Gravity-Water-Underground-bicharger',
      'Gravity-Water-Underground-store': 'Gravity-Water-Underground-store',
      'HighT-Molten-Salt-charger': 'HighT-Molten-Salt-charger',
      'HighT-Molten-Salt-discharger': 'HighT-Molten-Salt-discharger',
      'HighT-Molten-Salt-store': 'HighT-Molten-Salt-store',
      'Hydrogen-charger': 'Hydrogen-charger',
      'Hydrogen-discharger': 'Hydrogen-discharger',
      'Hydrogen-store': 'Hydrogen-store',
      'Lead-Acid-bicharger': 'Lead-Acid-bicharger',
      'Lead-Acid-store': 'Lead-Acid-store',
      'Liquid-Air-charger': 'Liquid-Air-charger',
      'Liquid-Air-discharger': 'Liquid-Air-discharger',
      'Liquid-Air-store': 'Liquid-Air-store',
      'Lithium-Ion-LFP-bicharger': 'Lithium-Ion-LFP-bicharger',
      'Lithium-Ion-LFP-store': 'Lithium-Ion-LFP-store',
      'Lithium-Ion-NMC-bicharger': 'Lithium-Ion-NMC-bicharger',
      'Lithium-Ion-NMC-store': 'Lithium-Ion-NMC-store',
      'LowT-Molten-Salt-charger': 'LowT-Molten-Salt-charger',
      'LowT-Molten-Salt-discharger': 'LowT-Molten-Salt-discharger',
      'LowT-Molten-Salt-store': 'LowT-Molten-Salt-store',
      'Ni-Zn-bicharger': 'Ni-Zn-bicharger',
      'Ni-Zn-store': 'Ni-Zn-store',
      'Pumped-Heat-charger': 'Pumped-Heat-charger',
      'Pumped-Heat-discharger': 'Pumped-Heat-discharger',
      'Pumped-Heat-store': 'Pumped-Heat-store',
      'Pumped-Storage-Hydro-bicharger': 'Pumped-Storage-Hydro-bicharger',
      'Pumped-Storage-Hydro-store': 'Pumped-Storage-Hydro-store',
      'Sand-charger': 'Sand-charger',
      'Sand-discharger': 'Sand-discharger',
      'Sand-store': 'Sand-store',
      'Vanadium-Redox-Flow-bicharger': 'Vanadium-Redox-Flow-bicharger',
      'Vanadium-Redox-Flow-store': 'Vanadium-Redox-Flow-store',
      'Zn-Air-bicharger': 'Zn-Air-bicharger',
      'Zn-Air-store': 'Zn-Air-store',
      'Zn-Br-Flow-bicharger': 'Zn-Br-Flow-bicharger',
      'Zn-Br-Flow-store': 'Zn-Br-Flow-store',
      'Zn-Br-Nonflow-bicharger': 'Zn-Br-Nonflow-bicharger',
      'Zn-Br-Nonflow-store': 'Zn-Br-Nonflow-store',
      }

dic_ref = {'Technology Data for Energy Plants for Electricity and District heating generation':'DEA_2019',
           'Impact of weighted average cost of capital, capital expenditure, and other parameters on future utility‐scale PV levelised cost of electricity': 'Vartiainen_2019',
           'European PV Technology and Innovation Platform' : 'Vartiainen_2017',
           'Lazard’s Levelized Cost of Energy Analysis - Version 13.0': 'Lazard_2019',
           'budischak2013':'Budischak_2013, DEA_2019',
           #'NREL http://www.nrel.gov/docs/fy09osti/45873.pdf;
           'IWES Interaktion':'Gerhardt_2015, DEA_2019',
           'Schaber thesis':'Schaber_2013',
           'Hagspiel et al. (2014): doi:10.1016/j.energy.2014.01.025 ': 'Hagspiel_2014',
           'Hagspiel':'Hagspiel_2014',
           #'Fasihi':'Fasihi_2017',
           'Fasihi et al 2017, table 1, https://www.mdpi.com/2071-1050/9/2/306':'Fasihi_2017',
           'HP' : ' ',
           'DIW DataDoc http://hdl.handle.net/10419/80348' : 'Schroeder_2013',
            888 : 'water tank charger',
           'BP 2019':'BP_2019',
           'https://www.eia.gov/environment/emissions/co2_vol_mass.php' : 'EIA_emission_coefficients',
           'DIW': 'Schroeder_2013',
           'IEA2011b' : 'BP_2019',
           'Is a 100% renewable European power system feasible by 2050?': 'Zappa_2019, JRC_biomass',
           'Entwicklung der spezifischen Kohlendioxid-Emissionen des deutschen Strommix in den Jahren 1990 - 2018': 'German_Environment_Agency',
           'IEA WEM2017 97USD/boe = http://www.iea.org/media/weowebsite/2017/WEM_Documentation_WEO2017.pdf':'IEA_WEO2017',
           'Danish Energy Agency': 'DEA_2019',
           'Danish Energy Agency, technology_data_for_el_and_dh.xlsx':'DEA_2019',
           'Danish Energy Agency, technology_data_for_el_and_dh_-_0009.xlsx':'DEA_2019',
           'Danish Energy Agency, technology_data_catalogue_for_energy_storage.xlsx':'DEA_2019',
           'Danish Energy Agency, technology_data_catalogue_for_energy_storage.xlsx, Note K.':'DEA_2019',
           'Danish Energy Agency, data_sheets_for_renewable_fuels.xlsx':'DEA_2019',
           'Danish Energy Agency, technology_data_for_industrial_process_heat_0002.xlsx':'DEA_2019',
           'Danish Energy Agency, technologydatafor_heating_installations_marts_2018.xlsx':'DEA_2019',
           'Lazard s Levelized Cost of Energy Analysis - Version 13.0':'Lazard_2019',
           'Global Energy System based on 100% Renewable Energy, Energywatchgroup/LTU University, 2019, Danish Energy Agency, technology_data_catalogue_for_energy_storage.xlsx' :'Ram_2019, DEA_2019',
           'Global Energy System based on 100% Renewable Energy, Energywatchgroup/LTU University, 2019, Danish Energy Agency, technology_data_catalogue_for_energy_storage.xlsx, Note K.' :'Ram_2019, DEA_2019',
           'TODO':'govUK',
           'Viswanathan_2022': 'Viswanathan_2022',
           'Georgiou_2018': 'Georgiou_2018',
}

# Solar thermal collector decentral & 270 & m$^{2}$ & 1.3 & 20 & variable & \cite{Henning20141003} \\
# Solar thermal collector central & 140 & m$^{2}$ & 1.4 & 20 & variable & \cite{Henning20141003} \\
# Building retrofitting\tnote{f} & see text &  & 1 & 50 & 1 & \cite{Henning20141003,PalzerThesis} \\
# High-density district heating network\tnote{f} & 220 & kW\th  & 1 & 40  & 1 & \cite{IEESWV} \\
# Gas distribution network\tnote{f} & 387 & kW\th & 2 & 40 & 1 & based on \cite{bnetza2017} \\

for technology in technologies:
    if idx[technology,'FOM'] in costs.index:
        FOM = str(round(costs.loc[idx[technology,'FOM'],'value'],1))
    else:
        FOM= ' '
    if idx[technology,'lifetime'] in costs.index:
        lifetime = str(int(costs.loc[idx[technology,'lifetime'],'value']))
    else:
        lifetime= ' '
    if idx[technology,'investment'] in costs.index:
        investment = str(int(costs.loc[idx[technology,'investment'],'value']/1000))
    else:
        investment= ' '
    if idx[technology,'efficiency'] in costs.index and technology not in ['onwind',
          'offwind', 'central gas CHP', 'biomass CHP', 'battery storage',
          'home battery storage', 'central coal CHP'
          'hydrogen storage underground', 'hydrogen storage tank type 1',
          'central water tank storage', 'decentral water tank storage',
          'decentral air-sourced heat pump', 'central ground-sourced heat pump',
          'decentral ground-sourced heat pump']:

        efficiency = str(round(costs.loc[idx[technology,'efficiency'],'value'],2))
    else:
        efficiency= ' '
    if technology not in ['water tank charger', 'hydro', 'ror', 'PHS',
                          'electrolysis', 'fuel cell', 'decentral water tank storage']:
        source = costs.loc[idx[technology,'lifetime'],'source']
    elif technology == 'decentral water tank storage':
        source = costs.loc[idx[technology,'investment'],'source']
    else:
        source = costs.loc[idx[technology,'efficiency'],'source']
    if technology == 'water tank charger':
       file.write(' ' + name[technology]
        + ' & ' + investment
        + ' & ' + FOM
        + ' & ' + lifetime
        + ' & ' + efficiency
        + ' & ' + ' \\' + ' ')
    else:        
        file.write(' ' + name[technology] 
        + ' & ' + investment
        + ' & ' + FOM
        + ' & ' + lifetime
        + ' & ' + efficiency
        + ' & ' + ' \\' + 'cite{' + dic_ref[source.split(sep=",")[0]] + '} ')

    file.write('\\')
    file.write('\\')
file.close()

#%%
"""
Table including costs as a function of years
"""
years=np.arange(2020,2055,5)
filename='table_costs.tex'
file = open(filename, 'w')
technologies=[t for t in technologies if t not in ['water tank charger']]
dic_units={'EUR/kWel':'\EUR/kW$_{el}$',
           'EUR/kWth':'\EUR/kW$_{th}$',
           'EUR/kWH2':'\EUR/kW$_{H_2}$',
           'EUR/kW_CH4':'\EUR/kW$_{CH4}$',
           'EUR/kWCH4':'\EUR/kW$_{CH4}$',
           'EUR/kWhth':'\EUR/kWh$_{th}$',
           'EUR/(tCO2/a)': '\EUR/(tCO$_2$/a)',
           'EUR/(tCO2/h)' :'\EUR/(tCO$_2$/h)',
           'EUR/m3':'\EUR/m$^3$',
           'EUR/MW/km':'\EUR/MWkm',
           'EUR/MW':'\EUR/MW',
           'USD/kWel':'USD/kW$_{el}$',
           'USD/kWh':'USD/kWh',
           'EUR/kWh': '\EUR/kWh',
           'EUR/kW': '\EUR/kW',
           'EUR/kW_e':'\EUR/kW$_{el}$',
           'EUR/kW_th - heat output':'\EUR/kW$_{th}$',
           'EUR/kW_th': '\EUR/kW$_{th}$',
           'EUR/kWhCapacity': '\EUR/kWh',
           'EUR/kW_th excluding drive energy': '\EUR/kW$_{th}$',
           'EUR/kW_FT/year':'\EUR/kW$_{FT}$/a',
           'EUR/kW_MeOH':'\EUR/kW$_{MeOH}$'
           }



for technology in technologies:
    file.write(' ' +name[technology] + ' & ')
    file.write(dic_units[costs.loc[idx[technology,'investment'],'unit']]+ ' & ' )

    for year in years:
        costs_year = pd.read_csv('../outputs/costs_' + str(year) +'.csv',index_col=list(range(2))).sort_index()
        if technology in ['hydrogen storage underground', 'central water tank storage']:
            file.write(str(round(costs_year.loc[idx[technology,'investment'],'value'],1))+ ' & ' )
        else:
            file.write(str(int(costs_year.loc[idx[technology,'investment'],'value']))+ ' & ' )

    if technology not in ['water tank charger', 'hydro', 'ror', 'PHS', 'decentral water tank storage']:
    # water tank charger has no lifetime, hydro reference for lifetime
    # is IEA2011, but for cost is DIW
        source = costs.loc[idx[technology,'lifetime'],'source']
    elif technology == 'decentral water tank storage':
        source = costs.loc[idx[technology,'investment'],'source']
    else:
        source = costs.loc[idx[technology,'efficiency'],'source']
    if technology == 'water tank charger':
        file.write( ' \\' + ' ')
    else:
        file.write( ' \\' + 'cite{' + dic_ref[source]+ '} ')
    file.write('\\')
    file.write('\\')
file.close()

#%%
"""
Table including fuel characteristics
"""

filename='table_fuels.tex'
file = open(filename, 'w')
for fuel in [ 'coal', 'lignite', 'gas', 'oil','nuclear', 'solid biomass']:
    if idx[fuel,'fuel'] in costs.index:
        cost = str(round(costs.loc[idx[fuel,'fuel'],'value'],1))
        source1 = costs.loc[idx[fuel,'fuel'],'source']
    else:
        cost = ' '

    if idx[fuel,'CO2 intensity'] in costs.index:
        emissions = str(round(costs.loc[idx[fuel,'CO2 intensity'],'value'],3))
        source2 = costs.loc[idx[fuel,'CO2 intensity'],'source']
    else:
        emissions = ' '
    if fuel not in ['nuclear', 'solid biomass','gas', 'oil','digestible biomass','biogas'] :
        file.write(' ' + fuel
                   + ' & ' +  cost
                   + ' & ' +
                   ' \\' + 'cite{' + dic_ref[source1]+ '} '
                   + ' & ' +  emissions
                   + ' & ' +
                   ' \\' + 'cite{' + dic_ref[source2]+ '} ')
    else:
       file.write(' ' + fuel
                   + ' & ' +  cost
                   + ' & ' +
                   ' \\' + 'cite{' + dic_ref[source1]+ '} '
                   + ' & ' +  str(0)
                   + ' & ' +
                   ' ')
    file.write('\\')
    file.write('\\')
file.close()
