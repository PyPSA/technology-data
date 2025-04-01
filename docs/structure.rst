..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

.. _structure:

##########################################
Structure of the repository
##########################################

This repository has the following structure:

-  **inputs**: input data in format .csv, .xlsx or .parquet

-  **outputs**: technology data saved as ``costs_{year}.csv`` format for defined years. In the output ``costs_{year}.csv`` are specified

    * technology (e.g. 'onwind')
    * parameter (e.g. FOM)
    * value (e.g. 1.18)
    * unit (e.g. %/year)
    * source (e.g. DEA, excel_file_name.xlsx)
    * further description (specific assumptions, sheet name if data from a multi-sheet Excel file)
    * currency_year (year used for adjusting economic values to reflect current purchasing power)

-  **outputs/US**: US-specific technology data saved as ``US/costs_{year}.csv`` format for defined years. In the outputs ``US/costs_{year}.csv`` are specified

    * technology (e.g. 'onwind')
    * parameter (e.g. FOM)
    * value (e.g. 1.18)
    * unit (e.g. %/year)
    * source (e.g. DEA, excel_file_name.xlsx)
    * further description (specific assumptions, sheet name if data from a multi-sheet Excel file)
    * currency_year (year used for adjusting economic values to reflect current purchasing power)
    * financial_case (financial assumptions for the definition of the cost of capital)
    * scenario (technology innovation scenario)

-  **config**:

	.. literalinclude:: ../config.yaml
	   :language: yaml
	   :lines: 4-41

 The following parameters can be set in the ``config.yaml``

    * years : numpy array of all the years of which an output costs csv should be created
    * nrel_atb_input_years : list of years that define the source files
    * nrel_atb_columns_to_keep : list of columns to use from the NREL/ATB source files
    * nrel_atb_core_metric_parameter_to_keep : list of parameters to use from the NREL/ATB source files
    * nrel_atb_technology_to_remove : list of technologies that should be excluded from NREL/ATB
    * nrel_atb_source_link : source url for the NREL/ATB data used
    * nrel_atb_further_description : Meaning of "scenario" and "financial case"
    * expectation : tech data uncertainty, possible options [None, "optimist", "pessimist"]
    * eur_year : year for EUR outputs
    * solar_utility_from_vartiaien : Bool (True/False) if solar utility data is taken from DEA or Vartiaien
    * solar_rooftop_from_etip : Bool (True/False) if solar rooftop data is taken from DEA or ETIP
    * h2_from_budischak : Bool (True/False) add fuel cell/electrolysis efficiencies from Budischak
    * ewg_home_battery: Bool (True/False) add home battery data derived from DEA data and EWG study
    * add_data: Bool (True/False) add storage data mainly from PNNL
    * approx_beyond_2030: ["geometric_series"] or ["same_as_2030"]
    * offwind_no_gridcosts : Bool (True/False) if offshore wind grid connection costs should be removed (they are calculated seperately in PyPSA-Eur)
    * salinity : estimated in PSU (Practical Salinity Unit) = kg/m^3
    * ndigits : number of significant digits

-  **scripts** :
    * :mod:`compile_cost_assumptions.py` converts input data from multiple sources to ``cost_{year}.csv`` for chosen year. Interpolates data for missing years or calculates the costs at a certain year based on the inflation rate. Technology data from the `Danish Energy Agency <https://github.com/PyPSA/technology-data>`_ are preferred. If data are missing from all sources, these are taken from the old PyPSA cost assumptions (with a printed warning).
    * :mod:`compile_cost_assumptions_usa.py` converts input data from NREL/ATB to ``US/cost_{year}.csv`` for chosen year. It starts from the cost assumptions files produced by `compile_cost_assumptions.py`. All technology-parameter pairs present in the NREL/ATB input data are updated. Those not present in NREL/ATB are left untouched.
    * :mod:`convert_pdf_fraunhofer_to_dataframe.py` converts table from Fraunhofer ISE report in pdf to csv format for input data. Script can be modified to convert other .pdf sources to .csv format
    * :mod:`retrieve_data_from_dea.py` downloads up-to-date technology data from DEA website and saves it in the **input** folder. Optional, also retrieves the documentation of the data into the folder **docu**

-  **docu**: reports, paper, additional information about the input data, format .pdf

-  **latex_tables**: .tex files with tables of the cost.csv and two python scripts

	* ``tables_in_latex.py`` to create .tex files with nice names
	* ``tables_in_csv.py`` to create csv files with nice name (which can be used in latex with csv autotabular)


The data licences and sources are given in the following table.


.. csv-table::
   :header-rows: 1
   :file: data.csv


