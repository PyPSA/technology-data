..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

.. _structure:

##########################################
Structure of the repository
##########################################

This repository has the following structure:

-  **inputs**: input data in format .csv or .xlsx

-  **outputs**: technology data saved as ``costs_{year}.csv`` format for defined years. In the output ``costs_{year}.csv`` are specified

	* technology (e.g. 'onwind')
	* parameter (e.g. FOM)
	* value (e.g. 1.18)
	* unit (e.g. %/year)
	* source (e.g. DEA, excel_file_name.xlsx)
	* further description (specific assumptions, sheet name if data from a multi-sheet Excel file)

-  **config**:

	.. literalinclude:: ../config.yaml
	   :language: yaml
	   :lines: 4-24

 the following parameters can be set in the ``config.yaml``

	* years : numpy array of all the years of which an output costs csv should be created
	* rate_inflation : inflation rate (currently: rate_inflation=0.02)
	* solar_utility_from_vartiaien : Bool (True/False) if solar utility data is taken from DEA or Vartiaien
	* solar_rooftop_from_etip : Bool (True/False) if solar rooftop data is taken from DEA or ETIP
	* h2_from_budischak : Bool (True/False) if fuel cell and electrolyzer efficiencies are taken from DEA or Budischak
	* offwind_no_gridcosts : Bool (True/False) if offshore wind grid connection costs should be removed (they are calculated seperately in PyPSA-Eur)

-  **scripts** :

	* :mod:`compile_cost_assumptions.py`  converts input data from multiple sources to ``cost_{year}.csv`` for chosen year. Interpolates data for missing years or calculates the costs at a certain year based on the inflation rate. Technology data from the `Danish Energy Agency <https://github.com/PyPSA/technology-data>`_ are preferred. If data are missing from all sources, these are taken from the old PyPSA cost assumptions (with a printed warning).
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


