# Technology data

This repository should sample and convert the technology data to be used in [PyPSA-EUR](https://github.com/PyPSA/pypsa-eur.git) and [PyPSA-EUR-Sec](https://github.com/PyPSA/pypsa-eur-sec.git).
Technology data for choosen years (e.g. [2020, 2030, 2050]) and different sources are merged in a standarized form with uniform cost years, technology names and units.

## Structure
this repository has the following structure:

* **inputs**: input data in format .csv or .xlsx
* **outputs**: technology data saved as costs_{year}.csv format for defined years.
               In the output costs_{year}.csv are specified:
  * technology (e.g. 'onwind')
  * parameter (e.g. FOM)
  * value (e.g. 1.18)
  * unit (e.g. %/year)
  * source (e.g. DEA, excel_file_name.xlsx)
  * further description (specific assumptions, sheet name if data from DEA)
* **docu**: reports, paper, additional information about the input data, format .pdf
* **scripts** :
  * compile_cost_assumptions.py
    converts input data from multiplte sources to cost_{year}.csv for chosen year. Interpolates data for missing years or calculates the costs at a certain year based on the inflation rate. Technology data from the [Danish Energy Agency Technology Database](https://ens.dk/en/our-services/projections-and-models/technology-data) are preferred.
If data are missing from all sources, these are taken from the old PyPSA cost
assumptions (with a printed warning).
The following parameters can be set at the beginning of the script (should be moved to a config.yaml):
      * years : numpy array of all the years of which an output costs csv should be created
      * rate_inflation : inflation rate (currently: rate_inflation=0.02)
      * path_in : str. of input folder (currently: path_in="../inputs/")
      * solar_utility_from_other : Bool (True/False) if solar utility data is taken from DEA or Vartiaien
      * solar_rooftop_from_other : Bool (True/False) if solar rooftop data is taken from DEA or ETIP
      * h2_from_budischak : Bool (True/False) if fuel cell and electrolyzer efficiencies are taken from DEA or Budischak
      * offwind_no_gridcosts : Bool (True/False) if offshore wind grid connection costs should be removed (they are calculated seperately in PyPSA-EUR)
  * convert_pdf_fraunhofer_to_dataframe.py
  converts table in pdf to csv format for input data. Script can be modified to convert other .pdf sources to .csv format
  * retrieve_data_from_dea.py
  takes up to date technology data from DEA website and saves it in the **input** folder. Optional, also retrieves the documentation of the data into the folder **docu**
* **latex_tables**: .tex files with tables of the cost.csv and 2 python scripts
  * tables_in_latex.py to create .tex files with nice names
  * tables_in_csv.py to create csv files with nice name (which can be used in latex with csv autotabular)


  ## Sources
  * most technologies
   [Danish Energy Agency Technology Database](https://ens.dk/en/our-services/projections-and-models/technology-data)
  * solar utility
  [Vartiaien et. al.](https://onlinelibrary.wiley.com/doi/full/10.1002/pip.3189)
  * solar rooftop
  [The European Technology and Innovation Platform for Photovoltaics](https://etip-pv.eu/)
  * conventional carriers (nuclear, coal, lignite)
  [Lazard](https://www.lazard.com/media/451086/lazards-levelized-cost-of-energy-version-130-vf.pdf)
  * fuel cost
  [Zappe et. al.](https://doi.org/10.1016/j.apenergy.2018.08.109)
  * CO2 intensity
  [Umweltbundesamt - Entwicklung Kohlendioxid](https://www.umweltbundesamt.de/publikationen/entwicklung-der-spezifischen-kohlendioxid-5)
  * gas pipeline costs (potenitally import costs, comparison to other data)
  [Fraunhofer ISE Studie](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/wege-zu-einem-klimaneutralen-energiesystem.html)
