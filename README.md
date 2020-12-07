# Energy System Technology Data

This script compiles assumptions on energy system technologies (such
as costs, efficiencies, lifetimes, etc.)  for chosen years
(e.g. [2020, 2030, 2050]) from a variety of sources into CSV files to
be read by energy system modelling software. The merged outputs have
standardized cost years, technology names, units and source information.


The outputs are used in
[PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) and
[PyPSA-Eur-Sec](https://github.com/PyPSA/pypsa-eur-sec).

## Structure

This repository has the following structure:

* **inputs**: input data in format .csv or .xlsx
* **outputs**: technology data saved as costs_{year}.csv format for defined years.
               In the output costs_{year}.csv are specified:
  * technology (e.g. 'onwind')
  * parameter (e.g. FOM)
  * value (e.g. 1.18)
  * unit (e.g. %/year)
  * source (e.g. DEA, excel_file_name.xlsx)
  * further description (specific assumptions, sheet name if data from a multi-sheet Excel file)
* **docu**: reports, paper, additional information about the input data, format .pdf
* **scripts** :
  * compile_cost_assumptions.py
    converts input data from multiple sources to cost_{year}.csv for chosen year. Interpolates data for missing years or calculates the costs at a certain year based on the inflation rate. Technology data from the [Danish Energy Agency Technology Database](https://ens.dk/en/our-services/projections-and-models/technology-data) are preferred.
If data are missing from all sources, these are taken from the old PyPSA cost
assumptions (with a printed warning).
The following parameters can be set at the beginning of the script (should be moved to a config.yaml):
      * years : numpy array of all the years of which an output costs csv should be created
      * rate_inflation : inflation rate (currently: rate_inflation=0.02)
      * solar_utility_from_other : Bool (True/False) if solar utility data is taken from DEA or Vartiaien
      * solar_rooftop_from_other : Bool (True/False) if solar rooftop data is taken from DEA or ETIP
      * h2_from_budischak : Bool (True/False) if fuel cell and electrolyzer efficiencies are taken from DEA or Budischak
      * offwind_no_gridcosts : Bool (True/False) if offshore wind grid connection costs should be removed (they are calculated seperately in PyPSA-Eur)
  * convert_pdf_fraunhofer_to_dataframe.py
  converts table from Fraunhofer ISE report in pdf to csv format for input data. Script can be modified to convert other .pdf sources to .csv format
  * retrieve_data_from_dea.py
  downloads up-to-date technology data from DEA website and saves it in the **input** folder. Optional, also retrieves the documentation of the data into the folder **docu**
* **latex_tables**: .tex files with tables of the cost.csv and 2 python scripts
  * tables_in_latex.py to create .tex files with nice names
  * tables_in_csv.py to create csv files with nice name (which can be used in latex with csv autotabular)


## Sources

* most technologies
   [Danish Energy Agency Technology Database](https://ens.dk/en/our-services/projections-and-models/technology-data)
* utility-scale solar PV
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
* Older data collected by the PyPSA team from a variety of sources


## Library dependencies

The current script has been tested with pandas up to version 1.1.0.

The PDF to CSV conversion requires the package `tabula`. For conda users, install with `conda install -c conda-forge tabula-py`.


## Release Notes

### technology-data 0.1.0 (21st August 2020)

This is the first release to coincide with the release of [PyPSA-Eur-Sec](https://github.com/PyPSA/pypsa-eur-sec) 0.1.0.

### Release process

* Update release notes.
* Update version number in `config.yaml`.
* Make a `git commit`.
* Tag the release with `git tag v0.x.x`, `git push`, `git push --tags`.
* Make a [GitHub release](https://github.com/PyPSA/technology-data/releases), which triggers archiving by [zenodo.org](https://zenodo.org/).

## Licence

Copyright 2019-2020 Marta Victoria (Aarhus University), Kun Zhu
(Aarhus University), Elisabeth Zeyen (KIT), Tom Brown (KIT)

The code in `scripts/` is released as free software under the
[GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html), see LICENSE.txt.
However, different licenses and terms of use may apply to the various
input data.
