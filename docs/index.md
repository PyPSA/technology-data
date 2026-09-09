<p align="center">
  <img src="assets/logo/technology_data_logo.png" alt="technologydata Header Logo" width="400"/>
</p>

[![PyPI version](https://img.shields.io/pypi/v/technologydata.svg)](https://pypi.python.org/pypi/technologydata)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/open-energy-transition/technology-data/blob/prototype-2/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/technologydata.svg)](https://pypi.python.org/pypi/technologydata)

# technologydata: techno-economic assumptions for energy models

`technologydata` is a Python package that supports the management of techno-economic assumptions for energy system models.
It provides a structured way to store, retrieve, and manipulate data related to various technologies used in energy systems, including unit-ful parameters, currency conversions, inflation adjustment, and temporal modelling.

The goal of this package is to make energy system modelling easier and more efficient, automating common tasks and transformations to reduce errors and allowing for easier data exchange between models.

Techno-economic catalogues are published in different currencies, price years, physical units, heating-value conventions and parameter naming schemes.
Harmonising to combine multiple sources into one energy system model therefore requires a sequence of conversions that is easy to get wrong and rarely recorded.

This package contains a data schema and represents the techno-economic assumpstions as Python objects.
The objects called `Parameter` carry not only values ("magnitude"), but also unit information, currency year, energy carrier, heating-value basis, provenance and bibliographic sources.
Commonly used conversions are available, parameters can be checked for consistency.
Provenance information for modifying parameters is automatically recorded, e.g. to keep track of currency conversion factors or formulas used for calculation.

## Installation

```bash
pip install technologydata
```

Or, using [uv](https://docs.astral.sh/uv/):

```bash
# Standalone
uv pip install technologydata

# As part of a project
uv add install technologydata
```

The package requires Python 3.12 or newer.
The bundled datasets are installed with it; no additional download step is needed.

## Features of the package

- **Unit-ful parameters.**
  Parameters come with units, including currencies an currency years, heating values and support for (energy) carriers.
  Unit conversion is handled automatically and unit-compatability is checked for all operations via [`pint`](https://pint.readthedocs.io/en/stable/).
- **Currency and inflation adjustment.**
  `to_currency()` combines exchange-rate conversion and deflation between price years using official inflation exchange rate data from the World Bank or alternatively the International Monetary Fund via [`pydeflate`](https://jm-rivera.github.io/pydeflate/).
- **Heating-value handling.**
  Parameters can record an `LHV` or `HHV` basis and be converted between them with `change_heating_value()`, using the speicifc energy content depending on the energy carrier.
- **Arithmetic with compatibility checks.**
  Adding parameters with incompatible units, carriers or heating-value bases raises `ValueError` rather than returning a misleading value.
- **Provenance and source tracking.**
  Every `Parameter` carries a provenance and source information.
  Operations and modifications on `Parameter` objects automatically track provenance information for improved reproducibility and insights.
- **Filtering and export.**
  `TechnologyCollection.get()` filters by case-insensitive regular expression on name, region, year, case and detailed technology;
  collections export to `pandas.DataFrame`, CSV and JSON, and reload from JSON.
- **Interpolation and extrapolation.**
  `GrowthModels` allow for easy gap filling when data is missing for particular years.
  Use model fitting and projections with commonly used growth models, including linear, exponential and logistic growth.

## Bundled datasets (Batteries included)

The package already ships with some parsed data catalogues for immediate use.
Prominent examples are:

- Energy Storage Catalogue from the Danish Energy Agency (automatically extracted)
- Annual Technology Baseline 2024 from NREL (manually collected in the old [`technologydata` repository](https://technology-data.readthedocs.io))

The raw files and parsing logic is also shipped alongside to provide the opportunity for verification.

!!! tip
    More packaged data to come!
    Your contributions are also very welcome!

## Example

The Danish Energy Agency storage catalogue provides investment cost of utility-scale lithium-ion battery storage.
NREL's ATB 2024 provides similar technology costs, but their currency units do not align.
The follow examples loads both datasets and harmonises them to the same units for comparison:

```python
import technologydata

# Load data provided by the DEA
dea = technologydata.DataAccessor(data_source="dea_energy_storage", version="v10").load()
# Load data from NREL's ATB2024
atb = DataAccessor(data_source="manual_input_usa", version="v0.13.4", data_path=data).load()

# We get the data specific to battery storage and adjust the currency year
# from 2022 to 2023. By default the inflation adjustment uses World Bank data
us = atb.technologies.get(
    name="battery storage", case="Moderate - Market",
    region=None, year=None, detailed_technology=None,
).to_currency("USD_2023")

# The DEA data is usually valid for the "EU" and the entries carry the region "EU"
# For currency conversion, this is an valid ISO 3166 alpha-3 code that can be used
# for currency conversion and inflation adjustment, so we need to specify a reference country
# like Germany explicitly
eu = dea.technologies.get(
    name="lithium ion battery", case="control",
    region=None, year=None, detailed_technology=None,
).to_currency("USD_2023", overwrite_country="DEU")

for t in eu.technologies:
    if t.year == 2030:
        p = t.parameters["specific investment"].to("USD_2023/kWh")
        print(f"DEA v10  {t.year}: {p.magnitude:7.2f} {p.units}")

for t in us.technologies:
    if t.year == 2030:
        p = t.parameters["investment"]
        print(f"NREL ATB {t.year}: {p.magnitude:7.2f} {p.units}")
```

```text
DEA v10  2030:  351.74 USD_2023 / kilowatt_hour
NREL ATB 2030:  264.23 USD_2023 / kilowatt_hour
```

Each converted `Parameter` retains the sources of the value it came from, so the origin of a number remains traceable after conversion.
The currency conversion uses exchange-rate and deflator series retrieved through [pydeflate](https://github.com/jm-rivera/pydeflate).
Data from the World Bank is used in the example above and is downloaded on first use automatically; the data is cached locally and later calls are much faster.

## Citing

If you use `technologydata` in your research, please cite it as:

```text
technologydata: Data for Energy Systems Models.

The package is available at: https://github.com/open-energy-transition/technology-data/tree/prototype-2.

Authors: Johannes Hampp, Fabrizio Finozzi
```

## License

`technologydata` is released under the [MIT license](home/license.md).

## Contacts

- For **bugs and feature requests**,
  use the [issue tracker](https://github.com/open-energy-transition/technology-data/issues).
- For **contributions**,
  open a pull request on [GitHub](https://github.com/open-energy-transition/technology-data). Ideas, suggestions and problem reports are equally welcome as issues.
