# Manual Input USA Parser Documentation

## Overview

The Manual Input USA data parser `manual_input_usa.py` demonstrates a data-cleaning and transformation pipeline for converting manually curated, USA-specific tabular data into the `technologydata` schema files `technologies.json` and `sources.json`. The parser is implemented in `src/technologydata/parsers/manual_input_usa/manual_input_usa.py`.

## Dataset Description

The original dataset is a manually curated CSV file containing USA-specific technology parameters available at this [link](https://github.com/PyPSA/technology-data/blob/v0.13.4/inputs/US/manual_input_usa.csv). The raw source file is included in the repository at `src/technologydata/parsers/raw/manual_input_usa.csv`.

The dataset is in CSV format and includes a flat table of technology parameters for various energy technologies relevant to the USA context. Columns include `technology`, `parameter`, `year`, `value`, `unit`, `currency_year`, `source`, `further_description`, `financial_case`, and `scenario`. Rows are individual parameter records (parameter value + unit + context) for technologies with different scenarios and financial cases.

## Parser description

The parser is articulated in the following steps.

### Command line argument parsing

Function `CommonsParser.parse_input_arguments()` defines and parses the command-line arguments:

- `--num_digits` (int, default 4) — number of decimals used when rounding numeric values. The default value is 4.
- `--store_source` (boolean flag) — whether to store the source on the Wayback Machine. The default value is `false`.

### Read the raw data

The script reads the raw data available at `src/technologydata/parsers/raw/manual_input_usa.csv` in a `pandas` dataframe. It uses `pandas.read_csv(..., dtype=str, na_values="None")`. All entries are handled as strings initially except for the `value` column which is converted to float.

### Data cleaning, validation and dealing with missing/null values

The data cleaning and validation happens with the following steps.

Function `extract_units_carriers_heating_value()` extracts standardized units, carriers, and heating values from input unit strings. This function maps complex unit representations to simplified unit, carrier, and heating value combinations using a predefined dictionary of special patterns. Examples include:

- `USD_2022/MW_FT` → unit: `USD_2022/MW`, carrier: `1/FT`, heating_value: `1/LHV`
- `MWh_H2/MWh_FT` → unit: `MWh/MWh`, carrier: `H2/FT`, heating_value: `LHV`
- `MWh_el/MWh_FT` → unit: `MWh/MWh`, carrier: `el/FT`, heating_value: `LHV`
- `t_CO2/MWh_FT` → unit: `t/MWh`, carrier: `CO2/FT`, heating_value: `LHV`
- `USD_2022/kWh_H2` → unit: `USD_2022/kWh`, carrier: `1/H2`, heating_value: `LHV`
- `USD_2023/t_CO2/h` → unit: `USD_2023/t/h`, carrier: `1/CO2`, heating_value: `None`
- `MWh_el/t_CO2` → unit: `MWh/t`, carrier: `el/CO2`, heating_value: `LHV`
- `MWh_th/t_CO2` → unit: `MWh/t`, carrier: `thermal/CO2`, heating_value: `LHV`

The parser also fills missing values in the `scenario` column with `"not_available"`.

The parser applies the following unit conversions:

- Convert `per unit` to `%` and multiply the corresponding `value` by 100.0, rounding to `num_digits` decimals.

Function `Commons.update_unit_with_currency_year(unit, currency_year)` appends `currency_year` information to currency units when present. This is because `technologydata` follows the currency pattern `\b(?P<cu_iso3>[A-Z]{3})_(?P<year>\d{4})\b`, as for example `USD_2022`.

### Populate and export the source and technology collections

Function `build_technology_collection()`:

- if `store_source` is set, constructs a `Source` object for the manual input USA dataset, calls `ensure_in_wayback()` and writes `sources.json`; otherwise reads an existing `sources.json`.
- groups the cleaned DataFrame by `scenario`, `year`, `technology`.
- for each group, builds a dictionary of `Parameter` objects (each with `magnitude`, `sources`, and optionally `carrier`, `heating_value`, `units`, `note`).
- captures the `financial_case` value from rows within each group to combine with `scenario`.
- creates a `case` value by combining `scenario` and `financial_case` in the format `"{scenario} - {financial_case}"` when `financial_case` is present; otherwise uses `scenario` alone.
- creates a `Technology` object for each group, with `name` = `technology`, `detailed_technology` = `technology`, `year` = `year`, `region` = `USA`, `case` = combined case value, and collects them into a `TechnologyCollection` object.
- writes the `TechnologyCollection` object to a `technologies.json`.

## Running the parser

### Execution instructions

From repository root:

- Basic run: `python src/technologydata/parsers/manual_input_usa/manual_input_usa.py`
- Example with options: `--num_digits 3 --store_source`

### Outputs

The parser generates the following outputs:

- `src/technologydata/parsers/manual_input_usa/technologies.json`.
- `src/technologydata/parsers/manual_input_usa/sources.json`.
- Optional schema files moved to `src/technologydata/parsers/schemas` when `--export_schema` is used.
