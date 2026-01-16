# Danish Energy Agency Parser Documentation

## Overview

The Danish Energy Agency (DEA) data parser `dea_energy_storage.py` demonstrates a full data-cleaning and transformation pipeline for converting raw tabular data into the `technologydata` schema files `technologies.json` and `sources.json`. The parser is implemented in `src/technologydata/package_data/dea_energy_storage/dea_energy_storage.py`.

## Dataset Description

The original dataset is available at this [link](https://ens.dk/media/6589/download). A full description of the dataset is available at this [link](https://ens.dk/media/6588/download). The raw source file is included in the repository at `src/technologydata/package_data/raw/Technology_datasheet_for_energy_storage.xlsx`.

The dataset is in Excel format, and it includes, under the data sheet `alldata_flat`, a flat table of technology parameters for a range of energy storage technologies. Columns include `Technology`, `ws`, `par` (parameter name), `val` (value), `unit`, `year`, `est` (case/estimate), `priceyear`, plus metadata columns such as `cat`, `ref`, `note`. Rows are individual parameter records (parameter value + unit + context) for technologies and estimation cases.

## Parser description

The parser is articulated in the following steps.

### Command line argument parsing

Function `parse_input_arguments()` defines and parses the command-line arguments:

- `--num_digits` (int, default 4) — number of decimals used when rounding numeric values. The default value is 4.
- `--store_source` (boolean flag) — whether to store the source on the Wayback Machine. The default value is `false`.
- `--filter_params` (boolean flag) — whether to limit exported parameters to a fixed allowed set. The default value is `false`.
- `--export_schema` (boolean flag) — export JSON schema files. The default value is `false`.

### Read the raw data

The script reads the raw data available at `src/technologydata/package_data/raw/Technology_datasheet_for_energy_storage.xlsx`, under sheet `alldata_flat`, in a `pandas` dataframe. It uses `pandas.read_excel(..., engine=calamine, dtype=str)`. All entries are handled as strings initially.

### Data cleaning, validation and dealing with missing/null values

The data cleaning and validation happens with the following steps.

Function `drop_invalid_rows(df)` validates whether required columns are present. It drops rows with missing/null or empty critical fields (`Technology`, `par`, `val`, `year`) and keeps rows where `year` contains a 4-digit year and `val` contains numeric characters and no comparator symbols (`<`, `>`, `≤`, `≥`).

Function `clean_technology_string()` normalizes text fields by removing leading 3-digit numeric codes, trims whitespace and lower-cases the string for consistent matching. It is applied to the columns `Technology` and `ws`. As an example, `clean_technology_string()` converts `151b Hydrogen Storage - LOHC` to `hydrogen storage - lohc`.

Function `extract_year()` extracts the first sequence of digits from the `year` column and converts it to an integer. The column contains in fact entries like `Uncertainty (2050)` (str) which are converted to `2050` (int).

Function `clean_parameter_string()` removes leading hyphens, removes text inside square brackets (units/notes), collapses extra spaces and lower-cases the parameter name. It is applied to the `par` column.

Function `standardize_units()` is applied to columns `par` and `unit`. It completes missing units based on parameter name (e.g., `energy storage capacity for one unit` is mapped to the unit `MWh`) via a parameter-to-unit map. Moreover, it replaces known incorrect unit strings as `⁰C` -> `C` or `m2` to `meter**2`. The unit substitutions are driven by the `pint` documentation available at this [link](https://github.com/hgrecco/pint/blob/master/pint/default_en.txt).

Function `Commons.update_unit_with_currency_year(unit, priceyear)`, if present, appends `priceyear` information to currency units. This is because `technologydata` follows the currency pattern `\b(?P<cu_iso3>[A-Z]{3})_(?P<year>\d{4})\b`, as for example `EUR_2021`.

Function `format_val_number(value, num_decimals)` parses numeric formats including comma decimal separators and scientific notation variants (e.g., `×10`) and converts them to float and rounds them to `num_decimals`.

The parser also applies the following corrections and substitutions:

- Convert `MEUR_2020` and `kEUR_2020`/`KEUR_2020` to `EUR_2020` and scale numeric `val` accordingly (×1e6 or ×1e3).
- Specific unit fixes (example: `mol/s/m/MPa1/2` → `mol/s/m/Pa` with value scaling).
- Certain `par` values (e.g., `energy storage capacity for one unit`, `tank volume of example`) are normalized to `capacity`.

Function `clean_est_string()` normalizes the `est` column by casefolding it and by replacing `ctrl` with `control`.

Function `filter_parameters(df, filter_flag)`, if `filter_flag` is true, keeps only an allowed set of parameters (e.g., `technical lifetime`, `fixed o&m`, `specific investment`, `variable o&m`, `charge efficiency`, `discharge efficiency`, `capacity`). Otherwise returns the full set.

### Populate and export the source and technology collections

Function `build_technology_collection()`:

- if `store_source` is set, constructs a `Source` object for the DEA dataset, calls `ensure_in_wayback()` and writes `sources.json`; otherwise reads an existing `sources.json`.
- groups the cleaned DataFrame by `est`, `year`, `ws`, `Technology`.
- for each group, builds a dictionary of `Parameter` objects (each with `magnitude`, `units`, `sources`, `provenance`).
- creates a `Technology` object for each group, with `name` = `ws`, `detailed_technology` = `Technology`, `year`=`year`, `region` = `EU`, `case` = `est` and collects them into a `TechnologyCollection` object.
- writes the `TechnologyCollection` object to a `technologies.json`.
- if `--export_schema` is used, schema files produced during export are moved to the sub-folder `src/technologydata/package_data/schemas`.

## Running the parser

### Execution instructions

From repository root:

- Basic run: `python src/technologydata/package_data/dea_energy_storage/dea_energy_storage.py`
- Example with options: `--num_digits 3 --store_source --filter_params --export_schema`

### Outputs

The parser generates the following outputs:

- `src/technologydata/package_data/dea_energy_storage/technologies.json`.
- `src/technologydata/package_data/dea_energy_storage/sources.json`.
- Optional schema files moved to `src/technologydata/package_data/schemas` when `--export_schema` is used.
