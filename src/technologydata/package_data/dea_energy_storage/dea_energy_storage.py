# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""
Data parser for the DEA energy storage data set.

How to run:
    From the repository root, execute:
        python src/technologydata/package_data/dea_energy_storage/dea_energy_storage.py

Configuration options (command-line arguments):
    --num_digits <int>         Number of significant digits to round the values. Default: 4
    --store_source             Store the source object on the Wayback Machine. Default: False
    --filter_params            Filter the parameters stored to technologies.json. Default: False

Example:
    python src/technologydata/package_data/dea_energy_storage/dea_energy_storage.py --num_digits 3 --store_source --filter_params

"""

import argparse
import logging
import pathlib
import re
import typing

import pandas as pd
import pydantic

from technologydata import (
    Commons,
    Parameter,
    Source,
    SourceCollection,
    Technology,
    TechnologyCollection,
)

path_cwd = pathlib.Path.cwd()

logger = logging.getLogger(__name__)


def drop_invalid_rows(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and filter a DataFrame by removing rows with invalid or incomplete data.

    This function performs multiple validation checks to ensure data quality:
    - Removes rows with None or NaN values in critical columns
    - Removes rows with empty or whitespace-only strings
    - Filters rows based on specific data integrity criteria
    - Discards rows where 'val' column contains comparator symbols or non-numeric values

    Parameters
    ----------
    dataframe : pd.DataFrame
        The input DataFrame to be cleaned and validated.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with invalid rows removed, maintaining data integrity.

    Notes
    -----
    Validation criteria include:
    - Non-empty 'Technology', 'par', and 'val' columns
    - 'year' column containing a valid 4-digit year
    - 'val' column containing only numeric values (no comparator symbols)

    """
    # Create a copy to avoid modifying the original DataFrame
    df_cleaned = dataframe.copy()

    # Validate column existence
    required_columns = ["Technology", "par", "val", "year"]
    missing_columns = [col for col in required_columns if col not in df_cleaned.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Remove rows with None or NaN values in critical columns
    df_cleaned = df_cleaned.dropna(subset=required_columns)

    # Remove rows with empty or whitespace-only strings
    for column in required_columns:
        df_cleaned = df_cleaned[df_cleaned[column].astype(str).str.strip() != ""]

    # Filter rows with valid year (4 consecutive digits)
    df_cleaned = df_cleaned[
        df_cleaned["year"].astype(str).str.contains(r"\d{4}", regex=True)
    ]

    # Remove rows with comparator symbols or without digits in 'val' column
    df_cleaned = df_cleaned[
        (~df_cleaned["val"].astype(str).str.contains(r"[<>≤≥]", regex=True))
        & (df_cleaned["val"].astype(str).str.contains(r"\d", regex=True))
    ]

    return df_cleaned


@pydantic.validate_call
def clean_parameter_string(text_string: str) -> str:
    """
    Remove any string between [] or [), any leading hyphen or double quotes from the input string. Lower-case all.

    Parameters
    ----------
    text_string : str
        input string to be cleaned.

    Returns
    -------
    str
        cleaned string with [] and [) and leading hyphen or double quotes removed.

    Examples
    --------
    >>> clean_parameter_string("- Charge efficiency [%]")
    charge efficiency
    >>> clean_parameter_string("Energy storage capacity for one unit [MWh)")
    energy storage capacity for one unit

    """
    # Remove leading hyphen
    text_string = text_string.lstrip("-")

    # Remove content inside square brackets including the brackets themselves
    result = re.sub(r"\[.*?\]", "", text_string)

    # Remove content inside square bracket and parenthesis including the brackets/parenthesis themselves
    result = re.sub(r"\[.*?\)", "", result)

    # Remove extra spaces resulting from the removal and set all to lower case
    result = re.sub(r"\s+", " ", result).strip().casefold()

    return result


@pydantic.validate_call
def clean_technology_string(tech_str: str) -> str:
    """
    Clean a technology string by removing numeric patterns and standardizing case.

    This function pre-processes technology-related strings by:
    - Removing three-digit numeric patterns (with optional letter)
    - Stripping leading and trailing whitespace
    - Converting to lowercase for case-insensitive comparison

    Parameters
    ----------
    tech_str : str
        Input technology string to be cleaned.

    Returns
    -------
    str
        Cleaned technology string with:
        - Numeric patterns (like '123' or '456a') removed
        - Whitespace stripped
        - Converted to lowercase

    Raises
    ------
    Exception
        If string conversion or processing fails, logs the error and returns the original input.

    Examples
    --------
    >>> clean_technology_string("143a Rock-based Carnot battery")
    rock-based carnot battery
    >>> clean_technology_string("Pit Thermal Energy Storage [PTES]")
    pit thermal energy storage [ptes]

    """
    try:
        # Remove three-digit patterns or three digits followed by a letter
        return re.sub(r"^(\d{3}[a-zA-Z]?)", "", tech_str.strip()).strip().casefold()
    except Exception as e:
        logger.error(f"Error cleaning technology '{tech_str}': {e}")
        return tech_str


@pydantic.validate_call
def format_val_number(input_value: str, num_decimals: int) -> float | None | typing.Any:
    """
    Parse various number formats into a float value.

    Parameters
    ----------
    input_value : str
        The input number in different formats, such as:
        - Scientific notation with "x10^": e.g., "2.84x10^23"
        - Numbers with commas as decimal separators: e.g., "1,1"
    num_decimals : int
        Number of decimals

    Returns
    -------
    float
        The parsed numerical value as a float.

    Raises
    ------
    ValueError
        If the input cannot be parsed into a float.

    Examples
    --------
    >>> format_val_number("1,1")
    1.1
    >>> format_val_numer("2.84×10-27")
    2.84e-27

    """
    s = str(input_value).strip()

    # Handle scientific notation like "2.84x10^23"
    match = re.match(r"([+-]?\d*\.?\d+)×10([+-]?\d+)", s)
    if match:
        base, exponent = match.groups()
        return round(float(base), num_decimals) * (10 ** int(exponent))

    # Replace comma with dot for decimal numbers
    s = s.replace(",", ".")
    try:
        return round(float(s), num_decimals)
    except ValueError:
        raise ValueError(f"Cannot parse number from input: {input_value}")


@pydantic.validate_call
def extract_year(year_str: str) -> int | None:
    """
    Extract the first year (integer) from a given input.

    Parameters
    ----------
    year_str : str
        Input value containing a potential year.

    Returns
    -------
    int, None
        Extracted first year.

    Examples
    --------
    >>> extract_year('uncertainty (2050)')
    2050

    """
    # Extract digits
    digits = re.findall(r"\d+", year_str)

    # Convert to integer
    return int(digits[0]) if digits else None


@pydantic.validate_call
def clean_est_string(est_str: str) -> str:
    """
    Casefold the 'est' string, trim whitespace and replace `ctrl` with `control`.

    Parameters
    ----------
    est_str : str
        The input 'est' string to be cleaned.

    Returns
    -------
    str
        The cleaned 'est' string.

    Examples
    --------
    >>> clean_est_string("Lower")
    lower
    >>> clean_est_string("ctrl")
    control

    """
    if est_str == "ctrl":
        cleaned_str = "control"
    else:
        cleaned_str = est_str.casefold().strip()
    return cleaned_str


def standardize_units(series: pd.Series) -> pd.Series:
    """
    Complete missing units based on parameter names and replace incorrect units.

    Parameters
    ----------
    series : pandas.Series
        A series containing two elements: [par, unit]

    Returns
    -------
    pandas.Series
        Updated series with completed and corrected unit.

    Notes
    -----
    The following substitutions are driven by the `pint` documentation available
     at https://github.com/hgrecco/pint/blob/master/pint/default_en.txt:
    - "pct.": "percent",
    - "m2": "meter**2",
    - "m3": "meter**3",

    """
    par, unit = series

    # Mapping of parameters to their default units
    param_unit_map = {
        "energy storage capacity for one unit": "MWh",
        "typical temperature difference in storage": "K",
        "fixed o&m": "pct./year",
        "lifetime in total number of cycles": "cycles",
        "cycle life": "cycles",
    }

    # Mapping of incorrect units to correct units
    unit_corrections = {
        "pct./period": "percent",
        "⁰C": "C",
        "°C": "C",
        "pct./30sec": "pct.",
        "m2": "meter**2",
        "m3": "meter**3",
        "MWhoutput": "MWh",
        "hot/cold,K": "K",
        "pct.investement": "percent",
        "pct.investment": "percent",
        "tank/": "",
        "pct.": "percent",
    }

    # Complete missing or empty units
    if (not isinstance(unit, str)) or (unit.strip() == ""):
        unit = param_unit_map.get(par, unit)

    # Replace wrong units
    for incorrect, correct in unit_corrections.items():
        if incorrect == unit or incorrect in unit:
            unit = unit.replace(incorrect, correct)

    return pd.Series([par, unit])


def filter_parameters(dataframe: pd.DataFrame, filter_flag: bool) -> pd.DataFrame:
    """
    Filter rows of a DataFrame by allowed technology parameters.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Input DataFrame containing at least a "Technology" column.
    filter_flag : Boolean
        If true, filter parameter `par` column

    Returns
    -------
    pandas.DataFrame
        The filtered DataFrame. The returned
        DataFrame contains only rows where the `par` column is one of:
        "technical lifetime", "fixed o&m", "specific investment", or
        "variable o&m".

    """
    allowed_set = {
        "technical lifetime",
        "fixed o&m",
        "specific investment",
        "variable o&m",
        "charge efficiency",
        "discharge efficiency",
        "capacity",
    }
    print("filter_flag", filter_flag)
    if filter_flag:
        # Filter the DataFrame based on the allowed set
        df_filtered = dataframe[dataframe["par"].isin(allowed_set)].reset_index(
            drop=True
        )
        logger.info(
            f"technologies.json contains a subset of the allowed parameters: {allowed_set}."
        )
    else:
        # Return the original DataFrame if filter_flag is False
        df_filtered = dataframe
        logger.info("All parameters are outputted to technologies.json")
    return df_filtered


def build_technology_collection(
    dataframe: pd.DataFrame,
    sources_path: pathlib.Path,
    store_source: bool = False,
    output_schema: bool = False,
) -> TechnologyCollection:
    """
    Compute a collection of technologies from a grouped DataFrame.

    Processes input DataFrame by grouping technologies and extracting their parameters,
    creating Technology instances for each unique group.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Input DataFrame containing technology parameters.
        Expected columns include:
        - 'est': Estimation or case identifier
        - 'year': Year of the technology
        - 'ws': Workspace or technology identifier
        - 'Technology': Detailed technology name
        - 'par': Parameter name
        - 'val': Parameter value
        - 'unit': Parameter units
    sources_path: pathlib.Path
        Output path for storing the SourceCollection object
    store_source: Optional[bool]
        Flag to decide whether to store the source object on the Wayback Machine. Default False.
    output_schema : Optional[bool]
        Flag to decide whether to export the source collection schema. Default False.

    Returns
    -------
    TechnologyCollection
        A collection of Technology instances, each representing a unique
        technology group with its associated parameters.

    Notes
    -----
    - The function groups the DataFrame by 'est', 'year', 'ws', and 'Technology'
    - For each group, it creates a dictionary of Parameters
    - Each Technology is instantiated with group-specific attributes

    """
    list_techs = []

    if store_source:
        source = Source(
            title="Technology Data for Energy storage (May 2025)",
            authors="Danish Energy Agency",
            url="https://ens.dk/media/6589/download",
            url_date="2025-10-08 09:24:00",
        )
        source.ensure_in_wayback()
        sources = SourceCollection(sources=[source])
        sources.to_json(sources_path, output_schema=output_schema)
    else:
        sources = SourceCollection.from_json(sources_path)

    for (est, year, ws, technology_name), group in dataframe.groupby(
        ["est", "year", "ws", "Technology"]
    ):
        parameters = {}
        for _, row in group.iterrows():
            parameters[row["par"]] = Parameter(
                magnitude=row["val"],
                units=row["unit"],
                sources=sources,
                provenance="Parsed from Excel file",
            )
        list_techs.append(
            Technology(
                name=ws,
                region="EU",
                year=year,
                parameters=parameters,
                case=est,
                detailed_technology=technology_name,
            )
        )
    return TechnologyCollection(technologies=list_techs)


@pydantic.validate_call
def parse_input_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments containing:
        - Number of significant digits
        - Store source flag

    """
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Parse the DEA technology storage dataset",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Define arguments
    parser.add_argument(
        "--num_digits",
        type=int,
        default=4,
        help="Name of significant digits to round the values. ",
    )

    parser.add_argument(
        "--store_source",
        action="store_true",
        help="store_source, store the source object on the wayback machine. Default: false",
    )

    parser.add_argument(
        "--filter_params",
        action="store_true",
        help="filter_params. Filter the parameters stored to technologies.json. Default: false",
    )

    parser.add_argument(
        "--export_schema",
        action="store_true",
        help="export_schema. Export the Source/TechnologyCollection schemas. Default: false",
    )

    # Parse arguments
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    # Parse input arguments
    input_args = parse_input_arguments()
    logger.info("Command line arguments parsed.")

    # Read the raw data
    dea_energy_storage_file_path = pathlib.Path(
        path_cwd,
        "src",
        "technologydata",
        "package_data",
        "raw",
        "Technology_datasheet_for_energy_storage.xlsx",
    )

    dea_energy_storage_df = pd.read_excel(
        dea_energy_storage_file_path,
        sheet_name="alldata_flat",
        engine="calamine",
        dtype=str,
    )
    logger.info("Input file read-in.")

    # Drop unnecessary rows
    cleaned_df = drop_invalid_rows(dea_energy_storage_df)
    logger.info("Unnecessary rows dropped.")

    # Clean technology (Technology) column
    cleaned_df["Technology"] = cleaned_df["Technology"].apply(clean_technology_string)

    # Clean ws column
    cleaned_df["ws"] = cleaned_df["ws"].apply(clean_technology_string)
    logger.info("`Technology` and `ws` cleaned.")

    # Clean year column
    cleaned_df["year"] = cleaned_df["year"].apply(extract_year)
    logger.info("`year` column cleaned.")

    # Clean parameter (par) column
    cleaned_df["par"] = cleaned_df["par"].apply(clean_parameter_string)
    logger.info("`par` column cleaned.")

    # Complete missing units based on parameter names and replace incorrect units.
    cleaned_df[["par", "unit"]] = cleaned_df[["par", "unit"]].apply(
        standardize_units, axis=1
    )
    logger.info("Missing units added and wrong units replaced.")

    # Include priceyear in unit if applicable
    cleaned_df["unit"] = cleaned_df.apply(
        lambda row: Commons.update_unit_with_currency_year(
            row["unit"], row["priceyear"]
        ),
        axis=1,
    )
    logger.info("`priceyear` included in `unit` column.")

    # Format value (val) column
    cleaned_df["val"] = cleaned_df["val"].apply(
        lambda x: format_val_number(x, input_args.num_digits)
    )
    logger.info("`val` column formatted.")

    # Replace "MEUR_2020" with "EUR_2020" and multiply val by 1_000_000
    mask_meur = cleaned_df["unit"].str.contains("MEUR_2020")
    cleaned_df.loc[mask_meur, "unit"] = cleaned_df.loc[mask_meur, "unit"].str.replace(
        "MEUR_2020", "EUR_2020"
    )
    cleaned_df.loc[mask_meur, "val"] = (
        cleaned_df.loc[mask_meur, "val"] * 1_000_000.0
    ).round(input_args.num_digits)

    # Replace "kEUR_2020" with "EUR_2020" and multiply val by 1_000
    mask_lower_keur = cleaned_df["unit"].str.contains("kEUR_2020")
    cleaned_df.loc[mask_lower_keur, "unit"] = cleaned_df.loc[
        mask_lower_keur, "unit"
    ].str.replace("kEUR_2020", "EUR_2020")
    cleaned_df.loc[mask_lower_keur, "val"] = (
        cleaned_df.loc[mask_lower_keur, "val"] * 1_000.0
    ).round(input_args.num_digits)

    # Replace "KEUR_2020" with "EUR_2020" and multiply val by 1_000
    mask_upper_keur = cleaned_df["unit"].str.contains("KEUR_2020")
    cleaned_df.loc[mask_upper_keur, "unit"] = cleaned_df.loc[
        mask_upper_keur, "unit"
    ].str.replace("KEUR_2020", "EUR_2020")
    cleaned_df.loc[mask_upper_keur, "val"] = (
        cleaned_df.loc[mask_upper_keur, "val"] * 1_000.0
    ).round(input_args.num_digits)

    # Replace "mol/s/m/MPa1/2" with "mol/s/m/Pa" and multiply val by 1_000_000
    mask_mols = cleaned_df["unit"].str.contains("mol/s/m/MPa1/2")
    cleaned_df.loc[mask_mols, "unit"] = cleaned_df.loc[mask_mols, "unit"].str.replace(
        "mol/s/m/MPa1/2", "mol/s/m/Pa"
    )
    cleaned_df.loc[mask_mols, "val"] = cleaned_df.loc[mask_mols, "val"] * 1_000_000.0

    # Replace, in column `par`, `energy storage capacity for one unit` and `tank volume of example` with `capacity`
    mask_capacity = cleaned_df["par"].isin(
        [
            "energy storage capacity for one unit",
            "tank volume of example",
        ]
    )
    cleaned_df.loc[mask_capacity, "par"] = "capacity"

    # Clean est column
    cleaned_df["est"] = cleaned_df["est"].apply(clean_est_string)
    logger.info("`est` column cleaned.")

    # Drop unnecessary columns
    columns_to_drop = ["cat", "priceyear", "ref", "note"]
    cleaned_df = cleaned_df.drop(columns=columns_to_drop, errors="ignore")
    logger.info("Unnecessary columns dropped.")

    filtered_df = filter_parameters(cleaned_df, input_args.filter_params)

    # Build TechnologyCollection
    dea_storage_path = pathlib.Path(
        path_cwd,
        "src",
        "technologydata",
        "package_data",
        "dea_energy_storage",
    )
    output_technologies_path = pathlib.Path(
        dea_storage_path,
        "technologies.json",
    )
    output_sources_path = pathlib.Path(
        dea_storage_path,
        "sources.json",
    )

    tech_col = build_technology_collection(
        filtered_df,
        output_sources_path,
        store_source=input_args.store_source,
        output_schema=input_args.export_schema,
    )
    logger.info("TechnologyCollection object instantiated.")
    tech_col.to_json(output_technologies_path, output_schema=input_args.export_schema)
    logger.info("TechnologyCollection object exported to json.")

    if input_args.export_schema:
        # Move schema files if they exist
        schema_folder = pathlib.Path(
            path_cwd, "src", "technologydata", "package_data", "schemas"
        )
        sources_schema = pathlib.Path(dea_storage_path, "sources.schema.json")
        technologies_schema = pathlib.Path(dea_storage_path, "technologies.schema.json")

        schema_folder.mkdir(parents=True, exist_ok=True)

        if sources_schema.exists():
            sources_schema.rename(schema_folder / "sources.schema.json")
            logger.info(f"Moved {sources_schema} to {schema_folder}")
        if technologies_schema.exists():
            technologies_schema.rename(schema_folder / "technologies.schema.json")
            logger.info(f"Moved {technologies_schema} to {schema_folder}")
