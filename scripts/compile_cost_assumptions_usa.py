# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

"""
Script creates cost csv for chosen years concatenating US-specific cost assumptions from NREL ATB.
The input files are in parquet format and can be downloaded from https://data.openei.org/s3_viewer?bucket=oedi-data-lake&prefix=ATB%2Felectricity%2Fparquet%2F
"""

import logging
import pathlib

import numpy as np
import pandas as pd
from _helpers import (
    adjust_for_inflation,
    configure_logging,
    mock_snakemake,
    prepare_inflation_rate,
)

logger = logging.getLogger(__name__)


def get_conversion_dictionary(flag: str) -> dict:
    """
    The function provides conversion dictionaries.

    Parameters
    ----------
    flag : str
        flag for setting the returned dictionary

    Returns
    -------
    Dictionary
        conversion dictionary such that
            flag == parameter: it returns a conversion dictionary that renames the parameter values to the PyPSA standard
            flag == pypsa_technology_name: it returns a conversion dictionary to rename the technology names to the PyPSA nomenclature
            flag == atb_technology_name: it returns a conversion dictionary to align the atb_e_2022 and atb_e_2024 nomenclatures
            flag == output_column: it returns a conversion dictionary that renames the column names of the cost dataframe
    """

    if flag.casefold() == "parameter":
        return {
            "CAPEX": "investment",
            "Fixed O&M": "FOM",
            "Variable O&M": "VOM",
            "Fuel": "fuel",
            "Additional OCC": "investment",
            "WACC Real": "discount rate",
        }
    elif flag.casefold() == "pypsa_technology_name":
        return {
            "Coal-new -> 2nd Gen Tech": "coal",
            "Coal-new": "coal",
            "NG F-Frame CT": "CCGT",
            "NG Combustion Turbine (F-Frame)": "CCGT",
            "Hydropower - NPD 1": "hydro",
            "Hydropower - NSD 1": "ror",
            "Pumped Storage Hydropower - National Class 1": "PHS",
            "Nuclear - Large": "nuclear",
            "Nuclear - AP1000": "nuclear",
            "Geothermal - Hydro / Flash": "geothermal",
            "Land-Based Wind - Class 4": "onwind",
            "Land-Based Wind - Class 4 - Technology 1": "onwind",
            "Offshore Wind - Class 3": "offwind",
            "Utility PV - Class 5": "solar-utility",
            "Commercial PV - Class 1": "solar-rooftop",
            "Utility-Scale Battery Storage - 6Hr": "battery storage",
            "Biopower": "biomass",
            "Biopower - Dedicated": "biomass",
            "CSP - Class 2": "csp-tower",
        }
    elif flag.casefold() == "atb_technology_name":
        return {
            "Land-Based Wind - Class 1": "Land-Based Wind - Class 1 - Technology 1",
            "Land-Based Wind - Class 2": "Land-Based Wind - Class 2 - Technology 1",
            "Land-Based Wind - Class 3": "Land-Based Wind - Class 3 - Technology 1",
            "Land-Based Wind - Class 5": "Land-Based Wind - Class 5 - Technology 1",
            "Land-Based Wind - Class 6": "Land-Based Wind - Class 6 - Technology 1",
            "Land-Based Wind - Class 7": "Land-Based Wind - Class 7 - Technology 1",
            "Land-Based Wind - Class 8": "Land-Based Wind - Class 8 - Technology 2",
            "Land-Based Wind - Class 9": "Land-Based Wind - Class 9 - Technology 3",
            "Land-Based Wind - Class 10": "Land-Based Wind - Class 10 - Technology 4",
            "NG F-Frame CC": "NG 2-on-1 Combined Cycle (F-Frame)",
            "NG H-Frame CC": "NG 2-on-1 Combined Cycle (H-Frame)",
            "NG combined cycle 95% CCS (F-frame basis -> 2nd Gen Tech)": "NG 2-on-1 Combined Cycle (F-Frame) 95% CCS",
            "NG combined cycle 95% CCS (H-frame basis -> 2nd Gen Tech)": "NG 2-on-1 Combined Cycle (H-Frame) 95% CCS",
            "NG combined cycle Max CCS (F-frame basis -> 2nd Gen Tech)": "NG 2-on-1 Combined Cycle (F-Frame) 97% CCS",
            "NG combined cycle Max CCS (H-frame basis -> 2nd Gen Tech)": "NG 2-on-1 Combined Cycle (H-Frame) 97% CCS",
            "Coal-CCS-95% -> 2nd Gen Tech": "Coal-95%-CCS",
            "Coal-Max-CCS -> 2nd Gen Tech": "Coal-99%-CCS",
            "Coal-IGCC": "Coal - IGCC",
            "CSP - Class 7": "CSP - Class 8",
            "Nuclear - Small Modular Reactor": "Nuclear - Small",
        }
    elif flag.casefold() == "output_column":
        return {
            "display_name": "technology",
            "core_metric_parameter": "parameter",
            "units": "unit",
            "atb_year": "currency_year",
            "core_metric_case": "financial_case",
        }
    else:
        raise Exception(
            "{} is not among the allowed choices: parameter, technology, output_column"
        )


def filter_atb_input_file(
    input_file_path: pathlib.Path,
    year: int,
    list_columns_to_keep: list,
    list_core_metric_parameter_to_keep: list,
    list_tech_to_remove: list,
) -> pd.DataFrame:
    """
    The function filters the input cost dataframe from NREL/ATB. Namely, it:
    - selects the necessary columns (the atb_e_2022 and atb_e_2024 have in fact a slightly different schema)
    - selects the rows corresponding to the necessary core_metric_parameter(s)
    - selects the rows corresponding to the year for the cost assumption (2020, 2025, 2030 etc.)
    - drops duplicated rows
    - drops rows corresponding to unnecessary technologies

    Parameters
    ----------
    input_file_path : pathlib.Path
        NREL/ATB file path
    year: int
        year for the cost assumption
    list_columns_to_keep: list
        columns from NREL/ATB dataset that are relevant
    list_core_metric_parameter_to_keep: list
        values of the core_metric_paramater that are relevant
    list_tech_to_remove: list
        technologies names that are should be excluded

    Returns
    -------
    pandas.DataFrame
        NREL/ATB cost dataframe
    """

    atb_file_df = pd.read_parquet(input_file_path)
    list_core_metric_parameter_to_keep = [
        str(x).casefold() for x in list_core_metric_parameter_to_keep
    ]
    list_tech_to_remove = [str(x).casefold() for x in list_tech_to_remove]
    year_string = str(year).casefold()

    # --> select columns
    for column_name in list_columns_to_keep:
        if column_name not in atb_file_df:
            print("missing column", column_name)
            atb_file_df[column_name] = pd.Series(dtype="str")
        else:
            pass
    atb_file_df = atb_file_df.loc[:, list_columns_to_keep]

    # --> select rows based on core_metric_parameter
    # Note: we do not apply any selection on core_metric_case and scenario.
    # Such selection can be done in the model config (e.g. PyPSA-Earth)
    atb_file_df = atb_file_df.loc[
        atb_file_df["core_metric_parameter"]
        .str.casefold()
        .isin(list_core_metric_parameter_to_keep)
    ]

    # --> select rows based on core_metric_variable
    if input_file_path.name == "atb_e_2022.parquet":
        # Note: 2020 data are fetched from the input file atb_e_2022.
        atb_file_df = atb_file_df.loc[
            atb_file_df["core_metric_variable"].astype(str).str.casefold()
            == year_string
        ]
    elif input_file_path.name == "atb_e_2024.parquet":
        # Note: 2025, 2030, 2035, 2040, 2045, 2050 data are fetched from the input file atb_e_2024.*
        atb_file_df = atb_file_df.loc[
            atb_file_df["core_metric_variable"].astype(str).str.casefold()
            == year_string
        ]
    else:
        raise Exception(
            f"{input_file_path.stem} - the input file considered is not among the needed ones: atb_e_2022.parquet, atb_e_2024.parquet"
        )

    # --> drop duplicated rows
    atb_file_df = atb_file_df.drop_duplicates(keep="first")

    # --> remove unnecessary technologies
    atb_file_df = atb_file_df.loc[
        ~atb_file_df["display_name"].str.casefold().isin(list_tech_to_remove)
    ]

    return atb_file_df


def get_query_string(
    column_list: list, column_to_exclude: list, fom_normalization_parameter: str
) -> str:
    """
    The Fixed O&M values from the NREL/ATB database ought to be normalized by Additional OCC
    (for retrofits technologies) or CAPEX (for any other technology). The function returns the
    query strings that links Fixed O&M values to the corresponding Additional OCC and CAPEX values

    Parameters
    ----------
    column_list : list
        columns to consider in the query
    column_to_exclude: list
        columns that should not be considered in the query
    fom_normalization_parameter: str
        either Additional OCC or CAPEX

    Returns
    -------
    str
        query string
    """

    if set(column_to_exclude).issubset(set(column_list)):
        column_to_use_list = list(set(column_list) - set(column_to_exclude))
        query_list = sorted(
            [
                f"{column_name}.str.casefold() == '{fom_normalization_parameter}'"
                if column_name == "core_metric_parameter"
                else f"{column_name} == @x.{column_name}"
                for column_name in column_to_use_list
            ]
        )
        query_string = " & ".join(query_list)
    else:
        exception_message = f"The following columns {list(set(column_to_exclude).difference(set(column_list)))} are not included in the original list"
        raise Exception(exception_message)
    return query_string.strip()


def calculate_fom_percentage(
    x: pd.Series, dataframe: pd.DataFrame, columns_list: list
) -> float:
    """
    The Fixed O&M values from the NREL/ATB database ought to be normalized by Additional OCC
    (for retrofits technologies) or CAPEX (for any other technology). The function returns the
    query strings that links Fixed O&M values to the corresponding Additional OCC and CAPEX values

    Parameters
    ----------
    x : row
        row of the cost dataframe
    dataframe : pandas.DataFrame
        cost DataFrame
    columns_list: list
        columns to consider in the query

    Returns
    -------
    float
        normalized value of Fixed O&M
    """

    if x["core_metric_parameter"].casefold() == "fixed o&m":
        if "retrofit" in x["technology"].casefold():
            query_string = get_query_string(
                columns_list, ["units", "value"], "additional occ"
            )
        else:
            query_string = get_query_string(columns_list, ["units", "value"], "capex")
        fom_perc_value = x.value / dataframe.query(query_string)["value"] * 100.0
        return round(fom_perc_value.values[0], 2)
    else:
        return x.value


def replace_value_name(
    dataframe: pd.DataFrame, conversion_dict: dict, column_name: str
) -> pd.DataFrame:
    """
    The function replaces values in a given dataframe column, using values from a conversion dictionary.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        cost dataframe
    conversion_dict : dict
        conversion dictionary
    column_name: str
        column name where values shall be replaced

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    dataframe[column_name] = dataframe[column_name].replace(conversion_dict)
    return dataframe


def pre_process_manual_input_usa(
    manual_input_usa_file_path: str,
    list_of_years: list,
    year: int,
) -> pd.DataFrame:
    """
    The function reads and modifies the manual_input_usa.csv file. Namely, it:
    - reads the input file
    - renames the column "further_description" to "further description"
    - starting from manual_input_usa.csv, it estimates the parameters for each technology for all the requested years
    - it selects the values for a given year
    - queries the necessary rows of the existing cost dataframe

    Parameters
    ----------
    manual_input_usa_file_path : str
        manual_input_usa.csv file path
    list_of_years : list
        years for which a cost assumption is provided
    year : int
        year from list_of_years

    Returns
    -------
    pandas.DataFrame
        updated manual input usa
    """

    # Read the input file
    manual_input_usa_file_df = pd.read_csv(
        manual_input_usa_file_path, quotechar='"', sep=",", keep_default_na=False
    )

    # Rename the column further_description
    manual_input_usa_file_df = manual_input_usa_file_df.rename(
        columns={"further_description": "further description"}
    )

    # Create cost estimates for all years
    list_dataframe_row = []
    for tech in manual_input_usa_file_df["technology"].unique():
        technology_filtered_df = manual_input_usa_file_df[
            manual_input_usa_file_df["technology"] == tech
        ]
        for param in technology_filtered_df["parameter"].unique():
            technology_parameter_filtered_df = manual_input_usa_file_df.query(
                "technology == @tech and parameter == @param"
            )

            # Consider differences among scenarios
            scenarios = technology_parameter_filtered_df["scenario"].dropna().unique()

            for scenario in scenarios:
                scenario_value = technology_parameter_filtered_df[
                    technology_parameter_filtered_df["scenario"] == scenario
                ]["value"].values  # Extract values for each scenario

                if scenario_value.size > 0:
                    scenario_years = technology_parameter_filtered_df[
                        technology_parameter_filtered_df["scenario"] == scenario
                    ]["year"].values
                    scenario_values = technology_parameter_filtered_df[
                        technology_parameter_filtered_df["scenario"] == scenario
                    ]["value"].values

                    interpolated_values = np.interp(
                        list_of_years, scenario_years, scenario_values
                    )

                    # Create a row for each scenario
                    s_copy = pd.Series(
                        index=list_of_years,
                        data=interpolated_values,  # values are now interpolated
                        name=param,
                    )

                    s_copy["parameter"] = param
                    s_copy["technology"] = tech
                    s_copy["scenario"] = scenario
                    try:
                        s_copy["currency_year"] = int(
                            technology_parameter_filtered_df["currency_year"].values[0]
                        )
                    except ValueError:
                        s_copy["currency_year"] = np.nan

                    # Add the other columns in the data file
                    for col in ["unit", "source", "further description"]:
                        s_copy[col] = technology_parameter_filtered_df[col].unique()[0]

                    # Add a separate row for each `financial_case`
                    for financial_case in technology_parameter_filtered_df[
                        "financial_case"
                    ].unique():
                        s_copy["financial_case"] = financial_case
                        list_dataframe_row.append(s_copy.copy())
    manual_input_usa_file_df = pd.DataFrame(list_dataframe_row).reset_index(drop=True)

    # Filter the information for a given year
    manual_input_usa_file_df = manual_input_usa_file_df[
        [
            "technology",
            "parameter",
            year,
            "unit",
            "source",
            "further description",
            "currency_year",
            "financial_case",
            "scenario",
        ]
    ].rename(columns={year: "value"})

    # Filter data to get technologies with scenario differentiation
    with_scenario_df = manual_input_usa_file_df[
        manual_input_usa_file_df["scenario"].notna()
    ]
    without_scenario_df = manual_input_usa_file_df[
        manual_input_usa_file_df["scenario"].isna()
    ]

    final_rows = []

    for tech in manual_input_usa_file_df["technology"].unique():
        tech_with_scenario = with_scenario_df[with_scenario_df["technology"] == tech]
        if len(tech_with_scenario) > 0:
            # Keep rows where a scenario exists
            final_rows.append(tech_with_scenario)
        else:
            # If a scenario is not defined, keep the row without scenario
            tech_without_scenario = without_scenario_df[
                without_scenario_df["technology"] == tech
            ]
            final_rows.append(tech_without_scenario)

    manual_input_usa_file_df = pd.concat(final_rows, ignore_index=True)

    # Cast the value column to float
    manual_input_usa_file_df["value"] = manual_input_usa_file_df["value"].astype(float)

    return manual_input_usa_file_df


def modify_cost_input_file(
    cost_dataframe: pd.DataFrame,
    manual_input_usa_dataframe: pd.DataFrame,
    list_of_years: list,
    year: int,
    n_digits: int,
) -> pd.DataFrame:
    """
    The function filters out from the existing cost dataframe the rows corresponding
    to the technology-parameter pairs from manual_input_usa.csv. It then concatenates manual_input_usa.csv and
    adjourns the estimates for "electrobiofuels". Namely, it:
    - creates a list of tuples (technology, parameter) from manual_input_usa.csv
    - filters out from the existing cost dataframe all rows with (technology, parameter) in manual_input_usa.csv
    - concatenates manual_input_usa.csv
    - updates the parameters for electrobiofuels

    Parameters
    ----------
    cost_dataframe : pandas.DataFrame
        existing cost dataframe
    manual_input_usa_dataframe : pandas.DataFrame
        manual_input_usa dataframe
    list_of_years : list
        years for which a cost assumption is provided
    year : int
        year from list_of_years
    n_digits : int
        number of significant digits

    Returns
    -------
    pandas.DataFrame
        updated cost dataframe
    """

    # Create a list of tuples (technology, parameter) from manual_input_usa.csv
    list_technology_parameter_tuples_manual_input_usa = [
        (str(x), str(y))
        for (x, y) in zip(
            manual_input_usa_dataframe.technology, manual_input_usa_dataframe.parameter
        )
    ]

    # Filter out rows from the existing cost dataframe corresponding
    # to the (technology, parameter) pairs from manual_input_usa.csv
    queried_cost_df = cost_dataframe[
        ~(
            pd.Series(
                list(zip(cost_dataframe["technology"], cost_dataframe["parameter"]))
            ).isin(list_technology_parameter_tuples_manual_input_usa)
        )
    ]

    # Concatenate manual_input_usa.csv to the filtered existing cost dataframe
    updated_cost_dataframe = pd.concat(
        [queried_cost_df, manual_input_usa_dataframe]
    ).reset_index(drop=True)

    # Update the parameters for the technology "electrobiofuels"
    btl_cost_data = np.interp(x=list_of_years, xp=[2020, 2050], fp=[3500, 2000])
    btl_cost = pd.Series(data=btl_cost_data, index=list_of_years)

    efuel_scaling_factor = (
        updated_cost_dataframe.query(
            "technology.str.casefold() == 'btl' & parameter.str.casefold() == 'c stored'"
        )["value"].values[0]
        * updated_cost_dataframe.query(
            "technology.str.casefold() == 'fischer-tropsch' & parameter.str.casefold() == 'capture rate'"
        )["value"].values[0]
    )

    investment_cost = np.round(
        (
            btl_cost[year]
            + updated_cost_dataframe.loc[
                (
                    updated_cost_dataframe["technology"].str.casefold()
                    == "fischer-tropsch"
                )
                & (updated_cost_dataframe["parameter"].str.casefold() == "investment"),
                "value",
            ].values[0]
            * efuel_scaling_factor
        ),
        n_digits,
    )

    updated_cost_dataframe.loc[
        (updated_cost_dataframe["technology"].str.casefold() == "electrobiofuels")
        & (updated_cost_dataframe["parameter"].str.casefold() == "efficiency-tot"),
        "value",
    ] = np.round(
        1.0
        / (
            1.0
            / updated_cost_dataframe.loc[
                (
                    updated_cost_dataframe["technology"].str.casefold()
                    == "electrobiofuels"
                )
                & (
                    updated_cost_dataframe["parameter"].str.casefold()
                    == "efficiency-hydrogen"
                ),
                "value",
            ].values[0]
            + 1.0
            / updated_cost_dataframe.loc[
                (
                    updated_cost_dataframe["technology"].str.casefold()
                    == "electrobiofuels"
                )
                & (updated_cost_dataframe["parameter"] == "efficiency-biomass"),
                "value",
            ].values[0]
        ),
        n_digits,
    )

    updated_cost_dataframe.loc[
        (updated_cost_dataframe["technology"] == "electrobiofuels")
        & (updated_cost_dataframe["parameter"] == "efficiency-hydrogen"),
        "value",
    ] = np.round(
        (
            updated_cost_dataframe.loc[
                (updated_cost_dataframe["technology"] == "Fischer-Tropsch")
                & (updated_cost_dataframe["parameter"] == "efficiency"),
                "value",
            ].values[0]
            / efuel_scaling_factor
        ),
        n_digits,
    )

    vom_value = np.round(
        (
            updated_cost_dataframe.loc[
                (updated_cost_dataframe["technology"] == "BtL")
                & (updated_cost_dataframe["parameter"] == "VOM"),
                "value",
            ].values[0]
            + updated_cost_dataframe.loc[
                (updated_cost_dataframe["technology"] == "Fischer-Tropsch")
                & (updated_cost_dataframe["parameter"] == "VOM"),
                "value",
            ].values[0]
            * efuel_scaling_factor
        ),
        n_digits,
    )

    if investment_cost > 0.0:
        updated_cost_dataframe.loc[
            (updated_cost_dataframe["technology"] == "electrobiofuels")
            & (updated_cost_dataframe["parameter"] == "investment"),
            "value",
        ] = investment_cost
        updated_cost_dataframe.loc[
            (updated_cost_dataframe["technology"] == "electrobiofuels")
            & (updated_cost_dataframe["parameter"] == "investment"),
            "currency_year",
        ] = 2022.0

    if vom_value > 0.0:
        updated_cost_dataframe.loc[
            (updated_cost_dataframe["technology"] == "electrobiofuels")
            & (updated_cost_dataframe["parameter"] == "VOM"),
            "value",
        ] = vom_value
        updated_cost_dataframe.loc[
            (updated_cost_dataframe["technology"] == "electrobiofuels")
            & (updated_cost_dataframe["parameter"] == "VOM"),
            "currency_year",
        ] = 2022.0

    return updated_cost_dataframe


def query_cost_dataframe(
    cost_dataframe: pd.DataFrame,
    technology_dictionary: dict,
    parameter_dictionary: dict,
) -> pd.DataFrame:
    """
    The function queries the rows of the existing cost dataframe.
    The selection is done by means of an OR operator. The two operands for this logical operation are returned
    by the following queries:
    - query_string_part_one: selects all the rows corresponding to the technologies NOT updated with NREL-ATB data
    - query_string_part_two: some of the techno-economic parameters (e.g., efficiency, capture rate) to be updated with NREL-ATB data are NOT present in the NREL-ATB dataset. They are instead added to the former cost csv files by means of the manual_input.csv. They should be kept in the final output. This query selects such rows

    Parameters
    ----------
    cost_dataframe: pandas.DataFrame
        existing cost dataframe
    technology_dictionary: dict
        dictionary of the technologies updated with NREL/ATB data
    parameter_dictionary: dict
        dictionary of the parameters for which NREL/ATB estimates are available

    Returns
    -------
    pandas.DataFrame
        queried cost dataframe
    """

    technologies_from_nrel = [
        str(x).casefold() for x in set(technology_dictionary.values())
    ]
    parameters_from_nrel = [
        str(x).casefold() for x in set(parameter_dictionary.values())
    ]

    # Query the rows from the former cost csv file.
    # --> The selection is done by means of an OR operator.
    # --> The two operands for this logical operation are returned by the queries below

    # --> QUERY 1: this query selects all the rows corresponding to the technologies NOT updated with NREL-ATB data
    query_string_part_one = "~technology.str.casefold().isin(@t)"

    # --> QUERY 2: some of the parameters of the technologies to be updated with NREL-ATB data are NOT present
    # --> in the NREL-ATB dataset. They are instead added to the former cost csv files by means of the
    # --> manual_input.csv. They should be kept in the final output. This query selects such rows
    query_string_part_two = (
        "technology.str.casefold().isin(@t) & ~parameter.str.casefold().isin(@p)"
    )

    query_string = f"{query_string_part_one} | {query_string_part_two}"

    queried_cost_dataframe = cost_dataframe.query(
        query_string,
        local_dict={"t": technologies_from_nrel, "p": parameters_from_nrel},
    ).reset_index(drop=True)

    return queried_cost_dataframe


def pre_process_cost_input_file(
    input_file_path: str, columns_to_add_list: list
) -> pd.DataFrame:
    """
    The function filters and cleans the existing cost file. Namely, it:
    - reads the input file
    - adds the columns from NREL/ATB not present in the existing cost dataframe
    - queries the necessary rows of the existing cost dataframe

    Parameters
    ----------
    input_file_path : str
        existing cost file path
    columns_to_add_list: list
        column names from NREL/ATB to be added to the existing cost dataframe

    Returns
    -------
    pandas.DataFrame
        updated NREL/ATB cost dataframe
    """

    cost_input_df = pd.read_csv(input_file_path)

    # The data coming from NREL/ATB contain the columns "financial_case", "scenario".
    # Such columns are NOT present in the former cost csv. They are added here
    for column_name in columns_to_add_list:
        cost_input_df[column_name] = pd.Series(dtype="str")

    cost_input_df = query_cost_dataframe(
        cost_input_df,
        get_conversion_dictionary("pypsa_technology_name"),
        get_conversion_dictionary("parameter"),
    )

    return cost_input_df.reset_index(drop=True)


def pre_process_atb_input_file(
    input_file_path: str,
    nrel_source: str,
    nrel_further_description: str,
    year: int,
    list_columns_to_keep: list,
    list_core_metric_parameter_to_keep: list,
    tech_to_remove: list,
) -> pd.DataFrame:
    """
    The function filters and cleans the input NREL/ATB cost file. Namely, it:
    - reads the input file
    - normalizes the Fixed O&M by Additional OCC (for retrofits technologies) or CAPEX (for any other technology)
    - changes the units
    - renames the technology names to the PyPSA nomenclature
    - aligns the atb_e_2022 nomenclature to the atb_e 2024 nomenclature

    Parameters
    ----------
    input_file_path : str
        NREL/ATB file path
    nrel_source: str
        link to the NREL/ATB source files. This information shall be used to populate the source column
    nrel_further_description: str
        text that details the further description field for the NREL/ATB sources
    year: int
        year for the cost assumption
    list_columns_to_keep: list
        columns from NREL/ATB dataset that are relevant
    list_core_metric_parameter_to_keep: list
        values of the core_metric_paramater that are relevant
    tech_to_remove: list
        technologies names that should be excluded from NREL/ATB

    Returns
    -------
    pandas.DataFrame
        updated NREL/ATB cost dataframe
    """

    # Read inputs and filter relevant columns and relevant core_metric_variables (i.e. the years to consider)
    atb_input_df = filter_atb_input_file(
        pathlib.Path(input_file_path),
        year,
        list_columns_to_keep,
        list_core_metric_parameter_to_keep,
        tech_to_remove,
    )

    # Normalize Fixed O&M by CAPEX (or Additional OCC for retrofit technologies)
    atb_input_df["value"] = atb_input_df.apply(
        lambda x: calculate_fom_percentage(x, atb_input_df, list_columns_to_keep),
        axis=1,
    )

    # Modify the unit of the normalized Fixed O&M to %/yr
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "%/year"
        if x["core_metric_parameter"].casefold() == "fixed o&m"
        else x["units"],
        axis=1,
    )

    # Modify the unit of CF to per unit
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "per unit"
        if x["core_metric_parameter"].casefold() == "cf"
        else x["units"],
        axis=1,
    )

    # Modify the unit of Additional OCC to USD/kW instead of $/kW
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "USD/kW"
        if x["core_metric_parameter"].casefold() == "additional occ"
        else x["units"],
        axis=1,
    )

    # Modify the unit of CAPEX to USD/kW instead of $/kW
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "USD/kW"
        if x["core_metric_parameter"].casefold() == "capex"
        else x["units"],
        axis=1,
    )

    # Modify the unit of Variable O&M to USD/MWh instead of $/MWh
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "USD/MWh"
        if x["core_metric_parameter"].casefold() == "variable o&m"
        else x["units"],
        axis=1,
    )

    # Modify the unit of Fuel cost O&M to USD/MWh instead of $/MWh
    atb_input_df["units"] = atb_input_df.apply(
        lambda x: "USD/MWh"
        if x["core_metric_parameter"].casefold() == "fuel"
        else x["units"],
        axis=1,
    )

    # Replace the display_name column values with PyPSA technology names
    technology_conversion_dict_pypsa = get_conversion_dictionary(
        "pypsa_technology_name"
    )
    atb_input_df = replace_value_name(
        atb_input_df, technology_conversion_dict_pypsa, "display_name"
    )

    # Uniform the display_name nomenclature of atb_e_2022 to the one of atb_e_2024
    technology_conversion_dict_atb = get_conversion_dictionary("atb_technology_name")
    atb_input_df = replace_value_name(
        atb_input_df, technology_conversion_dict_atb, "display_name"
    )

    # Add source column
    atb_input_df["source"] = nrel_source

    # Add further description column
    atb_input_df["further description"] = nrel_further_description

    # Rename columns and select just columns used in PyPSA
    column_rename_dict = get_conversion_dictionary("output_column")
    tuple_output_columns_to_keep = (
        "display_name",
        "core_metric_parameter",
        "value",
        "units",
        "source",
        "further description",
        "atb_year",
        "scenario",
        "core_metric_case",
    )
    atb_input_df = atb_input_df.loc[:, tuple_output_columns_to_keep].rename(
        columns=column_rename_dict
    )

    # Replace parameter with PyPSA cost parameter names
    parameter_conversion_dict = get_conversion_dictionary("parameter")
    atb_input_df = replace_value_name(
        atb_input_df, parameter_conversion_dict, "parameter"
    )

    # ATB currency year dates back to 2018 for ATB2020 and to 2022 for ATB2024
    atb_input_df["currency_year"] = atb_input_df["currency_year"] - 2

    # Cast currency_year from int to float
    atb_input_df["currency_year"] = atb_input_df["currency_year"].astype(float)

    return atb_input_df.reset_index(drop=True)


def duplicate_fuel_cost(input_file_path: str, list_of_years: list) -> pd.DataFrame:
    """
    The function reads-in the fuel cost file to a Pandas DataFrame and
    replicates the last available row for each technology. Namely, it
    - reads the fuel cost input file to a Pandas DataFrame
    - determines the list of the technologies
    - loops through the technology list and determines for each technology the last available year
    - creates a list with the missing years
    - replicates the estimation for the last available year for all the missing years

    Parameters
    ----------
    input_file_path : str
        fuel cost file path
    list_of_years: list
        years for which a cost assumption is provided

    Returns
    -------
    pandas.DataFrame
        updated fuel cost dataframe
    """

    # Read-in the fuel cost file for the US
    input_fuel_cost_df = pd.read_csv(input_file_path)

    # Determine a list of the technologies
    technology_list = list(input_fuel_cost_df["technology"].unique())

    # Create an empty dataframe
    replicated_fuel_cost_df = pd.DataFrame()

    # Loop through the available technologies
    for tech_value in technology_list:
        # For each technology, determine the last available year and extract the list of missing years
        max_year = np.max(
            input_fuel_cost_df[input_fuel_cost_df["technology"] == tech_value]["year"]
        )
        first_missing_year_index = list_of_years.index(max_year) + 1
        missing_year_list = list_of_years[first_missing_year_index:]

        # For each technology, loop through the list of missing years
        for i, val_year in enumerate(missing_year_list):
            # Extract the row corresponding to the last available year and replace the year with the missing year
            df_to_replicate = (
                input_fuel_cost_df.loc[
                    (input_fuel_cost_df["technology"] == tech_value)
                    & (input_fuel_cost_df["year"] == max_year)
                ]
                .copy(deep=True)
                .replace(max_year, val_year)
            )

            # Append the extracted and modified row to the (originally empty) dataframe
            replicated_fuel_cost_df = pd.concat(
                [replicated_fuel_cost_df, df_to_replicate]
            ).reset_index(drop=True)

    # Append the dataframe with the replicated rows to the original dataframe. Sort by technology and year
    input_fuel_cost_df = (
        pd.concat([replicated_fuel_cost_df, input_fuel_cost_df])
        .sort_values(by=["technology", "year"])
        .reset_index(drop=True)
    )
    return input_fuel_cost_df


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("compile_cost_assumptions_usa")

    configure_logging(snakemake)

    year_list = sorted(snakemake.config["years"])
    num_digits = snakemake.config["ndigits"]
    eur_reference_year = snakemake.config["eur_year"]
    input_file_list_atb = snakemake.input.nrel_atb_input_files
    input_file_discount_rate = snakemake.input.nrel_atb_input_discount_rate
    input_file_fuel_costs = snakemake.input.nrel_atb_input_fuel_costs
    input_file_manual_input_usa = snakemake.input.nrel_atb_manual_input_usa
    input_file_inflation_rate = snakemake.input.inflation_rate
    cost_file_list = snakemake.input.cost_files_to_modify
    nrel_atb_columns_to_keep = snakemake.config["nrel_atb"]["nrel_atb_columns_to_keep"]
    nrel_atb_core_metric_parameter_to_keep = snakemake.config["nrel_atb"][
        "nrel_atb_core_metric_parameter_to_keep"
    ]
    nrel_atb_source_link = snakemake.config["nrel_atb"]["nrel_atb_source_link"]
    nrel_atb_further_description = snakemake.config["nrel_atb"][
        "nrel_atb_further_description"
    ]
    nrel_atb_technology_to_remove = snakemake.config["nrel_atb"][
        "nrel_atb_technology_to_remove"
    ]

    logger.info("Configuration variables are set")

    if len(set(snakemake.config["years"])) < len(snakemake.config["years"]):
        raise Exception(
            "Please verify the list of cost files. It may contain duplicates."
        )

    if len(year_list) != len(cost_file_list):
        raise Exception(
            f"The cost files {year_list} are more than the considered years {cost_file_list}"
        )

    # get the discount rate values for the US
    discount_rate_df = pd.read_csv(input_file_discount_rate)

    logger.info("discount_rate file for the US has been read in")

    # get the fuel costs values for the US
    fuel_costs_df = duplicate_fuel_cost(input_file_fuel_costs, year_list)

    logger.info("fuel_cost file for the US has been read in")

    for year_val in year_list:
        # get the cost file to modify
        input_cost_path = [
            path
            for path in snakemake.input.cost_files_to_modify
            if str(year_val) in path
        ][0]

        # get the atb values for a given year
        if year_val == 2020:
            # choose atb_e_2022
            input_atb_path = input_file_list_atb[0]
        elif year_val in year_list[1:]:
            # choose atb_e_2024
            input_atb_path = input_file_list_atb[1]
        else:
            raise Exception(f"{year_val} is not a considered year")

        logger.info(f"The file {input_atb_path} is used for year {year_val}")

        manual_input_usa_df = pre_process_manual_input_usa(
            input_file_manual_input_usa,
            year_list,
            year_val,
        )

        cost_df = pre_process_cost_input_file(
            input_cost_path, ["financial_case", "scenario"]
        )

        cost_df = modify_cost_input_file(
            cost_df, manual_input_usa_df, year_list, year_val, num_digits
        )

        atb_e_df = pre_process_atb_input_file(
            input_atb_path,
            nrel_atb_source_link,
            nrel_atb_further_description,
            year_val,
            nrel_atb_columns_to_keep,
            nrel_atb_core_metric_parameter_to_keep,
            nrel_atb_technology_to_remove,
        )

        # get the discount rate file for the given year
        discount_rate_year_df = discount_rate_df.loc[
            discount_rate_df["year"] == year_val
        ]
        discount_rate_year_df = discount_rate_year_df.loc[
            :,
            (
                "technology",
                "parameter",
                "value",
                "unit",
                "source",
                "further description",
                "currency_year",
                "financial_case",
                "scenario",
            ),
        ].reset_index(drop=True)

        # get the fuel costs file for the given year
        fuel_costs_year_df = fuel_costs_df.loc[fuel_costs_df["year"] == year_val]
        fuel_costs_year_df = fuel_costs_year_df.loc[
            :,
            (
                "technology",
                "parameter",
                "value",
                "unit",
                "source",
                "further description",
                "currency_year",
                "financial_case",
                "scenario",
            ),
        ].reset_index(drop=True)

        # concatenate the existing and NREL/ATB cost dataframes
        updated_cost_df = pd.concat([cost_df, atb_e_df]).reset_index(drop=True)

        # add discount rate

        # --> remove possible entries which are going to be replaced with the data from discount_rates_usa.csv
        technology_discount_rate_list = list(
            discount_rate_year_df["technology"].unique()
        )
        query_string_discount_rate = "~(technology.str.casefold().isin(@technology_discount_rate_list) & parameter.str.casefold() == 'discount rate')"

        # --> concatenate the discount_rates_usa.csv data
        updated_cost_df = pd.concat(
            [updated_cost_df.query(query_string_discount_rate), discount_rate_year_df]
        ).reset_index(drop=True)

        # add fuel costs

        # --> remove possible entries which are going to be replaced with the data from fuel_costs_usa.csv
        technology_fuel_cost_list = list(fuel_costs_year_df["technology"].unique())
        query_string_fuel_cost = "~(technology.str.casefold().isin(@technology_fuel_cost_list) & parameter.str.casefold() == 'fuel')"

        # --> concatenate the fuel_costs_usa.csv data
        updated_cost_df = pd.concat(
            [updated_cost_df.query(query_string_fuel_cost), fuel_costs_year_df]
        ).reset_index(drop=True)

        # Cast "value" from float
        updated_cost_df["value"] = updated_cost_df["value"].astype(float)

        # Sort the modified cost dataframe by technology and parameter
        updated_cost_df = updated_cost_df.sort_values(
            by=["technology", "parameter"]
        ).reset_index(drop=True)

        # correct for inflation for technology-parameter pairs having units that contain USD
        inflation_rate_series_usd = prepare_inflation_rate(
            input_file_inflation_rate, "USD"
        )
        mask_usd = (
            updated_cost_df["unit"].str.casefold().str.startswith("usd", na=False)
        )
        inflation_adjusted_updated_cost_df = updated_cost_df.copy()
        inflation_adjusted_updated_cost_df.loc[mask_usd, "value"] = (
            adjust_for_inflation(
                inflation_rate_series_usd,
                inflation_adjusted_updated_cost_df.loc[mask_usd],
                inflation_adjusted_updated_cost_df.loc[mask_usd, "technology"].unique(),
                eur_reference_year,
                "value",
                usa_costs_flag=True,
            )["value"]
        )

        # round the value column
        rounded_df = inflation_adjusted_updated_cost_df.copy()
        rounded_df.loc[:, "value"] = round(
            inflation_adjusted_updated_cost_df.value.astype(float), num_digits
        )

        # output the modified cost dataframe
        output_cost_path_list = [
            path for path in snakemake.output if str(year_val) in path
        ]
        if len(output_cost_path_list) == 1:
            output_cost_path = output_cost_path_list[0]
            rounded_df.to_csv(output_cost_path, index=False)
            logger.info(
                f"The cost assumptions file for the US has been compiled for year {year_val}"
            )
        else:
            raise Exception(
                "Please verify the list of cost files. It may contain duplicates."
            )
