# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Parser for version 0.13.4 of the manual_input_usa.csv dataset."""

import logging
import pathlib
from typing import Any

import pandas

from technologydata import (
    Commons,
    Parameter,
    Source,
    SourceCollection,
    Technology,
    TechnologyCollection,
)
from technologydata.parsers.data_parser_base import ParserBase

path_cwd = pathlib.Path.cwd()

logger = logging.getLogger(__name__)


class ManualInputUSAV0134Parser(ParserBase):
    """Parser for v0.13.4 of the manual_input_usa.csv dataset."""

    @staticmethod
    def _extract_units_carriers_heating_value(
        input_unit: str,
    ) -> tuple[str, str | None, str | None]:
        """
        Extract standardized units and carriers from an input unit string. Add also heating_value.

        This function maps complex unit representations to simplified unit and carrier
        combinations using a predefined dictionary of special patterns.

        Parameters
        ----------
        input_unit : str
            A specialized unit string to be converted.

        Returns
        -------
        tuple[str, str | None, str | None]
            A tuple containing two elements:
            - The first element is the standardized unit
            - The second element is the corresponding carrier (or None if not found)
            - The third element is the corresponding heating value (or None if not found)

        """
        # Define conversion dictionary
        special_patterns = {
            "USD_2022/MW_FT": ("USD_2022/MW", "1/FT", "1/LHV"),
            "MWh_H2/MWh_FT": ("MWh/MWh", "H2/FT", "LHV"),
            "MWh_el/MWh_FT": ("MWh/MWh", "el/FT", "LHV"),
            "t_CO2/MWh_FT": ("t/MWh", "CO2/FT", "LHV"),
            "USD_2022/kWh_H2": ("USD_2022/kWh", "1/H2", "LHV"),
            "MWh_el/MWh_H2": ("MWh/MWh", "el/H2", "LHV"),
            "USD_2023/t_CO2/h": ("USD_2023/t/h", "1/CO2", None),
            "MWh_el/t_CO2": ("MWh/t", "el/CO2", "LHV"),
            "MWh_th/t_CO2": ("MWh/t", "thermal/CO2", "LHV"),
        }

        if isinstance(input_unit, str) and input_unit in special_patterns.keys():
            return special_patterns[input_unit]
        else:
            return input_unit, None, None

    @staticmethod
    def _build_technology_collection(
        dataframe: pandas.DataFrame,
        sources_path: pathlib.Path,
        archive_source: bool = False,
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
            - 'scenario': Estimation or case identifier
            - 'year': Year of the technology
            - 'technology': Detailed technology name
            - 'parameter': Parameter name
            - 'value': Parameter value
            - 'unit': Parameter units
            - 'further_description': Extra information about the technology
            - 'financial_case': Technology financial case
        sources_path: pathlib.Path
            Output path for storing the SourceCollection object
        archive_source: Optional[bool]
            Flag to decide whether to archive the source object on the Wayback Machine. Default False.
        output_schema : Optional[bool]
            Flag to decide whether to export the source collection schema. Default False.

        Returns
        -------
        TechnologyCollection
            A collection of Technology instances, each representing a unique
            technology group with its associated parameters.

        Notes
        -----
        - The function groups the DataFrame by ["scenario", "year", "technology"]
        - For each group, it creates a dictionary of Parameters
        - Each Technology is instantiated with group-specific attributes

        """
        list_techs = []

        if archive_source:
            source = Source(
                title="Energy system technology data for the US",
                authors="Contributors to technology-data. Data source: manual_input_usa.csv",
                url="https://github.com/PyPSA/technology-data/blob/master/inputs/US/manual_input_usa.csv",
            )
            source.ensure_in_wayback()
            sources = SourceCollection(sources=[source])
            sources.to_json(sources_path, output_schema=output_schema)
        else:
            sources = SourceCollection.from_json(sources_path)

        for (scenario, year, technology), group in dataframe.groupby(
            ["scenario", "year", "technology"]
        ):
            parameters = {}
            financial_case_for_tech = None
            for _, row in group.iterrows():
                unit, carrier, heating_value = (
                    ManualInputUSAV0134Parser._extract_units_carriers_heating_value(
                        row["unit"]
                    )
                )
                param_kwargs = {
                    "magnitude": row["value"],
                    "sources": sources,
                }
                if carrier is not None:
                    param_kwargs["carrier"] = carrier
                if heating_value is not None:
                    param_kwargs["heating_value"] = heating_value
                if unit is not None:
                    param_kwargs["units"] = unit
                if row["further_description"] is not None and isinstance(
                    row["further_description"], str
                ):
                    param_kwargs["note"] = row["further_description"]
                if row["financial_case"] is not None and isinstance(
                    row["financial_case"], str
                ):
                    financial_case_for_tech = str(row["financial_case"])
                parameters[row["parameter"]] = Parameter(**param_kwargs)

            # Combine scenario and financial_case for the case attribute
            case_value = str(scenario)
            if financial_case_for_tech is not None:
                case_value = f"{scenario} - {financial_case_for_tech}"

            list_techs.append(
                Technology(
                    name=technology,
                    region="USA",
                    year=year,
                    parameters=parameters,
                    case=case_value,
                    detailed_technology=technology,
                )
            )

        return TechnologyCollection(technologies=list_techs)

    def parse(
        self,
        input_path: pathlib.Path,
        num_digits: int,
        archive_source: bool,
        **kwargs: Any,
    ) -> None:
        """
        Parse and process version 0.13.4 of the manual_input_usa.csv dataset.

        This method reads the raw data from an Excel file, cleans and transforms
        it through a series of steps, and then builds a TechnologyCollection.
        The processed data is saved to JSON files.

        Parameters
        ----------
        input_path : pathlib.Path
            Path to the raw input data file (Excel).
        num_digits : int
            Number of significant digits to round numerical values.
        archive_source : bool
            If True, archives the source object on the Wayback Machine.
        **kwargs : bool
            export_schema : bool
                If True, exports the Pydantic schema for the data models.

        Returns
        -------
        TechnologyCollection
            A collection of parsed technology data.

        """
        export_schema = kwargs.get("export_schema", False)

        manual_input_usa_input_path = pathlib.Path(input_path)

        manual_input_usa_df = pandas.read_csv(
            manual_input_usa_input_path, dtype=str, na_values="None"
        )
        manual_input_usa_df["value"] = manual_input_usa_df["value"].astype(float)
        manual_input_usa_df["scenario"] = manual_input_usa_df["scenario"].fillna(
            "not_available"
        )

        # Replace "per unit" with "%" and multiply val by 100
        mask_per_unit = manual_input_usa_df["unit"].str.contains("per unit")
        manual_input_usa_df.loc[mask_per_unit, "unit"] = manual_input_usa_df.loc[
            mask_per_unit, "unit"
        ].str.replace("per unit", "%")
        manual_input_usa_df.loc[mask_per_unit, "value"] = (
            manual_input_usa_df.loc[mask_per_unit, "value"] * 100.0
        ).round(num_digits)
        logger.info(
            "`per unit` replaced by `%`. Corresponding value multiplied by 100."
        )

        # Include currency_year in unit if applicable
        manual_input_usa_df["unit"] = manual_input_usa_df.apply(
            lambda row: Commons.update_unit_with_currency_year(
                row["unit"], row["currency_year"]
            ),
            axis=1,
        )
        logger.info("`currency_year` included in `unit` column.")

        # Build TechnologyCollection
        manual_input_usa_base_path = pathlib.Path(
            path_cwd,
            "src",
            "technologydata",
            "parsers",
            "manual_input_usa",
        )
        output_technologies_path = pathlib.Path(
            manual_input_usa_base_path,
            "v0.13.4/technologies.json",
        )
        output_sources_path = pathlib.Path(
            manual_input_usa_base_path,
            "v0.13.4/sources.json",
        )

        tech_col = ManualInputUSAV0134Parser._build_technology_collection(
            manual_input_usa_df,
            output_sources_path,
            archive_source=archive_source,
            output_schema=export_schema,
        )

        logger.info("TechnologyCollection object instantiated.")
        tech_col.to_json(output_technologies_path, output_schema=export_schema)
        logger.info("TechnologyCollection object exported to json.")
