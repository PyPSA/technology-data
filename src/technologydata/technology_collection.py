# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""TechnologyCollection class for representing an iterable of Technology Objects."""

import csv
import json
import pathlib
import re
import typing
from collections.abc import Iterator
from typing import Annotated, Self

import pandas
import pydantic
import pydantic_core

from technologydata.technology import Technology


class TechnologyCollection(pydantic.BaseModel):  # type: ignore
    """
    Represent a collection of technologies.

    Attributes
    ----------
    technologies : List[Technology]
        List of Technology objects.

    """

    technologies: Annotated[
        list[Technology], pydantic.Field(description="List of Technology objects.")
    ]

    def __iter__(self) -> Iterator[Technology]:
        """
        Return an iterator over the list of Technology objects.

        Returns
        -------
        Iterator[Technology]
            An iterator over the Technology objects contained in the collection.

        """
        return iter(self.technologies)

    def __len__(self) -> int:
        """
        Return the number of technologies in this collection.

        Returns
        -------
        int
            The number of Technology objects in the technologies list.

        """
        return len(self.technologies)

    def get(
        self, name: str, region: str, year: int, case: str, detailed_technology: str
    ) -> Self:
        """
        Filter technologies based on regex patterns for non-optional attributes.

        Parameters
        ----------
        name : str
            Regex pattern to filter technology names.
        region : str
            Regex pattern to filter region identifiers.
        year : int
            Regex pattern to filter the year of the data.
        case : str
            Regex pattern to filter case or scenario identifiers.
        detailed_technology : str
            Regex pattern to filter detailed technology names.

        Returns
        -------
        TechnologyCollection
            A new TechnologyCollection with filtered technologies.

        """
        filtered_technologies = self.technologies

        if name is not None:
            pattern_name = re.compile(name, re.IGNORECASE)
            filtered_technologies = [
                t for t in filtered_technologies if pattern_name.search(t.name)
            ]

        if region is not None:
            pattern_region = re.compile(region, re.IGNORECASE)
            filtered_technologies = [
                t for t in filtered_technologies if pattern_region.search(t.region)
            ]

        if year is not None:
            pattern_year = re.compile(str(year), re.IGNORECASE)
            filtered_technologies = [
                t for t in filtered_technologies if pattern_year.search(str(t.year))
            ]

        if case is not None:
            pattern_case = re.compile(case, re.IGNORECASE)
            filtered_technologies = [
                t for t in filtered_technologies if pattern_case.search(t.case)
            ]

        if detailed_technology is not None:
            pattern_detailed_technology = re.compile(detailed_technology, re.IGNORECASE)
            filtered_technologies = [
                t
                for t in filtered_technologies
                if pattern_detailed_technology.search(t.detailed_technology)
            ]

        return TechnologyCollection(technologies=filtered_technologies)

    def to_dataframe(self) -> pandas.DataFrame:
        """
        Convert the TechnologyCollection to a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the technology data.

        """
        return pandas.DataFrame(
            [technology.model_dump() for technology in self.technologies]
        )

    def to_csv(self, **kwargs: pathlib.Path | str | bool) -> None:
        """
        Export the TechnologyCollection to a CSV file.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed to pandas.DataFrame.to_csv().
            Common options include:
            - path_or_buf : str or pathlib.Path or file-like object, optional
                File path or object, if None, the result is returned as a string.
                Default is None.
            - sep : str
                String of length 1. Field delimiter for the output file.
                Default is ','.
            - index : bool
                Write row names (index). Default is True.
            - encoding : str
                String representing the encoding to use in the output file.
                Default is 'utf-8'.

        Notes
        -----
        The method converts the collection to a pandas DataFrame using
        `self.to_dataframe()` and then writes it to a CSV file using the provided
        kwargs.

        """
        default_kwargs = {
            "sep": ",",
            "index": False,
            "encoding": "utf-8",
            "quoting": csv.QUOTE_ALL,
        }

        # Merge default_kwargs with user-provided kwargs, giving precedence to user kwargs
        merged_kwargs = {**default_kwargs, **kwargs}
        output_dataframe = self.to_dataframe()
        output_dataframe.to_csv(**merged_kwargs)

    def to_json(
        self, file_path: pathlib.Path, schema_path: pathlib.Path | None = None
    ) -> None:
        """
        Export the TechnologyCollection to a JSON file, together with a data schema.

        Parameters
        ----------
        file_path : pathlib.Path
            The path to the JSON file to be created.
        schema_path : pathlib.Path
            The path to the JSON schema file to be created. By default, created with a `schema` suffix next to `file_path`.

        """
        if schema_path is None:
            schema_path = file_path.with_suffix(".schema.json")

        # Export the model's schema with descriptions to a dict
        schema = self.model_json_schema()

        # Save the schema (which includes descriptions) to a JSON file
        with open(schema_path, "w") as f:
            json.dump(schema, f, indent=4)

        with open(file_path, mode="w", encoding="utf-8") as jsonfile:
            json_data = self.model_dump_json(indent=4)  # Convert to JSON string
            jsonfile.write(json_data)

    @classmethod
    def from_json(cls, file_path: pathlib.Path | str) -> Self:
        """
        Load a TechnologyCollection instance from a JSON file.

        Parameters
        ----------
        file_path : pathlib.Path or str
            Path to the JSON file containing the data. Can be a pathlib.Path object or a string path.

        Returns
        -------
        TechnologyCollection
            An instance of TechnologyCollection initialized with the data from the JSON file.

        Raises
        ------
        TypeError
            If `file_path` is not a pathlib.Path or str.

        """
        if isinstance(file_path, (pathlib.Path | str)):
            file_path = pathlib.Path(file_path)
        else:
            raise TypeError("file_path must be a pathlib.Path or str")

        with open(file_path, encoding="utf-8") as jsonfile:
            json_data = jsonfile.read()

        # pydantic_core.from_json return Any. Therefore, typing.cast makes sure that
        # the output is indeed a TechnologyCollection
        return typing.cast(
            TechnologyCollection,
            cls.model_validate(pydantic_core.from_json(json_data, allow_partial=True)),
        )

    def to_currency(
        self,
        target_currency: str,
        overwrite_country: None | str = None,
        source: str = "worldbank",
    ) -> Self:
        """
        Adjust the currency of all parameters of all contained Technology objects to the target currency.

        The conversion includes inflation and exchange rates based on each Technology objects's region.
        If a different country should be used for inflation adjustment, use `overwrite_country`.

        Parameters
        ----------
        target_currency : str
            The target currency (e.g., 'EUR_2020').
        overwrite_country : str, optional
            ISO 3166 alpha-3 country code to use for inflation adjustment instead of the object's region.
        source: str, optional
            The source of the inflation data, either "worldbank"/"wb" or "international_monetary_fund"/"imf".
            Defaults to "worldbank".
            Depending on the source, different years to adjust for inflation may be available.

        Returns
        -------
        TechnologyCollection
            A new TechnologyCollection object with all its parameters adjusted to the target currency.

        """
        new_techs = []

        for i, tech in enumerate(self.technologies):
            new_techs.append(
                tech.to_currency(
                    target_currency=target_currency,
                    overwrite_country=overwrite_country,
                    source=source,
                )
            )

        return TechnologyCollection(technologies=new_techs)
