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

import pandas
import pydantic

from technologydata.technology import Technology


class TechnologyCollection(pydantic.BaseModel):  # type: ignore
    """
    Represent a collection of technologies.

    Attributes
    ----------
    technologies : List[Technology]
        List of Technology objects.

    """

    technologies: typing.Annotated[
        list[Technology], pydantic.Field(description="List of Technology objects.")
    ]

    def __iter__(self) -> Iterator["Technology"]:
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
    ) -> "TechnologyCollection":
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
    def from_json(cls, file_path: pathlib.Path | str) -> "TechnologyCollection":
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

        Notes
        -----
        This method reads the JSON data from the specified file, creates `Technology` objects
        for each item in the JSON list using `Technology.from_dict()`, and returns a new
        `TechnologyCollection` containing these objects.

        """
        if isinstance(file_path, pathlib.Path | str):
            file_path = pathlib.Path(file_path)
        else:
            raise TypeError("file_path must be a pathlib.Path or str")
        with open(file_path, encoding="utf-8") as jsonfile:
            json_data = json.load(jsonfile)
        techs = []
        for item in json_data:
            techs.append(Technology.from_dict(item))
        return cls(technologies=techs)
