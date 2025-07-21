# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""SourceCollection class for representing an iterable of Source Objects."""

import csv
import json
import pathlib
import re

import pandas
import pydantic

from technologydata.source import Source


class SourceCollection(pydantic.BaseModel):  # type: ignore
    """
    Represent a collection of sources.

    Parameters
    ----------
    sources : List[Source]
        List of Source objects.

    Attributes
    ----------
    sources : List[Source]
        List of Source objects.

    """

    sources: list[Source] = pydantic.Field(..., description="List of Source objects.")

    def get(self, title: str, authors: str) -> "SourceCollection":
        """
        Filter sources based on regex patterns for non-optional attributes.

        Parameters
        ----------
        title : str
            Regex pattern to filter titles.
        authors : str
            Regex pattern to filter authors.

        Returns
        -------
        SourceCollection
            A new SourceCollection with filtered sources.

        """
        filtered_sources = self.sources

        if title is not None:
            pattern_title = re.compile(title, re.IGNORECASE)
            filtered_sources = [
                s for s in filtered_sources if pattern_title.search(s.title)
            ]

        if authors is not None:
            pattern_authors = re.compile(authors, re.IGNORECASE)
            filtered_sources = [
                s for s in filtered_sources if pattern_authors.search(s.authors)
            ]

        return SourceCollection(sources=filtered_sources)

    def __len__(self) -> int:
        """
        Return the number of sources in this collection.

        Returns
        -------
        int
            The number of Source objects in the sources list.

        """
        return len(self.sources)

    def retrieve_all_from_wayback(
        self, download_directory: pathlib.Path
    ) -> list[pathlib.Path | None]:
        """
        Download archived files for all sources in the collection using retrieve_from_wayback.

        Parameters
        ----------
        download_directory : pathlib.Path
            The base directory where all files will be saved.

        Returns
        -------
        list[pathlib.Path | None]
            List of paths where each file was stored, or None if download failed for a source.

        """
        return [
            source.retrieve_from_wayback(download_directory) for source in self.sources
        ]

    def to_dataframe(self) -> pandas.DataFrame:
        """
        Convert the SourceCollection to a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the source data.

        """
        return pandas.DataFrame([source.model_dump() for source in self.sources])

    def to_csv(self, **kwargs: pathlib.Path | str | bool) -> None:
        """
        Export the SourceCollection to a CSV file.

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
        Export the SourceCollection to a JSON file, together with a data schema.

        Parameters
        ----------
        file_path : pathlib.Path
            The path to the JSON file to be created.
        schema_path : pathlib.Path
            The path to the JSON schema file to be created. By default created with a `schema` suffix next to `file_path`.

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
    def from_json(cls, file_path: pathlib.Path | str) -> "SourceCollection":
        """
        Import the SourceCollection from a JSON file.

        Parameters
        ----------
        file_path : pathlib.Path
            The path to the JSON file to be imported.

        """
        if isinstance(file_path, pathlib.Path | str):
            file_path = pathlib.Path(file_path)
        else:
            raise TypeError("file_path must be a pathlib.Path or str")
        with open(file_path, encoding="utf-8") as jsonfile:
            json_data = json.load(jsonfile)
        return cls(**json_data)
