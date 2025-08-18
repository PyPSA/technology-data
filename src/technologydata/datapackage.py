# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""
DataPackage class for managing collections of Technology objects and batch operations.

Examples
--------
>>> dp = DataPackage.from_json("path/to/data_package.json")

"""

import logging
import pathlib
import typing

import pydantic

from technologydata.source_collection import SourceCollection
from technologydata.technology_collection import TechnologyCollection

logger = logging.getLogger(__name__)


class DataPackage(pydantic.BaseModel):  # type: ignore
    """
    Container for a collection of Technology objects and/or Source objects, with batch operations and loading utilities.

    Attributes
    ----------
    technologies : Optional[TechnologyCollection]
        List of Technology objects.
    sources : Optional[SourceCollection]
        List of Source objects.

    """

    technologies: typing.Annotated[
        TechnologyCollection | None,
        pydantic.Field(description="List of Technology objects."),
    ] = None
    sources: typing.Annotated[
        SourceCollection | None, pydantic.Field(description="List of Source objects.")
    ] = None

    def get_source_collection(self) -> None:
        """
        Get the SourceCollection associated with this DataPackage from the TechnologyCollection.

        Returns
        -------
        SourceCollection
            The SourceCollection instance.

        """
        sources_set = set()
        if self.sources is None:
            if self.technologies is not None:
                for technology in self.technologies:
                    for parameter in technology.parameters.values():
                        for source in parameter.sources:
                            sources_set.add(source)
            else:
                raise ValueError(
                    "No technologies available to extract a sources collection from."
                )
        else:
            logger.info("The data package already has a sources collection.")
        self.sources = SourceCollection(sources=list(sources_set))

    @classmethod
    def from_json(cls, path_to_folder: pathlib.Path | str) -> "DataPackage":
        """
        Load a DataPackage from a JSON file.

        Parameters
        ----------
        path_to_folder : pathlib.Path or str
            Path to the data package folder.

        Returns
        -------
        DataPackage
            The loaded DataPackage instance.

        """
        # Load technologies
        technologies_path = pathlib.Path(path_to_folder, "technologies.json")
        technologies = TechnologyCollection.from_json(technologies_path)

        # Create DataPackage instance
        data_package = cls(
            path=path_to_folder,
            technologies=technologies,
        )

        # Generate sources collection from technologies if sources is None
        data_package.get_source_collection()

        return data_package

    def to_json(self, folder_path: pathlib.Path) -> None:
        """
        Export the Datapackage to JSON files.

        The files are:
         - the 'technologies' attribute is exported to technologies.json, together with the corresponding data schema
         - the 'sources' attribute is exported to sources.json, together with the corresponding data schema

        Parameters
        ----------
        folder_path : pathlib.Path
            The path to the folder where the JSON files are created

        """
        if self.technologies is not None:
            technologies_path = pathlib.Path(folder_path, "technologies.json")
            self.technologies.to_json(technologies_path)

        if self.sources is not None:
            sources_path = pathlib.Path(folder_path, "sources.json")
            self.sources.to_json(sources_path)

    def to_csv(self, folder_path: pathlib.Path) -> None:
        """
        Export the Datapackage to CSV files.

        The files are:
         - the 'technologies' attribute is exported to technologies.csv
         - the 'sources' attribute is exported to sources.csv

        Parameters
        ----------
        folder_path : pathlib.Path
            The path to the folder where the CSV files are created

        """
        if self.technologies is not None:
            technologies_path = pathlib.Path(folder_path, "technologies.csv")
            self.technologies.to_csv(path_or_buf=technologies_path)

        if self.sources is not None:
            sources_path = pathlib.Path(folder_path, "sources.csv")
            self.sources.to_csv(path_or_buf=sources_path)
