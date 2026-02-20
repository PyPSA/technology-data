# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Provide a class to access data from a data source."""

import logging
import pathlib
import re
import sys
from enum import Enum
from typing import Annotated

import pydantic
from packaging.version import parse
from pydantic import field_validator

from technologydata import DataPackage
from technologydata.parsers.dea_energy_storage import DeaEnergyStorageParser
from technologydata.parsers.manual_input_usa import ManualInputUsaParser

path_cwd = pathlib.Path.cwd()

logger = logging.getLogger(__name__)


class DataSourceName(str, Enum):
    """An enumeration of available data sources."""

    DEA_ENERGY_STORAGE = "dea_energy_storage"
    MANUAL_INPUT_USA = "manual_input_usa"


class DataAccessor(pydantic.BaseModel):
    """
    Access data from a versioned data source.

    This class provides a standardized interface to locate and load technology
    datasets from predefined data sources. It can either load a specific version
    or automatically determine and load the latest available version.

    Attributes
    ----------
    data_source : str
        The name of the data source to access, as defined in the
        `DataSourceName` enumeration.
    version : str, optional
        The specific version string of the data to load (e.g., "v1.0.0").
        If not provided, the latest version will be automatically determined
        and used. Default is None.

    """

    data_source: Annotated[
        str, pydantic.Field(description="The name of the data source.")
    ]
    version: Annotated[
        str | None, pydantic.Field(description="The version of the data source.")
    ] = None

    @field_validator("data_source", mode="before")
    @classmethod
    def _validate_data_source_name(cls, v: str) -> DataSourceName:
        # Validate if the given string is a valid DataSourceName
        try:
            return DataSourceName(v)
        except ValueError:
            raise ValueError(
                f"{v} is not a valid DataSourceName. Available options: {[e for e in DataSourceName]}"
            )

    @staticmethod
    def get_latest_version_string(data_source_path_list: list[pathlib.Path]) -> str:
        """
        Find the latest version string for the data source.

        Returns
        -------
        str
            The string of the latest version (e.g., 'v10', 'v1.0.0').

        Raises
        ------
        FileNotFoundError
            If the data source directory or valid version directories are not found.

        """
        version_pattern = re.compile(r"^v(\d+(\.\d+)*)$")
        versions = []
        for item in data_source_path_list:
            if item.is_dir():
                match = version_pattern.match(item.name)
                if match:
                    versions.append(item.name)

        if not versions:
            raise FileNotFoundError("No valid version directories found.")

        latest_version_str = max(versions, key=lambda v: parse(v[1:]))
        return latest_version_str

    def load(self) -> DataPackage:
        """
        Load the default 'technologies.json' from the package data.

        Returns
        -------
        DataPackage
            An instance of DataPackage initialized with the requested data.

        Raises
        ------
        FileNotFoundError
            If the data source directory or the specified version directory is not found.
        ValueError
            If the specified version is not found. The user is notified of  the latest available version.

        """
        source_path = pathlib.Path(
            path_cwd, "src", "technologydata", "parsers", self.data_source
        )
        if not source_path.is_dir():
            raise FileNotFoundError(f"Data source directory not found: {source_path}")

        source_path_list = [p.name for p in source_path.iterdir() if p.is_dir()]

        if self.version and self.version in source_path_list:
            version = self.version
            logger.info(
                f"Data source directory corresponding to version {self.version} found."
            )
        else:
            version = self.get_latest_version_string(list(source_path.iterdir()))
            raise ValueError(
                f"Data source version '{self.version}' not found. The latest available version is {version}."
            )

        data_path = pathlib.Path(source_path, version)
        dp = DataPackage.from_json(self.data_source, self.version, data_path)
        return dp

    def parse(
        self,
        input_file_name: str,
        num_digits: int = 4,
        archive_source: bool = False,
        filter_params: bool = False,
        export_schema: bool = False,
    ) -> None:
        """
        Run the parser for the specified data source and version.

        This method locates the appropriate parser for the given data source
        and version, and executes it to generate the technology data package.

        Parameters
        ----------
        input_file_name : str
            The name of the input file in the 'raw' directory.
        num_digits : int, optional
            Number of significant digits to round the values. Default is 4.
        archive_source : bool, optional
            Store the source object on the Wayback Machine. Default is False.
        filter_params : bool, optional
            Filter the parameters stored to technologies.json. Default is False.
        export_schema : bool, optional
            Export the Source/TechnologyCollection schemas. Default is False.

        Raises
        ------
        ValueError
            If the specified data source or version is not supported.
        FileNotFoundError
            If the required input data file is not found.

        """
        parser: DeaEnergyStorageParser | ManualInputUsaParser

        if self.data_source == DataSourceName.DEA_ENERGY_STORAGE:
            parser = DeaEnergyStorageParser()
        elif self.data_source == DataSourceName.MANUAL_INPUT_USA:
            parser = ManualInputUsaParser()
        else:
            raise ValueError(
                f"Unsupported data source: {self.data_source}. "
                f"Supported data sources are: {[e for e in DataSourceName]}"
            )

        # Read the raw data
        input_data_path = pathlib.Path(
            path_cwd,
            "src",
            "technologydata",
            "parsers",
            "raw",
            input_file_name,
        )

        logger.info(f"Input data path set to: {input_data_path}")

        if self.version not in parser.get_supported_versions():
            logging.error(
                f"Version '{self.version}' is not supported. "
                f"Supported versions: {parser.get_supported_versions()}"
            )
            sys.exit(1)

        try:
            parser.parse(
                version=self.version,
                input_path=input_data_path,
                num_digits=num_digits,
                archive_source=archive_source,
                filter_params=filter_params,
                export_schema=export_schema,
            )

            logging.info(f"Successfully generated files for version {self.version} ")

        except (ValueError, FileNotFoundError, KeyError) as e:
            logging.error(f"An error occurred during parsing: {e}")
            sys.exit(1)
