# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the Datapackage class."""

import pathlib

import technologydata

path_cwd = pathlib.Path.cwd()


class TestDatapackage:
    """Test suite for the Datapackage class in the technologydata module."""

    def test_get_source_collection(self) -> None:
        """Test how the sources attributes is extracted from a TechnologyCollection instance."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "technologies.json",
        )

        technologies_collection = technologydata.TechnologyCollection.from_json(
            input_file
        )
        data_package = technologydata.DataPackage(
            name="solar_photovoltaics_example",
            version="v1.0",
            technologies=technologies_collection,
        )
        data_package.get_source_collection()
        assert isinstance(data_package.sources, technologydata.SourceCollection)
        assert len(data_package.sources) == 2
        assert data_package.name == "solar_photovoltaics_example"
        assert data_package.version == "v1.0"
