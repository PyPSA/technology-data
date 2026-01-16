# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the SourceCollection class."""

import pathlib

import pandas
import pytest

import technologydata

path_cwd = pathlib.Path.cwd()


class TestSourceCollection:
    """Test suite for the SourceCollection class in the technologydata module."""

    @pytest.mark.parametrize(
        "example_source_collection, expected_string",
        [
            (
                [
                    {
                        "source_title": "atb_nrel",
                        "source_authors": "NREL/ATB",
                        "source_url_date": "2025-05-22 15:08:02",
                    },
                    {
                        "source_title": "tech_data_generation",
                        "source_authors": "Danish Energy Agency",
                        "source_url_date_archive": "2025-05-06 16:02:04",
                    },
                ],
                "SourceCollection with 2 sources: "
                "'NREL/ATB': 'atb_nrel', last accessed on '2025-05-22 15:08:02', "
                "'Danish Energy Agency': 'tech_data_generation', on '2025-05-06 16:02:04'.",
            ),
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_str(
        self,
        example_source_collection: technologydata.SourceCollection,
        expected_string: str,
    ) -> None:
        """Check if the example source collection is cast to string as expected."""
        assert str(example_source_collection) == expected_string

    @pytest.mark.parametrize(
        "example_source_collection",
        [
            [
                {
                    "source_title": "Source 1",
                    "source_authors": "Author 1",
                },
                {
                    "source_title": "Source 2",
                    "source_authors": "Author 2",
                },
            ]
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_example_source_collection(
        self,
        example_source_collection: technologydata.SourceCollection,
    ) -> None:
        """Check if the example source collection is instantiated correctly."""
        # Check that the returned object is a SourceCollection
        assert isinstance(example_source_collection, technologydata.SourceCollection)

        # Check the number of sources
        assert len(example_source_collection.sources) == 2

        # Check the titles of the sources
        assert example_source_collection.sources[0].title == "Source 1"
        assert example_source_collection.sources[1].title == "Source 2"
        assert example_source_collection.sources[0].authors == "Author 1"
        assert example_source_collection.sources[1].authors == "Author 2"

    @pytest.mark.webarchive  # type: ignore
    @pytest.mark.parametrize(
        "example_source_collection",
        [
            [
                {
                    "source_title": "atb_nrel",
                    "source_authors": "NREL/ATB",
                    "source_url": "https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_archive": "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_date": "2025-05-22 15:08:02",
                    "source_url_date_archive": "2025-05-22 15:08:02",
                },
                {
                    "source_title": "tech_data_generation",
                    "source_authors": "Danish Energy Agency",
                    "source_url": "https://ens.dk/media/3273/download",
                    "source_url_archive": "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                    "source_url_date": "2025-05-06 16:02:04",
                    "source_url_date_archive": "2025-05-06 16:02:04",
                },
            ],
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_retrieve_all_from_wayback(
        self,
        example_source_collection: technologydata.SourceCollection,
    ) -> None:
        """Check if the example source collection is downloaded from the Internet Archive Wayback Machine."""
        storage_paths = example_source_collection.retrieve_all_from_wayback(path_cwd)

        # Check if storage_paths is not None and is a list
        assert storage_paths is not None, "Expected storage_paths to be not None."
        assert isinstance(storage_paths, list), "Expected storage_paths to be a list."

        # Filter out None values and check the remaining paths
        valid_paths = [path for path in storage_paths if path is not None]

        assert all(isinstance(path, pathlib.Path) for path in valid_paths), (
            f"Expected all elements in {valid_paths} to be instances of pathlib.Path."
        )
        assert all(path.is_file() for path in valid_paths), (
            f"Expected all elements in {valid_paths} to be a file, but it does not exist."
        )

        # Delete the downloaded files
        for path in valid_paths:
            path.unlink(missing_ok=True)

    @pytest.mark.parametrize(
        "example_source_collection",
        [
            [
                {
                    "source_title": "atb_nrel",
                    "source_authors": "NREL/ATB",
                    "source_url": "https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_archive": "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_date": "2025-05-22 15:08:02",
                    "source_url_date_archive": "2025-05-22 15:08:02",
                },
                {
                    "source_title": "tech_data_generation",
                    "source_authors": "Danish Energy Agency",
                    "source_url": "https://ens.dk/media/3273/download",
                    "source_url_archive": "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                    "source_url_date": "2025-05-06 16:02:04",
                    "source_url_date_archive": "2025-05-06 16:02:04",
                },
            ],
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_to_csv(
        self, example_source_collection: technologydata.SourceCollection
    ) -> None:
        """Check if the example source collection is exported to CSV."""
        output_file = pathlib.Path(path_cwd, "export.csv")
        example_source_collection.to_csv(path_or_buf=output_file)
        assert output_file.is_file()
        output_file.unlink(missing_ok=True)

    # python
    @pytest.mark.parametrize(
        "example_source_collection, output_schema",
        [
            (
                [
                    {
                        "source_title": "atb_nrel",
                        "source_authors": "NREL/ATB",
                        "source_url": "https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                        "source_url_archive": "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                        "source_url_date": "2025-05-22 15:08:02",
                        "source_url_date_archive": "2025-05-22 15:08:02",
                    },
                    {
                        "source_title": "tech_data_generation",
                        "source_authors": "Danish Energy Agency",
                        "source_url": "https://ens.dk/media/3273/download",
                        "source_url_archive": "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                        "source_url_date": "2025-05-06 16:02:04",
                        "source_url_date_archive": "2025-05-06 16:02:04",
                    },
                ],
                True,
            ),
            (
                [
                    {
                        "source_title": "atb_nrel",
                        "source_authors": "NREL/ATB",
                        "source_url": "https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                        "source_url_archive": "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                        "source_url_date": "2025-05-22 15:08:02",
                        "source_url_date_archive": "2025-05-22 15:08:02",
                    },
                    {
                        "source_title": "tech_data_generation",
                        "source_authors": "Danish Energy Agency",
                        "source_url": "https://ens.dk/media/3273/download",
                        "source_url_archive": "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                        "source_url_date": "2025-05-06 16:02:04",
                        "source_url_date_archive": "2025-05-06 16:02:04",
                    },
                ],
                False,
            ),
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_to_json(
        self,
        example_source_collection: technologydata.SourceCollection,
        output_schema: bool,
    ) -> None:
        """Check if the example source collection is exported to JSON."""
        output_file = pathlib.Path(path_cwd, "sources.json")
        schema_file = pathlib.Path(path_cwd, "sources.schema.json")

        example_source_collection.to_json(
            pathlib.Path(output_file), output_schema=output_schema
        )

        assert output_file.is_file()
        if output_schema:
            assert schema_file.is_file()
        else:
            assert not schema_file.is_file()

        output_file.unlink(missing_ok=True)
        schema_file.unlink(missing_ok=True)

    @pytest.mark.parametrize(
        "example_source_collection",
        [
            [
                {
                    "source_title": "atb_nrel",
                    "source_authors": "NREL/ATB",
                    "source_url": "https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_archive": "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                    "source_url_date": "2025-05-22 15:08:02",
                    "source_url_date_archive": "2025-05-22 15:08:02",
                },
                {
                    "source_title": "tech_data_generation",
                    "source_authors": "Danish Energy Agency",
                    "source_url": "https://ens.dk/media/3273/download",
                    "source_url_archive": "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                    "source_url_date": "2025-05-06 16:02:04",
                    "source_url_date_archive": "2025-05-06 16:02:04",
                },
            ],
        ],
        indirect=["example_source_collection"],
    )  # type: ignore
    def test_to_dataframe(
        self,
        example_source_collection: technologydata.SourceCollection,
    ) -> None:
        """Check if the example source collection is exported to pandas dataframe."""
        assert isinstance(example_source_collection.to_dataframe(), pandas.DataFrame)

    def test_from_json(self) -> None:
        """Check if the example source collection is imported from JSON."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "sources.json",
        )
        source_collection = technologydata.SourceCollection.from_json(
            file_path=input_file
        )
        assert isinstance(source_collection, technologydata.SourceCollection)
        assert len(source_collection) == 2

    def test_from_json_to_json(self) -> None:
        """Check whether reading with from_json and exporting with to_json yields the same file."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "sources.json",
        )
        source_collection = technologydata.SourceCollection.from_json(input_file)
        output_file = pathlib.Path("to_json_test.json")
        schema_file = pathlib.Path(path_cwd, "to_json_test.schema.json")
        source_collection.to_json(output_file, output_schema=True)
        # Read files and strip trailing whitespace/newlines before comparing
        with open(input_file) as f1, open(output_file) as f2:
            assert f1.read().rstrip() == f2.read().rstrip(), "Files are not identical"
        assert output_file.is_file()
        assert schema_file.is_file()
        output_file.unlink(missing_ok=True)
        schema_file.unlink(missing_ok=True)

    @pytest.mark.parametrize(
        "title_pattern, authors_pattern",
        [["ATB", "nrel"], ["TECH_DATA", "danish"]],
    )  # type: ignore
    def test_get(self, title_pattern: str, authors_pattern: str) -> None:
        """Check if the example source collection is filtered as requested."""
        input_file = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "solar_photovoltaics_example",
            "sources.json",
        )
        source_collection = technologydata.SourceCollection.from_json(
            file_path=input_file
        )
        result = source_collection.get(title=title_pattern, authors=authors_pattern)
        assert len(result.sources) == 1
