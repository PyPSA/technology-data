# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the Source class."""

import datetime
import pathlib

import pytest

import technologydata

path_cwd = pathlib.Path.cwd()


class TestSource:
    """Test suite for the Source class in the technologydata module."""

    @pytest.mark.parametrize(
        "example_source, expected_equal",
        [
            (
                {
                    "source_title": "atb_nrel",
                    "source_authors": "NREL/ATB",
                    "source_url": "https:download",
                    "source_url_archive": "http:/3273/download",
                    "source_url_date": "2025-05-22 15:08:02",
                    "source_url_date_archive": "2025-05-22 15:08:02",
                },
                False,  # Expect not equal
            ),
            (
                {
                    "source_title": "tech_data_generation",
                    "source_authors": "Danish Energy Agency",
                    "source_url": "https:download",
                    "source_url_archive": "http:/3273/download",
                    "source_url_date": "2025-05-06 16:02:04",
                    "source_url_date_archive": "2025-05-06 16:02:04",
                },
                True,  # Expect equal
            ),
            (
                {
                    "source_title": "tech_data_generation",
                    "source_authors": "Danish Energy Agency",
                },
                False,  # Expect not equal
            ),
        ],
        indirect=["example_source"],
    )  # type: ignore
    def test_eq(
        self, example_source: technologydata.Source, expected_equal: bool
    ) -> None:
        """Check if the override method __eq__ works as expected."""
        reference_source = technologydata.Source(
            title="tech_data_generation",
            authors="Danish Energy Agency",
            url="https:download",
            url_archive="http:/3273/download",
            url_date="2025-05-06 16:02:04",
            url_date_archive="2025-05-06 16:02:04",
        )
        assert (example_source == reference_source) == expected_equal

    @pytest.mark.parametrize(
        "example_source",
        [
            {
                "source_title": "atb_nrel",
                "source_authors": "NREL/ATB",
                "source_url_archive": "http:/3273/download",
                "source_url_date_archive": "2025-05-22 15:08:02",
            },
            {
                "source_title": "tech_data_generation",
                "source_authors": "Danish Energy Agency",
                "source_url": "https:download",
                "source_url_archive": "http:/3273/download",
                "source_url_date": "2025-05-06 16:02:04",
                "source_url_date_archive": "2025-05-06 16:02:04",
            },
        ],
        indirect=["example_source"],
    )  # type: ignore
    def test_hash(self, example_source: technologydata.Source) -> None:
        """Check if the override method __hash__ works as expected."""
        assert isinstance(hash(example_source), int)

    @pytest.mark.parametrize(
        "example_source, expected_string",
        [
            (
                {
                    "source_title": "OET project page",
                    "source_authors": "Open Energy Transition gGmbH",
                    "source_url": "https://outputs.html",
                    "source_url_archive": "https:archived.html",
                    "source_url_date": "2025-05-06 16:02:04",
                    "source_url_date_archive": "2025-05-06 16:02:04",
                },
                "'Open Energy Transition gGmbH': 'OET project page', from url 'https://outputs.html', last accessed on '2025-05-06 16:02:04', archived at 'https:archived.html', on '2025-05-06 16:02:04'.",
            )
        ],
        indirect=["example_source"],
    )  # type: ignore
    def test_str(
        self, example_source: technologydata.Source, expected_string: str
    ) -> None:
        """Check if the override method __str__ works as expected."""
        # Ensure the snapshot is created
        assert str(example_source) == expected_string

    @pytest.mark.parametrize(
        "example_source",
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
        indirect=["example_source"],
    )  # type: ignore
    def test_retrieve_from_wayback(self, example_source: technologydata.Source) -> None:
        """Check if the example source is downloaded from the Internet Archive Wayback Machine."""
        storage_path = example_source.retrieve_from_wayback(path_cwd)
        # Check if storage_paths is not None and is a list
        assert storage_path is not None, "Expected a valid storage path, but got None."
        assert isinstance(storage_path, pathlib.Path), (
            "Expected storage_path to be a pathlib.Path."
        )
        assert storage_path is not None, "Expected a valid storage path, but got None."
        assert storage_path.is_file(), (
            f"Expected {storage_path} to be a file, but it does not exist."
        )
        # Delete the downloaded file
        storage_path.unlink(missing_ok=True)

    def test_store_in_wayback(self) -> None:
        """Check if a given url is correctly stored as a snapshot on Internet Archive Wayback Machine."""
        url_to_archive = "https://www.openenergytransition.org/outputs.html"
        archived_info = technologydata.Source.store_in_wayback(url_to_archive)

        # Check if archived_info is None
        assert archived_info is not None, "archived_info should not be None"

        archived_url, new_capture, output_timestamp = archived_info

        assert "https://web.archive.org/web/" in archived_url
        assert url_to_archive in archived_url
        assert isinstance(new_capture, bool)

        assert output_timestamp is not None, "output_timestamp should not be None"
        try:
            datetime.datetime.strptime(
                output_timestamp, technologydata.DateFormatEnum.SOURCES_CSV
            )
        except ValueError:
            pytest.fail("Valid date-time string did not match the format")

    @pytest.mark.parametrize(
        "example_source",
        [
            {
                "source_title": "OET project page",
                "source_authors": "Open Energy Transition gGmbH",
                "source_url": "https://openenergytransition.org/outputs.html",
            },
        ],
        indirect=["example_source"],
    )  # type: ignore
    def test_ensure_in_wayback(self, example_source: technologydata.Source) -> None:
        """Check if the snapshot URL is created correctly."""
        # Ensure the snapshot is created
        example_source.ensure_in_wayback()
        assert example_source.url_date_archive is not None
        assert example_source.url_archive is not None

    @pytest.mark.parametrize(
        "url_archived, source_path, source_title, expected_path",
        [
            (
                "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                path_cwd,
                "title",
                pathlib.Path(path_cwd, "title.pdf"),
            ),
            (
                "https://web.archive.org/web/20250522150802/https://oedi-data-lake.s3.amazonaws.com/ATB/electricity/parquet/2024/v3.0.0/ATBe.parquet",
                path_cwd,
                "title",
                pathlib.Path(path_cwd, "title.parquet"),
            ),
        ],
    )  # type: ignore
    def test_get_save_path(
        self,
        url_archived: str,
        source_path: pathlib.Path,
        source_title: str,
        expected_path: pathlib.Path,
    ) -> None:
        """Check if the path where to store the file to download follows the requested pattern."""
        assert (
            technologydata.Source._get_save_path(
                url_archived, source_path, source_title
            )
            == expected_path
        )
