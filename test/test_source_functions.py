"""Test functionality related to the Source and Sources classes."""

import pytest

import technologydata as td
from technologydata import Source


def test_available_sources() -> None:
    """Check if the example source is correctly detected."""
    assert "example01" in td.AVAILABLE_SOURCES, (
        "example01 should be in the available sources"
    )


def test_source_validity() -> None:
    """Check if the example source is correctly detected as valid."""
    assert td.check_source_validity(
        td.AVAILABLE_SOURCES["example01"], "technologies", report_invalid=True
    ), "example01 should be a valid source"


def test_packaged_sources() -> None:
    """Check if by default the default constructors properly load all packaged sources."""
    technologies = td.Technologies(td.AVAILABLE_SOURCES)
    assert len(technologies.sources.sources) == len(td.AVAILABLE_SOURCES.keys()), (
        "The default constructor should load all packaged sources"
    )


def test_archive_file_on_internet_archive() -> None:
    url = "https://ens.dk/media/3273/download"
    returned_url = Source.archive_file_on_internet_archive(url)
    assert False


@pytest.mark.parametrize(
    "input_date, input_format, output_format, expected_date",
    [
        (
            "20250507105201",
            "%Y%m%d%H%M%S",
            "%Y-%m-%d %H:%M:%S",
            "2025-05-07 10:52:01",
        ),
        (
            "2025-05-07 10:52:01",
            "%Y-%m-%d %H:%M:%S",
            "%Y%m%d%H%M%S",
            "20250507105201",
        ),
        (
            "2025-05-07 10:52:01",
            "%Y-%m-%d",
            "%Y%m%d%H%M%S",
            None,
        ),
    ],
)
def test_change_datetime_format(
    input_date, input_format, output_format, expected_date
) -> None:
    output_date = Source.change_datetime_format(input_date, input_format, output_format)
    assert output_date == expected_date


@pytest.mark.parametrize(
    "url, expected",
    [
        (
            "https://ens.dk/media/3273/download",
            (
                "http://web.archive.org/web/20250506160204/https://ens.dk/media/3273/download",
                "2025-05-06 16:02:04",
                "200",
            ),
        ),
        ("https://ens.dk/media/1/download", None),
    ],
)
def test_get_wayback_snapshot(url, expected) -> None:
    output = Source.get_wayback_snapshot(url)
    assert output == expected
