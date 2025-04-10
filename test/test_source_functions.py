"""Test functionality related to the Source and Sources classes."""

import technologydata as td


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
    technologies = td.Technologies()
    assert technologies.sources == td.AVAILABLE_SOURCES, (
        "The default constructor should load all packaged sources"
    )

    technologies = td.Technologies(packaged_sources="all")
    assert technologies.sources == td.AVAILABLE_SOURCES, (
        "Using 'all' should load all packaged sources"
    )


def test_source_adding() -> None:
    """Check if the example source can be added correctly."""
    technologies = td.Technologies()
    technologies.add_source("example01", td.AVAILABLE_SOURCES["example01"])
    assert "example01" in technologies.sources, (
        "example01 should be in the sources after adding it"
    )


def test_additional_source() -> None:
    """Check if an additional source can be added correctly."""
    additional_source = "example01"
    technologies = td.Technologies(
        packaged_sources=[],
        additional_sources={additional_source: td.AVAILABLE_SOURCES[additional_source]},
    )
    assert additional_source in technologies.sources, (
        f"additional source '{additional_source}' should be in the sources after adding it"
    )
