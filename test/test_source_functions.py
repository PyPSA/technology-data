"""Test functionality related to the Source and Sources classes."""

from technologydata import AVAILABLE_SOURCES, Technologies, check_source_validity


def test_available_sources() -> None:
    """Check if the example source is correctly detected."""
    assert "example01" in AVAILABLE_SOURCES, (
        "example01 should be in the available sources"
    )


def test_source_validity() -> None:
    """Check if the example source is correctly detected as valid."""
    assert check_source_validity(
        AVAILABLE_SOURCES["example01"], "technologies", report_invalid=True
    ), "example01 should be a valid source"


def test_packaged_sources() -> None:
    """Check if by default the default constructors properly load all packaged sources."""
    technologies = Technologies(AVAILABLE_SOURCES)
    assert len(technologies.sources.sources) == len(AVAILABLE_SOURCES.keys()), (
        "The default constructor should load all packaged sources"
    )
