"""technologydata: A package for managing and analyzing technology data used for energy system models."""

from importlib.metadata import version

from technologydata.sources import AVAILABLE_SOURCES, Source, Sources
from technologydata.technologies import (
    Technologies,
    check_source_validity,
)

__version__ = version("technologydata")


__all__ = [
    "Source",
    "Sources",
    "Technologies",
    "AVAILABLE_SOURCES",
    "check_source_validity",
]
