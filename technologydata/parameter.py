# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""
Parameter class for encapsulating a value, its unit, provenance, notes, and sources.

Examples
--------
>>> from technologydata.unit_value import UnitValue
>>> from technologydata.source import Source
>>> uv = UnitValue(value=1000, unit="EUR_2020/kW")
>>> src = Source(name="Example Source", authors="some authors", url="http://example.com")
>>> param = Parameter(quantity=uv, provenance="literature", note="Estimated", sources=[src])

"""

import pydantic

from technologydata.source_collection import SourceCollection
from technologydata.unit_value import UnitValue


# TODO rework class logic
class Parameter(pydantic.BaseModel):  # type: ignore
    """
    Encapsulate a value with its unit, provenance, notes, and sources.

    Parameters
    ----------
    quantity : UnitValue
        The value and its unit.
    provenance : Optional[str]
        Description of the data's provenance.
    note : Optional[str]
        Additional notes about the parameter.
    sources : SourceCollection
        One or more sources for the parameter.

    Attributes
    ----------
    quantity : UnitValue
        The value and its unit.
    provenance : Optional[str]
        Description of the data's provenance.
    note : Optional[str]
        Additional notes about the parameter.
    sources : SourceCollection
        List of sources for the parameter.

    """

    quantity: UnitValue = pydantic.Field(..., description="The value and its unit.")
    provenance: str | None = pydantic.Field(None, description="Data provenance.")
    note: str | None = pydantic.Field(None, description="Additional notes.")
    sources: SourceCollection = pydantic.Field(
        ..., description="Collection of Sources."
    )

    @property
    def value(self) -> float:
        """
        The numerical value of the parameter.

        Returns
        -------
        float
            The value.

        """
        return self.quantity.value

    @property
    def unit(self) -> str:
        """
        The unit of the parameter.

        Returns
        -------
        str
            The unit.

        """
        return self.quantity.unit
