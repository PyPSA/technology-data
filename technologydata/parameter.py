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

import typing

import pydantic

from technologydata.source_collection import SourceCollection
from technologydata.unit_value import UnitValue


# TODO rework class logic
class Parameter(pydantic.BaseModel):  # type: ignore
    """
    Encapsulate a value with its unit, provenance, notes, and sources.

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

    quantity: typing.Annotated[
        UnitValue, pydantic.Field(description="The value and its unit.")
    ]
    provenance: typing.Annotated[
        str | None, pydantic.Field(description="Data provenance.")
    ] = None
    note: typing.Annotated[
        str | None, pydantic.Field(description="Additional notes.")
    ] = None
    sources: typing.Annotated[
        SourceCollection, pydantic.Field(description="Collection of Sources.")
    ]

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

    @classmethod
    def from_dict(cls, data: dict[str, typing.Any]) -> "Parameter":
        """
        Create an instance of the class from a dictionary.

        Parameters
        ----------
        cls : type
            The class to instantiate.
        data : dict
            A dictionary containing the data to initialize the class instance.
            Expected keys include:
                - "quantity" (dict): A dictionary representing a UnitValue, parsed via `UnitValue.parse_obj()`.
                - "provenance" (str or None): Optional provenance information.
                - "note" (str or None): Optional notes.
                - "sources" (list): A list of source data dictionaries, to be converted into a SourceCollection.

        Returns
        -------
        instance : cls
            An instance of the class initialized with the provided data.

        Notes
        -----
        This method converts the "sources" list into a `SourceCollection` using `SourceCollection.from_json()`.
        The "quantity" field is parsed into a `UnitValue` object using `UnitValue.parse_obj()`.

        """
        # Convert sources list into SourceCollection
        sources_data = data.get("sources", [])
        sources = SourceCollection.from_json(from_str=sources_data)
        return cls(
            quantity=UnitValue.model_validate(data["quantity"]),
            provenance=data.get("provenance"),
            note=data.get("note"),
            sources=sources,
        )
