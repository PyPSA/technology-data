# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT
# TODO implement (Johannes)
"""
UnitValue class for representing a value with an associated unit.

This class is designed to support flexible units, including energy carriers, currencies with years, and more.

Examples
--------
>>> uv = UnitValue(value=100, unit="EUR_2020")
>>> uv.value
100
>>> uv.unit
"EUR_2020"

"""

import pint
import pydantic

ureg = pint.UnitRegistry()


class UnitValue(pydantic.BaseModel):  # type: ignore
    """
    Represent a numerical value with an associated unit of measurement.

    Parameters
    ----------
    value : float
        The numerical value.
    unit : str
        The unit of measurement (e.g., 'EUR_2020', 'kWh_electricity', 'kWh_hydrogen_LHV').

    Attributes
    ----------
    value : float
        The numerical value.
    unit : str
        The unit of measurement.

    """

    value: float = pydantic.Field(..., description="The numerical value.")
    unit: str = pydantic.Field(..., description="The unit of measurement.")

    def to(self, new_unit: str) -> "UnitValue":
        """
        Convert the value to a new unit using pint.

        Parameters
        ----------
        new_unit : str
            The unit to convert to.

        Returns
        -------
        UnitValue
            A new UnitValue instance with the converted value and unit.

        Raises
        ------
        pint.errors.DimensionalityError
            If the units are not compatible.

        """
        q = self.value * ureg(self.unit)
        q_converted = q.to(new_unit)
        return UnitValue(value=q_converted.magnitude, unit=str(q_converted.units))
