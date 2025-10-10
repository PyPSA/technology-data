# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""
Parameter class for encapsulating a value, its unit, provenance, notes, and sources.

Examples
--------
>>> from technologydata.source import Source
>>> uv = pint.Quantity(1000, "EUR_2020/kW")
>>> src = Source(name="Example Source", authors="some authors", url="http://example.com")
>>> param = Parameter(quantity=uv, provenance="literature", note="Estimated", sources=[src])

"""

import logging
from typing import Annotated, Self

import pint
from pydantic import BaseModel, Field, PrivateAttr

import technologydata
from technologydata.source_collection import SourceCollection

logger = logging.getLogger(__name__)


class Parameter(BaseModel):  # type: ignore
    """
    Encapsulate a value with its unit, provenance, notes, sources, and more optional attributes required to describe technology parameters, like carrier, and heating value.

    Attributes
    ----------
    magnitude : int | float
        The numerical value of the parameter.
    units : Optional[str]
        The unit of the parameter.
    carrier : Optional[str]
        The energy carrier.
    heating_value : Optional[str]
        The heating value type.
    provenance : Optional[str]
        Description of the data's provenance.
    note : Optional[str]
        Additional notes about the parameter.
    sources : Optional[SourceCollection]
        List of sources for the parameter.

    """

    magnitude: Annotated[
        int | float, Field(description="The numerical value of the parameter.")
    ]
    units: Annotated[str | None, Field(description="The unit of the parameter.")] = None
    carrier: Annotated[
        str | None,
        Field(description="Carriers of the units, e.g. 'H2', 'el', 'H2O'."),
    ] = None
    heating_value: Annotated[
        str | None,
        Field(description="Heating value type for energy carriers ('LHV' or 'HHV')."),
    ] = None
    provenance: Annotated[str | None, Field(description="The data's provenance.")] = (
        None
    )
    note: Annotated[str | None, Field(description="Additional notes.")] = None
    sources: Annotated[
        SourceCollection,
        Field(description="List of sources for this parameter."),
    ] = SourceCollection(sources=[])

    # Private attributes for derived pint objects
    _pint_quantity: pint.Quantity = PrivateAttr(None)
    _pint_carrier: pint.Unit = PrivateAttr(None)
    _pint_heating_value: pint.Unit = PrivateAttr(None)

    def __init__(self, **data: float | str | SourceCollection | None) -> None:
        """Initialize Parameter and update pint attributes."""
        # pint uses canonical names for units, carriers, and heating values
        # Ensure the Parameter object is always created with these consistent names from pint
        if "units" in data and data["units"] is not None:
            technologydata.ureg.ensure_currency_is_unit(str(data["units"]))
            data["units"] = str(technologydata.ureg.Unit(data["units"]))
        if "carrier" in data and data["carrier"] is not None:
            data["carrier"] = str(technologydata.creg.Unit(data["carrier"]))
        if "heating_value" in data and data["heating_value"] is not None:
            data["heating_value"] = str(
                technologydata.hvreg.Unit(data["heating_value"])
            )

        super().__init__(**data)
        self._update_pint_attributes()

    def _update_pint_attributes(self) -> None:
        """
        Update internal pint attributes based on current object fields.

        This method initializes or updates the following attributes:
        - `_pint_quantity`: a pint Quantity created from `magnitude` and `units`.
        - `_pint_carrier`: a pint Unit created from `carrier`.
        - `_pint_heating_value`: a pint Unit created from `heating_value`, if applicable.

        Notes
        -----
        - Ensures that `units` are valid, especially for currency units.
        - Raises a ValueError if `heating_value` is set without a valid `carrier`.

        """
        # Create a pint quantity from magnitude and units
        if self.units:
            # `units` may contain an undefined currency unit - ensure the ureg can handle it
            technologydata.ureg.ensure_currency_is_unit(self.units)

            self._pint_quantity = technologydata.ureg.Quantity(
                self.magnitude, self.units
            )
        else:
            self._pint_quantity = technologydata.ureg.Quantity(self.magnitude)
        # Create the carrier as pint unit
        if self.carrier:
            self._pint_carrier = technologydata.creg.Unit(self.carrier)
        else:
            self._pint_carrier = None

        # Create the heating value as pint unit
        if self.heating_value and self.carrier:
            self._pint_heating_value = technologydata.hvreg.Unit(self.heating_value)
        elif self.heating_value and not self.carrier:
            raise ValueError(
                "Heating value cannot be set without a carrier. Please provide a valid carrier."
            )
        else:
            self._pint_heating_value = None

    def to(self, units: str) -> Self:
        """Convert the parameter's quantity to new units."""
        self._update_pint_attributes()

        # Do not allow for currency conversion here, as it requires additional information
        if technologydata.extract_currency_units(
            self._pint_quantity.units
        ) != technologydata.extract_currency_units(units):
            raise NotImplementedError(
                "Currency conversion is not supported in the `to` method. "
                "Use `change_currency` for currency conversions."
            )

        self._pint_quantity = self._pint_quantity.to(units)
        return Parameter(
            magnitude=self._pint_quantity.magnitude,
            units=str(self._pint_quantity.units),
            carrier=self.carrier,
            heating_value=self.heating_value,
            provenance=self.provenance,
            note=self.note,
            sources=self.sources,
        )

    def change_currency(
        self, to_currency: str, country: str, source: str = "worldbank"
    ) -> Self:
        """
        Change the currency of the parameter.

        This allows for conversion to a different currency as well as for inflation adjustments.
        To properly adjust for inflation, the function requires the `country` for which the inflation
        adjustment should be applied for.

        Note that this will harmonise all currencies used in the parameter's units,
        i.e. if the parameter `units` contains multiple different currencies,
        all of them will be converted to the target currency.

        Parameters
        ----------
        to_currency : str
            The target currency unit to convert to, e.g. "USD_2020", "EUR_2024", "CNY_2022".
        country : str
            The country for which the inflation adjustment should be made for.
            Must be the official ISO 3166-1 alpha-3 country code, e.g. "USA", "DEU", "CHN".
        source : str, optional
            The source of the inflation data, either "worldbank"/"wb" or "international_monetary_fund"/"imf".
            Defaults to "worldbank".
            Depending on the source, different years to adjust for inflation may be available.

        Returns
        -------
        Parameter
            A new Parameter object with the converted currency.

        Examples
        --------
        >>> param.change_currency("USD_2024", "USA")
        >>> param.change_currency("EUR_2020", "DEU", source="imf")
        >>> param.change_currency("EUR_2023", "USA", source="worldbank")

        """
        self._update_pint_attributes()

        # Ensure the target currency is a valid unit
        technologydata.ureg.ensure_currency_is_unit(to_currency)

        # Current unit and currency/currencies
        from_units = self._pint_quantity.units
        from_currencies = technologydata.extract_currency_units(from_units)
        # Replace all currency units in the from_units with the target currency
        to_units = technologydata.CURRENCY_UNIT_PATTERN.sub(
            to_currency, str(from_units)
        )

        # Create a temporary context to which we add the conversion rates
        # We use a temporary context to avoid polluting the global unit registry
        # with potentially invalid or incomplete conversion rates that do not
        # match the `country` and `source` parameters.
        context = technologydata.ureg.Context()

        # Conversion rates are all relative to the reference currency
        ref_currency = technologydata.ureg.get_reference_currency()
        ref_currency_p = technologydata.CURRENCY_UNIT_PATTERN.match(ref_currency)
        if ref_currency_p:
            ref_iso3 = technologydata.get_iso3_from_currency_code(
                ref_currency_p.group("cu_iso3")
            )
            ref_year = ref_currency_p.group("year")
        else:
            raise ValueError(
                f"Reference currency '{ref_currency}' does not match expected pattern."
            )

        # Get conversion rates for all involved currencies
        currencies = set(from_currencies).union({to_currency})
        # Avoid recursion error in pint definition by re-adding the reference currency
        currencies = currencies - {ref_currency}

        for currency in currencies:
            from_currency_p = technologydata.CURRENCY_UNIT_PATTERN.match(currency)
            if from_currency_p:
                from_iso3 = technologydata.get_iso3_from_currency_code(
                    from_currency_p.group("cu_iso3")
                )
                from_year = from_currency_p.group("year")
            else:
                raise ValueError(
                    f"Currency '{currency}' does not match expected pattern."
                )

            conversion_rate = technologydata.get_conversion_rate(
                from_iso3=from_iso3,
                from_year=from_year,
                to_iso3=ref_iso3,
                to_year=int(ref_year),
                country=country,
                source=source,
            )

            context.redefine(f"{currency} = {conversion_rate} * {ref_currency}")

        # Actual conversion using pint
        quantity = self._pint_quantity.to(to_units, context)

        return Parameter(
            magnitude=quantity.magnitude,
            units=str(quantity.units),
            carrier=self.carrier,
            heating_value=self.heating_value,
            provenance=self.provenance,
            note=self.note,
            sources=self.sources,
        )

    def change_heating_value(self, to_heating_value: str) -> Self:
        """
        Change the heating value of the parameter.

        This converts the parameter's heating value to another heating value,
        e.g. from "LHV" to "HHV", by taking into account the parameter's carrier.

        Parameters
        ----------
        to_heating_value : str
            The target heating value to convert to, e.g. "LHV", "HHV".

        Returns
        -------
        Parameter
            A new Parameter object with the converted heating value.

        Raises
        ------
        ValueError
            If the current parameter does not have a carrier or a heating value set,

        Examples
        --------
        >>> Parameter(magnitude=1, units="kWh", carrier="H2", heating_value="LHV").change_heating_value("HHV")
        <Parameter magnitude=1.5, units='kWh', carrier='H2', heating_value='HHV'>

        """
        if not self.carrier:
            raise ValueError(
                "Cannot change heating value without a carrier. Please provide a valid carrier."
            )
        if not self.heating_value:
            raise ValueError(
                "Cannot change heating value without a current heating value. "
                "Please provide a valid heating value."
            )
        if to_heating_value == self.heating_value:
            # No change needed, return the same parameter
            return self

        self._update_pint_attributes()

        from technologydata.constants import EnergyDensityHHV, EnergyDensityLHV

        # Create a dictionary of heating value ratios based on energy densities
        # The units of heating values are harmonized to "hv_units".
        # hv_units is the units attribute of the first element of EnergyDensityLHV
        hv_ratios = dict()

        # Access the key of the first element of the EnergyDensityLHV dictionary
        first_pair_key = next(iter(EnergyDensityLHV))

        # Get the units attribute of the first element of the EnergyDensityLHV dictionary
        hv_units = str(EnergyDensityLHV[first_pair_key].units)

        lhvs = {
            str(technologydata.creg.get_dimensionality(k)): v.to(hv_units)
            for k, v in EnergyDensityLHV.items()
        }
        hhvs = {
            str(technologydata.creg.get_dimensionality(k)): v.to(hv_units)
            for k, v in EnergyDensityHHV.items()
        }

        for dimension in self._pint_carrier.dimensionality.keys():
            if dimension in lhvs and dimension in hhvs:
                hv_ratios[dimension] = (
                    hhvs[dimension].magnitude / lhvs[dimension].magnitude
                )
            else:
                logger.error(
                    f"No heating values found for '{dimension}' in EnergyDensityLHV or EnergyDensityHHV. "
                    f"Assuming a ratio of 1."
                )
                hv_ratios[dimension] = 1.0

        # When converting from HHV -> LHV, we need to multiply by the ratios
        # When converting from LHV -> HHV, we need to divide by the ratios
        # We modify the hv_ratios dictionary to match the conversion direction
        if technologydata.hvreg.Unit(to_heating_value).is_compatible_with("HHV"):
            hv_ratios = hv_ratios
        elif technologydata.hvreg.Unit(to_heating_value).is_compatible_with("LHV"):
            hv_ratios = {k: 1 / v for k, v in hv_ratios.items()}

        multiplier = 1
        for dim, exponent in self._pint_carrier.dimensionality.items():
            if dim not in hv_ratios:
                raise NotImplementedError(
                    f"Heating value conversion not implemented for carrier dimension '{dim}'."
                )
            # Adjust the hv_ratios for the exponent of the carrier
            multiplier *= hv_ratios[dim] ** exponent

        return Parameter(
            magnitude=self.magnitude * multiplier,
            units=self.units,
            carrier=self.carrier,
            heating_value=to_heating_value,
            provenance=self.provenance,  # TODO implement for this function
            note=self.note,
            sources=self.sources,
        )

    def _check_parameter_compatibility(self, other: Self) -> None:
        """
        Check if two parameters are compatible in terms of units, carrier, and heating value.

        Parameters
        ----------
        other : Parameter
            The other Parameter instance to compare against.

        Raises
        ------
        ValueError
            If the carriers or heating values of the two parameters are not compatible.
            The error message specifies which attribute differs.

        """
        if self._pint_carrier != other._pint_carrier:
            raise ValueError(
                f"Operation not permitted on parameters with different carriers: "
                f"'{self._pint_carrier}' and '{other._pint_carrier}'."
            )
        if self._pint_heating_value != other._pint_heating_value:
            raise ValueError(
                f"Operation not permitted on parameters with different heating values: "
                f"'{self._pint_heating_value}' and '{other._pint_heating_value}'."
            )

    def __add__(self, other: Self) -> Self:
        """
        Add this Parameter to another Parameter.

        Parameters
        ----------
        other : Parameter
            The Parameter instance to add.

        Returns
        -------
        Parameter
            A new Parameter instance representing the sum of the two parameters.

        Notes
        -----
        This method checks for parameter compatibility before performing the addition.
        The resulting Parameter retains the carrier, heating value, and combines provenance,
        notes, and sources from both operands.

        """
        self._check_parameter_compatibility(other)
        new_quantity = self._pint_quantity + other._pint_quantity
        return Parameter(
            magnitude=new_quantity.magnitude,
            units=new_quantity.units,
            carrier=self.carrier,
            heating_value=self.heating_value,
            provenance=(self.provenance or "")
            + (other.provenance or ""),  # TODO make nicer
            note=(self.note or "") + (other.note or ""),  # TODO make nicer
            sources=SourceCollection(
                sources=(self.sources.sources + other.sources.sources)
            ),
        )

    def __sub__(self, other: Self) -> Self:
        """
        Subtract another Parameter from this Parameter.

        Parameters
        ----------
        other : Parameter
            The Parameter instance to subtract.

        Returns
        -------
        Parameter
            A new Parameter instance representing the result of the subtraction.

        Notes
        -----
        This method checks for parameter compatibility before performing the subtraction.
        The resulting Parameter retains the carrier, heating value, and combines provenance, notes, and sources.

        """
        self._check_parameter_compatibility(other)
        new_quantity = self._pint_quantity - other._pint_quantity
        return Parameter(
            magnitude=new_quantity.magnitude,
            units=str(new_quantity.units),
            carrier=self.carrier,
            heating_value=self.heating_value,
            provenance=(self.provenance or "")
            + (other.provenance or ""),  # TODO make nicer
            note=(self.note or "") + (other.note or ""),  # TODO make nicer
            sources=SourceCollection(
                sources=(self.sources.sources + other.sources.sources)
            ),
        )

    def __truediv__(self, other: int | float | Self) -> Self:
        """
        Divide this Parameter by another Parameter.

        Parameters
        ----------
        other : float | Parameter
            A scalar or a Parameter instance to divide by.

        Returns
        -------
        Parameter
            A new Parameter instance representing the division result.

        Raises
        ------
        ValueError
            If the heating values of the two parameters are different.

        Notes
        -----
        The method divides the quantities of the parameters and constructs a new Parameter.
        It also handles the division of carriers and heating values if present.

        """
        if isinstance(other, (int | float)):
            return Parameter(
                magnitude=self.magnitude / other,
                units=self.units,
                carrier=self.carrier,
                heating_value=self.heating_value,
                provenance=self.provenance,
                note=self.note,
                sources=self.sources,
            )

        # We don't check general compatibility here, as division is not a common operation for parameters.
        # Only ensure that the heating values are compatible.
        if self._pint_heating_value != other._pint_heating_value:
            raise ValueError(
                f"Cannot divide parameters with different heating values: "
                f"{self._pint_heating_value} and {other._pint_heating_value}."
            )

        new_quantity = self._pint_quantity / other._pint_quantity
        new_carrier = (
            self._pint_carrier / other._pint_carrier
            if self._pint_carrier and other._pint_carrier
            else None
        )
        new_heating_value = (
            self._pint_heating_value / other._pint_heating_value
            if self._pint_heating_value and other._pint_heating_value
            else None
        )

        return Parameter(
            magnitude=new_quantity.magnitude,
            units=str(new_quantity.units),
            carrier=new_carrier,
            heating_value=new_heating_value,
            provenance=(self.provenance or "")
            + (other.provenance or ""),  # TODO make nicer
            note=(self.note or "") + (other.note or ""),  # TODO make nicer
            sources=SourceCollection(
                sources=(self.sources.sources + other.sources.sources)
            ),
        )

    def __mul__(self, other: int | float | Self) -> Self:
        """
        Multiply two Parameter instances.

        Parameters
        ----------
        other : int | float | Parameter
            A scalar or a Parameter instance to multiply with.

        Returns
        -------
        Parameter
            A new Parameter instance representing the product of the two parameters.

        Raises
        ------
        ValueError
            If the heating values of the two parameters are not compatible (i.e., not equal).

        Notes
        -----
        - Multiplication is only performed if the heating values are compatible.
        - The method multiplies the underlying quantities and carriers (if present).
        - The heating value of the resulting parameter is the product of the input heating values.
        - Provenance, notes, and sources are combined from both parameters.
        - Compatibility checks beyond heating values are not performed.

        """
        if isinstance(other, int | float):
            return Parameter(
                magnitude=self.magnitude * other,
                units=self.units,
                carrier=self.carrier,
                heating_value=self.heating_value,
                provenance=self.provenance,
                note=self.note,
                sources=self.sources,
            )

        # We don't check general compatibility here, as multiplication is not a common operation for parameters.
        # Only ensure that the heating values are compatible.
        if self._pint_heating_value != other._pint_heating_value:
            raise ValueError(
                f"Cannot multiply parameters with different heating values: "
                f"{self._pint_heating_value} and {other._pint_heating_value}."
            )

        new_quantity = self._pint_quantity * other._pint_quantity
        new_carrier = (
            self._pint_carrier * other._pint_carrier
            if self._pint_carrier and other._pint_carrier
            else None
        )

        new_heating_value = self._pint_heating_value * other._pint_heating_value
        return Parameter(
            magnitude=new_quantity.magnitude,
            units=str(new_quantity.units),
            carrier=str(new_carrier),
            heating_value=str(new_heating_value),
            provenance=(self.provenance or "")
            + (other.provenance or ""),  # TODO make nicer
            note=(self.note or "") + (other.note or ""),  # TODO make nicer
            sources=SourceCollection(
                sources=(self.sources.sources + other.sources.sources)
            ),
        )

    def __eq__(self, other: object) -> bool:
        """
        Check for equality with another Parameter object.

        Compares all attributes of the current instance with those of the other object.

        Parameters
        ----------
        other : object
            The object to compare with. Expected to be an instance of Parameter.

        Returns
        -------
        bool
            True if all attributes are equal between self and other, False otherwise.
            Returns False if other is not a Parameter instance.

        """
        if not isinstance(other, Parameter):
            return NotImplemented

        self._update_pint_attributes()
        other._update_pint_attributes()

        for field in self.__class__.model_fields.keys():
            value_self = getattr(self, field)
            value_other = getattr(other, field)
            if value_self != value_other:
                return False
        return True

    def __pow__(self, exponent: float | int) -> Self:
        """
        Raise the parameter's value to a specified power.

        Parameters
        ----------
        exponent : float or int
            The exponent to raise the parameter's value to.

        Returns
        -------
        Parameter
            A new Parameter instance with the value raised to the specified power.

        Notes
        -----
        This method updates the internal pint attributes before applying the power operation.
        If the parameter has a carrier, it is also raised to the specified power.

        """
        self._update_pint_attributes()

        new_quantity = self._pint_quantity**exponent
        return Parameter(
            magnitude=new_quantity.magnitude,
            units=str(new_quantity.units),
            carrier=self._pint_carrier**exponent if self._pint_carrier else None,
            heating_value=self.heating_value,
            provenance=self.provenance,
            note=self.note,
            sources=self.sources,
        )
