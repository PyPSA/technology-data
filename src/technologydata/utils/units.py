# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT
"""Submodule containing pint.UnitRegistry subclasses and utility functions for handling units, conversions, and currency units."""

import json
import logging
import re
import typing
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import pint
import pydeflate
from frozendict import frozendict
from hdx.location.country import Country
from platformdirs import user_cache_dir

logger = logging.getLogger(__name__)

CURRENCY_UNIT_PATTERN = re.compile(r"\b(?P<cu_iso3>[A-Z]{3})_(?P<year>\d{4})\b")

# Set up cache directory and file for currency codes
CACHE_DIR = Path(user_cache_dir("technologydata"))  # TODO move to commons?
CACHE_DIR.mkdir(parents=True, exist_ok=True)  # TODO move to commons?
CURRENCY_CODES_CACHE = CACHE_DIR / "iso3_to_currency_codes.json"


SPECIAL_CASES_CURRENCY_CODE_TO_ISO3 = frozendict(
    {
        # Multi-country currencies (return codes directly)
        "EUR": "EUR",  # Eurozone
        # Single primary countries
        "AUD": "AUS",  # Australia
        "CHF": "CHE",  # Switzerland
        "DKK": "DNK",  # Denmark
        "GBP": "GBR",  # United Kingdom (largest GBP economy)
        "ILS": "ISR",  # Israel
        "MAD": "MAR",  # Morocco
        "NOK": "NOR",  # Norway
        "NZD": "NZL",  # New Zealand
        "USD": "USA",  # Dollar-ized economies
        # Special regional cases with proxy selection criteria
        "ANG": "CUW",  # Curaçao (GDP-weighted proxy)
        "XAF": "CAF",  # Central African Republic (lowest inflation differential)
        "XCD": "GRD",  # Grenada (lowest inflation differential)
        "XOF": "NER",  # Niger (lowest inflation differential 2015-2023)
        "XPF": "PYF",  # French Polynesia (data availability)
    }
)


def get_iso3_to_currency_codes(
    refresh: bool = False, ignore_cache: bool = False
) -> dict[str, str]:
    """
    Get all 3-letter currency codes from official UN data.

    Uses a persistent local cache to avoid unnecessary network requests.

    Parameters
    ----------
    refresh : bool, optional
        If True, the cache will be updated from the internet, by default False.
    ignore_cache : bool, optional
        If True, the cache will not be used and the data will always be fetched from the live feed.

    Returns
    -------
    dict[str, str]
        A dictionary mapping ISO3 country codes to their corresponding 3-letter currency codes.

    """
    currencies: dict[str, str] = {}

    if refresh:
        logger.debug("Deleting existing currency codes cache to refresh it.")
        CURRENCY_CODES_CACHE.unlink(missing_ok=True)

    if ignore_cache:
        logger.debug("Ignoring cache and fetching live currency codes.")
        currencies = Country.countriesdata()["currencies"]
    elif not CURRENCY_CODES_CACHE.exists():
        logger.debug(
            "Cache does not exist. Fetching live currency codes and creating cache."
        )
        currencies = Country.countriesdata()["currencies"]
        with open(CURRENCY_CODES_CACHE, "w") as f:
            json.dump(currencies, f)
    else:
        logger.debug("Reading currency codes from cache.")
        with open(CURRENCY_CODES_CACHE) as f:
            currencies = json.load(f)

    return currencies


def extract_currency_units(units: str | pint.Unit) -> list[str]:
    """
    Extract currency-like strings from a string or pint.Unit.

    Parameters
    ----------
    units : str or pint.Unit
        The units string or pint.Unit from which to extract currency-like strings.

    Returns
    -------
    list[str]
        A list of currency-like strings found in the input, formatted as "{3-letter currency code}_{year as YYYY}".
        If no matches are found, an empty list is returned.

    Examples
    --------
    >>> extract_currency_units("USD_2020/kW")
    ["USD_2020"]

    >>> extract_currency_units("EUR_2015/USD_2020")
    ["EUR_2015", "USD_2020"]

    """
    # Ensure that the input is a string
    units = str(units)

    # Get the 3-letter currency codes for all officially recognized currencies
    logger.debug("Retrieving all 3-letter currency codes from the `hdx-country`.")
    all_currency_codes = set(get_iso3_to_currency_codes().values())

    # Check if the units contain a currency-like string, defined as "{3-letter currency code}_{year as YYYY}"
    matches = CURRENCY_UNIT_PATTERN.findall(units)
    if len(matches) == 0:
        logger.debug("No currency-like string found in the units.")
        return []

    # Extract the currency codes from the matches using the regex groups
    logger.debug(f"Found currency-like strings in the units: {matches}")
    currency_codes = {code for code, year in matches}

    # Ensure that all currency codes are legitimate 3-letter currency codes
    invalid_codes = currency_codes - all_currency_codes
    if invalid_codes:
        invalid_currencies = [
            f"{code}_{year}" for code, year in matches if code in invalid_codes
        ]
        raise ValueError(
            f"The following unit(s) appear to be currency units, but have invalid 3-letter currency codes: {', '.join(invalid_currencies)}. "
        )

    # Reconstruct currency units from the matches
    matches = [f"{code}_{year}" for code, year in matches]

    return matches


@lru_cache
def get_conversion_rate(
    from_iso3: str,
    to_iso3: str,
    country: str,
    from_year: int,
    to_year: int,
    source: str = "worldbank",
) -> float:
    """
    Get the conversion rate from one currency (year, ISO3) to another currency (year, ISO3) from pydeflate.

    Parameters
    ----------
    from_iso3 : str
        The ISO3 code of the country of the currency to convert from (e.g., 'USA' if the source currency is USD).
    to_iso3 : str
        The ISO3 code of the country of the currency to convert to (e.g., 'DEU' if the target currency is EUR).
    country : str
        The ISO3 code of the country to adjust for inflation.
    from_year : int
        The julian year (YYYY) of the source currency.
    to_year : int
        The julian year (YYYY) of the target currency.
    source : str
        The source of the inflation data ('worldbank'/'wb' or 'international_monetary_fund'/'imf').

    """
    # Choose the deflation function based on the source
    deflation_function = {
        "worldbank": pydeflate.wb_gdp_deflate,
        "wb": pydeflate.wb_gdp_deflate,
        "international_monetary_fund": pydeflate.imf_gdp_deflate,
        "imf": pydeflate.imf_gdp_deflate,
    }[source]

    # Ensure that `country` is a valid ISO3 code; to_iso3 and from_iso3 should already have been parsed before the function was called
    if country not in get_iso3_to_currency_codes().keys():
        raise ValueError(f"Unknown ISO3 code for `country`: {country}.")

    # pydeflate only operates on pandas.DataFrame
    data = pd.DataFrame(
        {
            "iso3": [country],
            "from_year": [from_year],
            "value": [1],
        }
    )

    # Deflate values include currency conversion
    conversion_rates = deflation_function(
        data,
        source_currency=from_iso3,
        target_currency=to_iso3,
        id_column="iso3",
        year_column="from_year",
        base_year=to_year,
        value_column="value",
        target_value_column="new_value",
    )

    if conversion_rates.isna().any().any():
        raise ValueError(
            f"Conversion rate from {from_iso3} ({from_year}) to {to_iso3} ({to_year}) with inflation rate for {country} not found. "
        )

    return float(conversion_rates.loc[0, "new_value"])


@lru_cache
def get_iso3_from_currency_code(
    currency_code: str,
    special_cases: frozendict[str, str] = SPECIAL_CASES_CURRENCY_CODE_TO_ISO3,
) -> str:
    """
    Get the ISO3 country code from a 3-letter currency code using official UN data and opinionated assumptions.

    Parameters
    ----------
    currency_code : str
        The 3-letter currency code (e.g., 'USD', 'EUR').
    special_cases : dict[str, str], optional
        A dictionary mapping specific currency codes to their ISO3 country codes for special cases.
        Defaults to `technologydata.utils.units.SPECIAL_CASES_CURRENCY_CODE_TO_ISO3`.

    Returns
    -------
    str
        The ISO 3166 alpha 3 country code of the `currency_code`.

    Raises
    ------
    ValueError
        If the currency code is not found in the official list of currencies.

    """
    # Build reverse mapping: currency code -> list of ISO3 codes
    iso3_to_currency_codes = pd.DataFrame.from_dict(
        get_iso3_to_currency_codes(), orient="index", columns=["currency"]
    )
    iso3_to_currency_codes = iso3_to_currency_codes.reset_index(drop=False).rename(
        columns={"index": "iso3"}
    )
    currency_codes_to_iso3 = iso3_to_currency_codes.groupby("currency", as_index=False)[
        "iso3"
    ].agg(list)

    # Handle special cases for specific currency codes

    # Remove all currencies that are in the special cases from the mapping
    currency_codes_to_iso3 = currency_codes_to_iso3.loc[
        ~currency_codes_to_iso3["currency"].isin(special_cases.keys())
    ]

    # Mapping should now only contain unique currency codes to ISO3 codes mappings, safe to explode
    currency_codes_to_iso3 = currency_codes_to_iso3.explode("iso3")

    # Add the special cases to the mapping
    currency_codes_to_iso3 = pd.concat(
        [
            currency_codes_to_iso3,
            pd.DataFrame(
                list(special_cases.items()),
                columns=["currency", "iso3"],
            ),
        ],
        ignore_index=False,
    )

    # Special cases should handle all non-unique currency codes, check to make sure
    # and return the ones that are not handled in special cases

    if (
        duplicated_iso3 := currency_codes_to_iso3.explode("iso3")[
            "currency"
        ].duplicated()
    ).any():
        raise ValueError(
            "Some currency codes are used by multiple ISO3 codes but are not handled in `special_cases` "
            "and need to be added to the mapping: "
            f"{currency_codes_to_iso3[duplicated_iso3]}"
        )

    currency_codes_to_iso3 = currency_codes_to_iso3.set_index("currency")[
        "iso3"
    ].to_dict()

    try:
        return str(currency_codes_to_iso3[currency_code])
    except KeyError as e:
        raise ValueError(
            f"Currency code '{currency_code}' not found in the list of currencies. "
            "Please ensure it is a valid 3-letter currency code."
        ) from e


def patch_pint_registry_error_handling(registry: pint.registry.UnitRegistry) -> None:
    """
    Patch a Pint registry to use CustomUndefinedUnitError.

    Parameters
    ----------
    registry : pint.registry.UnitRegistry
        The Pint unit registry to patch.

    """
    # Store the original method
    original_get_name = registry.get_name

    def patched_get_name(
        self: pint.registry.UnitRegistry, name: str, *args: Any, **kwargs: Any
    ) -> Any:
        try:
            return original_get_name(name, *args, **kwargs)
        except pint.errors.UndefinedUnitError as e:
            # Raise the custom error with the same arguments
            raise CustomUndefinedUnitError(e.args[0]) from e

    # Replace the method
    registry.get_name = patched_get_name.__get__(registry)


class CustomUndefinedUnitError(pint.errors.UndefinedUnitError):  # type: ignore
    """
    Custom message for undefined unit errors.

    This custom error is raised when a unit is not defined in the unit registry.
    It provides more specific error messages, especially for currency units
    that are missing the currency year.

    Parameters
    ----------
    unit_names : list of str
        The names of the units that are not defined in the unit registry.

    Attributes
    ----------
    unit_names : list of str
        The names of the units that are not defined in the unit registry.

    Notes
    -----
    This error is a subclass of `pint.errors.UndefinedUnitError` and is designed
    to provide more specific error messages for currency units that are missing
    the currency year.

    """

    def __init__(self, unit_names: str | typing.Iterable[str]) -> None:
        """
        Initialize a CustomUndefinedUnitError instance.

        This constructor creates a new `CustomUndefinedUnitError` object,
        inheriting all default behaviors from `pint.errors.UndefinedUnitError`.

        Parameters
        ----------
        unit_names : str or iterable of str
            The name or names of the undefined units that caused the error.

        """
        # Use all defaults definitions from a standard pint.errors.UndefinedUnitError
        super().__init__(unit_names)

    def __str__(self) -> str:
        """
        Generate a custom error message string.

        This method generates a custom error message string based on the
        unit names that are not defined in the unit registry. It provides
        specific messages for currency units that are missing the currency year.

        Returns
        -------
        str
            The custom error message string.

        """
        # Retrieve unit names, defaulting to an empty list if not present
        unit_names = getattr(self, "unit_names", [])

        # Precompute valid currency codes for efficiency
        all_currency_codes = set(get_iso3_to_currency_codes().values())

        # Identify currency units without a year specification
        currency_errors = [
            code
            for unit in unit_names
            for code in re.findall(r"[A-Z]{3}", str(unit))
            if code in all_currency_codes
            and not CURRENCY_UNIT_PATTERN.search(str(unit))
        ]

        # Generate specific error message for currency units
        if currency_errors:
            missing_code = currency_errors[0]
            return f"Currency unit '{missing_code}' is missing the 4-digit currency year (e.g. {missing_code}_2020)."

        # Fallback to parent class error message if no specific currency error found
        return super().__str__()  # type: ignore


class SpecialUnitRegistry(pint.UnitRegistry):  # type: ignore
    """A special pint.UnitRegistry subclass that includes methods for handling currency units and conversion using pydeflate."""

    def __init__(self, *args: tuple[Any, ...], **kwargs: dict[str, Any]) -> None:
        """
        Initialize a SpecialUnitRegistry instance.

        This constructor creates a new `SpecialUnitRegistry` object,
        inheriting all default behaviors from `pint.UnitRegistry`.
        It also defines a reference currency unit (`USD_2020`) for
        handling currency conversions and related operations.

        Parameters
        ----------
        *args : tuple
            Positional arguments passed to the base `pint.UnitRegistry` constructor.
        **kwargs : dict
            Keyword arguments passed to the base `pint.UnitRegistry` constructor.

        """
        # Use all defaults definitions from a standard pint.UnitRegistry
        super().__init__(*args, **kwargs)

        # Define the reference currency unit
        # This is the base currency unit that all other currency units will be defined relative to
        # We use USD_2020 as the base currency unit because it is a currency that we can relate all other currencies to
        self.define("USD_2020 = [currency]")

    def get_reference_currency(self) -> str:
        """Get the reference currency from the unit registry."""
        reference_currency = [
            self._units[u].name
            for u in self._units
            if "[currency]" in self._units[u].reference
        ]
        if not reference_currency or len(reference_currency) != 1:
            raise ValueError(
                "The unit registry does not have a unique base currency defined as '[currency]'. Please define a base currency unit to proceed."
            )

        return str(reference_currency[0])

    def ensure_currency_is_unit(self, units: str) -> None:
        """
        Ensure that all currency units in the given string are valid units in the registry.

        Extracts all currency-like strings from the input and checks if they are defined
        in the unit registry. If they are not defined, they are added as valid units
        relative to the reference currency but without a conversion factor.

        Parameters
        ----------
        units : str
            The units string to check for currency units.

        Examples
        --------
        >>> ureg.ensure_currency_is_unit("USD_2020/kW")
        >>> ureg.ensure_currency_is_unit("EUR_2015/USD_2020")

        """
        logger.debug(f"Ensuring currency units of '{units}' are defined in `ureg`")
        currency_units = extract_currency_units(units)
        logger.debug(f"Found currency-like strings in the units: {currency_units}")

        if not currency_units:
            # Nothing to do
            return

        reference_currency = self.get_reference_currency()
        logger.debug(f"Reference currency is '{reference_currency}'.")

        # Check if the currency unit is already defined in the unit registry
        # if not, define it relative to the base currency USD_2015
        for currency_unit in currency_units:
            if currency_unit in self._units:
                logger.debug(
                    f"Currency unit '{currency_unit}' is already defined in the unit registry. Not redefining it."
                )
                continue

            logger.debug(
                f"Currency unit '{currency_unit}' not found in the unit registry. "
                f"Defining it without a conversion factor relative to the base currency '{reference_currency}'."
            )
            self.define(f"{currency_unit} = nan {reference_currency}")


# Unit registries used throughout the package for different purposes
ureg = SpecialUnitRegistry()  # For handling units, conversions, and currency units
creg = pint.UnitRegistry(
    filename=Path(__file__).parent / "carriers.txt"
)  # For tracking carriers and ensuring compatibility between them
hvreg = pint.UnitRegistry(
    filename=Path(__file__).parent / "heating_values.txt"
)  # For tracking heating values and ensuring compatibility between them

patch_pint_registry_error_handling(ureg)
