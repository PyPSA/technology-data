"""Classes for Currencies methods."""

import logging
import pathlib
import re
from collections.abc import Callable
from typing import Any

import pandas as pd
import pydeflate as pyd
from hdx.location.country import Country

logger = logging.getLogger(__name__)


deflation_function_registry: dict[Any, Any] = {}


def _register_deflator(name: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """
    Register a deflation function with a given name.

    This decorator function allows you to register a deflation function
    under a specified name in the deflation function registry.

    Parameters
    ----------
    name : str
        The name under which the deflation function will be registered.

    Returns
    -------
    Callable[[Callable], Callable]
        A decorator that registers the provided function as a deflation
        function in the registry.

    Examples
    --------
    >>> @_register_deflator('example_deflator')
    ... def example_function(data):
    ...     return data * 0.5

    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        if name.casefold() in deflation_function_registry.keys():
            raise KeyError(
                "Another deflator with the same name has already been registered."
            )
        deflation_function_registry[name.casefold()] = func
        return func

    return decorator


@_register_deflator("International Monetary Fund")
def _imf_gdp_deflate_wrapper(*args: Any, **kwargs: Any) -> pd.DataFrame:
    """
    Introduce wrapper function for pydeflate.imf_gdp_deflate.

    Uses GDP deflators and exchange rates from the IMF World Economic Outlook.

    This function acts as a registered deflator under the name "imf_gdp_deflate" and delegates
    all arguments to the underlying `imf_gdp_deflate` function from the pydeflate package.

    Parameters
    ----------
    *args : tuple
        Positional arguments to be passed to `pydeflate.imf_gdp_deflate`.
    **kwargs : dict
        Keyword arguments to be passed to `pydeflate.imf_gdp_deflate`.

    Returns
    -------
    pandas.DataFrame
        The DataFrame returned by `pydeflate.imf_gdp_deflate`, containing the
        deflated values according to the specified parameters.

    Notes
    -----
    This wrapper function is primarily used for registration purposes and does
    not modify the behavior or signature of the underlying `imf_gdp_deflate` function.

    """
    return pyd.imf_gdp_deflate(*args, **kwargs)


@_register_deflator("World Bank")
def _wb_gdp_deflate_wrapper(*args: Any, **kwargs: Any) -> pd.DataFrame:
    """
    Introduce wrapper function for pydeflate.wb_gdp_deflate.

    Uses GDP deflators and exchange rates from the World Bank.

    This function acts as a registered deflator under the name "wb_gdp_deflate"
    and delegates all arguments to the underlying `wb_gdp_deflate` function
    from the pydeflate package.

    Parameters
    ----------
    *args : tuple
        Positional arguments to be passed to `pydeflate.wb_gdp_deflate`.
    **kwargs : dict
        Keyword arguments to be passed to `pydeflate.wb_gdp_deflate`.

    Returns
    -------
    pandas.DataFrame
        The DataFrame returned by `pydeflate.wb_gdp_deflate`, containing the
        deflated values according to the specified parameters.

    Notes
    -----
    This wrapper function is primarily used for registration purposes and does
    not modify the behavior or signature of the underlying `wb_gdp_deflate` function.

    """
    return pyd.wb_gdp_deflate(*args, **kwargs)


@_register_deflator("World Bank (Linked)")
def _wb_gdp_linked_deflate_wrapper(*args: Any, **kwargs: Any) -> pd.DataFrame:
    """
    Introduce wrapper function for pydeflate.wb_gdp_linked_deflate.

    Uses the World Bank’s linked GDP deflator and exchange rates data.

    This function acts as a registered deflator under the name "wb_gdp_linked_deflate"
    and delegates all arguments to the underlying `wb_gdp_linked_deflate` function
    from the pydeflate package.

    Parameters
    ----------
    *args : tuple
        Positional arguments to be passed to `pydeflate.wb_gdp_linked_deflate`.
    **kwargs : dict
        Keyword arguments to be passed to `pydeflate.wb_gdp_linked_deflate`.

    Returns
    -------
    pandas.DataFrame
        The DataFrame returned by `pydeflate.wb_gdp_linked_deflate`, containing the
        deflated values according to the specified parameters.

    Notes
    -----
    This wrapper function is primarily used for registration purposes and does
    not modify the behavior or signature of the underlying `wb_gdp_linked_deflate` function.

    """
    return pyd.wb_gdp_linked_deflate(*args, **kwargs)


class Currencies:
    """
    A utility class for handling currency-related operations.

    This class provides static methods for ensuring currency units, deflating values based on a deflator function, and
    converting and adjusting currency values in a DataFrame.

    """

    PYDEFLATE_BASE_PATH = pathlib.Path(
        pathlib.Path(__file__).resolve().parent, "pydeflate"
    )
    CURRENCY_UNIT_DEFAULT_FORMAT = r"([A-Z]{3})_(\d{4})"
    _currency_countries_cache: dict[str, list[str]] = {}

    @staticmethod
    def get_country_from_currency(currency_code: str) -> str:
        """
        Retrieve a list of country ISO3 codes associated with a given currency code.

        This method accesses a predefined dictionary of countries from the hdx package and their corresponding
        currencies, and returns a list of countries that use the specified currency.

        Additionally, the method accommodates special cases, such as EUR/USD, and proxies for multi-country currencies based on established selection criteria derived from economic indicators. It is designed to be used alongside the 'pydeflate' functions, which accept currency codes instead of country codes for currencies of particular importance. Consequently, for the Euro and US Dollar, this method will return USD and EUR, respectively.

        Currency to Country Mapping Rules:
            - EUR/USD: Return currency code itself (special pydeflate handling)
            - GBP/NZD/etc: Return largest economy in currency zone
            - Regional currencies (XAF/XOF/XCD/XPF): Proxy countries based on inflation alignment
            - Others: First country found using the currency in HDX data

        Parameters
        ----------
        currency_code: str
            The ISO3 currency code for which to find associated countries.

        Returns
        -------
        str
            ISO3 country code that use the specified currency.

        Raises
        ------
        ValueError
            If currency code is not found in either special cases or HDX data

        Examples
        --------
        >>> Currencies.get_country_from_currency("AFN")
        "AFG"
        >>> Currencies.get_country_from_currency("USD")
        "USD"

        """
        # Response cache - builds only on first call
        if len(Currencies._currency_countries_cache) == 0:
            hdx_currency_dict = Country.countriesdata()["currencies"]
            for country, currency in hdx_currency_dict.items():
                if currency is not None:
                    if currency not in Currencies._currency_countries_cache.keys():
                        Currencies._currency_countries_cache[currency] = []
                    Currencies._currency_countries_cache[currency].append(country)

        # Handle special cases for specific currency codes
        special_cases = {
            # Multi-country currencies (return codes directly)
            "EUR": "EUR",  # Eurozone
            "USD": "USD",  # Dollarized economies
            # Single primary countries
            "GBP": "GBR",  # United Kingdom (largest GBP economy)
            "NZD": "NZL",  # New Zealand
            "NOK": "NOR",  # Norway
            "AUD": "AUS",  # Australia
            "ILS": "ISR",  # Israel
            "CHF": "CHE",  # Switzerland
            "MAD": "MAR",  # Morocco
            # Special regional cases with proxy selection criteria
            "ANG": "CUW",  # Curaçao (GDP-weighted proxy)
            "XPF": "PYF",  # French Polynesia (data availability)
            "XOF": "NER",  # Niger (lowest inflation differential 2015-2023)
            "XCD": "GRD",  # Grenada (lowest inflation differential)
            "XAF": "CAF",  # Central African Republic (lowest inflation differential)
        }

        # Return the special case if it exists
        if currency_code in special_cases:
            return special_cases[currency_code]
        elif currency_code not in Currencies._currency_countries_cache:
            raise KeyError(f"Unsupported currency code {currency_code}")
        else:
            if len(Currencies._currency_countries_cache[currency_code]) == 1:
                return Currencies._currency_countries_cache[currency_code][0]
            else:
                raise KeyError(
                    f"The currency {currency_code} corresponds to more than one country. Namely: {Currencies._currency_countries_cache[currency_code]}"
                )

    @staticmethod
    def extract_currency_unit(
        input_string: str, expected_format: str = CURRENCY_UNIT_DEFAULT_FORMAT
    ) -> str | None:
        r"""
        Check if the input string contains a currency unit.

        The method searches for a substring that matches the expected format and extracts it if found.

        Parameters
        ----------
        input_string : str
            The string to check.
        expected_format : str
            The string with the expected format.

        Returns
        -------
        str
            The matched currency unit if found, None otherwise.

        Raises
        ------
        ValueError
            If the input_string or the expected_format are not a string.

        Examples
        --------
        >>> Currencies.extract_currency_unit("The price is USD_2025/kW_el", r"[A-Z]{3}_\d{4}")
        'USD_2025'
        >>> Currencies.extract_currency_unit("No currency here", r"[A-Z]{3}_\d{4}")
        None
        >>> Currencies.extract_currency_unit("US-2025", r"[A-Z]{2}_\d{4}")
        None

        """
        if not isinstance(input_string, str) or not isinstance(expected_format, str):
            raise ValueError("Input must be a string.")

        currency_unit_pattern = re.compile(expected_format)
        match = currency_unit_pattern.search(input_string)
        return match.group(0) if match else None

    @staticmethod
    def update_currency_unit(
        input_string: str,
        new_currency_code: str,
        new_currency_year: str,
        expected_format: str = CURRENCY_UNIT_DEFAULT_FORMAT,
        separator: str = "_",
    ) -> str | None:
        """
        Replace the currency code and/or currency unit in the input string.

        Parameters
        ----------
        input_string : str
            The string containing the currency unit to be replaced.
        new_currency_code : str
            The new currency code to replace the existing one.
        new_currency_year : str
            The new currency year to replace the existing one.
        expected_format: str
            The string expected format.
        separator: str
            The currency unit separator.  #TODO make it automatically detectable

        Returns
        -------
        str
            The modified string with the currency code and/or currency_year replaced, None otherwise.

        Raises
        ------
        ValueError
            If the input_string or the expected_format or the new_currency_code or the new_currency_year are not a string.

        Examples
        --------
        >>> Currencies.update_currency_unit("The price is EUR_2025", "USD")
        'The price is USD_2025'
        >>> Currencies.update_currency_unit("The price is EUR_2025", "2024")
        'The price is EUR_2024'
        >>> Currencies.update_currency_unit("No currency here", "USD")
        None

        """
        if not isinstance(input_string, str) or not isinstance(expected_format, str):
            raise ValueError("Input must be a string.")
        if new_currency_code is not None and not isinstance(new_currency_code, str):
            raise ValueError("new_currency_code must be a string.")
        if new_currency_year is not None and not isinstance(new_currency_year, str):
            raise ValueError("new_currency_year must be a string.")

        currency_unit = Currencies.extract_currency_unit(input_string, expected_format)
        if currency_unit:
            currency_code = (
                new_currency_code
                if new_currency_code is not None
                else currency_unit.split(separator)[0]
            )
            currency_year = (
                new_currency_year
                if new_currency_year is not None
                else currency_unit.split(separator)[1]
            )
            new_currency_unit = f"{currency_code}{separator}{currency_year}"
            return input_string.replace(currency_unit, new_currency_unit)
        else:
            return None

    @staticmethod
    def get_deflate_row_function(
        base_year: int,
        deflator_name: str,
        year: str,
        iso_code: str,
        target_value_column: str,
        source_country: str,
        target_country: str,
    ) -> Callable[[pd.Series], pd.Series]:
        """
        Retrieve a function to deflate a row of data based on specified parameters.

        This method retrieves a deflation function from the registry using the provided deflator name
        and returns a closure that can be used to deflate a single row of data. The returned function
        takes a row as input and applies the deflation logic to it, returning the deflated row.

        Parameters
        ----------
        base_year : int
            The base year to which the values should be deflated.
        deflator_name : str
            The name of the deflation function to retrieve from the registry.
        year : str
            The name of the column that contains the year information for deflation.
        iso_code : str
            The name of the column that contains the ISO code for the country.
        target_value_column : str
            The name of the column that contains the target value to be adjusted.
        source_country : str
            The name of the column that contains the ISO3 country representing the source currency to convert from.
        target_country : str
            The ISO3 country representing the target currency to convert to.

        Returns
        -------
        Callable[[pd.Series], pd.Series]
            A function that takes a row (as a pandas Series) and returns the deflated row (as a pandas Series).

        Raises
        ------
        ValueError
            If the specified deflator function name is not found in the deflation function registry.

        Examples
        --------
        >>> deflator_func = Currencies.get_deflate_row_function(
        ...     base_year=2022,
        ...     deflator_name='cpi_deflator',
        ...     year='currency_year',
        ...     iso_code='region',
        ...     target_value_column='value',
        ...     source_country='currency_code',
        ...     target_country='USD',
        ... )

        """
        deflation_function = deflation_function_registry.get(deflator_name.casefold())
        if deflation_function is None:
            raise ValueError(
                f"Deflator function '{deflator_name}' not found in registry"
            )

        def deflate_row(row: pd.Series) -> pd.Series:
            row_df = pd.DataFrame([row])

            deflated_df = deflation_function(
                data=row_df,
                base_year=base_year,
                source_currency=row[source_country],
                target_currency=target_country,
                id_column=iso_code,
                year_column=year,
                value_column="value",
                target_value_column=target_value_column,
            )
            return deflated_df.iloc[0]

        return deflate_row

    @staticmethod
    def adjust_currency(
        base_year_val: int,
        target_currency: str,
        data: pd.DataFrame,
        deflator_function_name: str = "World Bank",
        pydeflate_path: pathlib.Path = PYDEFLATE_BASE_PATH,
        separator: str = "_",
    ) -> pd.DataFrame:
        """
        Convert and/or adjust currency values in a DataFrame to a target currency and base year using deflation.

        This function adjusts the currency values in the input DataFrame by deflating and/or converting
        them to a specified target currency and base year. It identifies rows with currency units
        matching the format `<CURRENCY_CODE>_<YEAR>`, applies a deflator function to adjust the values,
        and updates the currency codes in the 'unit' column accordingly.

        Parameters
        ----------
        base_year_val : int
            The base year to which the currency values should be adjusted.
        target_currency : str
            The ISO3 currency code representing the target currency to convert to.
        data : pandas.DataFrame
            The input DataFrame containing at least the columns 'unit', 'value', and 'region'.
            The 'unit' column must have currency codes in the format `<3-letter currency code>_<year>`.
        deflator_function_name : str
            The name of the deflation function to use from the deflation function registry. Default is "world_bank".
        pydeflate_path : pathlib.Path
            The file system path where deflator and exchange rate data will be saved or loaded from.
        separator : str
            Currency unit separator #TODO: try to make it automatically detectable

        Returns
        -------
        pandas.DataFrame
            A DataFrame with currency values deflated to the specified base year and converted to the
            target currency. The 'unit' column will have updated currency codes reflecting the target currency.

        Raises
        ------
        ValueError
            If `pydeflate_path` is None.
            If no rows in the DataFrame contain a valid currency unit matching the expected pattern.
            If the specified deflator function name is not found in the deflation function registry.
        KeyError
            If any of the required columns ('unit', 'value', 'region') are missing from the input DataFrame.

        Examples
        --------
        >>> data = pd.DataFrame({
        ...     'unit': ['WB_USD_2020', 'EUR_2021', 'JPY_2019'],
        ...     'value': [100, 200, 300],
        ...     'region': ['USA', 'FRA', 'JPN']
        ... })
        >>> adjusted_data = Currencies.adjust_currency(
        ...     base_year_val=2022,
        ...     target_currency='USD',
        ...     data=data,
        ...     deflator_function_name='International Monetary Fund',
        ... )
        >>> adjusted_data['unit']
        0    USD_2022
        1    USD_2022
        2    USD_2022
        Name: unit, dtype: object

        """
        # Specify the path where deflator and exchange data will be saved
        if pydeflate_path is not None:
            pyd.set_pydeflate_path(pydeflate_path)
            pydeflate_path.mkdir(parents=True, exist_ok=True)
        else:
            raise ValueError(
                "The path where the deflator and exchange data will be saved is None"
            )

        # Validate columns presence
        required_columns = {"unit", "value", "region"}
        if not required_columns.issubset(data.columns):
            missing = required_columns - set(data.columns)
            raise KeyError(f"Input dataFrame is missing required columns: {missing}")

        # Create a copy of the original data to avoid modifying input directly
        results = data.copy()

        # Select rows that correspond to a unit column that fulfills the format <3-letter currency code>-<currency year>
        has_currency_mask = (
            results["unit"].apply(Currencies.extract_currency_unit).notna()
        )
        currency_rows = results.loc[has_currency_mask].copy()

        if currency_rows.empty:
            logger.warning("No rows contain a valid currency unit.")
        else:
            # For each row, extract currency year and currency from the unit column
            currency_rows[["currency_code", "currency_year"]] = (
                currency_rows["unit"]
                .apply(Currencies.extract_currency_unit)
                .str.split(separator, expand=True)
            )

            # Cast 'currency_year' column to integer type
            currency_rows["currency_year"] = pd.to_numeric(
                currency_rows["currency_year"]
            ).astype(int)

        target_country = Currencies.get_country_from_currency(target_currency)
        currency_rows["source_country"] = currency_rows["currency_code"].apply(
            lambda c: Currencies.get_country_from_currency(c)
            if Currencies.get_country_from_currency(c)
            else None
        )

        deflate_row_func = Currencies.get_deflate_row_function(
            base_year=base_year_val,
            deflator_name=deflator_function_name,
            year="currency_year",
            iso_code="region",
            target_value_column="value",
            source_country="source_country",
            target_country=target_country,
        )

        adjusted_rows = currency_rows.apply(deflate_row_func, axis=1)

        adjusted_rows = adjusted_rows.drop(
            columns=["currency_code", "currency_year", "source_country"]
        )

        # Explicitly cast the data types of adjusted_rows to match results
        adjusted_rows["value"] = pd.to_numeric(
            adjusted_rows["value"], errors="coerce"
        ).astype(float)

        # Replace updated rows back into results DataFrame
        results.loc[has_currency_mask, adjusted_rows.columns] = adjusted_rows

        # Update the 'unit' column with the new currency code
        results.loc[has_currency_mask, "unit"] = currency_rows["unit"].apply(
            lambda unit: Currencies.update_currency_unit(
                unit,
                target_currency,
                str(base_year_val),
            )
        )

        return results
