# Introduction

## Currencies Class Documentation

### Overview

The `Currencies` class provides utility methods for handling currency-related operations, including deflating values based on specified deflator functions, converting currency values in a DataFrame, and ensuring currency units are correctly formatted. 
This class utilizes the `pydeflate` package for inflation adjustment and currency conversion calculations and the `hdx` package for relating country and currency ISO3 codes.

### Installation

To use the `Currencies` class, ensure you have the following packages installed:

`pip install pandas pydeflate hdx-python`

### Usage

#### Getting Country from Currency Code

To retrieve a list of country ISO3 codes associated with a given currency code, use the `get_country_from_currency` method. The method handles special cases (like EUR or USD) and proxies for multi-country currencies according to defined selection criteria based on economic indicators. Namely
```bibtex
"EUR": "EUR",  # special code in pydeflate
"USD": "USD",  # special code in pydeflate
"GBP": "GBR",  # Based on GDP, pop
"NZD": "NZL",  # Based on GDP, pop
"NOK": "NOR",  # Based on GDP, pop
"AUD": "AUS",  # Based on GDP, pop
"ILS": "ISR",  # Based on GDP, pop
"CHF": "CHE",  # Based on GDP, pop
"MAD": "MAR",  # Based on GDP, pop
"ANG": "CUW",  # Based on GDP, pop; Issue in hdx, since 04/2025 ANG and CUW use XCG (XCG is not supported by `pydeflate` and hdx yet)
"XPF": "PYF",  # Missing data for "NCL", "WLF" for a range of years in World Bank data / "pydeflate"
"XOF": "NER",  # Chosen as proxy, because the difference between NER inflation and unweighted-average inflation rate of all XOF member countries is the lowest for 2015-2023 (XOF average: 1.180816, NER: 1.173449)
"XCD": "GRD",  # Chosen as proxy, because the difference between GRD inflation and unweighted-average inflation rate of all XCD member countries is the lowest for 2015-2023 (XCD average: 1.147741, GRD: 1.156121)
"XAF": "CAF",  # Chosen as proxy, because the difference between CAF inflation and unweighted-average inflation rate of all XAF member countries is the lowest for 2015-2023 (XAF average: 1.300318, CAF: 1.277488)
```

#### Extracting Currency Unit

To check if a string contains a currency unit and extract it, use the `extract_currency_unit` method.

#### Updating Currency Unit

To replace the currency code and/or year in a string representing a currency unit, use the `update_currency_unit` method.

#### Adjusting Currency Values

To convert and adjust currency values in a DataFrame to a target currency and base year, use the `adjust_currency` method. The method uses deflator and exchange rates from a selected number of resource, namely the International Monetary Fund World Economic Outlook, GDP deflators and exchange rates from the World Bank and the World Bank’s linked GDP deflator and exchange rates. The currency adjustment is performed in-place and updates the input DataFrame.

### Limitations

The limitations of the class are:
- the `get_country_from_currency` method assumes that the currency code provided is valid and exists in the HDX country data. If the currency code is not found, it may return an empty list or raise a `KeyError`. 
- the `extract_currency_unit` method relies on a specific format for currency units (i.e., `<ISO3 CURRENCY_CODE>_<YEAR>`). If the input string does not match this format, it will return `None`. I.e. currency values need to provided with the correctly formatted currency unit, e.g. `USD_2020`, `EUR_2023` or `CNY_2025`. The currency symbols follow [ISO 4217](https://de.wikipedia.org/wiki/ISO_4217).
- the `adjust_currency` method requires the input DataFrame to contain specific columns: `unit`, `value`, and `region`. If any of these columns are missing, a `KeyError` will be raised.
- Currency conversion and adjustment methods only follow World Bank and International Monetary Fund publications. Furthermore, due to data availability, they provide results only for YYYY-2 (WB) or YYYY+2 (IMF).

### Assumptions

The assumptions are:
- the class assumes that the currency codes provided are valid ISO3 currency codes.
- the `get_country_from_currency `method assumes that the HDX package's country data is up-to-date and accurately reflects the current currency usage by countries.
- the `adjust_currency` method assumes that the `pydeflate` package is correctly configured and that the necessary deflation functions are registered in the `deflation_function_registry`.

### Conclusion

The `Currencies` class is a powerful utility for managing currency operations in Python, leveraging the capabilities of the pydeflate and hdx packages. By following the usage examples and being aware of the limitations and assumptions, users can effectively integrate currency handling into their data processing workflows.

## Technologies Class Documentation

The `Technologies` class provides a convenient interface for working with collections of technology data, including loading datasets and performing various adjustments such as currency conversion.

### Creating a Technologies object

You can create a `Technologies` object by providing a list of sources. Sources can be specified as strings (referring to available datasets), `Source` objects, or dictionaries mapping source names to paths.

```bibtex
from technologydata import Technologies

# Example: Load all available sources
techs = Technologies(sources=["source1", "source2"])

# Or, if you want to load a single source
techs = Technologies(sources="source1")

# You can also use a dictionary mapping
# techs = Technologies(sources={"source1": Path("/path/to/source1")})
```

### Adjusting the currencies

The `adjust_currency` method allows you to convert all currency values in the dataset to a specified target currency and year. It uses the `adjust_currency` method from the `Currencies` class. The currency adjustment is performed in-place and updates the `data` attribute of the `Technologies` object.