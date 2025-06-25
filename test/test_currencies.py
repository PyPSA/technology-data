"""Test the utility methods."""

import pandas as pd
import pytest

import technologydata as td


@pytest.mark.parametrize(
    "input_string, expected_format, expected_result",
    [
        ("EUR_2025", r"^[A-Z]{3}_\d{4}$", "EUR_2025"),
        ("EU_2025", r"^[A-Z]{2}_\d{4}$", "EU_2025"),
        ("The currency unit is EU_2025", r"[A-Z]{3}_\d{4}", None),
        ("The currency unit is EUR_202", r"[A-Z]{3}_\d{4}", None),
        ("The currency unit is EUR_2025/kW_el", r"[A-Z]{3}_\d{4}", "EUR_2025"),
        ("The currency unit is US_2025", r"[A-Z]{2}_\d{4}", "US_2025"),
        ("USD_2025", r"[A-Z]{3}_\d{4}", "USD_2025"),
        ("USD_2025", r"[A-Z]{3}-\d{4}", None),
        (123, r"[A-Z]{3}_\d{4}", ValueError),
        (r"[A-Z]{3}_\d{4}", 123, ValueError),
    ],
)  # type: ignore
def test_extract_currency_unit(
    input_string: str, expected_format: str, expected_result: str | None | ValueError
) -> None:
    """Check if a currency unit follows the wished format."""
    if isinstance(expected_result, type) and expected_result is ValueError:
        with pytest.raises(ValueError, match="Input must be a string."):
            td.Currencies.extract_currency_unit(input_string, expected_format)
    else:
        assert (
            td.Currencies.extract_currency_unit(input_string, expected_format)
            == expected_result
        )


@pytest.mark.parametrize(
    "input_string, new_currency_code, new_currency_year, expected_format, expected_result, expected_exception_message",
    [
        ("EUR_2025", "USD", None, r"^[A-Z]{3}_\d{4}$", "USD_2025", None),
        (
            "EUR_2025",
            "USD",
            "2023",
            td.Currencies.CURRENCY_UNIT_DEFAULT_FORMAT,
            "USD_2023",
            None,
        ),
        ("EUR_2025", None, "2021", r"^[A-Z]{3}_\d{4}$", "EUR_2021", None),
        (
            "The currency unit is EUR_2025",
            "GPD",
            "2021",
            r"[A-Z]{3}_\d{4}",
            "The currency unit is GPD_2021",
            None,
        ),
        ("The currency unit", "USD", "2019", r"[A-Z]{3}_\d{4}", None, None),
        (
            12345,
            "USD",
            "2019",
            r"[A-Z]{3}_\d{4}",
            ValueError,
            "Input must be a string.",
        ),
        (
            "The currency unit",
            123,
            "2019",
            r"[A-Z]{3}_\d{4}",
            ValueError,
            "new_currency_code must be a string.",
        ),
        (
            "The currency unit",
            "USD",
            123,
            r"[A-Z]{3}_\d{4}",
            ValueError,
            "new_currency_year must be a string.",
        ),
        (
            "The currency unit",
            "USD",
            "2019",
            123,
            ValueError,
            "Input must be a string.",
        ),
    ],
)  # type: ignore
def test_update_currency_unit(
    input_string: str,
    new_currency_code: str,
    new_currency_year: str,
    expected_format: str,
    expected_result: str | None | ValueError,
    expected_exception_message: str | None,
) -> None:
    """Check if a currency unit is correctly replaced."""
    if isinstance(expected_result, type) and expected_result is ValueError:
        with pytest.raises(ValueError, match=expected_exception_message):
            td.Currencies.update_currency_unit(
                input_string, new_currency_code, new_currency_year, expected_format
            )
    else:
        result = td.Currencies.update_currency_unit(
            input_string, new_currency_code, new_currency_year, expected_format
        )
        assert result == expected_result


@pytest.mark.parametrize(
    "base_year_val, deflator_function_name, input_dataframe, target_currency, expected_result, "
    "expected_exception_message",
    [
        (
            2021,
            "internaTiOnAl Monetary fUnD",
            pd.DataFrame(
                {
                    "region": ["FRA", "USA", "CAN", "ITA"],
                    "unit": ["EUR_2020/MWh_el", "USD_2020", "CAD_2020", "MWh"],
                    "value": [50.0, 100.0, 200.0, 300.0],
                }
            ),
            "USD",
            pd.DataFrame(
                {
                    "region": ["FRA", "USA", "CAN", "ITA"],
                    "unit": ["USD_2021/MWh_el", "USD_2021", "USD_2021", "MWh"],
                    "value": [59.93, 104.57, 171.92, 300],
                }
            ),
            None,
        ),
        (
            2020,
            "International Monetary Fund",
            pd.DataFrame(
                {
                    "unit": ["EUR_2015/MWh_el", "USD_2015", "CAD_2015", "MWh"],
                    "value": [50.0, 100.0, 200.0, 300.0],
                }
            ),
            "USD",
            KeyError,
            "Input dataFrame is missing required columns:",
        ),
        (
            2020,
            "random_deflate",
            pd.DataFrame(
                {
                    "region": ["FRA", "USA", "CAN", "ITA"],
                    "unit": ["EUR_2020/MWh_el", "USD_2020", "CAD_2020", "MWh"],
                    "value": [50.0, 100.0, 200.0, 300.0],
                }
            ),
            "USD",
            ValueError,
            "Deflator function 'random_deflate' not found in registry",
        ),
    ],
)  # type: ignore
def test_adjust_currency(
    base_year_val: int,
    deflator_function_name: str,
    input_dataframe: pd.DataFrame,
    target_currency: str,
    expected_result: pd.DataFrame | ValueError | KeyError,
    expected_exception_message: str | None,
) -> None:
    """Check if currency conversion and inflation adjustment work correctly."""
    if isinstance(expected_result, type) and expected_result is ValueError:
        with pytest.raises(ValueError, match=expected_exception_message):
            td.Currencies.adjust_currency(
                base_year_val,
                target_currency,
                input_dataframe,
                deflator_function_name,
            )
    elif isinstance(expected_result, type) and expected_result is KeyError:
        with pytest.raises(KeyError, match=expected_exception_message):
            td.Currencies.adjust_currency(
                base_year_val,
                target_currency,
                input_dataframe,
                deflator_function_name,
            )
    else:
        # Assume td.CurrencyUtils is imported in the test context
        new_dataframe = td.Currencies.adjust_currency(
            base_year_val,
            target_currency,
            input_dataframe,
            deflator_function_name,
        )

        new_dataframe["value"] = new_dataframe["value"].astype(float).round(2)

        pd.testing.assert_frame_equal(new_dataframe, expected_result)


@pytest.mark.parametrize(
    "input_currency, expected_country, expected_exception_message",
    [
        ("AFN", "AFG", None),
        ("CNY", "CHN", None),
        ("EUR", "EUR", None),
        ("USD", "USD", None),
        ("GBP", "GBR", None),
        ("NZD", "NZL", None),
        ("NOK", "NOR", None),
        ("AUD", "AUS", None),
        ("ILS", "ISR", None),
        ("CHF", "CHE", None),
        ("MAD", "MAR", None),
        ("ANG", "CUW", None),
        ("XPF", "PYF", None),
        ("XOF", "NER", None),
        ("XCD", "GRD", None),
        ("XAF", "CAF", None),
        ("YEA", KeyError, "Unsupported currency code"),
    ],
)  # type: ignore
def test_get_country_from_currency(
    input_currency: str,
    expected_country: str | KeyError,
    expected_exception_message: str | None,
) -> None:
    """Verify that the country(ies) ISO3 code(s) are correctly returned for a given currency ISO3 code."""
    if isinstance(expected_country, type) and expected_country is KeyError:
        with pytest.raises(KeyError, match=expected_exception_message):
            td.Currencies.get_country_from_currency(input_currency)
    else:
        result = td.Currencies.get_country_from_currency(input_currency)
        assert result == expected_country, (
            f"Expected {expected_country} but got {result} for currency {input_currency}"
        )
