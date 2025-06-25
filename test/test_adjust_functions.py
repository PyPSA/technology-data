"""Test functionality related to adjust_X functions in the Technologies class."""

import pathlib

import pandas as pd
import pytest

import technologydata as td


@pytest.mark.parametrize(
    "example_technologies",
    [
        {
            "technologies_name": "example01",
            "technologies_path": pathlib.Path(
                "technologydata", "datasources", "example01"
            ),
        },
        {
            "technologies_name": "example02",
            "technologies_path": pathlib.Path(
                "technologydata", "datasources", "example02"
            ),
        },
    ],
    indirect=True,
)  # type: ignore
def test_no_economies_of_scale(example_technologies: td.Technologies) -> None:
    """Without economies of scale the value should not change."""
    org_data = example_technologies.data.copy()

    # Get the adjusted data
    example_technologies.adjust_scale(
        new_scale=999,
        unit="MW",
        scaling_exponent=1,
    )
    adjusted_data = example_technologies.data["value"].copy()

    # Check if all values in adjusted_data are equal to the original values
    assert all(
        adjusted_value == original_value
        for adjusted_value, original_value in zip(adjusted_data, org_data["value"])
    ), "Scaling with exponent 1 should not change the value"


@pytest.mark.parametrize(
    "example_technologies",
    [
        {
            "technologies_name": "example01",
            "technologies_path": pathlib.Path(
                "technologydata", "datasources", "example01"
            ),
        },
        {
            "technologies_name": "example02",
            "technologies_path": pathlib.Path(
                "technologydata", "datasources", "example02"
            ),
        },
    ],
    indirect=True,
)  # type: ignore
def test_economies_of_scale(example_technologies: td.Technologies) -> None:
    """Test with common economies of scale with exponent 0.5."""
    org_data = example_technologies.data.copy()

    # Get the adjusted data
    example_technologies.adjust_scale(
        new_scale=2,
        unit="MW",
        scaling_exponent=0.5,
    )
    adjusted_data = example_technologies.data["value"].copy()

    # Calculate the expected values based on the scaling formula
    expected_values = org_data["value"] * (2 / org_data["scale"]) ** (0.5 - 1)

    # Check if all adjusted values are approximately equal to the expected values
    assert all(
        abs(adjusted_value - expected_value)
        < 1e-6  # Use a small tolerance for floating-point comparison
        for adjusted_value, expected_value in zip(adjusted_data, expected_values)
    ), "Scaling with exponent 0.5 should change the value to approx 0.7"


@pytest.mark.parametrize(
    "example_technologies",
    [
        {
            "technologies_name": "forecast01",
            "technologies_path": pathlib.Path(
                "test", "test_adjust_functions", "forecast01"
            ),
        }
    ],
    indirect=True,
)  # type: ignore
def test_adjust_year_linear_interpolation(
    example_technologies: td.Technologies,
) -> None:
    """Test linear forecasting through a middle value in 2025."""
    forecast = example_technologies.adjust_year(year=2025, model={"method": "linear"})

    assert forecast.shape[0] == 3, "Forecasted data missing additional row"
    assert forecast.iloc[1]["value"] == 150, (
        "Linear forecasted value for 2025 is incorrect"
    )


@pytest.mark.parametrize(
    "example_technologies",
    [
        {
            "technologies_name": "forecast01",
            "technologies_path": pathlib.Path(
                "test", "test_adjust_functions", "forecast01"
            ),
        }
    ],
    indirect=True,
)  # type: ignore
def test_adjust_year_linear_extrapolation(
    example_technologies: td.Technologies,
) -> None:
    """Test linear forecasting for a value outside the range of entries provided."""
    forecast = example_technologies.adjust_year(
        year=2040,
        model={
            "method": "linear",
            "order": 1,
            "limit_area": None,
            "limit_direction": "forward",
        },
    )

    assert forecast.shape[0] == 3, "Forecasted data missing additional row"
    assert forecast.iloc[2]["value"] == 200, (
        "Linear forecasted value for 2040 is incorrect and should be equal to the last support year value, i.e. 2030"
    )


@pytest.mark.parametrize(
    "example_technologies, to_currency, source, reference_technologies",
    [
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2020",
            "World Bank",
            {
                "output_USD_2020": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "currency_conversion01",
                    "output",
                    "WB_USD_2020",
                ),
            },
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "EUR_2020",
            "World Bank",
            {
                "output_EUR_2020": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "currency_conversion01",
                    "output",
                    "WB_EUR_2020",
                ),
            },
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "CNY_2020",
            "World Bank",
            {
                "output_CNY_2020": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "currency_conversion01",
                    "output",
                    "WB_CNY_2020",
                ),
            },
        ),
    ],
    indirect=["example_technologies"],
)  # type: ignore
def test_adjust_currency(
    example_technologies: td.Technologies,
    to_currency: str,
    source: str,
    reference_technologies: dict[str, pathlib.Path],
) -> None:
    """Test currency conversion and inflation adjustments."""
    # Load the reference technologies data
    references_techs = td.Technologies(reference_technologies)

    # Adjust the currency using the specified method
    example_technologies.adjust_currency(to_currency=to_currency, source=source)

    # Fill known NaN values in comments
    references_techs.data = references_techs.data.fillna({"comments": ""})
    example_technologies.data = example_technologies.data.fillna({"comments": ""})

    pd.testing.assert_frame_equal(references_techs.data, example_technologies.data)


@pytest.mark.parametrize(
    "example_technologies, to_currency, source, expected_result",
    [
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2020",
            "International Monetary Fund",
            "no failure",
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2023",
            "International Monetary Fund",
            "no failure",
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2031",
            "International Monetary Fund",
            ValueError,
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2020",
            "World Bank",
            "no failure",
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2023",
            "World Bank",
            "no failure",
        ),
        (
            {
                "technologies_name": "input",
                "technologies_path": pathlib.Path(
                    "test", "test_adjust_functions", "currency_conversion01", "input"
                ),
            },
            "USD_2031",
            "World Bank",
            ValueError,
        ),
    ],
    indirect=["example_technologies"],
)  # type: ignore
def test_adjust_currency_fails(
    example_technologies: td.Technologies,
    to_currency: str,
    source: str,
    expected_result: str | ValueError,
) -> None:
    """Test if currency conversion and inflation adjustments work for different years."""
    if isinstance(expected_result, type) and expected_result is ValueError:
        with pytest.raises(ValueError, match="No data found for base year"):
            example_technologies.adjust_currency(to_currency=to_currency, source=source)
    else:
        example_technologies.adjust_currency(to_currency=to_currency, source=source)
