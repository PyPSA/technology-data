"""Test functionality related to adjust_X functions in the Technologies class."""

import pathlib

import pandas as pd
import pytest

import technologydata as td


@pytest.mark.parametrize(
    "input_technologies, to_capacity, to_capacity_unit, scaling_exponent, scaled_parameters, absolute_parameters, expected_technologies",
    [
        (
            {
                "input": pathlib.Path(
                    "test", "test_adjust_functions", "adjust_scale", "input"
                ),
            },
            1,  # Original capacity equals the target capacity mean no change in scale.
            "MW",
            0.5,
            None,  # default
            None,  # default
            {
                "output": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "adjust_scale",
                    "output",
                    "no-change",
                ),
            },
        ),
        (
            {
                "input": pathlib.Path(
                    "test", "test_adjust_functions", "adjust_scale", "input"
                ),
            },
            2,  # This should double absolute values, but not change specific values.
            "MW",
            1.0,
            None,  # default
            None,  # default
            {
                "output": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "adjust_scale",
                    "output",
                    "absolute-2_scale-1",
                ),
            },
        ),
        (
            # This is the interesting case where the scaling exponent is not 1.0
            # and the capacity is not equal to the target capacity.
            {
                "input": pathlib.Path(
                    "test", "test_adjust_functions", "adjust_scale", "input"
                ),
            },
            10,
            "MW",
            0.7,
            None,  # default
            None,  # default
            {
                "output": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "adjust_scale",
                    "output",
                    "absolute-10_scale-0.7",
                ),
            },
        ),
        (
            # Test to not scale certain parameters
            {
                "input": pathlib.Path(
                    "test", "test_adjust_functions", "adjust_scale", "input"
                ),
            },
            10,
            "MW",
            0.7,
            None,  # default
            ["capacity"],  # Only scale the capacity parameter
            {
                "output": pathlib.Path(
                    "test",
                    "test_adjust_functions",
                    "adjust_scale",
                    "output",
                    "absolute-10_scale-0.7_unmodified-investment",
                ),
            },
        ),
        (
            {
                "input": pathlib.Path(
                    "test", "test_adjust_functions", "adjust_scale", "input"
                ),
            },
            1,
            "kg",  # Incompatible unit for scaling
            1.0,
            None,  # default
            None,  # default
            NotImplementedError,
        ),
    ],
)  # type: ignore
def test_adjust_scale(
    input_technologies: dict[str, pathlib.Path],
    to_capacity: float,
    to_capacity_unit: str,
    scaling_exponent: float,
    scaled_parameters: list[str] | None,
    absolute_parameters: list[str] | None,
    expected_technologies: dict[str, pathlib.Path] | type[Exception],
) -> None:
    """Test adjust_scale using input/output files and Technologies object comparison."""
    tech = td.Technologies(input_technologies)
    if isinstance(expected_technologies, dict):
        tech.adjust_scale(
            to_capacity=to_capacity,
            to_capacity_unit=to_capacity_unit,
            scaling_exponent=scaling_exponent,
            scaled_parameters=scaled_parameters,
            absolute_parameters=absolute_parameters,
        )
        expected = td.Technologies(expected_technologies)
        pd.testing.assert_frame_equal(
            tech.data, expected.data, check_like=True, atol=1e-4
        )
    else:
        with pytest.raises(expected_technologies):
            tech.adjust_scale(
                to_capacity=to_capacity,
                to_capacity_unit=to_capacity_unit,
                scaling_exponent=scaling_exponent,
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
