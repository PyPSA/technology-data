import pandas as pd
import pytest

import technologydata as td


@pytest.fixture
def example_source():
    """Fixture to provide the example01 source."""
    return td.Technologies(packaged_sources=["example01"], load=True)


@pytest.fixture
def forecast_source():
    """Fixture to provide an example dataset for time-related forecasting."""
    df = pd.DataFrame(
        {
            "source": ["test", "test"],
            "technology": ["example tech", "example tech"],
            "detailed_technology": [
                "example tech - detailed",
                "example tech - detailed",
            ],
            "case": ["example case", "example case"],
            "region": ["IE", "IE"],
            "parameter": ["investment", "investment"],
            "year": [2020, 2030],
            "value": [100, 200],
            "unit": ["EUR_2020/MW_electric", "EUR_2020/MW_electric"],
            "scale": [1, 1],
            "scale_unit": ["MW", "MW"],
            "comment": ["", ""],
        }
    )
    return td.Technologies().from_pandas(df)


def test_no_economies_of_scale(example_source) -> None:
    """Without economies of scale the value should not change."""
    org_data = example_source.data.copy()
    assert all(
        example_source.adjust_scale(
            new_scale=999,
            unit="MW",
            scaling_exponent=1,
        ).data["value"]
        == org_data["value"]
    ), "Scaling with exponent 1 should not change the value"


def test_economies_of_scale(example_source) -> None:
    """Test with common economies of scale with exponent 0.5."""
    org_data = example_source.data.copy()
    assert all(
        example_source.adjust_scale(
            new_scale=2,
            unit="MW",
            scaling_exponent=0.5,
        ).data["value"]
        == org_data["value"] * (2 / org_data["scale"]) ** (0.5 - 1)
    ), "Scaling with exponent 0.5 should change the value to approx 0.7"


def test_adjust_year_linear_interpolation(forecast_source) -> None:
    """Test linear forecasting through a middle value in 2025."""
    forecast = forecast_source.adjust_year(year=2025, model={"method": "linear"})

    assert forecast.shape[0] == 3, "Forecasted data missing additional row"
    assert forecast.iloc[1]["value"] == 150, (
        "Linear forecasted value for 2025 is incorrect"
    )


def test_adjust_year_linear_extrapolation(forecast_source) -> None:
    """Test linear forecasting for a value outside the range of entries provided."""
    forecast = forecast_source.adjust_year(
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
