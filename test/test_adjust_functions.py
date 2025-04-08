import pytest

import technologydata as td


@pytest.fixture
def example_source():
    """Fixture to provide the example01 source."""
    return td.Technologies(packaged_sources=["example01"], load=True)


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
