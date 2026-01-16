# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the initialization and methods of the Parameter class."""

import pathlib

import numpy as np
import pandas as pd
import pint
import pytest

import technologydata
from technologydata.constants import EnergyDensityLHV
from technologydata.utils.units import CustomUndefinedUnitError

path_cwd = pathlib.Path.cwd()


class TestParameter:
    """Test suite for the Parameter class in the technologydata module."""

    def test_parameter_creation(self) -> None:
        """Test the creation of a Parameter instance with various units."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
            sources=technologydata.SourceCollection(
                sources=[
                    technologydata.Source(
                        authors="some authors",
                        title="Example Title",
                    )
                ]
            ),
        )
        assert param.magnitude == 1000
        assert param.units == "USD_2020 / kilowatt"
        assert param.provenance == "literature"
        assert param.note == "Estimated"
        assert param.sources is not None

        # Internal pint quantities
        # pint autoconverts to the canonical name for the carrier and heating value
        # so comparing with == "H2" would fail; instead check that the units are compatible
        assert param._pint_quantity.units.is_compatible_with("USD_2020 / kW")
        assert param._pint_carrier.is_compatible_with("H2")
        assert param._pint_heating_value.is_compatible_with("LHV")

    def test_parameter_invalid_units(self) -> None:
        """Test that an error is raised when invalid units are provided."""
        with pytest.raises(CustomUndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="USD",
            )
        assert (
            str(excinfo.value)
            == "Currency unit 'USD' is missing the 4-digit currency year (e.g. USD_2020)."
        )
        assert excinfo.type == CustomUndefinedUnitError

        with pytest.raises(CustomUndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="USD_20/kW",
            )
        assert (
            str(excinfo.value)
            == "Currency unit 'USD' is missing the 4-digit currency year (e.g. USD_2020)."
        )
        assert excinfo.type == CustomUndefinedUnitError

        with pytest.raises(CustomUndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="abc_USD_20/kW",
            )
        assert (
            str(excinfo.value)
            == "Currency unit 'USD' is missing the 4-digit currency year (e.g. USD_2020)."
        )
        assert excinfo.type == CustomUndefinedUnitError

        with pytest.raises(pint.errors.UndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="INVALID_UNIT",
            )
        assert "INVALID_UNIT" in str(excinfo.value)

        with pytest.raises(pint.errors.UndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="USD_2020/kW",
                carrier="H2",
                heating_value="INVALID_HEATING_VALUE",
            )
        assert "INVALID_HEATING_VALUE" in str(excinfo.value)

        with pytest.raises(pint.errors.UndefinedUnitError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="USD_2020/kW",
                carrier="INVALID_CARRIER",
            )
        assert "INVALID_CARRIER" in str(excinfo.value)

        with pytest.raises(ValueError) as excinfo:
            technologydata.Parameter(
                magnitude=1000,
                units="USD_2020/kW",
                heating_value="LHV",
            )
        assert "Heating value cannot be set without a carrier" in str(excinfo.value)

    def test_parameter_to_conversion_fail_on_currency_conversion(self) -> None:
        """Test the unit conversion of a Parameter instance to fail on currency conversion."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
        )
        with pytest.raises(NotImplementedError) as excinfo:
            param.to("EUR_2025 / kilowatt")
            assert "Currency conversion is not supported" in str(excinfo.value)

    def test_parameter_to_conversion(self) -> None:
        """Test the conversion of a Parameter instance to different units."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
        )

        ref = technologydata.Parameter(
            magnitude=param.magnitude * 1000,
            units="USD_2020 / megawatt",
        )

        converted = param.to("USD_2020 / megawatt")

        assert isinstance(converted, technologydata.Parameter)
        assert converted.units == ref.units
        assert converted.magnitude == ref.magnitude

    def test_pint_attributes_update(self) -> None:
        """Test that pint attributes are updated correctly when attributes change and a method that uses the pint fields is called."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
        )
        # Change magnitude and units
        param.magnitude = 2000
        param.units = "EUR_2025 / kilowatt"
        param._update_pint_attributes()
        assert param._pint_quantity.magnitude == 2000
        assert str(param._pint_quantity.units) == "EUR_2025 / kilowatt"

        # Change magnitude and units again and check if they are updated
        # when calling a method that uses them
        param.magnitude = 3000
        param.units = "USD_2020 / kWh"
        param = param.to("USD_2020 / kWh")
        assert param._pint_quantity.magnitude == 3000
        assert str(param._pint_quantity.units) == str(
            technologydata.ureg.Unit("USD_2020 / kWh")
        )

    def test_parameter_to_currency(self) -> None:
        """Test currency conversion with inflation adjustment."""
        param = technologydata.Parameter(
            magnitude=1,
            units="USD_2020/kW",
        )

        # Convert to EUR with inflation adjustment for Germany
        converted = param.to_currency("EUR_2023", "DEU")
        assert isinstance(converted, technologydata.Parameter)
        assert converted.units is not None
        assert "EUR_2023" in converted.units
        assert converted._pint_quantity is not None
        assert converted._pint_quantity.is_compatible_with("EUR_2023 / kW")

        # Check that magnitude changed due to currency conversion
        assert not np.isnan(converted.magnitude)
        assert converted.magnitude != param.magnitude

    def test_parameter_change_currency_explicit_source(self) -> None:
        """Test currency conversion with explicit inflation data source."""
        param = technologydata.Parameter(
            magnitude=1,
            units="EUR_2019/MWh",
        )

        # Convert using IMF data source
        converted = param.to_currency("USD_2022", "USA", source="worldbank")
        assert isinstance(converted, technologydata.Parameter)
        assert converted.units is not None
        assert "USD_2022" in converted.units
        assert converted._pint_quantity is not None
        assert not np.isnan(converted.magnitude)
        assert converted._pint_quantity.is_compatible_with("USD_2022 / MWh")

    def test_parameter_change_currency_different_source(self) -> None:
        """Test currency conversion with different inflation data source."""
        param = technologydata.Parameter(
            magnitude=1,
            units="EUR_2019/MWh",
        )

        # Convert using IMF data source
        converted = param.to_currency("USD_2022", "USA", source="imf")
        assert isinstance(converted, technologydata.Parameter)
        assert converted.units is not None
        assert "USD_2022" in converted.units
        assert converted._pint_quantity.is_compatible_with("USD_2022 / MWh")

    def test_parameter_change_currency_multiple_currencies(self) -> None:
        """Test currency conversion when units contain multiple currencies."""
        param = technologydata.Parameter(
            magnitude=1,
            units="USD_2020 * EUR_2021 / kW",
        )

        # Convert all currencies to CNY_2023
        converted = param.to_currency("CNY_2023", "CHN")
        assert isinstance(converted, technologydata.Parameter)
        # Both USD_2020 and EUR_2021 should be replaced with CNY_2023
        assert converted.units is not None
        assert "CNY_2023" in converted.units
        assert "USD_2020" not in converted.units
        assert "EUR_2021" not in converted.units

    def test_parameter_change_currency_same_currency(self) -> None:
        """Test currency conversion to the same currency (inflation adjustment only)."""
        param = technologydata.Parameter(
            magnitude=1,
            units="USD_2019/kW",
        )

        # Convert to USD but different year (inflation adjustment)
        converted = param.to_currency("USD_2023", "USA")
        assert isinstance(converted, technologydata.Parameter)
        assert converted.units == "USD_2023 / kilowatt"
        # Magnitude should change due to inflation adjustment
        assert not np.isnan(converted.magnitude)
        assert converted.magnitude != param.magnitude

    def test_parameter_no_currency_change(self) -> None:
        """Test that no currency change occurs when the target currency is the same as the current one."""
        param = technologydata.Parameter(
            magnitude=1,
            units="USD_2020/kW",
        )

        # Convert to the same currency and year
        converted = param.to_currency("USD_2020", "USA")
        assert isinstance(converted, technologydata.Parameter)
        assert converted._pint_quantity.is_compatible_with("USD_2020 / kW")
        # Magnitude should remain unchanged
        assert not np.isnan(converted.magnitude)
        assert converted.magnitude == param.magnitude

    def test_parameter_change_currency_invalid_country(self) -> None:
        """Test that invalid country codes raise appropriate errors."""
        param = technologydata.Parameter(
            magnitude=1,
            units="USD_2020/kW",
        )

        # Invalid country code should raise an error
        with pytest.raises((ValueError, KeyError)):
            param.to_currency("EUR_2023", "USB")

    def test_parameter_change_currency_invalid_source(self) -> None:
        """Test that invalid inflation data sources raise appropriate errors."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
        )

        # Invalid source should raise an error
        with pytest.raises(KeyError):
            param.to_currency("EUR_2023", "DEU", source="invalid_source")

    def test_parameter_change_currency_no_units(self) -> None:
        """Test currency conversion with parameter that has no units."""
        param = technologydata.Parameter(
            magnitude=42,
        )

        # Should handle parameters without currency units gracefully
        converted = param.to_currency("EUR_2023", "DEU")
        assert isinstance(converted, technologydata.Parameter)
        assert converted.magnitude == 42
        assert converted.units is None or "EUR_2023" not in str(converted.units)

    def test_parameter_unchanged_other_attributes(self) -> None:
        """Test that other attributes remain unchanged after currency conversion."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
        )

        # Convert to EUR_2023
        converted = param.to_currency("EUR_2023", "DEU")

        # Check that other attributes remain unchanged
        assert converted.carrier == param.carrier
        assert converted.heating_value == param.heating_value
        assert converted.provenance == param.provenance
        assert converted.note == param.note

    def test_parameter_heating_value_compatibility(self) -> None:
        """Test that we do not permit operations on mixed heating values."""
        param_lhv = technologydata.Parameter(
            magnitude=1,
            carrier="H2",
            heating_value="LHV",
        )
        param_hhv = technologydata.Parameter(
            magnitude=1,
            carrier="H2",
            heating_value="HHV",
        )

        # Different error messages for + and -
        with pytest.raises(ValueError) as excinfo:
            param_lhv + param_hhv
            assert (
                "Operation not permitted on parameters with different heating values"
                in str(excinfo.value)
            )
        with pytest.raises(ValueError) as excinfo:
            param_lhv - param_hhv
            assert (
                "Operation not permitted on parameters with different heating values"
                in str(excinfo.value)
            )

        # Different error messages for * and /
        with pytest.raises(ValueError) as excinfo:
            param_lhv * param_hhv
            assert "Cannot multiply parameters with different heating values" in str(
                excinfo.value
            )

        with pytest.raises(ValueError) as excinfo:
            param_lhv / param_hhv
            assert "Cannot divide parameters with different heating values" in str(
                excinfo.value
            )

    def test_parameter_carrier_compatibility(self) -> None:
        """Test that we do only permit certain operations on mixed carriers."""
        param_h2 = technologydata.Parameter(
            magnitude=1,
            carrier="H2",
            heating_value="LHV",
        )
        param_ch4 = technologydata.Parameter(
            magnitude=1,
            carrier="CH4",
            heating_value="LHV",
        )

        # Different error messages for + and -
        with pytest.raises(ValueError) as excinfo:
            param_h2 + param_ch4
            assert (
                "Operation not permitted on parameters with different carriers"
                in str(excinfo.value)
            )
        with pytest.raises(ValueError) as excinfo:
            param_h2 - param_ch4
            assert (
                "Operation not permitted on parameters with different carriers"
                in str(excinfo.value)
            )

        # * and / are permitted with different carriers and should yield mixed carriers
        assert (
            param_h2 * param_ch4
        ).carrier == param_h2._pint_carrier * param_ch4._pint_carrier
        assert (
            param_h2 / param_ch4
        ).carrier == param_h2._pint_carrier / param_ch4._pint_carrier

    def test_parameter_equality(self) -> None:
        """Test equality comparison of Parameter objects."""
        # Create two identical parameters
        param1 = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
            sources=technologydata.SourceCollection(
                sources=[
                    technologydata.Source(
                        authors="some authors",
                        title="Example Title",
                    )
                ]
            ),
        )
        param2 = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
            sources=technologydata.SourceCollection(
                sources=[
                    technologydata.Source(
                        authors="some authors",
                        title="Example Title",
                    )
                ]
            ),
        )

        # Should be equal
        assert param1 == param2
        assert param2 == param1

    def test_parameter_equality_different_magnitude(self) -> None:
        """Test that parameters with different magnitudes are not equal."""
        param1 = technologydata.Parameter(magnitude=1000, units="USD_2020/kW")
        param2 = technologydata.Parameter(magnitude=2000, units="USD_2020/kW")

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_units(self) -> None:
        """Test that parameters with different units are not equal."""
        param1 = technologydata.Parameter(magnitude=1000, units="USD_2020/kW")
        param2 = technologydata.Parameter(magnitude=1000, units="EUR_2020/kW")

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_carrier(self) -> None:
        """Test that parameters with different carriers are not equal."""
        param1 = technologydata.Parameter(magnitude=1000, carrier="H2")
        param2 = technologydata.Parameter(magnitude=1000, carrier="CH4")

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_heating_value(self) -> None:
        """Test that parameters with different heating values are not equal."""
        param1 = technologydata.Parameter(
            magnitude=1000, carrier="H2", heating_value="LHV"
        )
        param2 = technologydata.Parameter(
            magnitude=1000, carrier="H2", heating_value="HHV"
        )

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_provenance(self) -> None:
        """Test that parameters with different provenance are not equal."""
        param1 = technologydata.Parameter(magnitude=1000, provenance="literature")
        param2 = technologydata.Parameter(magnitude=1000, provenance="expert_estimate")

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_note(self) -> None:
        """Test that parameters with different notes are not equal."""
        param1 = technologydata.Parameter(magnitude=1000, note="Estimated")
        param2 = technologydata.Parameter(magnitude=1000, note="Measured")

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_different_sources(self) -> None:
        """Test that parameters with different sources are not equal."""
        source1 = technologydata.Source(authors="Author A", title="Title A")
        source2 = technologydata.Source(authors="Author B", title="Title B")

        param1 = technologydata.Parameter(
            magnitude=1000, sources=technologydata.SourceCollection(sources=[source1])
        )
        param2 = technologydata.Parameter(
            magnitude=1000, sources=technologydata.SourceCollection(sources=[source2])
        )

        assert param1 != param2
        assert param2 != param1

    def test_parameter_equality_none_values(self) -> None:
        """Test equality with None values for optional fields."""
        param1 = technologydata.Parameter(magnitude=1000)
        param2 = technologydata.Parameter(magnitude=1000)

        # Both have None for optional fields, should be equal
        assert param1 == param2

        # One has None, the other has a value
        param3 = technologydata.Parameter(magnitude=1000, units="USD_2020/kW")
        assert param1 != param3
        assert param3 != param1

    def test_parameter_equality_with_non_parameter(self) -> None:
        """Test equality comparison with non-Parameter objects."""
        param = technologydata.Parameter(magnitude=1000)

        # Should return NotImplemented for non-Parameter objects
        assert param.__eq__("not a parameter") == NotImplemented
        assert param.__eq__(42) == NotImplemented
        assert param.__eq__(None) == NotImplemented

    def test_parameter_equality_canonical_units(self) -> None:
        """Test that parameters with equivalent but differently formatted units are equal."""
        # Units should be canonicalized during initialization
        param1 = technologydata.Parameter(magnitude=1000, units="USD_2020/kW")
        param2 = technologydata.Parameter(magnitude=1000, units="USD_2020/kilowatt")

        # Should be equal because units are canonicalized
        assert param1 == param2
        assert param2 == param1

    def test_parameter_equality_self_reference(self) -> None:
        """Test that a parameter is equal to itself."""
        param = technologydata.Parameter(
            magnitude=1000,
            units="USD_2020/kW",
            carrier="H2",
            heating_value="LHV",
            provenance="literature",
            note="Estimated",
        )

        assert param == param

    @pytest.mark.parametrize(
        "example_parameter",
        [
            {
                "parameter_magnitude": 1000,
                "parameter_units": "USD_2020/kW",
                "parameter_carrier": "H2",
                "parameter_heating_value": "LHV",
                "parameter_provenance": "literature",
                "parameter_note": "Estimated",
                "parameter_sources": [
                    technologydata.Source(title="title", authors="authors")
                ],
            }
        ],
        indirect=["example_parameter"],
    )  # type: ignore
    def test_example_parameter(
        self, example_parameter: technologydata.Parameter
    ) -> None:
        """Test that the fixture example_parameter yields a parameter object."""
        assert isinstance(example_parameter, technologydata.Parameter)

    def test_parameter_mul_scalar(self) -> None:
        """Test multiplication of a Parameter by a scalar."""
        factor = 3.0
        param = technologydata.Parameter(magnitude=2, units="kW")
        result = param * factor
        assert isinstance(result, technologydata.Parameter)
        assert result.magnitude == param.magnitude * factor
        assert result.units == param.units

        # Test multiplication by an integer
        factor = int(3.0)
        param = technologydata.Parameter(magnitude=2, units="kW")
        result = param * factor
        assert isinstance(result, technologydata.Parameter)
        assert result.magnitude == param.magnitude * factor
        assert result.units == param.units

    def test_parameter_div_scalar(self) -> None:
        """Test division of a Parameter by a scalar."""
        # Test division by an float
        divisor = 3.0
        param = technologydata.Parameter(magnitude=6, units="kW")
        result = param / divisor
        assert isinstance(result, technologydata.Parameter)
        assert result.magnitude == param.magnitude / divisor
        assert result.units == param.units

        # Test division by an integer
        divisor = int(3.0)
        param = technologydata.Parameter(magnitude=6, units="kW")
        result = param / divisor
        assert isinstance(result, technologydata.Parameter)
        assert result.magnitude == param.magnitude / divisor
        assert result.units == param.units

    def test_parameter_pow_basic(self) -> None:
        """Test integer exponentiation."""
        param = technologydata.Parameter(magnitude=2, units="kW")
        result = param**3
        assert isinstance(result, technologydata.Parameter)
        assert result.magnitude == 8
        assert result.units == "kilowatt ** 3"

    def test_parameter_pow_fractional(self) -> None:
        """Test fractional exponentiation."""
        param = technologydata.Parameter(magnitude=9, units="m**2")
        result = param**0.5
        assert pytest.approx(result.magnitude) == 3
        assert result.units == "meter"

    def test_parameter_pow_zero(self) -> None:
        """Test zero exponent returns dimensionless."""
        param = technologydata.Parameter(magnitude=5, units="J")
        result = param**0
        assert result.magnitude == 1
        assert result.units == "dimensionless"

    def test_parameter_pow_negative(self) -> None:
        """Test negative exponentiation and unit handling."""
        param = technologydata.Parameter(magnitude=2, units="W")
        result = param**-2
        assert pytest.approx(result.magnitude) == 0.25
        # Compare units using pint, not string equality
        assert technologydata.ureg.Unit(result.units) == technologydata.ureg.Unit(
            "watt ** -2"
        )

    def test_parameter_pow_carrier(self) -> None:
        """Test that the carrier attribute is also affected."""
        param = technologydata.Parameter(
            magnitude=3,
            units="kg",
            carrier="H2",
        )
        result = param**2
        assert result.carrier == f"{param.carrier} ** 2"

    def test_parameter_pow_preserves_metadata(self) -> None:
        """Test that metadata is preserved after exponentiation."""
        param = technologydata.Parameter(
            magnitude=3,
            units="kg",
            carrier="H2",
            heating_value="LHV",
            provenance="test",
            note="note",
        )
        result = param**2
        assert result.carrier == f"{param.carrier} ** 2"
        assert result.heating_value == param.heating_value
        assert result.provenance == param.provenance
        assert result.note == param.note

    def test_change_heating_value_h2_lhv_to_hhv(self) -> None:
        """Test LHV to HHV conversion for H2."""
        p = technologydata.Parameter(
            magnitude=119.6,
            units="kilowatt_hour",
            carrier="hydrogen",
            heating_value="lower_heating_value",
        )
        p2 = p.change_heating_value("higher_heating_value")
        assert pytest.approx(p2.magnitude) == 141.2278
        assert p2.heating_value == "higher_heating_value"
        assert p2.carrier == "hydrogen"
        assert p2.units == "kilowatt_hour"

    def test_change_heating_value_h2_hhv_to_lhv(self) -> None:
        """Test HHV to LHV conversion for H2."""
        p = technologydata.Parameter(
            magnitude=141.8,
            units="kilowatt_hour",
            carrier="hydrogen",
            heating_value="higher_heating_value",
        )
        p2 = p.change_heating_value("lower_heating_value")
        assert pytest.approx(p2.magnitude) == 120.0848
        assert p2.heating_value == "lower_heating_value"
        assert p2.carrier == "hydrogen"
        assert p2.units == "kilowatt_hour"

    def test_change_heating_value_ch4_lhv_to_hhv(self) -> None:
        """Test LHV to HHV conversion for CH4."""
        p = technologydata.Parameter(
            magnitude=10,
            units="kilowatt_hour",
            carrier="methane",
            heating_value="lower_heating_value",
        )
        p2 = p.change_heating_value("higher_heating_value")
        assert pytest.approx(p2.magnitude) == 11.1
        assert p2.heating_value == "higher_heating_value"
        assert p2.carrier == "methane"
        assert p2.units == "kilowatt_hour"

    def test_change_heating_value_ch4_hhv_to_lhv(self) -> None:
        """Test HHV to LHV conversion for CH4."""
        p = technologydata.Parameter(
            magnitude=11.1,
            units="kilowatt_hour",
            carrier="methane",
            heating_value="higher_heating_value",
        )
        p2 = p.change_heating_value("lower_heating_value")
        assert pytest.approx(p2.magnitude) == 10.0
        assert p2.heating_value == "lower_heating_value"
        assert p2.carrier == "methane"
        assert p2.units == "kilowatt_hour"

    def test_change_heating_value_ch4_hhv_to_lhv_adapt_units(self) -> None:
        """Test HHV to LHV conversion for CH4, where HHV has a different unit than LHV."""
        p = technologydata.Parameter(
            magnitude=11.1,
            units="kilowatt_hour",
            carrier="methane",
            heating_value="higher_heating_value",
        )
        # Note: metric_ton = 1e3 * kilogram = t = tonne
        # Note: ton = 2e3 * pound = _ = short_ton
        EnergyDensityLHV["methane"].units = "MJ/metric_ton"
        EnergyDensityLHV["methane"].magnitude = 50000
        p2 = p.change_heating_value("lower_heating_value")
        assert p2.units == "kilowatt_hour"
        assert pytest.approx(p2.magnitude) == 10.0
        assert p2.heating_value == "lower_heating_value"
        assert p2.carrier == "methane"

    def test_change_heating_value_no_carrier_in_units(self) -> None:
        """Test conversion when carrier does not appear in units (should treat as 1 appearance)."""
        p = technologydata.Parameter(
            magnitude=1,
            units="kilowatt_hour",
            carrier="electricity",
            heating_value="lower_heating_value",
        )
        p2 = p.change_heating_value("higher_heating_value")
        assert p2.magnitude == p.magnitude
        assert p2.heating_value == "higher_heating_value"
        assert p2.units == "kilowatt_hour"

    def test_change_heating_value_same_hv(self) -> None:
        """Test that no conversion occurs if heating value is unchanged."""
        p = technologydata.Parameter(
            magnitude=1,
            units="kilowatt_hour",
            carrier="hydrogen",
            heating_value="lower_heating_value",
        )
        p2 = p.change_heating_value("lower_heating_value")
        assert p2.magnitude == 1
        assert p2.heating_value == "lower_heating_value"

    @pytest.mark.parametrize(
        "folder_id",
        ["WB_CNY_2020", "WB_EUR_2020", "WB_USD_2020"],
    )  # type: ignore
    def test_to_currency(self, folder_id: str) -> None:
        """Validate the currency conversion rates."""
        input_path = pathlib.Path(
            path_cwd,
            "test",
            "test_data",
            "currency_conversion",
            folder_id,
            "parameters.csv",
        )
        input_dataframe = pd.read_csv(input_path).reset_index()
        for _, row in input_dataframe.iterrows():
            reference_param = technologydata.Parameter(
                magnitude=np.round(row["reference_magnitude"], 2),
                units=row["reference_units"],
            )

            output_param = technologydata.Parameter(
                magnitude=row["input_magnitude"],
                units=row["input_units"],
            ).to_currency(
                target_currency=technologydata.extract_currency_units(
                    row["reference_units"]
                )[0],
                country=row["country"],
                source=row["source"],
            )
            assert output_param.magnitude == pytest.approx(
                reference_param.magnitude, rel=1e-2
            )
            assert output_param.units == reference_param.units
