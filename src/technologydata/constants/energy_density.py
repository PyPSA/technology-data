# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""
Energy density parameters for various carriers.

This module defines dictionaries containing energy density parameters
for different fuel carriers, specifically focusing on lower heating value (LHV)
and higher heating value (HHV). Each entry maps a carrier name to a `Parameter`
object that includes magnitude, unit, and carrier type.

Attributes
----------
EnergyDensityLHV : dict[str, Parameter]
    Dictionary mapping carrier names to their lower heating value (LHV)
    parameters.
    - Key: carrier name (e.g., 'hydrogen', 'methane')
    - Value: Parameter object with magnitude, unit, and carrier info.

EnergyDensityHHV : dict[str, Parameter]
    Dictionary mapping carrier names to their higher heating value (HHV)
    parameters.
    - Key: carrier name (e.g., 'hydrogen', 'methane')
    - Value: Parameter object with magnitude, unit, and carrier info.

"""

from technologydata import Parameter

EnergyDensityLHV: dict[str, Parameter] = dict(
    hydrogen=Parameter(
        magnitude=119.6,
        unit="GJ/t",
        carrier="hydrogen",
        # source=  # TODO
    ),
    methane=Parameter(
        magnitude=50.0,
        unit="GJ/t",
        carrier="methane",
        # source=  # TODO
    ),
    # Add more energy densities as needed
)

EnergyDensityHHV: dict[str, Parameter] = dict(
    hydrogen=Parameter(
        magnitude=141.8,
        unit="GJ/t",
        carrier="hydrogen",
        # source=  # TODO
    ),
    methane=Parameter(
        magnitude=55.5,
        unit="GJ/t",
        carrier="methane",
        # source=  # TODO
    ),
    # Add more energy densities as needed
)
