# SPDX-FileCopyrightText: technologydata contributors
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

from technologydata import Parameter, Source, SourceCollection

EnergyDensityLHV: dict[str, Parameter] = dict(
    hydrogen=Parameter(
        magnitude=120,
        units="MJ/kg",
        carrier="hydrogen",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    methane=Parameter(
        magnitude=50.0,
        units="MJ/kg",
        carrier="methane",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    natural_gas=Parameter(
        magnitude=47.1,
        units="MJ/kg",
        carrier="natural_gas",
        note="Value for natural gas in the US market",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    ammonia=Parameter(
        magnitude=18.646,
        units="MJ/kg",
        carrier="ammonia",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Heat of combustion",
                    authors="Wikipedia contributors",
                    url="https://en.wikipedia.org/w/index.php?title=Heat_of_combustion&oldid=1307512525",
                    url_date="2025-09-09",
                ),
            ],
        ),
    ),
    wood=Parameter(
        magnitude=15.4,
        units="MJ/kg",
        carrier="wood",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    carbon=Parameter(
        magnitude=32.8,
        units="MJ/kg",
        carrier="carbon",
        note="For pure carbon we assume LHV=HHV",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    methanol=Parameter(
        magnitude=19.9,
        units="MJ/kg",
        carrier="methanol",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    gasoline=Parameter(
        magnitude=43.4,
        units="MJ/kg",
        carrier="gasoline",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    # Add more energy densities as needed
)

EnergyDensityHHV: dict[str, Parameter] = dict(
    hydrogen=Parameter(
        magnitude=141.7,
        units="MJ/kg",
        carrier="hydrogen",
        sources=SourceCollection(
            sources=[
                Source(
                    title="The Engineering ToolBox",
                    authors="The Engineering ToolBox (2003). Higher Calorific Values of Common Fuels: Reference & Data. [online] Available at:https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html [Accessed 6 September 2025].",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                ),
            ],
        ),
    ),
    methane=Parameter(
        magnitude=55.5,
        units="MJ/kg",
        carrier="methane",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    natural_gas=Parameter(
        magnitude=52.2,
        units="MJ/kg",
        carrier="natural_gas",
        note="Value for natural gas in the US market",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    ammonia=Parameter(
        magnitude=22.5,
        units="MJ/kg",
        carrier="ammonia",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    wood=Parameter(
        magnitude=16.2,
        units="MJ/kg",
        carrier="wood",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    carbon=Parameter(
        magnitude=32.8,
        units="MJ/kg",
        carrier="carbon",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    lignite=Parameter(
        magnitude=14.0,
        units="MJ/kg",
        carrier="lignite",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    coal=Parameter(
        magnitude=32.6,
        units="MJ/kg",
        carrier="coal",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    methanol=Parameter(
        magnitude=23.0,
        units="MJ/kg",
        carrier="methanol",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    jet_fuel_a1=Parameter(
        magnitude=46.2,
        units="MJ/kg",
        carrier="jet_fuel_a1",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Jet fuel",
                    authors="Wikipedia contributors",
                    url="https://en.wikipedia.org/w/index.php?title=Jet_fuel&oldid=1310114781",
                    url_date="2025-09-07",
                ),
            ],
        ),
    ),
    gasoline=Parameter(
        magnitude=46.4,
        units="MJ/kg",
        carrier="gasoline",
        sources=SourceCollection(
            sources=[
                Source(
                    title="Higher Calorific Values of Common Fuels: Reference & Data",
                    authors="The Engineering ToolBox (2003)",
                    url="https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html",
                    url_date="2025-09-06",
                ),
            ],
        ),
    ),
    # Add more energy densities as needed
)
