<!--
SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
SPDX-License-Identifier: GPL-3.0-only
-->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/pypsa/technology-data?include_prereleases)
[![Documentation](https://readthedocs.org/projects/technology-data/badge/?version=latest)](https://technology-data.readthedocs.io/en/latest/?badge=latest)
![Licence](https://img.shields.io/github/license/pypsa/technology-data)
![Size](https://img.shields.io/github/repo-size/pypsa/technology-data)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3994163.svg)](https://doi.org/10.5281/zenodo.3994163)
[![Gitter](https://badges.gitter.im/PyPSA/community.svg)](https://gitter.im/PyPSA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)


# Energy System Technology Data

This script compiles assumptions on energy system technologies (such
as costs, efficiencies, lifetimes, etc.)  for chosen years
(e.g. [2020, 2030, 2050]) from a variety of sources into CSV files to
be read by energy system modelling software. The merged outputs have
standardized cost years, technology names, units and source information. For further information about the structure and how to add new technologies, see the [documentation](https://technology-data.readthedocs.io/en/latest/).


The outputs are used in
[PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) and [PyPSA-Earth](https://github.com/pypsa-meets-earth/pypsa-earth).

## Contributions welcome!

We are happy to receive contributions, either as feedback, corrections, modifications or addition of new technology entries.
You can either open an *Issue* on GitHub for reporting potentials problems, providing feedback or suggesting additional entries.
Alternatively you can also directly open a *Pull Request* already containing suggested modifications.

All sizes and qualities of contributions are welcome. After an the *Issue* or *Pull Request* we'll assist with comments and leads to get your contribution included into the repository.

## Unit conventions

The following conventions are recommended for new additions to the repository.
Care must be taken with legacy entries, of which not all follow these conventions.

* Energy units: Use MWh and MW or kWh and kW
  Thermal energy content: When referring to the thermal energy content of a mass or volume, use the Lower Heating Value (LHV)
* Currency-values are specified in EUR
* Ambigiuous units are avoided by specifying whether the unit applies to the input or the output of a process, e.g. the capacity of hydrogen electrolysis "MW_e" refers to "MW" of electricity input capacity whilst "MW_H2" refers to "MW" of output hydrogen
* Specification is always done using a subscript, e.g. "MWh_e" (MWh of electricity), "MWh_th" (MWh of thermal energy), "t_CO2" (t of CO2)
* Multiple units are concatenated using an asterisk an space before and after the asterisk " * ", e.g. "m^3 * h"
* Combined units are written as a single unit, e.g. "MWh" (instead of "MW*h")
* Multiple units in the demoninator are encapsulated in brackets, e.g. "EUR/(m^3_water * h)" instead of "EUR/m^3_water/h


## Licence

Copyright 2019-2025 [Contributors](https://github.com/PyPSA/technology-data/graphs/contributors) to technology-data

The code in `scripts/` is released as free software under the
[GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html), see LICENSE.txt.
However, different licenses and terms of use may apply to the various
input data.
