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


## Licence

Copyright 2019-2025 [Contributors](https://github.com/PyPSA/technology-data/graphs/contributors) to technology-data

The code in `scripts/` is released as free software under the
[GPLv3](http://www.gnu.org/licenses/gpl-3.0.en.html), see LICENSE.txt.
However, different licenses and terms of use may apply to the various
input data.
