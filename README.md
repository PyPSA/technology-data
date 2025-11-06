# technologydata

<!--
SPDX-FileCopyrightText: The technology-data authors

SPDX-License-Identifier: MIT

-->

<!--

TODO: Add badges
**Suggestions:**

- Use http://shields.io or a similar service to create and host the images.
- Add the [Standard Readme badge](https://github.com/RichardLitt/standard-readme#badge).

-->

A Python package to manage techno-economic assumptions for energy system models.

## Overview

`technologydata` is a Python package that supports the management of techno-economic assumptions for energy system models.
It provides a structured way to store, retrieve, and manipulate data related to various technologies used in energy systems,
including unit-ful parameters, currency conversions, inflation adjustment, and temporal modelling.

In the future it will include pre-parsed data from common public data sources such NREL's ATB or DEA's Technology Catalogue.
For these datasources it will also include parsers that can be modified to extract data in a custom way.

The goal of this package is to make energy system modelling easier and more efficient,
automating common tasks and transformations to reduce errors and allowing for easier data exchange between models.

## Table of Contents

1. [Background](#background)
2. [Install](#install)
3. [Usage](#usage)
4. [Maintainers](#maintainers)
5. [Thanks](#thanks)
6. [Contributing](#contributing)
7. [License](#license)

## Background

> Modelling is 10% science, 10% art, and 80% finding the right data and getting it into the right format.
> — Every energy modeller ever

Modelling energy systems requires a lot of data.
Techno-economic data, i.e. data about the costs for building operating technologies and their technical
characteristics, is a key input to many energy system models today.

Techno-economic data is usually collected from a variety of scattered sources and then manually processed
into the format required by a specific model.
This is repeated for every new modelling project, leading to a lot of duplicated effort.
The manual processing also carries a high risk of errors, which can lead to misleading results.

When projects are finished, the processed data is often discarded, leading to a loss of valuable information.
In better cases, the processed data is also published along with the model results, but usually in a
non-standardised and non-machine-readable format and without information about the data provenance and processing steps.

It sounds abstract, but if someone used cost assumptions for the US in 2015 USD, then there are many wrong ways
and a few right ways to convert these to e.g. EUR and adjust it for inflation to 2023.

## Install

The package is currently under development and not yet published to `PyPI` or `conda-forge`.

To install the package locally from GitHub, first clone the package and then use `uv` to install it in editable mode:

```bash
git clone https://github.com/open-energy-transition/technology-data/tree/prototype-2
cd technology-data
git checkout prototype-2
uv sync --group dev --group docs
```

## Usage

Detailed usage instructions and examples can be found in the [documentation](https://technology-data--240.org.readthedocs.build/en/240/).

To create a `Technology` object with some parameters and then adjust for inflation and convert to a different currency, you can use the following code:

```python
from technologydata import Technology, Parameter

tech = Technology(
    name="Solar PV",
    detailed_technology="Utility-scale PV",
    region="USA",
    year=2020,
    parameters={
        "capital_cost": Parameter(value=600, unit="USD_2020/kW"),
        "efficiency": Parameter(value=15, unit="percent"),
    }
)

tech = tech.to_currency("EUR_2023")
```

## Maintainers

This repository is currently maintained by [Open Energy Transition](https://openenergytransition.org/) with the maintainers and developers being:

- [euronion](https://github.com/euronion)
- [finozzifa](https://github.com/finozzifa)

## Thanks

Development of this prototype package would not have been possible without the funding from [Breakthrough Energy](https://www.breakthroughenergy.org/).

## Contributing

For contributing instructions, guidelines and our code of conduct, please refer to the [contributing section](https://technology-data--240.org.readthedocs.build/en/240/contributing/instructions/) in the documentation.

## License

This project is licensed under the [MIT License](LICENSES/MIT.txt).

Primary data included in the repository may be licensed under specific terms.
Processed data included in the project is licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSES/CC-BY-4.0.txt).

To make it easier to identify which data is licensed under which terms, this repository follows the [REUSE](https://reuse.software/) specification.
This means you can find the license information for each file either located in its header or in [REUSE.toml](REUSE.toml).
