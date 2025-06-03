
# Technology Data for Energy System Models

[![CI](https://github.com/PyPSA/technology-data/actions/workflows/ci.yaml/badge.svg)](https://github.com/PyPSA/technology-data/actions/workflows/ci.yaml)
![GitHub Release](https://img.shields.io/github/v/release/pypsa/technology-data)
![GitHub last commit](https://img.shields.io/github/last-commit/pypsa/technology-data)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fopen-energy-transition%2Ftechnology-data%2Frefs%2Fheads%2Fprototype%2Fpyproject.toml)
![PyPI - Downloads](https://img.shields.io/pypi/dm/technologydata)

A Python package providing and harmonising techno-economic assumptions for Energy System Models.

*TODOs*:
* Check badges if correctly working (all sections)

**Table of Contents:**
1. [Features](#features)
2. [Background](#background)
3. [Installation](#installation)
4. [Examples](#examples)
5. [Documentation](#documentation)
6. [FAQ](#faq)
7. [Support and Feedback](#support-and-feedback)
8. [Contributing](#contributing)
9. [Used By](#used-by)
10. [Acknowledgements](#acknowledgements)
11. [Citing](#citing)
12. [License](#license)


## Features

Simplify integration of techno-econmic assumptions into your Energy System Model.
This package ships with:

- A pre-existing database of assumptions ready to got for your model
- Functions to pre-process and harmonise multiple datasets
- Common models for gap-filling and projection of parameters
- A data format specification to use custom and new datasets with this package

## Background

Energy System Models require a range of assumptions that need to be processed and harmonised in order to generate meaningful results.


*TODOs*:
* Expand introduction paragraph a bit
* Add info on motivation / vision
* Check badges if correctly working (all sections)

## Installation

![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fopen-energy-transition%2Ftechnology-data%2Frefs%2Fheads%2Fprototype%2Fpyproject.toml)
[![PyPI version](https://img.shields.io/pypi/v/technologydata.svg)](https://pypi.python.org/pypi/technologydata)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/technologydata.svg)](https://anaconda.org/conda-forge/technologydata)

Using `pip` or `uv` from `pypi`

```bash
# pip
pip install technologydata

# uv
uv pip install technologydata
```

Using `conda` from `conda-forge`

```bash
conda install -c conda-forge technologydata
```

For a development installation, clone the repository first and then install the package in `-e/--editable`:

```bash
git clone https://github.com/PyPSA/technology-data.git  <target_directory>
cd <target_directory>

# pip
pip install -e .

# uv
uv pip install -e .
```

*TODOs*:
* Fix link in badge "Python Version from PEP 621" to point to the right directory, shiled generator here: https://shields.io/badges/python-version-from-pep-621-toml
    
## Examples

*TODOs*:
* add one or two examples
* add a reference to the full documentation for more examples

```python
import technologydata as td

# TODO 
```


## Documentation

[![Documentation Status](https://readthedocs.org/projects/technology-data/badge/?version=latest)](https://technology-data.readthedocs.io/en/latest/?badge=latest)

The full documentation for the package [can be found here](http://technologydata.readthedocs.io/).


## FAQ

*TODOs*:
* Evolve the FAQ section, but focus on a high-level possible FAQs
* For all other FAQs, refer to the main documentation

## Support and Feedback

[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)

To get in touch with the community please join our [Discord Server](https://discord.gg/AnuJBk23FU).
The Discord server is also the right place for feedback.

If you have found a bug, want to suggest a new data source or feature, please open a [new issue](https://github.com/PyPSA/technology-data/issues/new/choose) in the GitHub repository.
## Contributing

Any contributions, however minor, are welcome.
If you plan to make code or data contributions, we will try to assist you as good as we can.


To ensure a high quality of the project, please also have a look at the instructions on contributing that we outline in the documentation.

*TODOs*:
* Create section in the documentation about contributing guidelines
* Add link to contribution guidelines and instructions in the documentation
## Used By
![GitHub Repo stars](https://img.shields.io/github/stars/pypsa/technology-data)
![Pepy Total Downloads](https://img.shields.io/pepy/dt/technologydata)

An extensive list of users can be found in the [documentation here](TODO://add.the.correct.link.here).

If you use this repository, we would love to hear from your and about your use case!
Drop us a message on Discord or add yourselves to the user list in the documentation via a Pull Request.


The following is a selection of projects and users of this package:

- [PyPSA-Eur](https://pypsa-eur.readthedocs.io): An open-source sector-coupled high-resolution Energy System Model for Europe
- [PyPSA-Earth](https://pypsa-earth.readthedocs.io): An open-source sector-coupled high-resolution Energy System Model with global coverage

*TODOs:*
* Add correct link to documentation above
## Acknowledgements

*TODOs*:
* Acknowledgements to funders that helped create and evolve the project
* Acknowledgements to major code contributors or all code authors and data integrators (?)
## Citing

To cite this package in your work, we suggest you use an attribution similar to the following one:

> \<Authors>, Technology Data for Energy System Models, 2025, https://github.com/PyPSA/technology-data, Version 0.12.1 .

Make sure to adjust the "Version" to match the actual version that you used for your work.

*TODOs*:
* Update \<Authors>
* Optional: Add CITATIONS.CFF file to repo for easier citing (note that the version in the CFF file should also be updated based on releases ... somehow)
* Include authors of previous TD versions as authors as well?

## License

![License](https://img.shields.io/pypi/l/technologydata.svg)
![REUSE Compliance](https://img.shields.io/reuse/compliance/github.com%2FPyPSA%2Ftechnology-data)


This repository is licensed under an [MIT license](TODO://add.the.right.link.here).

This repository is [REUSE compatible](https://reuse.software/) for more clarity and easier resuability.

*TODOs*
* Add the correct link to license in the repository in the text above
* Check licensing, if MIT is enough, or also other license for data, e.g. CC0 or CC BY 4.0
* Check that badges are working correctly