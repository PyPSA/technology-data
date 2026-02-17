# `DataAccessor` Class Documentation

<!--
SPDX-FileCopyrightText: technologydata contributors

SPDX-License-Identifier: MIT

-->

## Overview

The `DataAccessor` class provides a standardized way to access and load technology datasets from versioned data sources. It simplifies the process of selecting a specific data source and version, automatically handling the resolution of file paths and loading the data into a `DataPackage` object.

It is designed to work with a specific directory structure where datasets are organized by name and version.

## Features

- **Data Source Selection**: Easily specify which dataset to load using the `DataSourceName` enumeration.
- **Version Management**: Load a specific version of a dataset or automatically detect and use the latest available version.
- **Standardized Loading**: Loads a `DataPackage` from a versioned folder, which is expected to contain `technologies.json` and optionally `sources.json`.
- **Error Handling**: Raises `FileNotFoundError` if the specified data source or version directories cannot be found.
- **Seamless Integration**: Returns a `DataPackage` object, ready for use with other components of the `technologydata` library.

## Usage Examples

### Creating a DataAccessor

To use the `DataAccessor`, you first create an instance, specifying the `data_source_name`. You can optionally provide a `data_version`. If no version is specified, the latest available version for this dataset is used.

```python
from technologydata.parsers.data_accessor import DataAccessor

# Create an accessor for a specific version
accessor_v1 = DataAccessor(
    data_source="manual_input_usa",
    version="v1.0.0"
)

# Create an accessor that will use the latest version
accessor_latest = DataAccessor(
    data_source_name="dea_energy_storage"
)
```

### Accessing Data

Once the `DataAccessor` is instantiated, call the `access_data()` method to load the dataset. The method handles finding the correct directory and then calls `DataPackage.from_json()` on that directory.

The directory structure is expected to be: `src/technologydata/parsers/<data_source_name>/<version>/`.

#### Loading a Specific Version

If `data_version` was specified during instantiation, `access_data()` will look for that exact version.

```python
# Assuming the path .../parsers/manual_input_usa/v1.0.0/ exists
dp_v1 = accessor_v1.access_data()

# dp_v1 is now a DataPackage object containing the data from v1.0.0
print(type(dp_v1))
# <class 'technologydata.datapackage.DataPackage'>
```

#### Loading the Latest Version

If `data_version` is not provided or not found, `access_data()` automatically identifies the latest version based on semantic versioning rules (e.g., `v2.1.0` is newer than `v2.0.0`, and `v10` is newer than `v2`).

```python
# Assuming .../parsers/dea_energy_storage/ contains version directories like 'v1', 'v2'
# The accessor will load data from the latest version found (e.g., 'v2')
dp_latest = accessor_latest.access_data()

print(type(dp_latest))
# <class 'technologydata.datapackage.DataPackage'>
```

## API Reference

Please refer to the [API documentation](../api/data_accessor.md) for detailed information on the `DataAccessor` class methods and attributes.

## Limitations & Notes

-   **Directory Structure**: The `DataAccessor` expects a specific directory structure within the project: `src/technologydata/parsers/<data_source_name>/<version>/`. It will not find data located elsewhere.
-   **Version Naming**: Version directories must be prefixed with a `v` and follow a pattern that can be parsed by `packaging.version` (e.g., `v1`, `v2.0`, `v1.0.1-alpha`). Directories that do not match this pattern will be ignored when searching for the latest version.
-   **Target Data**: The class is designed to load a `DataPackage` from a folder. See the `DataPackage` documentation for more details on the expected contents of that folder (i.e., `technologies.json`).
