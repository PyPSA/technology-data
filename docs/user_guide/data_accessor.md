# `DataAccessor` Class Documentation

<!--
SPDX-FileCopyrightText: technologydata contributors

SPDX-License-Identifier: MIT

-->

## Overview

The `DataAccessor` class provides a standardized way to access and load technology datasets from versioned data sources. It simplifies the process of selecting a specific data source and version, automatically handling the resolution of file paths and loading the data into a `DataPackage` object. It also provides an interface to run the parsers that generate this data.

It is designed to work with a specific directory structure where datasets are organized by name and version.

## Features

- **Data Source Selection**: Easily specify which dataset to load using the `DataSourceName` enumeration.
- **Version Management**: Load a specific version of a dataset or automatically detect and use the latest available version.
- **Standardized Loading**: Loads a `DataPackage` from a versioned folder, which is expected to contain `technologies.json` and optionally `sources.json`.
- **Remote Data Download**: Download and load technology data directly from remote URLs using the `download()` method.
- **Parsing Interface**: Provides a `parse()` method to run the appropriate parser for a given data source and version to generate the data files.
- **Error Handling**: Raises `FileNotFoundError` if the specified data source or version directories cannot be found, and `ValueError` for unsupported sources or versions.
- **Seamless Integration**: The `load()` and `download()` methods return a `DataPackage` object, ready for use with other components of the `technologydata` library.
- **Path Safety Helper**: Provides `ensure_path_exists()` to create missing directories (including parents) safely.

## Usage Examples

### Creating a DataAccessor

To use `DataAccessor`, create an instance with `data_source`. If the `version` is not provided, the `DataAccessor` will automatically determine the latest available data version and use it. Best practice is to always specify the version to ensure reproducibility, otherwise package updates may change the data returned when newer versions for this package become available.
Finally, you can optionally override `data_path`.

```python
from technologydata import DataAccessor

# Create an accessor for a specific version
accessor_v1 = DataAccessor(
    data_source="manual_input_usa",
    version="v1.0.0"
)
```

### Loading Data

Once the `DataAccessor` is instantiated, call the `load()` method to load the dataset. The method handles finding the correct directory and then calls `DataPackage.from_json()` on that directory.

The directory structure is expected to be: `src/technologydata/parsers/<data_source_name>/<version>/`.

The `load()` method will look for the exact version specified during instantiation. If the version is not provided, it will log a warning and use the latest available version. If the version is provided but not found, it will raise a `ValueError` and inform you of the latest available version.

```python
# Assuming the path .../parsers/manual_input_usa/v1.0.0/ exists
dp_v1 = accessor_v1.load()

# dp_v1 is now a DataPackage object containing the data from v1.0.0
print(type(dp_v1))
# <class 'technologydata.datapackage.DataPackage'>
```

### Downloading Data from Remote URLs

The `download()` method enables you to download and load technology data directly from a remote URL. This is useful when you want to access datasets hosted on external servers without manually downloading files.

```python
from technologydata import DataAccessor

# Create an accessor for remote data
remote_accessor = DataAccessor(
    data_source="dea_energy_storage",
    version="v10"
)

# Download and load data from a remote URL
base_url = "https://example.com/data/dea_energy_storage/v10/"
dp_remote = remote_accessor.download(base_url)

# dp_remote is now a DataPackage object containing the downloaded data
print(type(dp_remote))
# <class 'technologydata.datapackage.DataPackage'>
```

The `download()` method will:

1. Download `technologies.json` from the specified base URL
2. Attempt to download `sources.json` (optional - if not found, sources will be extracted from technologies)
3. Save the files to the configured data path
4. Load and return a `DataPackage` object

**Note**: The base URL should point to the directory containing the JSON files. The method will automatically append the file names to construct the full URLs.

#### Important: URL and Version Consistency

!!! warning "Caller Responsibility"
    The `data_source` and `version` attributes of the `DataAccessor` instance determine the **local storage location** for downloaded files. The `base_url` parameter determines **what data is actually downloaded**.

```text
**The method does not validate that the URL content matches the specified version.** It is the caller's responsibility to ensure these are consistent.
```

**What this means in practice:**

- The `data_source` and `version` attributes control where files are saved on disk
- The `base_url` controls which remote files are downloaded
- If these don't match, you'll download one version's data and label it as another version

**Example of a problematic mismatch:**

```python
# ⚠️ BAD: This will download v10 data but store it as v9!
accessor = DataAccessor(data_source="dea_energy_storage", version="v9")
base_url = "https://example.com/data/dea_energy_storage/v10/"  # URL points to v10
dp = accessor.download(base_url)  # Downloads v10 data, stores at .../v9/
```

**Correct usage:**

```python
# ✅ GOOD: Version in accessor matches version in URL
accessor = DataAccessor(data_source="dea_energy_storage", version="v10")
base_url = "https://example.com/data/dea_energy_storage/v10/"
dp = accessor.download(base_url)  # Downloads v10 data, stores at .../v10/
```

**Storage location:**

Downloaded files are saved to:

```text
{data_path}/{data_source}/{version}/technologies.json
{data_path}/{data_source}/{version}/sources.json
```

Where `data_path` defaults to `src/technologydata/parsers/`. Once downloaded, the data can be accessed later using the `load()` method without re-downloading.

### Parsing Raw Data

The `parse()` method is used to execute the data processing pipeline for a specific data source and version. It takes the raw data file as input and generates the structured `technologies.json` and `sources.json` files.

```python
from technologydata import DataAccessor

# Create an accessor for the version to be parsed
parser_accessor = DataAccessor(
    data_source="dea_energy_storage",
    version="v10"
)

# Run the parser
parser_accessor.parse(
    input_file_name="dea_energy_storage_v10.xlsx",
    num_digits=2,
    archive_source=True
)
```

## API Reference

Please refer to the [API documentation](../api/data_accessor.md) for detailed information on the `DataAccessor` class methods and attributes.

## Limitations & Notes

-   **Directory Structure**: The `DataAccessor` expects a specific directory structure within the project: `src/technologydata/parsers/<data_source_name>/<version>/ unless `data_path` is overridden.
-   **Version Naming**: Version directories must be prefixed with a `v` and follow a pattern that can be parsed by `packaging.version` (e.g., `v1`, `v2.0`, `v1.0.1-alpha`). Directories that do not match this pattern will be ignored when searching for the latest version.
-   **Target Data**: The `load()` method is designed to load a `DataPackage` from a folder. See the [DataPackage](./datapackage.md) documentation for more details on the expected contents of that folder (i.e.,`technologies.json`).
-   **Input Location for Parse**: `parse()` expects `input_file_name` under `data_path/raw/`.
-   **Download Method Validation**: The `download()` method does **not** validate that the remote URL content matches the `data_source` and `version` attributes specified in the `DataAccessor` instance. Users must ensure consistency between these parameters to avoid downloading incorrect data or mislabeling versions. See the "Important: URL and Version Consistency" section above for details.
-   **Version Required for Download**: The `version` attribute must be set when using the `download()` method. A `ValueError` will be raised if the version is `None`.
