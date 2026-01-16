# `DataPackage` Class Documentation

<!--
SPDX-FileCopyrightText: technologydata contributors

SPDX-License-Identifier: MIT

-->

## Overview

The `DataPackage` class in `technologydata` provides a container for managing collections of `Technology` and `Source` objects, supporting batch operations and import/export utilities. It is designed to facilitate the organization, sharing, and processing of technology datasets, including provenance tracking and source management.

## Features

- **Technology Collection**: Stores a collection of `Technology` objects via the `TechnologyCollection` class.
- **Source Collection**: Stores a collection of `Source` objects via the `SourceCollection` class.
- **Batch Operations**: Supports batch export to JSON and CSV formats.
- **Source Extraction**: Automatically extracts and aggregates sources from all parameters in the technology collection.
- **Loading Utilities**: Provides methods to load a data package from JSON files.

## Usage Examples

### Creating a DataPackage

You can create a `DataPackage` by instantiating it directly or by loading from JSON files.

```python
from technologydata import DataPackage, TechnologyCollection, SourceCollection

# Create a DataPackage with existing collections
dp = DataPackage(
    technologies=TechnologyCollection(...),
    sources=SourceCollection(...),
)
```

### Loading from JSON

To load a `DataPackage` from a folder containing `technologies.json` and (optionally) `sources.json`:

```python
from technologydata import DataPackage
dp = DataPackage.from_json("path/to/data_package_folder")
```

This will automatically extract sources from the technologies if not already present.

### Exporting to JSON

Export the data package to JSON files in a specified folder:

```python
from technologydata import DataPackage, TechnologyCollection, SourceCollection

# Create a DataPackage with existing collections
dp = DataPackage(
    technologies=TechnologyCollection(...),
    sources=SourceCollection(...),
)
dp.to_json("path/to/output_folder")
```

### Exporting to CSV

Export the data package to CSV files:

```python
from technologydata import DataPackage, TechnologyCollection, SourceCollection

# Create a DataPackage with existing collections
dp = DataPackage(
    technologies=TechnologyCollection(...),
    sources=SourceCollection(...),
)

dp.to_csv("path/to/output_folder")
# Creates technologies.csv and sources.csv in the output folder
```

### Extracting Source Collection

The `sources` attribute of the `DataPackage` can be automatically populated by extracting the sources from the `TechnologyCollection`.

In this context, `extracting` means scanning the `TechnologyCollection` for all `Source` references that appear in the technology parameters, and aggregating them into a single `SourceCollection`. The extraction process yields a collection of unique sources, by removing duplicates based on all source attributes.

```python
from technologydata.datapackage import DataPackage
from technologydata.technology_collection import TechnologyCollection

# Create a DataPackage with existing collections
dp = DataPackage(
    technologies=TechnologyCollection(...),
)

# Populate dp.sources with all unique sources from the technology collection
dp.get_source_collection()
```

Extracting the source collection can be useful in scenarios such as:

- When loading a data package that does not include a `sources.json` file, to ensure that all sources referenced in the technologies are captured.
- Before exporting the data package (to `sources.json`, CSV, or for sharing) so the package includes a consistent, central catalog of sources.
- When you need to produce provenance, citation lists, or run validations that require an explicit `SourceCollection`.

## API Reference

Please refer to the [API documentation](../api/datapackage.md) for detailed information on the `DataPackage` class methods and attributes.

## Limitations & Notes

- **Error Handling**: If neither technologies nor sources are available, source extraction will raise a `ValueError`.
- **No Data Validation**: The class assumes that the underlying `TechnologyCollection` and `SourceCollection` are valid and compatible.
