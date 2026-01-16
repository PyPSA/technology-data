# `SourceCollection` Class Documentation

## Overview

The `SourceCollection` class in `technologydata` represents a collection of `Source` objects, providing tools for filtering, exporting, archiving, and data conversion. It is designed to manage multiple bibliographic or web sources, supporting reproducibility and provenance tracking for technology datasets.

## Features

- **Collection Management**: Stores and iterates over multiple `Source` objects.
- **Filtering**: Supports regex-based filtering by title and authors.
- **Archiving**: Downloads archived files for all sources using the Wayback Machine.
- **Data Export**: Converts the collection to pandas DataFrame, CSV, and JSON formats.
- **Schema Export**: Exports a JSON schema describing the collection structure.
- **Integration**: Designed for use with technology parameters and datasets.

## Usage Examples

### Creating a SourceCollection

```python
from technologydata.source import Source
from technologydata.source_collection import SourceCollection

src1 = Source(title="Source One", authors="Author A", url="http://example.com/1")
src2 = Source(title="Source Two", authors="Author B", url="http://example.com/2")
collection = SourceCollection(sources=[src1, src2])
```

### Filtering Sources

```python
from technologydata.source import Source
from technologydata.source_collection import SourceCollection

src1 = Source(title="Source One", authors="Author A", url="http://example.com/1")
src2 = Source(title="Source Two", authors="Author B", url="http://example.com/2")
collection = SourceCollection(sources=[src1, src2])

filtered = collection.get(title="Source One", authors="Author A")
print(filtered)  # SourceCollection with sources matching the patterns
```

### Exporting to CSV

```python
from technologydata.source import Source
from technologydata.source_collection import SourceCollection

src1 = Source(title="Source One", authors="Author A", url="http://example.com/1")
src2 = Source(title="Source Two", authors="Author B", url="http://example.com/2")
collection = SourceCollection(sources=[src1, src2])

collection.to_csv(path_or_buf="sources.csv")
```

### Exporting to JSON and Schema

```python
import pathlib

from technologydata.source import Source
from technologydata.source_collection import SourceCollection

src1 = Source(title="Source One", authors="Author A", url="http://example.com/1")
src2 = Source(title="Source Two", authors="Author B", url="http://example.com/2")
collection = SourceCollection(sources=[src1, src2])

collection.to_json(file_path=pathlib.Path("sources.json"))
```

### Downloading Archived Files

```python
import pathlib

from technologydata.source import Source
from technologydata.source_collection import SourceCollection

src1 = Source(title="Source One", authors="Author A", url="http://example.com/1")
src2 = Source(title="Source Two", authors="Author B", url="http://example.com/2")
collection = SourceCollection(sources=[src1, src2])

downloaded_paths = collection.retrieve_all_from_wayback(pathlib.Path("downloads/"))
print(downloaded_paths)
```

## API Reference

Please refer to the [API documentation](../api/source_collection.md) for detailed information on the `SourceCollection` class methods and attributes.

## Notes

- **Filtering**: Regex patterns are case-insensitive and applied to non-optional attributes.
- **Archiving**: Download failures return `None` for the corresponding source.
- **Export**: Default CSV export uses UTF-8 encoding and quotes all fields.
- **Schema**: JSON schema is generated automatically and includes field descriptions.
- **Type Checking**: The class uses Pydantic for validation and type enforcement.
