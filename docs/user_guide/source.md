# `Source` Class Documentation

## Overview

The `Source` class in `technologydata` represents bibliographic and web sources, supporting metadata, archiving, and retrieval from the Wayback Machine. It is designed to track provenance, ensure reproducibility, and facilitate the management of references for technology parameters and datasets.

## Features

- **Bibliographic Metadata**: Stores title, authors, and optional URL, access date, archive URL, and archive date.
- **Equality and Hashing**: Implements equality and hashing for use in sets and as dictionary keys.
- **String Representation**: Provides a readable string summary of the source.
- **Wayback Machine Archiving**: Ensures URLs are archived and retrieves archive URLs and timestamps from the Wayback Machine.
- **File Retrieval**: Downloads archived files from the Wayback Machine to a specified directory.
- **Automatic File Naming**: Determines file extension and save path based on content type or URL.

## Usage Examples

### Creating a Source

```python
from technologydata.source import Source
src = Source(title="Example Source", authors="The Authors", url="http://example.com")
```

A `Source` object can also be used to archive and retrieve PDFs or other files format as Excel.

### Archiving a URL

```python
from technologydata.source import Source
src = Source(title="Example Source", authors="The Authors", url="http://example.com")

src.ensure_in_wayback()
print(src.url_archive)  # Archived URL
print(src.url_date_archive)  # Archive timestamp
```

### Downloading an Archived File

```python
import pathlib
from technologydata.source import Source
src = Source(title="Example Source", authors="The Authors", url="http://example.com")

output_path = src.retrieve_from_wayback(pathlib.Path("downloads/"))
print(output_path)  # Path to downloaded file
```

### Example usage for an Excel file

```python
import pathlib
import pandas as pd
from technologydata.source import Source

# Create a Source for an Excel file
excel_src = Source(title="Example Spreadsheet", authors="The Authors", url="http://example.com/data.xlsx")

# Ensure the URL is archived (creates an archive if missing)
excel_src.ensure_in_wayback()
print("Archive URL:", excel_src.url_archive)
print("Archive timestamp:", excel_src.url_date_archive)

# Retrieve the archived file to `downloads/`
out_path = excel_src.retrieve_from_wayback(pathlib.Path("downloads/"))
print("Downloaded to:", out_path)

# Read the downloaded Excel file with pandas
df = pd.read_excel(out_path)
print(df.head())
```

## API Reference

Please refer to the [API documentation](../api/source.md) for detailed information on the `Source` class methods and attributes.

## Notes

- **Archiving**: If the URL is not set, archiving will raise a `ValueError`.
- **File Extensions**: File extension is inferred from content type or URL; unsupported types raise a `ValueError`.
- **HTTP Errors**: Download and content type retrieval may raise `requests.exceptions.RequestException`.
- **Duplicates**: Equality and hashing are based on all attributes; sources with identical metadata are considered equal.
