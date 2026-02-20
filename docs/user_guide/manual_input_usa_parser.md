# Manual Input USA Parser

The `ManualInputUsaParser` is responsible for parsing data from the `manual_input_usa.csv` dataset. It is designed to handle different versions of the dataset, with a specific implementation for `v0.13.4`.

## Main Parser: `ManualInputUsaParser`

The main parser, `ManualInputUsaParser`, acts as a dispatcher that selects the appropriate version-specific parser based on the user's request.

### Key Features:

-   **Version Dispatching**: Dynamically selects the correct parser for a given dataset version (e.g., `v0.13.4`).
-   **Supported Versions**: Provides a method to get a list of all supported dataset versions.

### Usage:

While it is possible to use the parser directly, the recommended way to access the data is through the `DataAccessor` class, which provides a higher-level interface and handles the parsing process internally. See the [Data Accessor](./data_accessor.md) documentation for more details.

If you need to use the parser directly, you can instantiate `ManualInputUsaParser` and call the `parse` method with the desired version and other parameters.

```python
from technologydata.parsers.manual_input_usa import ManualInputUsaParser
import pathlib

# Instantiate the main parser
manual_input_parser = ManualInputUsaParser()

# Define parameters
version = "v0.13.4"
input_file = pathlib.Path("path/to/your/manual_input_usa.csv")
num_digits = 3
archive_source = False
filter_params = True
export_schema = False

# Parse the data
manual_input_parser.parse(
    version=version,
    input_path=input_file,
    num_digits=num_digits,
    archive_source=archive_source,
    filter_params=filter_params,
    export_schema=export_schema,
)
```

## Version-Specific Parser: `ManualInputUSAV0134Parser`

The `ManualInputUSAV0134Parser` is a concrete implementation that handles version `v0.13.4` of the `manual_input_usa.csv` dataset. It inherits from `ParserBase` and contains the logic for reading, cleaning, and transforming the raw data.

### Key Responsibilities:

-   **Data Loading**: Reads data from the specified CSV file.
-   **Data Cleaning**:
  -   Handles missing values in the `scenario` column.
  -   Converts `per unit` to `%` and adjusts the corresponding values.
  -   Includes `currency_year` in the `unit` column where applicable.
-   **Data Transformation**:
  -   Extracts standardized units, carriers, and heating values from complex unit strings using the `_extract_units_carriers_heating_value` method.
-   **Object Creation**: Builds a `TechnologyCollection` from the processed data using the `_build_technology_collection` method.
-   **Output Generation**: Exports the final `TechnologyCollection` and `SourceCollection` to JSON files.

### `parse` Method

The `parse` method orchestrates the entire parsing process for the `v0.13.4` dataset.

**Parameters:**

-   `input_path` (`pathlib.Path`): Path to the raw input CSV file.
-   `num_digits` (`int`): Number of significant digits for rounding numerical values.
-   `archive_source` (`bool`): If `True`, archives the data source on the Wayback Machine.
-   `**kwargs`:
  -   `export_schema` (`bool`): If `True`, exports the Pydantic schema for the data models.

The processed data is saved to `technologies.json` and `sources.json` in the `src/technologydata/parsers/manual_input_usa/v0.13.4/` directory.

## API Reference

Please refer to the [API documentation](../api/manual_input_usa_parser.md) for detailed information on the class methods and attributes.
