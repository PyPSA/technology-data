# DEA Energy Storage Parser

The `DeaEnergyStorageParser` is responsible for parsing data from the Danish Energy Agency (DEA) Energy Storage dataset. It is designed to handle different versions of the dataset, with a specific implementation for `v10`.

## Main Parser: `DeaEnergyStorageParser`

The main parser, `DeaEnergyStorageParser`, acts as a dispatcher that selects the appropriate version-specific parser based on the user's request.

### Key Features:

-   **Version Dispatching**: Dynamically selects the correct parser for a given dataset version (e.g., `v10`).
-   **Supported Versions**: Provides a method to get a list of all supported dataset versions.

### Usage:

While it is possible to use the parser directly, the recommended way to access the data is through the `DataAccessor` class, which provides a higher-level interface and handles the parsing process internally. See the [Data Accessor](./data_accessor.md) documentation for more details.

If you need to use the parser directly, you can instantiate `DeaEnergyStorageParser` and call the `parse` method with the desired version and other parameters.

```python
from technologydata.parsers.dea_energy_storage import DeaEnergyStorageParser
import pathlib

# Instantiate the main parser
dea_parser = DeaEnergyStorageParser()

# Define parameters
version = "v10"
input_file = pathlib.Path("path/to/your/data.xlsx")
num_digits = 3
archive_source = False
filter_params = True
export_schema = False

# Parse the data
dea_parser.parse(
    version=version,
    input_path=input_file,
    num_digits=num_digits,
    archive_source=archive_source,
    filter_params=filter_params,
    export_schema=export_schema,
)
```

## Version-Specific Parser: `DeaEnergyStorageV10Parser`

The `DeaEnergyStorageV10Parser` is a concrete implementation that handles version 10 of the DEA Energy Storage dataset. It inherits from `ParserBase` and contains the logic for reading, cleaning, and transforming the raw data.

### Key Responsibilities:

-   **Data Loading**: Reads data from the specified Excel file (`alldata_flat` sheet).
-   **Data Cleaning**:
  -   Removes invalid or incomplete rows.
  -   Cleans and standardizes technology names, parameters, and units.
  -   Extracts and formats year values.
-   **Data Transformation**:
  -   Converts currency units (e.g., `MEUR_2020` to `EUR_2020`) and adjusts values accordingly.
  -   Standardizes parameter names (e.g., maps `energy storage capacity for one unit` to `capacity`).
-   **Object Creation**: Builds a `TechnologyCollection` from the processed data.
-   **Output Generation**: Exports the final `TechnologyCollection` and `SourceCollection` to JSON files.

### `parse` Method

The `parse` method orchestrates the entire parsing process for the `v10` dataset.

**Parameters:**

-   `input_path` (`pathlib.Path`): Path to the raw input Excel file.
-   `num_digits` (`int`): Number of significant digits for rounding numerical values.
-   `archive_source` (`bool`): If `True`, archives the data source on the Wayback Machine.
-   `**kwargs`:
  -   `filter_params` (`bool`): If `True`, filters parameters to a predefined allowed set.
  -   `export_schema` (`bool`): If `True`, exports the Pydantic schema for the data models.

The processed data is saved to `technologies.json` and `sources.json` in the `src/technologydata/parsers/dea_energy_storage/v10/` directory.

## API Reference

Please refer to the [API documentation](../api/dea_energy_storage_parser.md) for detailed information on the class methods and attributes.
