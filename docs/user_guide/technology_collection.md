# `TechnologyCollection` Class Documentation

## Overview

The `TechnologyCollection` class in `technologydata` represents a collection of `Technology` objects, providing tools for filtering, exporting, currency adjustment, model fitting, and projection. It is designed to manage multiple technology datasets, supporting reproducibility, scenario analysis, and future projections.

## Features

- **Collection Management**: Stores and iterates over multiple `Technology` objects.
- **Filtering**: Supports regex-based filtering by attributes `name`, `region`, `year`, `case`, and `detailed technology`.
- **Data Export**: Converts the collection to pandas DataFrame, CSV, and JSON formats, with schema export.
- **Currency Adjustment**: Harmonizes all technology parameters to a target currency, including inflation and exchange rates.
- **Model Fitting**: Fits growth models to technology parameters across the collection.
- **Projection**: Projects parameters to future years using growth models or statistical options.
- **Integration**: Designed for use with energy system modeling and technology parameter analysis.

## Usage Examples

### Creating a TechnologyCollection

```python
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])
```

### Filtering Technologies

```python
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])

filtered = collection.get(name="Tech", region="DEU", year=2020, case="Base", detailed_technology="Solar")
print(filtered)  # TechnologyCollection with matching technologies
```

### Exporting to CSV

```python
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])

collection.to_csv(path_or_buf="technologies.csv")
```

### Exporting to JSON and Schema

```python
import pathlib
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])

collection.to_json(file_path=pathlib.Path("technologies.json"))
```

### Currency Adjustment

```python
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])
converted = collection.to_currency("USD_2025", source="worldbank")
```

### Fitting a Growth Model

```python
from technologydata.technologies.growth_models import LinearGrowth
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])

fitted_model = collection.fit(parameter="installed capacity", model=LinearGrowth(m=0.5, A=10))
```

### Projecting Parameters

```python
from technologydata.technologies.growth_models import LinearGrowth
from technologydata.technology import Technology
from technologydata.technology_collection import TechnologyCollection

tech1 = Technology(name="Tech1", region="DEU", year=2020, case="Base", detailed_technology="Solar PV", parameters={})
tech2 = Technology(name="Tech2", region="DEU", year=2021, case="Base", detailed_technology="Wind", parameters={})
collection = TechnologyCollection(technologies=[tech1, tech2])

projected = collection.project(
    to_years=[2030, 2040],
    parameters={
        "installed capacity": LinearGrowth(m=0.5, A=10),
        "lifetime": "mean",
        "efficiency": "NaN"
    }
)
```

## API Reference

Please refer to the [API documentation](../api/technology_collection.md) for detailed information on the `TechnologyCollection` class methods and attributes.

## Notes

- **Filtering**: Regex patterns are case-insensitive and applied to non-optional attributes.
- **Export**: Default CSV export uses UTF-8 encoding and quotes all fields.
- **Schema**: JSON schema is generated automatically and includes field descriptions.
- **Type Checking**: The class uses Pydantic for validation and type enforcement.
- **Projection**: The 'closest' option for parameter projection is not yet implemented.
