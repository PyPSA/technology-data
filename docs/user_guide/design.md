# Use cases we have designed `technology-data` for

<!--
SPDX-FileCopyrightText: The technology-data authors
SPDX-License-Identifier: MIT

-->

## Design

`technologydata` is designed for energy system modellers in mind.
It intendeds to serve common use cases encountered during the design, development and execution of energy system model experiments, such as:

- Screening techno-economic inputs for energy systems
- Comparing techno-economics indicators
- Harmonizing inputs
- Modelling of input assumptions for gap filling
- Transformation of assumptions into model-ready input formats

The project was started for PyPSA-Eur and has since then been expanded to serve more models and purposes.
As such it has been leaning towards serving the purposes of PyPSA-Eur and related models.
We hope to expand it to serve the wider energy system modeling community and a such welcome suggestions and contributions.

## Target audience

The users we target are Energy System Modellers and Analysts.
We assume they are all:

- Familiar with the concept of techno-economics on how to describe technical and economic parameters of technologies for energy system models
- Familiar with the methodology of transforming techno-economic parameters for alignment and harmonization

The users differ in their experience in these fields, but are generally aware of the methodological background behind
the transformations that we make available, like inflation adjustments, currency conversions or unit conversions.

- They prefer simplicity and robustness in use over being able to customize the transformations.
- They would like to be offered a range of options to choose from, but not too many.
- They would like to be able to use the package without having to read too much documentation, but require clear documentation on the transformations that are applied.
- Data provenance and reproducibility are important to them, so they need to be able to trace data back to its source and understand all individual steps that were applied in the transformation process to the data.

The users differ in their experience with Python and programming in general, we aim to serve three main user types

### Programmers and Data Engineers

These users are:

- Well familiar with using Python and object-oriented programming languages, data processing and exchange formats.
- Interacts with the package through Python scripts, Python modules and Jupyter notebooks.

#### Energy System Modeller

These users are:

- Only basic Python programming skills.
- Interacts with the package through a Jupyter notebook or a Python script; may want to simply access and inspect the data without writing and executing code.

#### Energy Analyst

These users are:

- No programming skills or only very basic Python skills like using pandas or DataFrames.
- Interacts with the package either through a Jupyter notebook or wants to be able to use csv / Spreadsheet files for inspection and use of the data.

## Use Cases

Below follow central use cases that we want to serve with `technologydata`.

### 🧾 UC-001: Screen techno-economic inputs for validity and completeness

#### 🧑‍💻 Actor(s)
- **Primary**: Energy System Modeller, Programmer
- **Secondary**: Energy Analyst (if reviewing pre-screened data)

#### 🎯 Goal

Detect and correct inconsistencies or omissions in techno-economic input data, ensuring it adheres to the package's schema and parameter constraints.

#### 📚 Pre-conditions
- Python environment with `technologydata` installed
- Input data provided as:
  - JSON (conforming to data schema)
  - package-provided JSON files (conforming to schema)
  - user-provided e.g. through the `technologydata` class-interface or through a special pandas DataFrame
- Users are familiar with the structure of `Technology`, `Parameter`, and `Source` classes

#### 🚦 Trigger
- Instantiation of one or more `Technology` object with user-provided or pre-compiled data
- Manual invocation of validation or consistency-check method

#### 🧵 Main Flow

1. User gets/loads or defines technology input (via JSON, DataFrame, or the packaged data).
2. `Technology` object is instantiated, triggering automatic validation.
3. Schema checks enforce:

- Presence of required fields (e.g., name, parameter value, unit, source)
- Consistency between parameters (e.g., energy unit alignment)

4. User manually runs `.check_consistency()` or similar method to detect conflicting or incomplete parameters.
5. User manually runs `.calculate_parameters()` to derive

- specific missing parameters based on known rules (e.g., specific investment, EAC), or
- all missing parameters that can be derived from the existing parameters

6. The User can manually update parameters, either overwriting existing or adding missing ones.
7. Validated and completed `Technology` object is now ready for transformation or analysis.

#### 🔁 Alternate Flows

- **Invalid schema**:
  - System raises a validation exception and rejects the data.
- **Inconsistent parameters**:
  - Warnings are logged; object remains instantiable but marked as incomplete.
- **Partial data**:
  - User is able to complete missing fields to the `Technology` object

#### ✅ Post-conditions
- One or more `Technology` objects validated and completed with all parameters that can be derived, or
- Errors on schema-violations and Warning about inconsistencies are logged.

#### 🧪 Sample Input/Output

```python
from technologydata import DataPackage, Parameter, Technology, Source
dp = DataPackage.from_json("path/to/data_package") # triggers validation
techs = dp.technologies # Access the TechnologyCollection

techs.check_consistency() # Checks all technologies for consistency
techs = techs.calculate_parameters(parameters=["specific-investment", "eac"])

tech = techs[0] # Access a specific Technology object

tech.check_consistency()  # Check consistency of a single Technology object
tech = tech.calculate_parameters(parameters="<missing>")  # Calculate missing parameters

# Manually created Technology object
src = Source(title="A source", authors="example authors", url="http://example.com/source")
src2 = Source(title="Another source", authors="example authors", url="http://example.com/source2")
params: dict[str, Parameter] = {
    "efficiency": Parameter(value=0.85, unit="fraction"),
    "cost": Parameter(value=1500.0, unit="USD/kW"),
}
tech = Technology(
    name="Solar Photovoltaic",
    region="North America",
    year=2023,
    parameters=params,
    case="Best Case",
    detailed_technology="Monocrystalline Solar Panels"
)
```

#### 📊 Importance & Frequency

- Importance: High
- Usage Frequency: Frequent, core workflow entry point

#### 📌 Notes

- Consistency logic includes checks for units, and dependency constraints between parameters (e.g. one parameter may be derived from one or more other parameters).
- Schema-based validation is extensible to new parameter types and sources.

### 🧾 UC-002: Harmonize multiple input datasets

#### 🧑‍💻 Actor(s)
- **Primary**: Energy System Modeller, Energy Analyst
- **Secondary**: Data Engineer

#### 🎯 Goal

Enable the user to bring multiple techno-economic datasets to a common basis (currency, year, units) using explicit, user-invoked transformation methods, so that they can be compared or combined.

#### 📚 Pre-conditions
- Python environment with `technologydata` installed
- One or more Technology objects loaded as `DataPackage` or `TechnologyCollection` objects
- User is familiar with available transformation methods (e.g., `adjust_currency`, `adjust_scale`, `adjust_region`)

#### 🚦 Trigger
- User loads multiple datasets and wishes to harmonize them for comparison or integration

#### 🧵 Main Flow

1. User loads datasets (e.g., from JSON, CSV, or DataFrame) into separate `DataPackage` or `TechnologyCollection` objects or creates them programmatically through the package's class interface.
2. User inspects the datasets to identify differences in currency, year, units, or other conventions.
3. User applies transformation methods as needed:

- `.adjust_currency(target_currency)`
- `.adjust_scale(target_capacity, scaling_exponent)`
- `.adjust_region(target_region)`
- Unit conversions per parameter via Technology or TechnologyCollection level methods

4. User repeats or chains transformations as required or desired for each dataset.
5. User verifies harmonization by inspecting key parameters and units.
6. Harmonized datasets are now ready for comparison, merging, or further analysis.

#### 🔁 Alternate Flows

- **Unsupported transformation**: System raises an error if a requested transformation is not supported due to missing parameters in one or more of the Technology objects.
- **Partial harmonization**: User can harmonize only a subset of parameters.

#### ✅ Post-conditions
- All datasets are harmonized to the user-specified conventions, e.g. currency, currency year, units.

#### 🧪 Sample Input/Output

```python
from technologydata import DataPackage

dp1 = DataPackage.from_json("dataset1.json")
dp2 = DataPackage.from_json("dataset2.json")

dp1.technologies = dp1.technologies.adjust_currency(to_currency="EUR_2020")
dp2.technologies = dp2.technologies.adjust_currency(to_currency="EUR_2020")
dp1.technologies = dp1.technologies.adjust_scale(to_capacity=100, scaling_exponent=0.5)
dp2.technologies = dp2.technologies.adjust_scale(to_capacity=100, scaling_exponent=0.5)

dp1.technologies = dp1.technologies.adjust_region(to_region="EUR")
dp2.technologies = dp2.technologies.adjust_region(to_region="EUR")

dp1.technologies = dp1.technologies.adjust_units(parameter="specific-investment", to_unit="EUR/kW")
# ... further harmonization as needed
```

#### 📊 Importance & Frequency

- Importance: High
- Usage Frequency: Frequent, especially when integrating or comparing datasets

#### 📌 Notes

- All harmonization steps are explicit and user-driven; no automatic harmonization is performed.
- The user is responsible for the order and combination of transformations.
- Optionally, each transformation could be logged as data provenance, allowing users to trace back the steps taken and record them in e.g. a output file for documentation.

### 🧾 UC-003: Transform assumptions into model-ready formats

#### 🧑‍💻 Actor(s)
- **Primary**: Energy System Modeller, Programmer

#### 🎯 Goal

Allow the user to derive and access all model-relevant parameters (e.g., EAC, specific investment) from harmonized data, ready for direct use in energy system models such as PyPSA-Eur.

#### 📚 Pre-conditions
- One `Technology` object or multiple in a `TechnologyCollection` or `DataPackage` available
- User knows which parameters are required for the target model

#### 🚦 Trigger
- User wants to prepare data for model input, e.g., for PyPSA-Eur

#### 🧵 Main Flow

1. User ensures all required base parameters (e.g., WACC, lifetime, investment) are present and harmonized.
2. User invokes calculation methods to derive model-ready parameters:

- `.calculate_parameters(parameters="EAC")`
- `.calculate_parameters(parameters="specific-investment")`

3. System computes and adds the derived parameters to the relevant `Technology` objects.
4. User accesses the parameters directly from the `Technology` objects for export or further use.

#### 🔁 Alternate Flows

- **Missing base parameters**: System raises an error if parameters required for the calculation are missing.
- **Calculation error**: System logs the error and aborts the calculation.

#### ✅ Post-conditions
- All required model parameters are present and accessible in the `Technology` objects, ready for export or direct use.

#### 🧪 Sample Input/Output

```python
techs = dp.technologies
techs = techs.calculate_parameters(parameters=["specific-investment"])

tech = techs[0]  # Access a specific Technology object
tech.calculate_parameters(parameters="EAC")
tech["EAC"].value  # Access the calculated EAC parameter value
```

#### 📊 Importance & Frequency

- Importance: High
- Usage Frequency: Frequent, especially before running model scenarios

#### 📌 Notes

-

### 🧾 UC-004: Compare techno-economic indicators across datasets

#### 🧑‍💻 Actor(s)
- **Primary**: Energy Analyst, Energy System Modeller

#### 🎯 Goal

Enable the user to systematically compare key techno-economic parameters (e.g., CAPEX, OPEX, efficiency) across multiple harmonized datasets in a tabular format.

#### 📚 Pre-conditions
- Two or more `Technology` objects available as `TechnologyCollection` or `DataPackage` objects
- User knows which parameters and technologies to compare

#### 🚦 Trigger
- User wants to compare indicators across datasets for quality control, reporting, or analysis

#### 🧵 Main Flow

1. User selects datasets for comparison and technologies based on similar features, e.g. the same technology name.
2. User aligns datasets on specified parameters.
3. User generates a comparison table (e.g., pandas DataFrame) showing values from each dataset side by side.
4. User reviews the DataFrame, exports it for further analysis and optionally creates visualizations using external tools.

#### 🔁 Alternate Flows

- **Manually created comparison**: User can manually create a comparison table by selecting specific parameters and technologies individually through their `.values` attributes.

#### ✅ Post-conditions
- Tabular comparison of selected parameters across datasets is available for review and export.

#### 🧪 Sample Input/Output

```python
from technologydata import DataPackage

dp = DataPackage([dp1, dp2])

techs = dp.technologies

techs = techs.get(technology="Solar PV", region="EUR")
comparison_df.to_dataframe()

techs["lifetime"].values  # Access lifetime values across technologies
```

#### 📊 Importance & Frequency

- Importance: Medium
- Usage Frequency: Regular, especially for data quality checks and exploring values to be included into a model

#### 📌 Notes

- Tabular comparison is the core feature; visualizations can be build on top of the DataFrame by the user themselves.
- Optional: Outlier detection and summary statistics could be a nice feature, but are also part of `pandas` already, so we can put this into the documentation as a suggestion for the user to explore themselves.

### 🧾 UC-005: Audit data provenance and transformation trace

#### 🧑‍💻 Actor(s)
- **Primary**: Energy System Modeller
- **Secondary**: Energy Analyst

#### 🎯 Goal

Allow the user to trace the origin and transformation history of each data point, enabling transparency and reproducibility.

#### 📚 Pre-conditions
- Data loaded and/or transformed using `technologydata`

#### 🚦 Trigger
- User requests provenance or transformation history for a specific parameter, object or Collection.

#### 🧵 Main Flow

1. User selects a `parameter` of a `Technology` object.
2. The user requests the provenance information for that parameter.
3. System provides the information about:

- Original source(s) of the data (from `Source` and `SourceCollection`) of the parameter
- All transformations applied by our package (e.g., currency adjustment, scaling, calculations, including the values that were used for these transformations)

4. User can export or document the provenance trace for reporting or documentation.

#### 🔁 Alternate Flows

- **Manually changed information**: If the user made manual changes at some point, then the system only  notifies and provides trace over what was done by the package; User changes are not tracked.
- **Requesting provenance for a Technology or TechnologyCollection**: The user can access the provenance information for all parameters of a Technology or TechnologyCollection by exporting it to a DataFrame or JSON.

#### ✅ Post-conditions
- User has access to a detailed provenance and transformation trace for any data point.
- Reports or logs can be exported for documentation.

#### 🧪 Sample Input/Output

```python
print(tech["EAC"].provenance)

techs = dp.technologies
techs.to_dataframe() # Includes provenance information in the DataFrame as a dedicated column

dp.to_json("<datapackage_folder>") # Includes provenance information in the JSON output
```

#### 📊 Importance & Frequency

- Importance: Medium (for transparency and reproducibility)
- Usage Frequency: Occasional, but critical for future-us, report writing and rapport

#### 📌 Notes

- Only transformations performed by the package are tracked; user-made changes must be recorded manually.
- Provenance tracking should be automatic and cover all package-driven transformations.
