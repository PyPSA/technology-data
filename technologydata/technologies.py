"""Base class for collection of technologies and a single technology."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import frictionless as ftl
import numpy as np
import pandas as pd

import technologydata as td

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

# TODO Is there a better way to change the behaviour of frictionless?
# By default it does not allow reading from absolute or parent-relative paths
ftl.system.trusted = True

DATASOURCES_PATH = Path(__file__).parent / "datasources"
SPECIFICATIONS_PATH = DATASOURCES_PATH / "specification"


def check_source_validity(
    source_path: Path,
    schema: str,
    report_invalid: bool = False,
) -> bool:
    """
    Check whether the provided source is valid, i.e. it exists with the properly named CSV file that matches the schema.

    If the source is invalid, it optionally provides a report containing details on why the source is invalid.

    Params
    -------
    source_path: Path
        The path to the source folder.
    schema: str
        The name of the schema to validate against, i.e. the filename without the '.schema.json' suffix.
    report: bool
        If `True`, a report is provided for invalid sources containing details on why the source is invalid.
    """
    valid = True
    report = None
    if not source_path.is_dir():
        valid = False
        report = f"Source {source_path} is not a directory."

    if valid and not (source_path / (schema + ".csv")).exists():
        valid = False
        report = f"Source {source_path} does not contain a valid CSV file."

    if valid:
        # Validate the CSV file against the schema using frictionless
        validation_report = ftl.Resource(
            path=str(source_path / (schema + ".csv")),
            schema=str(SPECIFICATIONS_PATH / (schema + ".schema.json")),
        ).validate()

        valid = validation_report.valid
        report = validation_report if not valid else None

    if not valid and report_invalid:
        logger.info(f"Source {source_path} is invalid. Report: {report}")

    return valid


# Provide a list of all available sources currently available in the package's datasources folder.
def _get_available_sources(schema: str) -> dict[str, Path]:
    """
    Determine all available sources based on the folders in datasources.

    Params
    -------
    source_type: str
        The type of source to filter by, must match the file stem, e.g. 'technologies' for all sources providing 'technologies.csv'.
    """
    # Valid sources are folders and contain a file with the same name as the source_type + '.csv' ending
    sources = [
        folder
        for folder in DATASOURCES_PATH.iterdir()
        if check_source_validity(folder, schema, report_invalid=False)
    ]

    return {folder.stem: folder for folder in sources}


class Technologies:
    """Class to hold a collection of technologies."""

    feature = "Technologies"
    schema_name = "technologies"
    # Load the datapackage schema to be able to validate against it
    schema = ftl.Schema(str(SPECIFICATIONS_PATH / (schema_name + ".schema.json")))

    def __init__(
        self,
        sources: str | td.Source | td.Sources | list[str | td.Source] | dict[str, Path],
        sort_data: bool = True,
    ) -> None:
        """
        Initialize the Technologies class.

        Params
        -------
        sources: list[str | Source] | dict[str | Path] | Sources
            The sources to load, either as a list of source names or as a dictionary with source names and paths.
            See technologydata.Source and technologydata.Sources for more detailed options.
        sort_data: bool
            Automatically sort the data after loading following the sort order defined in `self.default_sort_by`.
        """
        # Make the sources available
        self.sources = (
            sources if isinstance(sources, td.Sources) else td.Sources(sources)
        )

        # Only keep sources that provide 'Technologies' as a feature
        sources_without_feature = [
            source
            for source in self.sources.sources
            if self.feature not in source.available_features
        ]
        if sources_without_feature:
            logger.warning(
                f"The following sources do not provide the feature '{self.feature}' and will not be used: {', '.join([source.name for source in sources_without_feature])}"
            )
        self.sources = td.Sources(
            [
                source
                for source in self.sources.sources
                if source not in sources_without_feature
            ]
        )

        self.data = pd.DataFrame()
        self.default_sort_by = [
            "source",
            "technology",
            "detailed_technology",
            "case",
            "region",
            "parameter",
            "unit",
            "year",
            "value",
        ]

        self.data = self._load()
        # Optional automatic sorting of data after loading
        if sort_data:
            self.sort_data()

    def sort_data(self) -> Technologies:
        """Sort the data by the default sort order as defined in `self.default_sort_by`."""
        self.data = self.data.sort_values(
            self.default_sort_by,
            ignore_index=True,
        )
        return self

    def _load(self) -> pd.DataFrame:
        """Load the declared sources."""
        if not self.sources:
            raise ValueError("No sources to load.")

        # Load all sources using frictionless and then combine them to a single pandas DataFrame
        logger.debug(
            f"Loading {len(self.sources.sources)} sources: {', '.join(source.name for source in self.sources.sources)}"
        )
        resources = []
        for source in self.sources.sources:
            if source.path is None:
                logger.error(f"Source path is None for source: {source.name}")
                raise ValueError(f"Source path is None for source: {source.name}")

            with ftl.Resource(
                path=str(source.path / (self.schema_name + ".csv")),
                schema=self.schema,
            ) as resource:
                resources.append(resource.to_pandas())

        loaded_data = pd.concat(resources, ignore_index=True, verify_integrity=True)

        return loaded_data

    def __repr__(self) -> pd.DataFrame:
        """Return a string representation of the data."""
        # TODO the __repr__ could be improved
        if self.data.empty:
            return "No data loaded."
        else:
            # Return the pd.DataFrame representation
            return self.data.__repr__()

    def __getattr__(self, name: str) -> Any:
        """Delegate attribute access to the `data` pandas.DataFrame if the attribute is not found in the class."""
        if not self.data.empty and hasattr(self.data, name):
            return_value = getattr(self.data, name)

            # Evaluate the function and return the object with the modified data (pd.DataFrame) for operations like .query()
            # or return the return value of the function if it is not a DataFrame, e.g. for operations like .plot()
            # TODO fix: the forwarding to pandas works, but the return value is still a DataFrame, not the object
            if isinstance(return_value, pd.DataFrame):
                self.data = return_value
            else:
                return return_value
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def adjust_currency(
        self,
        to_currency: str,
        source: str = "World Bank",
    ) -> None:
        """
        Adjust currency values in the instance's data based on the specified parameters.

        This method modifies the instance's data by converting and/or adjusting currency values
        to a target currency and base year using a specified deflation source. It determines
        the appropriate conversion function based on the provided source and applies the
        adjustment to the data.

        Parameters
        ----------
        to_currency : str
            The target currency units matching the format `<CURRENCY_CODE>_<YEAR>`,
            with the currency code specified with the ISO 4217 standard (e.g. 'USD', 'EUR').

        source : str, optional
            The source of the deflation data to be used for the adjustment. Supported values
            are "World Bank", "World Bank (Linked)", and "International Monetary Fund".
            The default is "World Bank".

        Raises
        ------
        ValueError
            If the specified source is not recognized. Supported sources are
            'World Bank', 'World Bank (Linked)', and 'International Monetary Fund'.

        Notes
        -----
        This method directly modifies the instance's `data` attribute, which is expected
        to be a DataFrame containing currency values that need to be adjusted.

        Examples
        --------
        >>> techs.adjust_currency("USD-2022", "International Monetary Fund")

        """
        match = re.match(td.Currencies.CURRENCY_UNIT_DEFAULT_FORMAT, to_currency)
        if match:
            currency = match.group(1)
            base_year = int(match.group(2))
            self.data = td.Currencies.adjust_currency(
                base_year, currency, self.data, source
            )
        else:
            raise ValueError(
                "The target currency unit does not match the requested format `<3-letter currency code>_<year>`."
            )

    def adjust_scale(
        self,
        to_capacity: float,
        to_capacity_unit: str,
        scaling_exponent: float,
        scaled_parameters: list[str] | None = None,
        absolute_parameters: list[str] | None = None,
    ) -> None:
        """
        Adjust parameters for all technologies to account for economies of scale using a scaling exponent.

        A scaling factor is calculated with each technology's capacity and unit, and the scaling exponent:

            scaling factor = (to_capacity / original capacity) ** (scaling_exponent - 1)

        The scaling factor is applied to all parameters specified in by `scaled_parameters`.
        The parameters affected by the scaling factor are the ones that are typically related to specific values, such as 'specific-investment' or 'specific-capex';
        these parameters are defined through the `scaling_parameters` argument.

        Some parameters are always scale with a scaling exponent of 1, such as 'capacity' or 'investment';
        these parameters are defined through the `absolute_parameters` argument.

        Parameters
        ----------
        to_capacity: float
            The new capacity value to adjust to.
        to_capacity_unit: str
            The unit associated with the capacity value, will be used to calculate the scaling factor.
        scaling_exponent: float
            The scaling exponent used to calculate the scaling factor, typical scaling exponents are 0.5 or 0.7.
        scaled_parameters: list[str] | None
            Parameters that are affected by the scaling factor. Defaults to ['specific-investment', 'specific-capex'].
        absolute_parameters: list[str] | None
            Parameters that are affected by the absolute capacity value changes, i.e. a scaling exponent of 1. Defaults to ['capacity', 'investment']

        Example
        -------
        >>> from technologydata import Technologies
        >>> tech = Technologies()
        >>> tech.adjust_scale(to_capacity=1000, to_capacity_unit='MW', scaling_exponent=0.7)
        >>> tech.data.head()

        """
        # Default parameters that are affected
        absolute_parameters = absolute_parameters or [
            "capacity",
            "investment",
        ]
        scaled_parameters = scaled_parameters or [
            "specific-investment",
            "specific-capex",
        ]

        # Find all unique groups / technologies based on these cols (excluding: parameter, year, value, unit, comment)
        group_cols = ["source", "technology", "detailed_technology", "case", "region"]
        grouped = self.data.groupby(group_cols, observed=True)
        updated_rows = []

        # Iterate over each technology and rescale if possible
        for group_keys, group_df in grouped:
            # Find the capacity row for this group
            cap_row = group_df.loc[(group_df["parameter"] == "capacity")]
            if cap_row.empty:
                logger.warning(
                    f"No row with parameter='capacity' found for technology {group_keys}. Skipping scaling for this technology."
                )
            elif len(cap_row) == 1:
                if cap_row.iloc[0]["unit"] != to_capacity_unit:
                    # TODO consider to_capacity_unit conversion here, if needed, and pass over technology if the units are not compatible
                    raise NotImplementedError(
                        f"Unit conversion from {cap_row.iloc[0]['unit']} to {to_capacity_unit} is not implemented yet."
                    )

                original_capacity = cap_row.iloc[0]["value"]

                # Update parameters that are scaled by the absolute change
                scaling_base = to_capacity / original_capacity
                mask = group_df["parameter"].isin(absolute_parameters)
                group_df.loc[mask, "value"] *= scaling_base

                # Update parameters that are scaled by the scaling exponent
                scaling_factor = scaling_base ** (scaling_exponent - 1)
                mask = group_df["parameter"].isin(scaled_parameters)
                group_df.loc[mask, "value"] *= scaling_factor

            else:
                raise ValueError(
                    f"Multiple rows with parameter='capacity' found for technology {group_keys}. Skipping scaling for this technology."
                )

            # technology processed and rows updated, append to the list
            updated_rows.append(group_df)

        self.data = pd.concat(updated_rows, ignore_index=True)
        self.sort_data()

    def adjust_year(
        self,
        year: int,
        model: dict[str, Any],
        parameters: list[str] | None = None,
    ) -> Technologies:
        """
        Create a forecast of the data into the future using a forecasting model.

        Forecasting models can be simple interpolation or extrapolation,
        or more complex models such as learning curves or s-curves.

        Parameters
        ----------
        year: int
            The year to adjust to.
        model: dict
            The model and details to use for the adjustment. For example:
            - 'learning curve': Includes details like the learning rate, starting and target capacities.
            - 'linear': Specifies linear interpolation or extrapolation.
        parameters: list[str]
            Parameters which are affected by the scaling, by default the function will only rows related to the indicators
            ['investment', 'capex', 'opex'].

        Example
        -------
        >>> from technologydata import Technologies
        >>> tech = Technologies()
        >>> tech.adjust_year(year=2030, model={'learning curve': {'learning_rate': 0.1, 'annual_growth': 0.05}})
        >>> tech.data.head()

        >>> tech.adjust_year(year=2030, model={"method": "linear"})

        """
        # Default parameters
        parameters = (
            parameters
            if parameters
            else [
                "investment",
                "capex",
                "opex",
            ]
        )

        # Filling of NaN values within this group of data
        grouping_columns = [
            "source",
            "technology",
            "detailed_technology",
            "case",
            "region",
            "parameter",
            "unit",
        ]
        # i.e. these are x and f(x), where we fill missing values of f(x) at a requested x
        x_col = "year"
        y_col = "value"

        # Each model uses the underlying data differently and must indicate
        # the data that is affected and the data that is not affected by the model
        unchanged_data = None
        changed_data = None

        # TODO models could be outsourced into separate classes, KISS for now for prototype
        if model["method"] == "linear":
            # TODO instead of the linear model, we can pass kwargs to pandas.interpolate
            # to use different interpolation methods, but KISS for now
            # also challenging to do extrapolation with default pandas because of bug here: https://github.com/pandas-dev/pandas/issues/31949
            model["method"] = "index"  # different name in pandas

            # Linear interpolation uses all data points for each parameter
            unchanged_data = self.data.query(f"parameter not in {parameters}")
            changed_data = self.data.query(f"parameter in {parameters}")

            # Add a duplicate entry for the target year for each group with NaN value that will be filled
            nan_data = changed_data.assign(
                **{x_col: year, y_col: np.nan},
            ).drop_duplicates(
                subset=grouping_columns + [x_col],
                keep="last",
            )
            changed_data = pd.concat(
                [changed_data, nan_data],
                ignore_index=True,
            )

            # Sorting by x_col required for interpolation to work properly
            changed_data = changed_data.sort_values(
                grouping_columns + [x_col],
                ignore_index=True,
            )

            # Interpolate group-wise the missing values using the 'year' index and scipy interpolation
            groups = []
            for label, group in changed_data.groupby(
                by=grouping_columns, observed=True
            ):
                group[y_col] = (
                    group.set_index(x_col)[y_col].interpolate(**model).to_numpy()
                )
                groups.append(group)

            # Recombine the individual groups
            changed_data = pd.concat(groups, ignore_index=True)

        elif model["method"] == "learning curve":
            # Learning curve model uses for each parameter only the data point that is closest and earlier than the target year,
            # we get this value by first restricting to the eligible parameters and years, then choosing the maximum year from
            # the remaining data for each group
            changed_data = self.data.loc[
                self.data.query(f"parameter in {parameters} and year < {year}")
                .groupby(grouping_columns, observed=True)[x_col]
                .idxmax()
            ]
            unchanged_data = self.data.loc[~self.data.index.isin(changed_data.index)]

            # TODO implement learning curve model
            raise NotImplementedError("Learning curve model is not implemented yet.")
        else:
            # TODO implement other models
            raise NotImplementedError(
                f"Model {model['method']} not implemented yet. Supported models are 'linear' and 'learning curve'."
            )

        # Recombine all data and restore the default order
        self.data = pd.concat(
            [unchanged_data, changed_data],
            ignore_index=True,
        )
        self.sort_data()

        return self

    def to_excel(self, path: str | Path) -> None:
        """
        Export the data and sources to an Excel file with sheet names based on the class names.

        Params
        -------
        path : str|Path
            The path to save the Excel file to.
        """
        with pd.ExcelWriter(path=path) as writer:
            self.data.to_excel(writer, sheet_name=type(self).__name__)
            self.sources.details.to_excel(
                writer, sheet_name=type(self.sources).__name__
            )

    def to_csv(self, path: str | Path) -> None:
        """
        Export the data to a CSV file, does not include the sources; for source exports, use Sources.to_csv(...).

        Params
        -------
        path : str|Path
            The path to save the CSV file to.
        """
        self.data.to_csv(path, index=False)

    def to_datapackage(self, path: str | Path, overwrite: bool = False) -> None:
        """
        Export the data to a folder following the datapackage specification.

        Params
        -------
        path : str|Path
            The path to save the data including sources and schema to.
            Must be a non-existing or empty folder.
            To overwrite existing files in a folder, use the `overwrite` parameter.
        overwrite : bool
            Existing files with the same name in the target path will be overwritten, default is False.
        """
        path = Path(path)

        # Check if the path exists and is empty
        if path.exists():
            if not path.is_dir():
                raise ValueError(f"Path {path} is not a directory.")
            if len(list(path.iterdir())) > 0 and not overwrite:
                raise ValueError(
                    f"Path {path} is not empty. To overwrite existing files, set `overwrite=True`."
                )

        # Safe to write beyond this point

        # Write the object data and its schema
        ftl.Resource(self.data).write(path=str(path / f"{self.schema_name}.csv"))
        self.schema.to_json(path=str(path / f"{self.schema_name}.schema.json"))

        # Write the sources and its schema
        self.sources.to_datapackage(
            path=str(path),
            overwrite=overwrite,
        )
