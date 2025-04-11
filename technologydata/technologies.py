"""Base class for collection of technologies and a single technology."""

from __future__ import annotations

import logging
from pathlib import Path

import frictionless as ftl
import numpy as np
import pandas as pd

from technologydata import Source, Sources

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

    def __init__(
        self,
        sources: str | Source | Sources | list[str | Source] | dict[str | Path],
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
        self.sources = sources if isinstance(sources, Sources) else Sources(sources)

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
        self.sources = Sources(
            [
                source
                for source in self.sources.sources
                if source not in sources_without_feature
            ]
        )

        self.data = None
        self.schema = None
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

        # Load the datapackage schema to be able to validate against it
        self.schema = ftl.Schema(
            str(SPECIFICATIONS_PATH / (self.schema_name + ".schema.json"))
        )

        self.load()
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

    def load(self) -> Technologies:
        """Load the declared sources."""
        if not self.sources:
            raise ValueError("No sources to load.")

        # Load all sources using frictionless and then combine them to a single pandas DataFrame
        logger.debug(
            f"Loading {len(self.sources.sources)} sources: {', '.join(source.name for source in self.sources.sources)}"
        )
        resources = []
        for source in self.sources.sources:
            with ftl.Resource(
                path=str(source.path / (self.schema_name + ".csv")),
                schema=self.schema,
            ) as resource:
                resources.append(resource.to_pandas())

        self.data = pd.concat(resources, ignore_index=True, verify_integrity=True)

        return self

    def from_pandas(self, data: pd.DataFrame) -> Technologies:
        """Load the data from a pandas DataFrame."""
        if not isinstance(data, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame.")

        # TODO add validation of the DataFrame against the schema
        resource = ftl.Resource(
            data=data,
            schema=self.schema,
        )

        report = resource.validate()
        if report.valid:
            self.data = resource.to_pandas()
            self.sort_data()
            self.sources = {"from pandas.DataFrame": None}
        else:
            raise ValueError(
                f"Data does not match the schema, see the report for details.{report}"
            )

        return self

    def __repr__(self) -> pd.DataFrame:
        """Return a string representation of the data."""
        # TODO the __repr__ could be improved
        if self.data is None:
            return "No data loaded."
        else:
            # Return the pd.DataFrame representation
            return self.data.__repr__()

    def __getattr__(self, name: str) -> None:
        """Delegate attribute access to the `data` pandas.DataFrame if the attribute is not found in the class."""
        if self.data is not None and hasattr(self.data, name):
            return_value = getattr(self.data, name)

            # Evaluate the function and return the object with the modified data (pd.DataFrame) for operations like .query()
            # or return the return value of the function if it is not a DataFrame, e.g. for operations like .plot()
            # TODO fix: the forwarding to pandas works, but the return value is still a DataFrame, not the object
            if isinstance(return_value, pd.DataFrame):
                self.data = return_value
                return self
            else:
                return return_value
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def adjust_currency(
        self,
        currency: str,
        currency_year: int,
        source: str = "World Bank",
        parameters: list[str] | None = None,
    ) -> Technologies:
        """
        Change the currency of the data and adjust for inflation based on the region the data is from.

        Params
        -------
        currency: str
            The currency to convert to.
        currency_year: int
            The year to adjust for inflation.
        region: str
            The region to adjust for inflation.
        source: str
            Source of the data to use for inflation adjustments, supports 'World Bank' and 'International Monetary Fund'.
            'World Bank' is the default based on historic data, but is limited to the past (current year -2).
            'International Monetary Fund' includes past values and projections into the future.
        parameters: list[str]
            Parameters which are affected by the currency conversion, by default the function will only rows related to the indicators
            ['investment', 'capex', 'opex'].
        """
        # TODO implement
        raise NotImplementedError("Currency conversion is not implemented yet.")

        # import pydeflate

        # if source == "World Bank":
        #     conversion_function = pydeflate.wb_exchange
        # elif source == "International Monetary Fund":
        #     conversion_function = pydeflate.imf_exchange
        # else:
        #     raise ValueError(
        #         f"Unknown source '{source}'. Supported sources are 'World Bank' and 'International Monetary Fund'."
        #     )

        # # Default parameters
        # if parameters is None:
        #     parameters = [
        #         "investment",
        #         "capex",
        #         "opex",
        #     ]

    def adjust_scale(
        self,
        new_scale: float,
        unit: str,
        scaling_exponent: float,
        parameters: list[str] | None = None,
    ) -> None:
        """
        Adjust the scale of the data to account for economies of scale, by applying a calculated scaling factor.

        The scaling factor is calculated based on the provided scale and unit, and the scaling exponent.
            scaling_factor = (new_scale / original scale of data) ** (scaling_exponent - 1)

        where it is assumed that all values are specific values and are subject to the scaling factor accordingly.

        Parameters
        ----------
        new_scale: float
            The scale value to adjust to, used to calculate the scaling factor.
        unit: str
            The unit associated with the scale value, needed to correctly calculate the scaling factor.
        scaling_exponent: float
            The scaling exponent used to calculate the scaling factor, typical scaling exponents are 0.5 or 0.7.
        parameters: list[str]
            Parameters which are affected by the scaling, by default will affect rows of the indicators ['investment', 'capex', 'opex'].

        Example
        -------
        >>> from technologydata import Technologies
        >>> tech = Technologies()
        >>> tech.adjust_scale(new_scale=1000, unit='MW', scaling_exponent=0.7)
        >>> tech.data.head()

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

        changed_data = None
        unchanged_data = None

        changed_data = self.data.query(f"parameter in {parameters}")
        unchanged_data = self.data.loc[~self.data.index.isin(changed_data.index)]

        if changed_data["scale"].isna().any():
            logger.warning(
                "Some data entries have no associated `scale` and will yield NaN values after scaling. "
                "Set their `scale` or remove them from the collection before using this function."
            )

        # TODO create a copy with converted units that align with the requested scale
        # and use this data for calculating the scaling factors, but don't use the changed values for returning, rather the original units
        if any(changed_data["scale_unit"] != unit):
            raise NotImplementedError(
                f"Unit conversion for different units to {unit} is not implemented yet."
            )

        # Calculate the scaling factor, where we assume that all entries are specific values
        changed_data["scaling_factor"] = (new_scale / changed_data["scale"]) ** (
            scaling_exponent - 1
        )

        # Apply the scaling factor
        changed_data["value"] = changed_data["value"] * changed_data["scaling_factor"]

        # Remove obsolete scaling_factor column
        changed_data = changed_data.drop(columns=["scaling_factor"])

        # Add information about new scale + unit
        changed_data["scale"] = new_scale
        changed_data["scale_unit"] = unit

        # Recombine all data and restore the default order
        self.data = pd.concat(
            [unchanged_data, changed_data],
            ignore_index=True,
        )
        self.sort_data()

        return self

    def adjust_year(
        self,
        year: int,
        model: dict,
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
        self.data.to_csv(path=path, index=False)

