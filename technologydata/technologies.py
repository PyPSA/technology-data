"""Base class for collection of technologies and a single technology."""

from __future__ import annotations

import logging
from pathlib import Path

import frictionless as ftl
import pandas as pd

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
def _get_available_sources(schema: str) -> list:
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


AVAILABLE_SOURCES = _get_available_sources("technologies")


class Technologies:
    """Class to hold a collection of technologies."""

    SCHEMA = "technologies"

    def __init__(
        self,
        packaged_sources: list[str] | str = "all",
        additional_sources: dict[str, Path] = None,
        load: bool = True,
        sort_data: bool = True,
    ) -> None:
        """
        Initialize the Technologies class.

        Params
        -------
        packaged_sources: list[str] | str
            A list of sources to include. If 'all', all sources are included. Otherwise, a list of source names.
        additional_sources: dict[str, Path]
            A dictionary of additional sources to include. The keys are the source names and the values are the paths to the sources.
        load: bool
            If `True`, the sources are loaded immediately. Otherwise, they can be loaded later using the `load` method.
        sort_data: bool
            Automatically sort the data after loading following the sort order defined in `self.default_sort_by`.
        """
        self.data = None
        self.schema = None
        self.sources = {}
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
            str(SPECIFICATIONS_PATH / (self.SCHEMA + ".schema.json"))
        )

        packaged_sources = (
            [packaged_sources]
            if isinstance(packaged_sources, str)
            else packaged_sources
        )
        for source in packaged_sources:
            if source == "all":
                add_source = AVAILABLE_SOURCES
            elif source in AVAILABLE_SOURCES.keys():
                add_source = {source: AVAILABLE_SOURCES[source]}
            else:
                raise ValueError(
                    f"Unknown packaged source: {source}. Check `AVAILABLE_SOURCES` for valid sources."
                )

            self.sources.update(add_source)

        # Add additional sources if provided
        if additional_sources:
            for source, path in additional_sources.items():
                self.add_source(source, path)

        # Load the data automatically if requested
        if load:
            self.load()
            # When requested, also sort the data
            if sort_data:
                self.sort_data()

    def add_source(self, name: str, path: Path, overwrite: bool = True) -> None:
        """
        Add another source to the collection from a local folder.

        If another source with the same name already exists, it will be replaced.

        Params
        -------
        name: str
            The name of the source to add.
        path: Path
            The path to the source folder.
        overwrite: bool
            Sources with the same name will be replaced if `True`, otherwise an error is raised.
        """
        if name in self.sources and not overwrite:
            raise ValueError(
                f"Source {name} already exists. Use `overwrite=True` to replace it."
            )

        if not check_source_validity(path, self.SCHEMA):
            raise ValueError(
                f"Source {path} is invalid. For details on why the source is invalid, use `check_source_validity` with `report_invalid=True`."
            )

        self.sources[name] = path
        logger.debug(f"Added source {name} to the collection.")

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
            f"Loading {len(self.sources)} sources: {', '.join(self.sources.keys())}"
        )
        resources = []
        for source in self.sources.values():
            with ftl.Resource(
                path=str(source / (self.SCHEMA + ".csv")),
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

    def __getattr__(self, name) -> None:
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
