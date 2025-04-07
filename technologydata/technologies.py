"""Base class for collection of technologies and a single technology."""

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
        """
        self.data = None
        self.sources = {}

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

    def load(self) -> None:
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
                schema=str(SPECIFICATIONS_PATH / (self.SCHEMA + ".schema.json")),
            ) as resource:
                resources.append(resource.to_pandas())

        self.data = pd.concat(resources, ignore_index=True, verify_integrity=True)

        return self

    def __getattr__(self, name):
        """Delegate attribute access to the `data` pandas.DataFrame if the attribute is not found in the class."""
        if self.data is not None and hasattr(self.data, name):
            return getattr(self.data, name)
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )
