from pathlib import Path

import frictionless
import pytest

dataset_folders = sorted(
    [
        folder
        for folder in Path("technologydata/datasources").iterdir()
        if folder.is_dir() and any(folder.glob("*.csv"))
    ]
)


@pytest.mark.parametrize(
    "folder",
    dataset_folders,
    ids=[folder.name for folder in dataset_folders],
)
def test_validate_datasets(folder: Path) -> None:
    """Validate all tables in the given dataset folder against their schemas."""
    # Locations where specifications for the CSVs are stored
    specifications_dir = Path("technologydata/datasources/specification")

    for csv_path in folder.glob("*.csv"):
        schema_name = csv_path.stem + ".schema.json"
        schema_path = specifications_dir / schema_name

        if schema_path.exists():
            report = frictionless.validate(str(csv_path), schema=str(schema_path))

            assert report.valid, (
                f"Validation failed for {Path(csv_path.parents[0].name) / csv_path.name}"
            )
        else:
            raise FileNotFoundError(f"No matching schema found for {csv_path}")
