# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the DataAccessor class."""

import pathlib

import pytest

from technologydata.parsers.data_accessor import DataAccessor, DataSourceName


class TestDataAccessor:
    """Test suite for the DataAccessor class in the technologydata module."""

    @pytest.mark.parametrize(
        ("versions", "expected"),
        [
            (["v1", "v2", "v3", "v4"], "v4"),
            (["v0.1.0", "v0.1.1"], "v0.1.1"),
            (["v1", "v10", "v2"], "v10"),
            (["v1.0.0", "v0.2.1", "v0.1.0"], "v1.0.0"),
        ],
    )  # type: ignore
    def test_get_latest_version_string(
        self, tmp_path: pathlib.Path, versions: list[pathlib.Path], expected: str
    ) -> None:
        """Test get_latest_version_string."""
        versions_dir = [pathlib.Path(tmp_path, version) for version in versions]
        for version in versions_dir:
            version.mkdir()

        assert DataAccessor.get_latest_version_string(versions_dir) == expected

    def test_get_latest_version_string_raises_error(
        self, tmp_path: pathlib.Path
    ) -> None:
        """Test if get_latest_version_string raises FileNotFoundError for no valid versions."""
        (tmp_path / "invalid1").mkdir()
        (tmp_path / "another_invalid").mkdir()

        path_list = [p for p in tmp_path.iterdir() if p.is_dir()]

        with pytest.raises(
            FileNotFoundError, match="No valid version directories found."
        ):
            DataAccessor.get_latest_version_string(path_list)

    def test_access_data_dea_energy_storage(self) -> None:
        """Test access_data."""
        data_accessor = DataAccessor(
            data_source_name="dea_energy_storage", data_version="v10"
        )
        data_package = data_accessor.load()

        assert data_accessor.data_source_name == DataSourceName.DEA_ENERGY_STORAGE
        assert data_accessor.data_version == "v10"
        assert data_package is not None
        assert data_package.technologies is not None
        assert data_package.sources is not None
        assert len(data_package.technologies) == 136

    def test_load_raises_value_error_for_invalid_version(self) -> None:
        """Test if load raises ValueError for an invalid version."""
        with pytest.raises(
            ValueError,
            match="Data source version 'v11' not found. The latest available version is v10.",
        ):
            DataAccessor(
                data_source_name="dea_energy_storage", data_version="v11"
            ).load()
