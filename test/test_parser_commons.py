# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Test the utility methods."""

import typing

from technologydata.parsers.commons import CommonsParser


class TestCommonsUtils:
    """Test suite for the Commons utility functions in the technologydata module."""

    def test_defaults_and_flag(self, monkeypatch: typing.Any) -> None:
        """Check if parse_input_arguments works when both default arg and the store_true flag are provided."""
        monkeypatch.setattr(
            "sys.argv",
            [
                "prog",
                "--num_digits",
                "2",
                "--archive_source",
                "--input_file_name",
                "file_name",
            ],
        )
        args = CommonsParser.parse_input_arguments()
        assert args.num_digits == 2
        assert args.archive_source is True
        assert args.input_file_name == "file_name"

    def test_default_values_when_not_provided(self, monkeypatch: typing.Any) -> None:
        """Check if parse_input_arguments works when no args are provided (so defaults are used)."""
        # no args -> defaults apply
        monkeypatch.setattr("sys.argv", ["prog", "--input_file_name", "file_name"])
        args = CommonsParser.parse_input_arguments()
        assert args.input_file_name == "file_name"
        assert args.num_digits == 4  # default from the function
        assert args.archive_source is False
