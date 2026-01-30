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
        monkeypatch.setattr("sys.argv", ["prog", "--num_digits", "2", "--store_source"])
        args = CommonsParser.parse_input_arguments()
        assert args.num_digits == 2
        assert args.store_source is True

    def test_default_values_when_not_provided(self, monkeypatch: typing.Any) -> None:
        """Check if parse_input_arguments works when no args are provided (so defaults are used)."""
        # no args -> defaults apply
        monkeypatch.setattr("sys.argv", ["prog"])
        args = CommonsParser.parse_input_arguments()
        assert args.num_digits == 4  # default from the function
        assert args.store_source is False
