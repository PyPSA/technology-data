# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Abstract base class for a specific version of the data parser."""

import abc
import pathlib
from typing import Any


class ParserBase(abc.ABC):
    """Abstract base class for a specific version of the data parser."""

    @abc.abstractmethod
    def parse(
        self,
        input_path: pathlib.Path,
        num_digits: int,
        archive_source: bool,
        **kwargs: Any,
    ) -> None:
        """
        Parse a specific version of a dataset.

        Parameters
        ----------
        input_path : pathlib.Path
            Path to the raw input data file.
        num_digits : int
            Number of significant digits to round the values.
        archive_source : bool
            If True, archive the source object on the Wayback Machine.
        **kwargs : Any
            filter_params : bool, optional
                If True, filter the parameters stored in the output.
            export_schema : bool, optional
                If True, export the Pydantic schema for the data models.

        """
        raise NotImplementedError
