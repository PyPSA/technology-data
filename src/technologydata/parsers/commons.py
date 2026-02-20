# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""Classes for Commons methods for the data parsers."""

import argparse
from typing import Annotated, Any

import pydantic
from pydantic import BaseModel, ConfigDict


class ArgumentConfig(BaseModel):
    """
    Pydantic model for defining argument configurations.

    Allows flexible configuration of command-line arguments with type checking
    and validation.
    """

    name: Annotated[str, pydantic.Field(description="Name of the argument config")]
    arg_type: Annotated[
        type | None,
        pydantic.Field(
            description="The type to which the command-line argument should be converted."
        ),
    ] = None
    default: Annotated[
        Any | None, pydantic.Field(description="Default value of the argument config")
    ] = None
    help: Annotated[
        str | None,
        pydantic.Field(description="A brief description of what the argument does."),
    ] = None
    action: Annotated[
        str | None,
        pydantic.Field(
            description="Specification of how the command-line arguments should be handled"
        ),
    ] = None
    required: Annotated[
        bool, pydantic.Field(description="Flag to check whether field is mondatory")
    ] = False

    # Allow extra fields for maximum flexibility
    model_config = ConfigDict(extra="allow")


class CommonsParser:
    """Commons methods for the data parsers."""

    @staticmethod
    @pydantic.validate_call
    def parse_input_arguments(
        additional_arguments: list[ArgumentConfig] | None = None,
        description: str = "Flexible command line argument parser",
    ) -> argparse.Namespace:
        """
        Parse command line arguments with robust configuration.

        Parameters
        ----------
        additional_arguments : Optional[List[ArgumentConfig]]
            A list of ArgumentConfig objects defining extra arguments.
        description : str
            Description for the argument parser. Defaults to a generic message.

        Returns
        -------
        argparse.Namespace
            Parsed command line arguments

        Examples
        --------
        >>> extra_args = [
        ...     ArgumentConfig(
        ...         name="--input_file",
        ...         arg_type=str,
        ...         required=True,
        ...         help="Path to input CSV file"
        ...     ),
        ...     ArgumentConfig(
        ...         name="--verbose",
        ...         action="store_true",
        ...         help="Enable verbose output"
        ...     )
        ... ]
        >>> args = CommonsParser.parse_input_arguments(additional_arguments=extra_args)

        """
        # Create parser with provided or default description
        parser = argparse.ArgumentParser(
            description=description,
            formatter_class=argparse.RawTextHelpFormatter,
        )

        # Default arguments
        default_args = [
            ArgumentConfig(
                name="--num_digits",
                arg_type=int,
                default=4,
                help="Number of significant digits to round the values.",
            ),
            ArgumentConfig(
                name="--archive_source",
                action="store_true",
                help="Archive_source, store the source object on the wayback machine. Default: false",
            ),
            ArgumentConfig(
                name="--version",
                arg_type=str,
                help="Version of the dataset to parse.",
            ),
            ArgumentConfig(
                name="--input_file_name",
                arg_type=str,
                help="Name of the dataset file to parse. Default: None",
                required=True,
            ),
        ]

        # Combine default and additional arguments
        all_arguments = default_args + (additional_arguments or [])

        # Add arguments to parser (Option 1)
        for arg_config in all_arguments:
            # Convert Pydantic model to argparse-compatible dictionary
            arg_dict = {
                k: v
                for k, v in arg_config.model_dump().items()
                if v is not None and k != "name"
            }

            if arg_dict.get("arg_type") is not None:
                arg_dict["type"] = arg_dict.pop("arg_type")

            print("arg_dict", arg_dict)

            # Add argument to parser
            parser.add_argument(arg_config.name, **arg_dict)

        # Parse arguments
        args = parser.parse_args()

        return args
