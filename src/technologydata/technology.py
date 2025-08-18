# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT


# TODO replaceholder

"""Technology class for representing a technology with parameters and transformation methods."""

import typing

import pydantic

from technologydata.parameter import Parameter


class Technology(pydantic.BaseModel):  # type: ignore
    """
    Represent a technology with region, year, and a flexible set of parameters.

    Attributes
    ----------
    name : str
        Name of the technology.
    region : str
        Region identifier.
    year : int
        Year of the data.
    parameters : Dict[str, Parameter]
        Dictionary of parameter names to Parameter objects.
    case : str
        Case or scenario identifier.
    detailed_technology : str
        More detailed technology name.

    """

    name: typing.Annotated[str, pydantic.Field(description="Name of the technology.")]
    region: typing.Annotated[str, pydantic.Field(description="Region identifier.")]
    year: typing.Annotated[int, pydantic.Field(description="Year of the data.")]
    parameters: typing.Annotated[
        dict[str, Parameter],
        pydantic.Field(default_factory=dict, description="Parameters."),
    ]
    case: typing.Annotated[
        str, pydantic.Field(description="Case or scenario identifier.")
    ]
    detailed_technology: typing.Annotated[
        str, pydantic.Field(description="Detailed technology name.")
    ]

    def __getitem__(self, key: str) -> Parameter:
        """
        Access a parameter by name.

        Parameters
        ----------
        key : str
            Parameter name.

        Returns
        -------
        Parameter
            The requested parameter.

        """
        return self.parameters[key]

    def __setitem__(self, key: str, value: Parameter) -> None:
        """
        Set a parameter by name.

        Parameters
        ----------
        key : str
            Parameter name.
        value : Parameter
            The parameter to set.

        """
        self.parameters[key] = value

    def check_consistency(self) -> bool:
        """
        Check for consistency and completeness of parameters.

        Returns
        -------
        bool
            True if consistent, False otherwise.

        """
        # Example: check required parameters
        required = ["specific_investment", "investment", "lifetime"]
        missing = [p for p in required if p not in self.parameters]
        return len(missing) == 0

    def calculate_parameters(
        self, parameters: typing.Any | None = None
    ) -> "Technology":
        """
        Calculate missing or derived parameters.

        Parameters
        ----------
        parameters : Optional[Any]
            List of parameter names to calculate, or "<missing>" for all missing.

        Returns
        -------
        Technology
            A new Technology object with calculated parameters.

        """
        # Placeholder: implement calculation logic as needed
        return self

    def adjust_currency(self, target_currency: str) -> "Technology":
        """
        Adjust all currency parameters to a target currency.

        Parameters
        ----------
        target_currency : str
            The target currency (e.g., 'EUR_2020').

        Returns
        -------
        Technology
            A new Technology object with adjusted currency.

        """
        # Placeholder: implement currency adjustment logic
        return self

    def adjust_region(self, target_region: str) -> "Technology":
        """
        Adjust technology parameters to match a different region.

        Parameters
        ----------
        target_region : str
            The target region.

        Returns
        -------
        Technology
            A new Technology object with adjusted region.

        """
        # Placeholder: implement region adjustment logic
        return self

    def adjust_scale(self, scaling_factor: float) -> "Technology":
        """
        Scale parameter values by a scaling factor.

        Parameters
        ----------
        scaling_factor : float
            The scaling factor to apply.

        Returns
        -------
        Technology
            A new Technology object with scaled parameters.

        """
        # Placeholder: implement scaling logic
        return self

    @classmethod
    def from_dict(cls, data: dict[str, typing.Any]) -> "Technology":
        """
        Create an instance of the class from a dictionary.

        Parameters
        ----------
        cls : type
            The class to instantiate.
        data : dict
            A dictionary containing the data to initialize the class instance.
            Expected keys include:
                - "region" (str): The region associated with the instance.
                - "case" (str): The case identifier.
                - "year" (int): The year value.
                - "name" (str): The name of the instance.
                - "detailed_technology" (str): Details about the technology.
                - "parameters" (dict): A dictionary where each key maps to a parameter data
                  dictionary, which will be converted to a Parameter object.

        Returns
        -------
        instance : cls
            An instance of the class initialized with the provided data.

        Notes
        -----
        This method processes the "parameters" field in the input data by converting each
        parameter dictionary into a Parameter object using `Parameter.from_dict()`. It then
        constructs and returns an instance of the class with all the provided attributes.

        """
        params = {}
        for key, param_data in data.get("parameters", {}).items():
            params[key] = Parameter.from_dict(param_data)
        return cls(
            region=data["region"],
            case=data["case"],
            year=data["year"],
            name=data["name"],
            detailed_technology=data["detailed_technology"],
            parameters=params,
        )
