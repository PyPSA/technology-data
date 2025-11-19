# SPDX-FileCopyrightText: The technology-data authors
#
# SPDX-License-Identifier: MIT

"""Technology class for representing a technology with parameters and transformation methods."""

from typing import Annotated, Any, Self

import pydantic

from technologydata.parameter import Parameter


class Technology(pydantic.BaseModel):
    """
    Represent a technology with region, year, and a flexible set of parameters.

    Attributes
    ----------
    name : str
        Name of the technology.
    detailed_technology : str
        More detailed technology name.
    case : str
        Case or scenario identifier.
    region : str
        Region identifier.
    year : int
        Year of the data.
    parameters : Dict[str, Parameter]
        Dictionary of parameter names to Parameter objects.

    """

    name: Annotated[str, pydantic.Field(description="Name of the technology.")]
    detailed_technology: Annotated[
        str, pydantic.Field(description="Detailed technology name.")
    ]
    case: Annotated[str, pydantic.Field(description="Case or scenario identifier.")]
    region: Annotated[str, pydantic.Field(description="Region identifier.")]
    year: Annotated[int, pydantic.Field(description="Year of the data.")]
    parameters: Annotated[
        dict[str, Parameter],
        pydantic.Field(default_factory=dict, description="Parameters."),
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

    def calculate_parameters(self, parameters: Any | None = None) -> Self:
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

    def to_currency(
        self,
        target_currency: str,
        overwrite_country: None | str = None,
        source: str = "worldbank",
    ) -> Self:
        """
        Adjust the currency of all parameters of the technology to the target currency.

        The conversion includes inflation and exchange rates based on the object's region.
        If a different country should be used for inflation adjustment, use `overwrite_country`.

        Parameters
        ----------
        target_currency : str
            The target currency (e.g., 'EUR_2020').
        overwrite_country : str, optional
            ISO 3166 alpha-3 country code to use for inflation adjustment instead of the object's region.
        source: str, optional
            The source of the inflation data, either "worldbank"/"wb" or "international_monetary_fund"/"imf".
            Defaults to "worldbank".
            Depending on the source, different years to adjust for inflation may be available.

        Returns
        -------
        Technology
            A new Technology object with all its parameters adjusted to the target currency.

        """
        country = self.region
        if overwrite_country:
            country = overwrite_country

        # Copy the Technology object
        new_tech: Self = self.model_copy(deep=True)

        # Iterate over parameters and convert their currency
        for name, param in new_tech.parameters.items():
            new_tech.parameters[name] = param.to_currency(
                target_currency=target_currency,
                country=country,
                source=source,
            )

        return new_tech

    def adjust_region(self, target_region: str) -> Self:
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

    def adjust_scale(self, scaling_factor: float) -> Self:
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
