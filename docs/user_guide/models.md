# Growth Models Documentation

## Models

Different models can be used to modify assumptions to fit specific scenarios.

## Supported Model Types

`technologydata` currently supports the following model types:

- **Growth models**: For projecting technology parameters forward in time using mathematical models.

These are implemented as Python classes and can be used for fitting to data and making projections.

## Growth Models

Growth models are mathematical models for projecting technology parameters over time. They are implemented as Pydantic classes in `technologydata.technologies.growth_models` and can be used for both fitting to data and making projections. The following growth models are available:

- `LinearGrowth`: Linear model, $f(x) = m \cdot (x - x_0) + c$
- `ExponentialGrowth`: Exponential model, $f(x) = A \cdot \exp(k \cdot (x - x_0))$
- `LogisticGrowth`: Logistic (sigmoid) model, $f(x) = \frac{L}{1 + \exp(-k \cdot (x - x_0))}$
- `GeneralLogisticGrowth`: Generalized logistic model, $f(x) = A + \frac{K - A}{(C + Q \cdot \exp(-B \cdot (x - x_0)))^{1/\nu}}$
- `GompertzGrowth`: Gompertz model, $f(x) = A \cdot \exp(-b \cdot \exp(-k \cdot (x - x_0)))$

Each model exposes:

- A `function(x, ...)` method for the mathematical form
- A `fit()` method to fit parameters to the data points registered with the model
- A `project(to_year)` method to project to a given year. Requires all parameters of the models to be set, either manually or via fitting.
- Data points can be added via each models constructor or later via the `add_data((x, y))` function

### Model Parameters

Each model has its own parameters, e.g. for `LinearGrowth` those are `x0`, `A` and `m`.
These parameters can either be set when instantiating the model or later by setting the attributes directly.
Not all parameters need to be set, e.g. if you are planning on fitting one parameter based on the data of the model, it makes sense to only set the other parameters.
Omitting parameters can be done by either not providing them at all or setting them to `None`.

```python
from technologydata.technologies.growth_models import LinearGrowth
model = LinearGrowth(x0=2020, A=None)  # Set only x0; A and m are not set

# Set A after instantiation of the model
model.A = 5_000_000

# All principal model parameters are available as
print(model.model_parameters)
> ['x0', 'A', 'm']

# All already set model parameters are available as
print(model.provided_parameters)
> ['x0', 'A']

# All parameters that are still missing (None) are available as
print(model.missing_parameters)
> ['m']
```

The meaning of the parameters is documented in the Pydantic field descriptions:

```python
print(LinearGrowth.model_fields["A"].description)
> Starting value of the linear function.
```

and the function docstring:

```python
::: technologydata.technologies.growth_models.LinearGrowth.function
```

### Creating a model and projecting a value

To create a model, e.g. for the growth of electric vehicles over time, instantiate the model with the known parameters:

```python
from technologydata.technologies.growth_models import LinearGrowth

# start with 5 million in 2020, and grow linearly by 2 million per year
model = LinearGrowth(x0=2020, A=5_000_000, m=2_000_000)

# Project to 2030, should be 5 million + 2 million * (2030 - 2020) = 25 million
model.project(2030)
> 25000000

```

The underlying function for each model can also be inspected and called directly, e.g.:

```python
LinearGrowth().function(2030, x0=2020, A=5_000_000, m=2_000_000)
> 25000000
```

Note that when using the `function` method directly, all parameters must be explicitly provided.
(TODO: maybe function should be a classmethod or read the parameters from the instance if they are set?)

### Fitting a model to data

To fit a model to data points, provide the data points to the model via the `data_points` argument when instantiating the model or add them later via the `add_data((x, y))` method:

```python
import numpy as np
from technologydata.technologies.growth_models import LinearGrowth

# Dummy data to fit to
x0 = 2020
m = 2_000_000
A = 5_000_000
x = np.arange(2020, 2025)
y = A + m * (x - x0)

# Create the model with data points
# note that a linear growth is not a great model for this type data
model = LinearGrowth(data_points=[*zip(x, y)])

# More data points can be added later as (x, y) tuples
model.add_data((2030 - x0, A + m * (2030 - x0)))

# The model can now be fitted to the data points
model.fit()

# After fitting, the model parameters are set and can be inspected ...
print(model.x0)
> -12173.437639145683
print(model.m)
> 121803685.25127366
print(model.A)
> 121803685.25127366

# ... or used for projection
print(f"Projected: {model.project(2050)}, expected: {A + m * (2050 - x0)}")
> Projected: 8779498.59552154, expected: 65000000
```

In the above example, the model yielded a projected value and the fitted parameters were not a good what was expected compared to the original values for `x0`, `A` and `m`.
This can have multiple reasons, common problems are:

1. Insufficient or poor quality data points
2. The model being not suitable for matching the data trend
3. Parameters of the model competing against each other during fitting

In the example above, the third case applies.
The model has three parameters, `x0`, `A` and `m`, of which `x0` and `A` are competing against each other.
To get a better fit of the model to the data the model parameters can either be partially fixed such that only some of the parameters need to be fitted, or initial guesses for the parameters can be provided to guide the fitting process.

To fix parameters, simply set them when instantiating the model or later by setting the attributes directly:

```python
# Fix x0 for the model
model = LinearGrowth(x0=2020, data_points=[*zip(x, y)])

# Fix m for the model
model.m = 2_000_000

# Now only A is a free and missing parameter
print(model.missing_parameters)
> ['A']

# Fitting the model again now yields a better fit
model.fit()
print(model.A)
> 5000000.0
```

Instead of fixing parameters, initial guesses can be provided via the `p0` argument to the `fit` method.
Initial guesses are ways of providing the model with some context information about the data and the expected parameter values.
This method does not guarantee that the parameters will be fitted to the provided values, but it can help to guide the fitting process and find a better fit.
By default, if no initial guesses are provided, the fitting algorithm will start with `1.0` for parameters for which no starting value is provided.
In the example `x0` should be around the current year or slightly in the past and `m` should be a large positive number.
This knowledge is provided as initial guess to the model fitting:

```python
model = LinearGrowth(data_points=[*zip(x, y)])

# Guess x0 to be around the 2010s and m to be a large number
# nothing is known about A, so it is not provided
model.fit(p0={"x0": 2015, "m": 100_000})

# Resulting in parameters
print(model.x0, model.A, model.m)
> 2016.759002811191 -1481994.3776179291 2000000.0000000002
```

As we can see, the fit for `m` now is the expected value.
The value for `x0` is still a bit off, and because `A` and `x0` are still competing against each other, the fit for `A` is also not as expected.
The model does however yield projected values that are corresponding to the expected values from the original function:

```python
print(f"Projected: {model.project(2050)}, expected: {A + m * (2050 - x0)}")
> Projected: 64999999.99999994, expected: 65000000
```

Best results for the quality of a fitted model can be achieved by combining the approaches and partially fixing a model as well as providing initial guesses.

## Using models to create consistent datasets

The main intention for models is for them to be used on combination with a `TechnologyCollection` to create consistent scenario for parameters over time.

Given a `TechnologyCollection` with multiple entries for a technology, e.g. solar photovoltaics and one or more parameters, e.g. electricity supply and lifetime, models can be used to create consistent scenarios for this technology and parameters.

Take for example the following `TechnologyCollection` for utility scale solar photovoltaics:

```python
from technologydata import TechnologyCollection, Technology, Parameter

tc = TechnologyCollection(
    technologies=[
        Technology(
            case="historic",
            name="Solar Photovoltaics",
            detailed_technology="Any",
            region="Global",
            year=2010,
            parameters={
                "electricity supply": Parameter(magnitude=1_000, units="PJ", carrier="electricity"),
                "lifetime": Parameter(magnitude=20, units="years"),
                "efficiency": Parameter(magnitude=0.22, units="%"),
                "specific investment cost": Parameter(magnitude=600, units="USD_2024/kW"),
            },
        ),
        Technology(
            case="historic",
            name="Solar Photovoltaics",
            detailed_technology="Any",
            region="Global",
            year=2022,
            parameters={
                "electricity supply": Parameter(magnitude=6_000, units="PJ", carrier="electricity"),
                "lifetime": Parameter(magnitude=24, units="years"),
                "efficiency": Parameter(magnitude=0.20, units="%"),
                "specific investment cost": Parameter(magnitude=400, units="USD_2024/kW"),
            },
        ),
        Technology(
            case="IEA STEPS 2024",
            name="Solar Photovoltaics",
            detailed_technology="Any",
            region="Global",
            year=2030,
            parameters={
                "electricity supply": Parameter(magnitude=26_000, units="PJ", carrier="electricity"),
                "lifetime": Parameter(magnitude=30, units="years"),
                "efficiency": Parameter(magnitude=0.21, units="%"),
                "specific investment cost": Parameter(magnitude=300, units="USD_2024/kW"),
            },
        ),
    ]
)
```

By using models it is possible to create consistent scenarios for the parameters over time, e.g. for 2040 and 2050.
Using the `.project()` method of the `TechnologyCollection`, models can be specified for each parameter to be projected to specified years:

```python
from technologydata.technologies.growth_models import ExponentialGrowth, LinearGrowth
tc.project(
    to_years=[2010, 2020, 2030, 2040, 2050],
    parameters={
        "electricity supply": ExponentialGrowth(x0=2010, A=1_000),
        "lifetime": LinearGrowth(x0=2010),
        "efficiency": "mean",
        "specific investment cost": "NaN",
    },
)
```

This returns a new `TechnologyCollection` with the entries for the years 2010, 2020, 2030, 2040 and 2050 and the four parameters projected according to the specified models:

- for `electricity supply` an exponential growth model is assumed and fitted to the three data points from the original collection
- for `lifetime` a linear growth model is assumed, with the starting year fixed to 2010 and the remaining model parameters fitted to the data points
- for `efficiency` the mean of the data points is taken for all projected years, representing high quality and expensive solar panels being deployed in earlier years, then cheaper and less efficient panels being deployed in later years and then efficiency increasing again as technology improves and production methods get better
- for `specific investment cost` no model is provided and instead the value is set to `NaN` for all projected years; this carries over the `Parameter` into the projected collection, and only keeps it as a placeholder; this can be useful for cases where other methods will be used to fill the parameter values later, e.g. via an external model or manual input.

The returning `TechnologyCollection` now has new values for all specified years.
The included values do not represent the original values, but modelled values using the original values.
That is also the reason why values from the original collection are different from the new collection, e.g. for 2010:
With the provided data and the requested model (`ExponentialGrowth`), the modelled yields a value of `1426.5 PJ` for 2010 instead of the original `1000 PJ`.

An already existing and initialised model can also be provided to the `project` method.
This allows models to be reused across multiple projections, e.g. for different regions or technologies.

A model can also be fit to data points before being provided to the `project` method, allowing more control over the fitting process like providing initial guesses `p0`:

```python
from technologydata.technologies.growth_models import ExponentialGrowth
model = tc.fit(
    parameter="electricity supply",
    model=ExponentialGrowth(x0=2010),
    p0={"A": 1_000, "k": 0.2},
)
print(model.model_dump(include=model.provided_parameters))
> {'x0': 2010.0, 'A': 434.75658280114783, 'm': 565.2434171988522, 'k': 0.19058662818514885}
```

TODO: Might be helpful to allow passing `p0` directly to the `project` method or while instantiating a model
TODO: Would a method to allow keeping the original values in the `TechnologyCollection` be useful, or should this rather be done by the user merging the original and the returned collection? I believe the latter would be better design.

### Model API

- `add_data((x, y))`: Add a data point for fitting
- `fit()`: Fit model parameters to data
- `project(to_year)`: Project to a given year (parameters must be set)

See the Python docstrings for each model for parameter details.
