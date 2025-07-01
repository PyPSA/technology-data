## General

* Description of use cases in `docs/design.md`
* Class diagramm in `docs/class-diagram.puml`

## Tools

* `ruff` for formatting and linting
* `uv` for package management with `pyproject.toml` for configuration and package metadata
* `pytest` for testing all classes and methods
    * a `tests/conftest.py` for all the fixtures
    * a `tests/<test_name>/input` and `tests/<test_name>/output` for input and output files
* `mypy` for type checking
* `pre-commit` for setting up linting, formatting and checking hooks, configured via `.pre-commit-config.yaml`
* `pydantic` for data validation
* `mkdocs` for documentation generation, configured via `docs/mkdocs.yaml`

## Packages to depend on

* `pydeflate` for currency conversions and inflation adjustments
* `pint` for unit conversions
* `savepagenow` for archiving web pages using the Wayback Machine