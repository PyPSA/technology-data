# Copilot Coding Agent Instructions for `technologydata`

<!--
SPDX-FileCopyrightText: The technology-data authors

SPDX-License-Identifier: MIT

-->

## Project Overview
- `technologydata` is a Python package for energy system modellers to screen, harmonize, and transform techno-economic input data for energy system models.
- Major use cases: data validation, harmonization, transformation to model-ready formats, and provenance/audit tracing.
- Core classes: `DataPackage`, `Technology`, `Parameter`, `Source`, and `SourceCollection` (see `technologydata/`).
- Data flows: Input data (JSON or DataFrame) → `DataPackage`/`Technology` objects → validation/transformation → output for modeling/analysis.

## Key Workflows
- **Install dependencies:** Use `uv` (see `pyproject.toml`, `uv.lock`). Example: `uv sync`.
- **Activate virtual environment:** Use `source .venv/bin/activate` before running commands.
- **Linting/formatting:** Use `ruff` (`ruff check .` and `ruff format .`).
- **Type checking:** Use `mypy`.
- **Testing:** Use `pytest` (tests in `test/`). Fixtures in `test/conftest.py`. Run: `pytest`.
- **Pre-commit hooks:** Set up with `pre-commit` using `.pre-commit-config.yaml`, use this preferably over using `mypy` and `ruff` manually.
- **Docs:** Built with `mkdocs` (see `docs/`).

## Project-Specific Patterns
- **Validation:** Instantiating `Technology` or `DataPackage` triggers schema validation. Use `.check_consistency()` and `.calculate_parameters()` for further checks and derivations.
- **Parameter Derivation:** Use `.calculate_parameters()` to fill in missing values based on rules.
- **Unit Handling:** Uses `pint` for units and `pydeflate` for currency/inflation adjustments. Custom `UnitRegistry` in `technologydata/utils/units.py`.
- **Data Provenance:** Each `Parameter` and `Technology` tracks its source and transformation history.
- **Data Input:** Prefer JSON conforming to the internal schema. See `test/test_data/` for examples.

## External Integrations
- `pydeflate` for currency/inflation
- `pint` for units
- `hdx.location.country` for country/currency codes
- `savepagenow` for web archiving

## References
- Design: `docs/design.md`
- Class diagram: `docs/class-diagram.puml`
- Example data: `test/test_data/`

---
If any conventions or workflows are unclear, please ask for clarification or check `docs/design.md` for rationale and details.
