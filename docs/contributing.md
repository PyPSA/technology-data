<!--
SPDX-FileCopyrightText: The technology-data authors
SPDX-License-Identifier: MIT
-->
# Contributing

## Development dependencies

To contribute to this project, you will need to install the development dependencies.
These dependencies include everything needed for development, testing and building the documentation of the project.

To install the development dependencies, run the following command inside the project directory:

```bash
uv sync
```

The development dependencies include `ipykernel` to support Jupyter notebooks and integration into e.g. VS Code.

To view the full list of development dependencies, you can check the `pyproject.toml` file under the `[dependency-groups]` section as `dev` and `docs`, which are both included in the `default-groups`.

## Testing

To run all tests in parallel, you can use the following command:

```bash
pytest -n auto
```

Some tests may not be thread-safe and therefore fail unexpectedly when run in parallel.
If you encounter such issues, you can run the tests sequentially by omitting the `-n auto` option:

```bash
pytest
```
