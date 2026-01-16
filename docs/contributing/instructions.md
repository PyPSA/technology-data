# Instructions for contributing to the project

First of all, thank you for your contributions to `technologydata`!

We enthusiastically invite anyone interested in `technologydata` to share new ideas, provide suggestions, submit bug reports, or contribute code changes.

## How to contribute

- [Code contributions](#code-contributions): Implement new features, fix bugs, or improve the performance.
- [Documentation contributions](#documentation-contributions): Improve the documentation by adding new sections, fixing typos, or clarifying existing content.

## Where to go for help

- To **discuss** with other `technologydata` users, organise projects, share news, and get in touch with the community, please refer to the [Contacts](../home/contacts.md) page.
- For **guidelines to contribute**, stay right here.

## Code contributions

## Contribution workflow in a nutshell

- Fork the repository on GitHub. Please make sure to check the box `Copy the master branch only`
- Clone your fork: `git clone https://github.com/<your-username>/technology-data.git`
- Set up the `upstream` repository (see [Set up the upstream repository](#set-up-the-upstream-repository))
- Fetch the upstream tags `git fetch --tags upstream`
- Install with dependencies in editable mode (see [Install dependencies](#install-dependencies))
- Setup linter and formatter, e.g `pre-commit install` (see [Code linting and formatting](#code-linting-and-formatting))
- Write your code (preferably on a new branch)
- Run tests: `pytest` (see [Testing](#testing))
- Push your changes to your fork and create a pull request on GitHub

### Set up the upstream repository

To keep your fork up to date with the original repository, you need to set up the `upstream` remote repository:

```bash
git remote add upstream https://github.com/PyPSA/technology-data.git
```

To verify that the remotes have been set up correctly, you can run:

```bash
git remote -v
```

The output should look like this:

```plaintext
origin  https://github.com/<your-username>/technology-data.git (fetch)
origin  https://github.com/<your-username>/technology-data.git (push)
upstream  https://github.com/PyPSA/technology-data.git (fetch)
upstream  https://github.com/PyPSA/technology-data.git (push)
```

### Install dependencies

To contribute to this project, you will need to install the development dependencies.
These dependencies are listed in the `pyproject.toml` file and include everything needed for development, testing and building the documentation of the project.

To install the development dependencies, run the following command inside the project directory:

```bash
uv sync
```

This will create a virtual environment at `.venv`.

The development dependencies include `ipykernel` to support Jupyter notebooks and integration into e.g. VS Code.

To view the full list of development dependencies, you can check the `pyproject.toml` file under the `[dependency-groups]` section as `dev` and `docs`, which are both included in the `default-groups`.

The virtual environment can be "activated" to make its packages available:

=== "macOS and Linux"

    ```bash
    source .venv/bin/activate
    ```

=== "Windows"

    ```pwsh-session
    PS> .venv\Scripts\activate
    ```

To exit a virtual environment, use the `deactivate` command:

```bash
deactivate
```

#### Add new development dependencies

To add new dependencies to a specific dependency group in the `pyproject.toml`:

```bash
uv add package_name --group group_name
```

To install the new dependency run `uv sync`.

### Code linting and formatting

#### pre-commit

[pre-commit](https://pre-commit.com) is used to ensure a consistent code style and to catch common programming errors or bad practices before they are committed.

The tool is installed with development dependencies. It can be also installed manually via `pip install pre-commit` or `conda install -c conda-forge pre-commit`.

To use it automatically before every commit (recommended), just run once:

```bash
pre-commit install
```

This will automatically check the changes which are staged before you commit them.

To manually run it, use:

```bash
pre-commit run --all
```

This will check all files in the repository.

#### Ruff

[Ruff](https://docs.astral.sh/ruff) is used as the project linter and formatter. It combines common tools like Flake8, Black, etc.
Besides pre-commit, you can also run it via your CLI (see [Ruff installation](https://docs.astral.sh/ruff/installation/)) or IDE (e.g. VSCode [plugin](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff)).

Ruff is also already installed with the development dependencies, but you can also install it
manually using `pip install ruff`.

To use the linter in your CLI, run:

```bash
ruff check . --fix
```

This will check all files in the repository and gives you hints on what to improve. The
`--fix` flag will also automatically fix some of the issues, if possible. Some
issues need to be fixed manually.

And to run the formatter, use:

```bash
ruff format .
```

This will format all the files in the repository and immediately apply the changes to
them. It is basically [the same](https://docs.astral.sh/ruff/faq/#how-does-ruffs-formatter-compare-to-black)
as Black.

> **Note**: It is not mandatory to use either Ruff or pre-commit. We will also be running it in
> our CI/CD pipeline. But it's highly recommended, to make everyone's life easier.

### Testing

To run all tests in parallel, you can use the following command:

```bash
pytest -n auto
```

Some tests may not be thread-safe and therefore fail unexpectedly when run in parallel.
If you encounter such issues, you can run the tests sequentially by omitting the `-n auto` option:

```bash
pytest
```

## Documentation contributions

The documentation is generated with [MkDocs](https://www.mkdocs.org/). The documentation source files are written in Markdown and are available under the `/docs` sub-folder. The documentation is configured with the `mkdocs.yaml` file.

!!! note
    If you are not familiar with Markdown, consult the following [quick guide](https://www.markdownguide.org/basic-syntax/).

MkDocs offers the possibility to start a built-in development server to preview the documentation as you work on it.  To start the development server run:

### Building the documentation locally

```bash
mkdocs serve
```
