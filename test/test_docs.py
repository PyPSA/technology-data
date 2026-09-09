# SPDX-FileCopyrightText: technologydata contributors
#
# SPDX-License-Identifier: MIT

"""
Run the Python code blocks in the documentation as doctests.

Every fenced ` ``` py ` block in a documented page is written in doctest
syntax (`>>>`/`...` prompts, expected output below). All blocks in a page
are concatenated and run in one shared namespace, in order, so later
blocks can reuse names defined by earlier ones -- matching how a reader
would paste them into a single session. This mirrors the approach used
by PyPSA (see PyPSA/PyPSA's test/test_docs.py).
"""

import doctest
import re
from pathlib import Path

import pytest

DOCS = [Path(__file__).parent.parent / "docs" / "tutorial" / "index.md"]


@pytest.mark.parametrize("fpath", DOCS, ids=str)  # type: ignore
def test_doctest_docs(fpath: Path, test_docs_flag: bool) -> None:
    """Test Python code blocks in a documentation page using doctest."""
    if not test_docs_flag:
        pytest.skip("Need --test-docs option to run documentation tests")

    content = fpath.read_text()

    # Extract ` ``` py ` fenced code blocks.
    python_blocks = re.findall(r"``` py\n(.*?)\n```", content, re.DOTALL)
    assert python_blocks, f"No doctest-style ` ``` py ` blocks found in {fpath}"

    # Combine all blocks into one docstring-like content so later blocks
    # can reuse names defined by earlier ones.
    combined_content = "\n\n".join(python_blocks)

    parser = doctest.DocTestParser()
    test = parser.get_doctest(
        combined_content, globs={}, name=str(fpath), filename=str(fpath), lineno=0
    )

    runner = doctest.DocTestRunner(optionflags=doctest.NORMALIZE_WHITESPACE)
    failures = runner.run(test).failed

    assert failures == 0, f"{failures} doctest(s) failed in {fpath}"
