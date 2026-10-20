..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

.. _review:

##########################################
Data review checklist
##########################################

This guide outlines how to review a pull request that adds or changes data in the repository.
The goal is to give

* Reviewers clear instructions that are as complete as possible, and
* Contributors a clear idea of which criteria a contribution needs to fulfill.

The main thing to keep in mind is that technology-data offers a reference dataset for various energy system models, including PyPSA-Eur and PyPSA-Earth.
These models depend on the data in this repository as well as on the naming of technologies; changes in this repository can lead to breaking changes or changing results in the models that use this data.
Special care must therefore be taken for all changes. A review is about the data, where it is coming from, and making sure that others can track its provenance or understand specific details in the future; the review is less about code quality.
For every number in this repository we want to be able to say: I can find it in the source, I understand how it was converted, and I can find more information about it if I need to.

Work through the list in order:

#. **Understand the intent:** Find out which technologies are new or changed, where the data comes from, and whether existing values are corrected, updated or replaced.
#. **Only data:** Make sure that only data was changed. Code, documentation and CI changes should live in separate PRs.
#. **CI and outputs:** CI must be passing. The branch should be up to date with the default branch, or at least show no merge conflicts. In the output files ``outputs/**.csv``, only the technologies mentioned in the PR should change. Any update to existing values needs an explanation of what was wrong with the old values and why the new values are more appropriate.
#. **No unintended changes:** Make sure that only rows referenced in the PR description are changed. Make sure no extra columns or stray files are added to the repository.
#. **Traceable sources:** Each new row in the output needs to cite the author, title, year, URL or DOI of the original data source, the exact page/table/sheet, and an accessed date (YYYY-MM-DD). Check whether you can open the source and find the original number in it.
#. **Calculations written down:** The ``further_description`` column gives the original value and unit and all conversion steps taken to get to the number included in the repository, e.g. exchange rates, heating value, scaling exponent, reference size. Check the arithmetic.
#. **Units:** Whenever possible, use the same units as other technologies use for the same parameters. Use SI units when possible and avoid SI prefixes. Inputs and outputs use specific units (e.g. ``MWh_el/MWh_H2``). Always assume lower heating values (LHV, net calorific value, NCV). Currency values are in EUR.
#. **year** is the year the data refers to, not the publication year.
#. **currency_year** is the price base and must be the same across all years of an entry (identified by ``technology``, ``parameter`` and ``year``).
#. **Parameter names** are exact (``VOM``, not ``Variable O&M``); use existing parameter names where possible.
#. **Technology names** go from general to specific and follow existing naming.
#. **Sorting order** in ``inputs/manual_input.csv`` should not change; add new technologies to the end of the file, with each ``technology`` grouped and ``year`` ascending.
#. **Backward compatibility** should be kept whenever possible. Renames and removals are breaking changes that need an appropriate release note. Updates to existing values should be made with care, accounting for the potential impact on models using the data. Adding new entries is comparably low risk.
#. **Plausibility:** Check values for plausibility by comparing them with similar technologies in the existing data and with other publications. Make sure that system boundaries are comparable.
#. **Data over model assumptions:** Prefer publications about data and avoid adding data that is just an assumption of another energy system model.
#. **Trajectories:** Include several years for a technology if the source gives a trajectory, even if only historic. The workflow automatically inter- and extrapolates linearly until 2050.
#. **US data** in ``inputs/US/manual_input_usa.csv`` gets the same checks, plus its ``financial_case`` and ``scenario`` columns.
#. **Release note:** One concise sentence under *Upcoming Release* in ``docs/release_notes.rst`` that names the technologies and links the PR. People get details from the PR description and the release notes help them understand whether they are affected by a change or not.
#. **Write the review:** Sort comments into *must fix*, *should fix*, *questions* and *nits*. It is usually easiest to attach comments to specific lines using GitHub's review functionality. Be as concrete as possible and keep the tone short and friendly.

Thank you for contributing! :)