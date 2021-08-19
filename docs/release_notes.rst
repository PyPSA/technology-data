##########################################
Release Notes
##########################################

Technology-Data 0.3.0 (xxth August 2021)
===========================================

This release includes several new technologies (see list below), the possibility
to easily add a new technology via a manual input and an update of the H2
Electrolysis assumptions.

It is released to coincide with `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ Version 0.6.0, and is known to work with this release.

Features in more detail:

**New**:

* new technologies:

  - seawater desalination (SWRO)
  - clean water tank storage
  - industrial heat pump for medium process temperatures
  - H2 and CH4 pipelines and compressors
  - shipping of CH4 (l), NH3 (l), LOHC, MeOH and H2 (l), Fischer-Tropsch
  - H2 liquefaction and evaporation
  - LOHC liquefication, hydrogenation and dehydrogenation
  - NH3 production (Haber-Bosch synthesis and air separation unit)
  - Fischer-Tropsch synthesis
  - costs for SMR (methane and methanol) and ammonia cracking
  - home battery storage and
  - CO2 pipeline
  - costs for retrofitting CH4 pipelines to H2 pipelines
* new function to adjust the investment costs according to the inflation. This is based on in the ``config.yaml`` specified rate of inflation and considered year
* new option to allow manual input via an additional csv file ``inputs/manual_inputs.csv``
* update of the H2 electrolyser assumptions based on new DEA release
* rudimentary CI and templates for pull requests and issues
* update of the latex tables for displaying the technology data


**Bugfixes**:

* adjust battery inverter lifetime to DEA footnote
* unit consistency, typos

Technology-Data 0.2.0 (11th December 2020)
===========================================

This release allows to include uncertainty bounds from the Danish Energy Agency (DEA), fixes inconsistencies with the handling of combined heat and power plants, and includes the latest data from the DEA on carbon capture technologies.

It is released to coincide with `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ Version 0.3.0 and `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ Version 0.4.0, and is known to work with these releases.

Features in more detail:

* Using the ``expectation`` parameter in ``config.yaml`` you can control whether the upper and lower uncertainty bounds on technology parameters are read in from the DEA datasets.
* The biomass and gas combined heat and power (CHP) parameters ``c_v`` and ``c_b`` were read in assuming they were extraction plants rather than back pressure plants. The data is now corrected and they are implemented in PyPSA-Eur-Sec Version 0.4.0 with a fixed ratio of electricity to heat output. The old assumptions underestimated the heat output.
* The updated assumptions from the DEA for carbon capture technologies have been incorporated, including direct air capture and post-combustion capture on CHPs, cement kilns and other industrial facilities. These are then used in PyPSA-Eur-Sec Version 0.4.0.


Technology-Data 0.1.0 (21st August 2020)
========================================

This is the first release to coincide with the release of `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ Version 0.2.0.


Release Process
===============

* Finalise release notes at ``docs/release_notes.rst``.

* Update version number in ``docs/conf.py`` and ``config.yaml``.

* Make a ``git commit``.

* Tag a release by running ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Make a `GitHub release <https://github.com/PyPSA/pypsa-eur-sec/releases>`_, which automatically triggers archiving by `zenodo <https://doi.org/10.5281/zenodo.3994163>`_.

* Send announcement on the `PyPSA mailing list <https://groups.google.com/forum/#!forum/pypsa>`_.
