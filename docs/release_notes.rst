##########################################
Release Notes
##########################################

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
