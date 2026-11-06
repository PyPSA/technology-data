..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

####################
Technology-data base
####################

.. image:: https://img.shields.io/github/v/release/pypsa/technology-data?include_prereleases
    :alt: GitHub release (latest by date including pre-releases)

.. image:: https://readthedocs.org/projects/technology-data/badge/?version=latest
    :target: https://technology-data.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/license/pypsa/technology-data
    :alt: GitHub

.. image:: https://img.shields.io/github/repo-size/pypsa/technology-data
    :alt: GitHub repo size

.. image:: https://badges.gitter.im/PyPSA/community.svg
    :target: https://gitter.im/PyPSA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Chat on Gitter


This script compiles assumptions on energy system technologies (such as costs, efficiencies, lifetimes, etc.) for chosen years (e.g. [2020, 2030, 2050]) from a variety
of sources into CSV files to be read by energy system modelling software. The merged outputs have standardized cost years, technology names, units and source information.


This project is maintained by the  `Photovoltaic Solar Energy group <https://eng.au.dk/en/research/mechanical-engineering/energy-systems-and-thermodynamics/photovoltaic-solar-energy/>`_ at the Department of Engineering at
`Aarhus University <https://international.au.dk/>`_ and the `Department of Digital Transformation in Energy Systems <https://www.ensys.tu-berlin.de>`_ at the `Technische Universität Berlin <https://www.tu.berlin>`_.


Documentation
=============

**Getting Started**

* :doc:`installation`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Started

   installation

**Structure**

* :doc:`structure`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: structure

   structure

**How to add a new technology**

* :doc:`addnew`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: How to add a new technology

   addnew

**References**

* :doc:`release_notes`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: References

   release_notes



Licence
=======

The code in technology-data is released as free software under the `GPLv3
<http://www.gnu.org/licenses/gpl-3.0.en.html>`_, see
`LICENSE <https://github.com/PyPSA/technology-data/blob/master/LICENSE.txt>`_.
However, different licenses and terms of use may apply to the various input data.
