..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Install Technology-data base
==============================

Clone the repository `technology-data <https://github.com/PyPSA/technology-data>`_


.. code:: bash

    projects % git clone https://github.com/PyPSA/technology-data.git

If you want to use the technology-data with `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ or `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ the repository should be cloned in a parallel directory.


Environment/package requirements
================================

Preferred method: pixi
----------------------

We recommend using the package and environment manager `pixi <https://pixi.sh>`_ to install the python dependencies.
Follow the `pixi installation guide <https://pixi.sh/latest/installation/>`_ for your operating system.

The dependencies are declared in `pixi.toml <https://github.com/PyPSA/technology-data/blob/master/pixi.toml>`_ and locked for all supported platforms in ``pixi.lock``.
Pixi installs the environment automatically when you first run a command in it.
Either activate the environment in your shell

.. code:: bash

    .../technology-data % pixi shell

or prefix individual commands with ``pixi run``, e.g.

.. code:: bash

    .../technology-data % pixi run snakemake -j1

Legacy method: conda
--------------------

Alternatively, the environment can be installed with ``conda``.
The `envs <https://github.com/PyPSA/technology-data/tree/master/envs>`_ directory contains a pinned environment file for each platform (``linux-64``, ``osx-64``, ``osx-arm64``, ``win-64``), exported from ``pixi.lock``.
These files contain the runtime dependencies only, without test and development tools.
Create and activate the environment for your platform using

.. code:: bash

    .../technology-data % conda create -n technology-data --file envs/default_<platform>.pin.txt

    .../technology-data % conda activate technology-data

If no pinned file fits your platform, use the unpinned environment file instead:

.. code:: bash

    .../technology-data % conda env create -f envs/environment.yaml

.. note::
    Note that activation is local to the currently open shell!
    After opening a new terminal window, one needs to reissue the activation command!



Data requirements
=================

Currently all `required data <https://github.com/PyPSA/technology-data/tree/master/inputs>`_ and the corresponding `documentation <https://github.com/PyPSA/technology-data/tree/master/docu>`_ is part of the github repository and directly cloned.




Getting started
===============


In ``config.yaml`` you can control the settings for your technology data, such as the years, the expectation concerning data uncertainty (None/optimist/pessimist),
the rate of inflation, the data source for e.g. solar PV or hydrogen electrolyers and if offshore wind farm grid connection costs should be included.

To generate the technology data with your settings run:

.. code:: bash

    projects/technology-data % pixi run snakemake -j1


