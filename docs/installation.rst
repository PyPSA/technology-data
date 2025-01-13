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

We recommend using the package manager and environment management system ``conda`` to install python dependencies.
Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are created in the `environment.yaml <https://github.com/PyPSA/technology-data/blob/master/environment.yaml>`_ file.
The environment can be installed and activated using

.. code:: bash

    .../technology-data % conda env create -f environment.yaml

    .../technology-data % conda activate technology-data

.. note::
    Note that activation is local to the currently open shell!
    After opening a new terminal window, one needs to reissue the second command!

.. note::
    If you have troubles with a slow ``conda`` installation, we recommend to install
    `mamba <https://github.com/QuantStack/mamba>`_ as a fast drop-in replacement via

    .. code:: bash
        
        conda install -c conda-forge mamba

    and then install the environment with

    .. code:: bash

        mamba env create -f environment.yaml



Data requirements
=================

Currently all `required data <https://github.com/PyPSA/technology-data/tree/master/inputs>`_ and the corresponding `documentation <https://github.com/PyPSA/technology-data/tree/master/docu>`_ is part of the github repository and directly cloned.




Getting started
===============


In ``config.yaml`` you can control the settings for your technology data, such as the years, the expectation concerning data uncertainty (None/optimist/pessimist),
the rate of inflation, the data source for e.g. solar PV or hydrogen electrolyers and if offshore wind farm grid connection costs should be included.

To generate the technology data with your settings run:

.. code:: bash

    projects/technology-data % snakemake -j1


