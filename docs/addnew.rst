..
  SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>

  SPDX-License-Identifier: GPL-3.0-only

.. _addnew:

##########################################
How to add a new technology
##########################################

Add a new technology from the Danish Energy Agency (DEA)
========================================================

You can add a new technology from the Danish Energy Agency database which are
saved in ``inputs/*.xlsx``. One can add the wished technology name (e.g. "onwind")
as a key and the excel sheet name of the technology in DEA (e.g. '20 Onshore turbines')
as a value to the dictionary ``sheet_names`` in the ``compile_cost_assumptions.py``.
For using the technology data uncertainty of DEA, add to the dictionary
``uncrtnty_lookup`` the columns in the excel sheet which contain the uncertainty
estimations (e.g. ``'onwind':'J:K'``)
It is recommended to still check in the output csv if all new technology data is
compiled correctly.

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [DEA-sheet-names]
   :end-before: [DEA-sheet-names]

If you want to extract further parameters from the DEA excel sheets you can
add them to the ``parameter`` list in the function ``get_data_DEA()``.

Add a new technology from another source
========================================================

A new technology can be added to the technology-data with an additional function in the ``compile_cost_assumptions.py``.
In this script a pandas Dataframe called ``costs`` is created for every year
specified in the ``config.yaml``. ``costs`` has a MultiIndex (``['technology', 'parameter']``)
and four columns (``['value', 'unit', 'source', 'further description']``).

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [RTD-target-multiindex-df]
   :lines: 1-4

The function for the newly added technology should extend this pandas Dataframe
``costs``. This is done for example for solar PV costs in the function `add_solar_from_other(costs) <https://github.com/PyPSA/technology-data/blob/2f4da6d75f07ef9457f4070b26d690f1e5e932a5/scripts/compile_cost_assumptions.py#L396>`_

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [add-solar-from-others]
   :lines: 1-4

If a new technology is added, existing parameter names (e.g. "investment") and
units (e.g. EUR/MW) should be used. If the energy output is not distinct, clarify
e.g. EUR/MW_el or EUR/MW_th. Efficiencies for converting/using hydrogen or
methane are given for Lower Heating Values (LHV).
All parameters and corresponding units are
summarised in the table below. It is not necessary to add for one technology all
parameters, e.g. if there are only investment costs for a technology only those
can be added without specifying e.g. VOM.

.. csv-table::
   :header-rows: 1
   :file: parameter.csv

Try to make the source of the data as clear as possible (e.g. including page
number, DOI) and add the source to the `source_dict <https://github.com/PyPSA/technology-data/blob/2f4da6d75f07ef9457f4070b26d690f1e5e932a5/scripts/compile_cost_assumptions.py#L35>`_
in the ``compile_cost_assumptions.py``. Convert the units within the script (this can be done `in part (2) <https://github.com/PyPSA/technology-data/blob/2f4da6d75f07ef9457f4070b26d690f1e5e932a5/scripts/compile_cost_assumptions.py#L1016>`_ in the script)
to be as transparent as possible.

The inflation rate should be considered for the investment costs,
as it is done for example for the cost assumptions from DIW from 2010

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [unify-diw-inflation]
   :lines: 1-5

The output cost assumptions are given for different years. So either the added
cost assumptions have to be interpolated for different in the ``config.yaml``
specified years, as done e.g. with the technology data from DEA `here <https://github.com/PyPSA/technology-data/blob/2f4da6d75f07ef9457f4070b26d690f1e5e932a5/scripts/compile_cost_assumptions.py#L259>`_

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [RTD-interpolation-example]
   :lines: 1-3

or the technology data is assumed to be constant for the different years, as e.g.
done if other electrolyzer data is assumed in this `function <https://github.com/PyPSA/technology-data/blob/2f4da6d75f07ef9457f4070b26d690f1e5e932a5/scripts/compile_cost_assumptions.py#L443>`_.

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :start-after: [add-h2-from-other]
   :lines: 1-10
