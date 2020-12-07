.. _addnew:

##########################################
How to add a new technology
##########################################

Add a new technology from the Danish Energy Agency (DEA)
========================================================
You can add a new technology from the Danish Energy Agency database which are saved in inputs/*.xlsx . One can add the wished technology name (e.g. "onwind") as a key and the excel sheet name of the technology in DEA (e.g. '20 Onshore turbines') as a value to the dictionary ``sheet_names`` in the ``compile_cost_assumptions.py``. It is recommended to still check in the output csv if all new technology data is compiled correctly.

.. literalinclude:: ../scripts/compile_cost_assumptions.py
   :language: python
   :lines: 54-60, 94

Add a new technology from another source
========================================================
