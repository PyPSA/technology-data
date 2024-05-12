##########################################
Release Notes
##########################################

.. Upcoming Release
.. ================

Technology-Data 0.9.0 (12 May 2024)
===================================

* update electrolyser investment costs based on latest communications (https://github.com/PyPSA/technology-data/pull/129)

* add heavy duty and shipping technology assumptions from DEA (https://github.com/PyPSA/technology-data/pull/128)

* add data for Organic Rankine Cycles (ORC) and geothermal energy (https://github.com/PyPSA/technology-data/pull/111)

* bugfix for retrieving optimistic and pessimistic value ranges from DEA (https://github.com/PyPSA/technology-data/pull/130)

* update ``mock_snakemake()`` to work with snakemake v8 (https://github.com/PyPSA/technology-data/pull/127)

* compatibility with newer pandas versions (https://github.com/PyPSA/technology-data/pull/126)

Technology-Data 0.8.1 (28 February 2024)
========================================

* adjust currency year in some DEA input data

Technology-Data 0.8.0 (19 February 2024)
=======================================

* Update currency year from 2015 to 2020. Add a currency year for each data
  input. The inflation rate is taken `Eurostat HICP - annual data (average index
  and rate of change)
  <https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/dataflow/ESTAT/prc_hicp_aind/1.0?references=descendants&detail=referencepartial&format=sdmx_2.1_generic&compressed=true>`_.

Technology-Data 0.7.0 (7 February 2024)
=======================================

* Updated to latest release of DEA renewable fuels (released January 2024). With the following changes
  * The following technologies have updated assumptions: ['BioSNG', 'BtL', 'Fischer-Tropsch', 'Haber-Bosch', 'air separation unit', 'biogas', 'biogas CC', 'biogas plus hydrogen', 'biogas upgrading', 'biomass-to-methanol', 'electrobiofuels', 'electrolysis', 'methanolisation']
  * biogas upgrading and biogas plant are differentiated in new data set between small and large plant, we assume small plant here
  * methanol from power changed to methanol from hydrogen, VOM are zero in new data set
  * CAREFUL: biogas upgrading units changed for VOM and investment from per input to per output units

* Add floating wind cost data

* Add biomass-to-methanol route from DEA.

* Add electric compression losses for hydrogen and gas pipelines from DEA.

* Add methanol-to-kerosene from `Concawe report <https://www.concawe.eu/wp-content/uploads/Rpt_22-17.pdf>`_.

* Add methanol-to-olefins/aromatics and steam cracker mostly from `DECHEMA report <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry.pdf>`_.

* Added FOM for enhanced geothermal systems.

* Added data for Organic Rankine Cycles.

* Moved efficiency for electricity generation from geothermal to ORC.

* Moved addition of geothermal data from `compile_cost_assumptions.py` to `manual_input.csv`.

* Costs for 'fuel' provided in the manual_inputs.csv are now also adjusted for inflation.

* Updated cost assumptions for 'nuclear', 'coal', and 'lignite' to Lazard's LCoE V16 (2023).

* Updated source for 'fuel' costs of 'gas', 'uranium', 'coal', and 'lignite' to DIW (2013) data.

* Updated hydrogen pipeline costs based on most recent `EHB report <https://ehb.eu/files/downloads/EHB-2023-20-Nov-FINAL-design.pdf>`_.

Technology-Data 0.6.2 (7 August 2023)
=====================================

* Use DEA electrolysis and fuel cell assumptions by default.

* Add steam generation of methanolisation process.

* Use HVDC submarine cable cost from Härtel et al. (2017).

Technology-Data 0.6.1 (4 August 2023)
===========================================

* New technologies
  - direct iron ore reduction (cost assumptions, conversion efficiencies)
  - dry bulk carrier Capesize (cost assumptions)
  - electric arc furnace (cost assumptions, conversion efficiencies)
  - shipping fuel methanol (cost assumptions, emission intensity)
  - iron ore DRI-ready (cost assumptions)

Technology-Data 0.6.0 (24 May 2023)
===========================================

* General:
  - Fix 'further_description' column from 'manual_inputs.csv' not being correctly parsed by the workflow
  - Adjust electrolysis currency year to 2015

* Updated technologies
  - updated cost assumptions for 'digestible biomass to hydrogen' and "solid biomass to hydrogen"
  - Fix: Unit for methanation investment costs is now correctly displayed as "EUR/kW_CH4" (`#82 <https://github.com/PyPSA/technology-data/issues/82#event-8638160137>`_)
  - Fix source and description for 'solar' and 'solar-rooftop' to correctly indicate how they are calculated

* New technologies
  - 17 new energy storage technologies
  - new biomass technologies ('biogas CC', 'central gas CHP CC', 'central hydrogen CHP', 'central solid biomass CHP CC', 'central solid biomass CHP powerboost CC',
'direct firing gas', 'direct firing gas CC', 'direct firing solid fuels', 'direct firing solid fuels CC', 'electrobiofuels', 'solid biomass boiler steam CC', 'waste CHP', 'waste CHP CC',
pelletizing cost for pellets from solid biomass residues)
  - "utility-scale single-axis tracking" PV (cost assumptions)
  - H2 liquefaction (Conversion efficiency)
  - CH4 liquefaction (Conversion efficiency)
  - CO2 liquefaction (Conversion efficiency)
  - LOHC hydrogenation (Conversion efficiency)
  - Ammonia crackier (Conversion efficiency)
  - Steam methane reforming (Conversion efficiency)
  - Methanol steam reforming (Conversion efficiency)
  - Fischer-Tropsch (Conversion efficiency)
  - seawater RO desalination (Conversion efficiency)
  - Haber-Bosch (Conversion efficiency)
  - air separation unit (Conversion efficiency)
  - direct air capture (Conversion efficiency)
  - Added data for Enhanced Geothermal Systems, as given by Aghahosseini, Breyer 2020

* Changed technologies
  - methanation (Conversion efficiency)
  - methanolisation (Conversion efficiency)

* Features
 - energy penalties for biomass usages, biogas and gas boilers with carbon capture (calculations will be provided in an upcoming paper
* Bug fixes
  - Fixed a bug that ommited 'further description' from manually added data

Technology-Data 0.5.0 (08 Februrary 2023)
===========================================

* Updated technologies
  - biomass CHP: changed from Wood Pellets to Wood Chips which are generally used in biomass CHPs and more expensive.
  - solid biomass fuel to 12 EUR/MWh_th based on JRC ENSPRESO datasets

* New technologies
  - new biomass technologies (BioSNG, BtL, biogas, biogas plus hydrogen, digestible biomass,digestible biomass to hydrogen, electric boiler steam, gas boiler steam, industrial heat pump high temperature, solid biomass boiler steam, solid bioass to hydrogen, biomass boiler for decentral heating
  - hydrogen storage tank type 1: Low pressure hydrogen tank storage (up to 200 bar)
  - hydrogen storage compressor: Compressor for filling hydrogen storage tanks (compression from 30 to 250 bar)
  - 18 new energy storage technologies from PNNL "Energy Storage Grand Challenge Cost and Performance Assessment 2022"

* Enable easy debugging of scripts by allowing python execution/ debugging in scripts

* Breaking changes
  - Renamed "hydrogen storage tank incl. compressor" to "hydrogen storage tank type 1 including compressor" for more clarity on the technology and consistency
  - Removed "hydrogen storage tank" technology assumption from old PyPSA assumptions which is superceeded by the "hydrogen storage tank type 1" assumptions

Technology-Data 0.4.0 (22 July 2022)
===========================================

* **WARNING**: For some technologies the units used were changed. Check for correct usage in automatic workflows.
* **WARNING**: The technology name "Haber-Bosch synthesis" was changed to "Haber-Bosch" for consistency.

* Updated technology data datasheets from DEA:
  - Industrial Process Heat (Version 11/2021)
  - Carbon Capture, Transport and Storage (Version 11/2021)
  - Renewable Fuels (Version 04/2022)

* Updated technologies (based on reviewer comments and subsequent investigation): (cf.`Pull Request #57 <https://github.com/PyPSA/technology-data/pull/57>`_)
  - Methanation:
    + Less optimistic number from report comparing multiple sources (incl. the source of the original number)
  - Fischer-Tropsch:
    + Mature technology (Hydrogen + Syngas to FTFs)
    + Account for economies of scale (previous numbers for very small installations)
    + Do not take value from DEA which is more focues on integrated Power-To-Liquid plant with low integration TRL
    + Use same value for Fischer-Tropsch and Methanolisation based on source report
    + Remove VOM for FTF, not reported in many sources and DEA numbers not reproduceable with original source
  - Methanolisation:
    + Mature technology (Hydrogen + CO2 to MeOH)
    + Account for economies of scale (previous numbers for very small installations)
    + Do not take value from DEA which is more focues on integrated Power-To-Liquid plant with low integration TRL
    + Use same value for Fischer-Tropsch and Methanolisation based on source report
  - Ammonia cracker:
    + Mixed existing/new technology with existing large plants (for different purpose)
    + Consider plant size: Higher scale up based on previously considered reference with expected economies of scale
  - H2 liquefaction:
    + Consider larger plant sizes based on recent IRENA report leading to economies of scale
    + added: lower 2050 value
    + Match plant size to other similar facility sizes (LOHC hydrogenation) in repository
  - H2 evaporation:
    + Previous value for very small-scale dispensing station
    + Consider larger plant sizes based on recent IRENA report leading to economies of scale
    + added: lower 2050 value
    + Match plant size to other similar facility sizes (LOHC dehydrogenation) in repository
  - LOHC hydrogenation:
    + Small change in investment value due to change in caluclation method
  - LOHC dehydrogenation:
    + Same calulcation method as LOHC hydrogenation applied
    + Larger facility considered with resulting economies of scale
    + Distinguishing between "LOHC dehydrogenation (small scale)" e.g. a hydrogen refueling station,
      and "LOHC dehydrogenation" for large scale applications like large scale hydrogen imports
  - Haber-Bosch:
    + Use numbers based on DEA
  - air separation unit:
    + Use numbers based on DEA from Haber-Bosch ammonia plant for consistency
  - CH4 liquefaction:
    + Fix cost, similar to issue already reported in issue #54 and PR #55
  - HVAC overhead
    + Add correct source attribution
  - HVDC overhead:
    + Add correct source attribution
  - HVDC inverter pair:
    + Add correct source attribution

Technology-Data 0.3.0 (1 October 2021)
===========================================

This release includes several new technologies (see list below), the possibility
to easily add a new technology via a manual input and an update of the H2
Electrolysis assumptions.

It is released to coincide with `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ Version 0.6.0, and is known to work with this release.

Features in more detail:

**New**:
  - update offshore wind assumptions according to DEA release in March 2022
  - update solar PV assumptions according to DEA release in Februrary 2022

* new technologies:

  - solar-rooftop residential
  - solar-rooftop commercial
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
