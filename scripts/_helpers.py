# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import os
import re
from pathlib import Path

import snakemake as sm
from snakemake.script import Snakemake

dea_sheet_names = {
    "onwind": "20 Onshore turbines",
    "offwind": "21 Offshore turbines",
    "solar-utility": "22 Utility-scale PV",
    "solar-utility single-axis tracking": "22 Utility-scale PV tracker",
    "solar-rooftop residential": "22 Rooftop PV residential",
    "solar-rooftop commercial": "22 Rooftop PV commercial",
    "OCGT": "52 OCGT - Natural gas",
    "CCGT": "05 Gas turb. CC, steam extract.",
    "oil": "50 Diesel engine farm",
    "biomass CHP": "09c Straw, Large, 40 degree",
    "biomass EOP": "09c Straw, Large, 40 degree",
    "biomass HOP": "09c Straw HOP",
    "central coal CHP": "01 Coal CHP",
    "central gas CHP": "04 Gas turb. simple cycle, L",
    "central gas CHP CC": "04 Gas turb. simple cycle, L",
    "central solid biomass CHP": "09a Wood Chips, Large 50 degree",
    "central solid biomass CHP CC": "09a Wood Chips, Large 50 degree",
    "central solid biomass CHP powerboost CC": "09a Wood Chips, Large 50 degree",
    "central air-sourced heat pump": "40 Comp. hp, airsource 3 MW",
    "central geothermal-sourced heat pump": "45.1.a Geothermal DH, 1200m, E",
    "central geothermal heat source": "45.1.a Geothermal DH, 1200m, E",
    "central excess-heat-sourced heat pump": "40 Comp. hp, excess heat 10 MW",
    "central water-sourced heat pump": "40 Comp. hp, seawater 20 MW",
    "central ground-sourced heat pump": "40 Absorption heat pump, DH",
    "central resistive heater": "41 Electric Boilers",
    "central gas boiler": "44 Natural Gas DH Only",
    "decentral gas boiler": "202 Natural gas boiler",
    "direct firing gas": "312.a Direct firing Natural Gas",
    "direct firing gas CC": "312.a Direct firing Natural Gas",
    "direct firing solid fuels": "312.b Direct firing Sold Fuels",
    "direct firing solid fuels CC": "312.b Direct firing Sold Fuels",
    "decentral ground-sourced heat pump": "207.7 Ground source existing",
    "decentral air-sourced heat pump": "207.3 Air to water existing",
    "central water pit storage": "140 PTES seasonal",
    "central water tank storage": "141 Large hot water tank",
    "decentral water tank storage": "142 Small scale hot water tank",
    "fuel cell": "12 LT-PEMFC CHP",
    "hydrogen storage underground": "151c Hydrogen Storage - Caverns",
    "hydrogen storage tank type 1 including compressor": "151a Hydrogen Storage - Tanks",
    "micro CHP": "219 LT-PEMFC mCHP - natural gas",
    "biogas": "81 Biogas, Basic plant, small",
    "biogas CC": "81 Biogas, Basic plant, small",
    "biogas upgrading": "82 Upgrading 3,000 Nm3 per h",
    "battery": "180 Lithium Ion Battery",
    "industrial heat pump medium temperature": "302.a High temp. hp Up to 125 C",
    "industrial heat pump high temperature": "302.b High temp. hp Up to 150",
    "electric boiler steam": "310.1 Electric boiler steam  ",
    "gas boiler steam": "311.1c Steam boiler Gas",
    "solid biomass boiler steam": "311.1e Steam boiler Wood",
    "solid biomass boiler steam CC": "311.1e Steam boiler Wood",
    "biomass boiler": "204 Biomass boiler, automatic",
    "electrolysis": "86 AEC 100 MW",
    "direct air capture": "403.a Direct air capture",
    "biomass CHP capture": "401.a Post comb - small CHP",
    "cement capture": "401.c Post comb - Cement kiln",
    "BioSNG": "84 Gasif. CFB, Bio-SNG",
    "BtL": "85 Gasif. Ent. Flow FT, liq fu ",
    "biomass-to-methanol": "97 Methanol from biomass gasif.",
    "biogas plus hydrogen": "99 SNG from methan. of biogas",
    "methanolisation": "98 Methanol from hydrogen",
    "Fischer-Tropsch": "102 Hydrogen to Jet",
    "central hydrogen CHP": "12 LT-PEMFC CHP",
    "Haber-Bosch": "103 Hydrogen to Ammonia",
    "air separation unit": "103 Hydrogen to Ammonia",
    "waste CHP": "08 WtE CHP, Large, 50 degree",
    "waste CHP CC": "08 WtE CHP, Large, 50 degree",
    "biochar pyrolysis": "105 Slow pyrolysis, Straw",
    "electrolysis small": "86 AEC 10 MW",
}


class Dict(dict):
    """
    Dict is a subclass of dict, which allows you to get AND SET items in the
    dict using the attribute syntax!

    Stripped down from addict https://github.com/mewwts/addict/ used in from pypsa.descriptor import Dict.
    """

    def __setattr__(self, name, value):
        """
        Setattr is called when the syntax a.b = 2 is used to set a value.
        """
        if hasattr(Dict, name):
            raise AttributeError(f"'Dict' object attribute '{name}' is read-only")
        self[name] = value

    def __getattr__(self, item):
        try:
            return self.__getitem__(item)
        except KeyError as e:
            raise AttributeError(e.args[0])

    def __delattr__(self, name):
        """
        Is invoked when del some_addict.b is called.
        """
        del self[name]

    _re_pattern = re.compile("[a-zA-Z_][a-zA-Z0-9_]*")

    def __dir__(self):
        """
        Return a list of object attributes.

        This includes key names of any dict entries, filtered to the
        subset of valid attribute names (e.g. alphanumeric strings
        beginning with a letter or underscore).  Also includes
        attributes of parent dict class.
        """
        dict_keys = []
        for k in self.keys():
            if isinstance(k, str):
                m = self._re_pattern.match(k)
                if m:
                    dict_keys.append(m.string)

        obj_attrs = list(dir(Dict))

        return dict_keys + obj_attrs


def mock_snakemake(
    rulename,
    root_dir=None,
    configfiles=None,
    submodule_dir="workflow/submodules/pypsa-eur",
    **wildcards,
):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    root_dir: str/path-like
        path to the root directory of the snakemake project
    configfiles: list, str
        list of configfiles to be used to update the config
    submodule_dir: str, Path
        in case PyPSA-Eur is used as a submodule, submodule_dir is
        the path of pypsa-eur relative to the project directory.
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """

    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    user_in_script_dir = Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in sm.SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = sm.ResourceSettings()
        config_settings = sm.ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = sm.WorkflowSettings()
        storage_settings = sm.StorageSettings()
        dag_settings = sm.DAGSettings(rerun_triggers=[])
        workflow = sm.Workflow(
            config_settings,
            resource_settings,
            workflow_settings,
            storage_settings,
            dag_settings,
            storage_provider_settings=dict(),
        )
        workflow.include(snakefile)

        if configfiles:
            for f in configfiles:
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Config file {f} does not exist.")
                workflow.configfile(f)

        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i, _ in enumerate(io):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake
