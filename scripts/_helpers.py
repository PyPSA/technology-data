# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

import logging
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


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


def get_relative_fn(fn):
    if isinstance(fn, str):
        fn = Path(fn).resolve()
    return fn.relative_to(os.path.commonpath([fn, os.getcwd()]))


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
    import os

    import snakemake as sm
    from snakemake.api import Workflow
    from snakemake.common import SNAKEFILE_CHOICES
    from snakemake.script import Snakemake
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

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
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = ResourceSettings()
        config_settings = ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = WorkflowSettings()
        storage_settings = StorageSettings()
        dag_settings = DAGSettings(rerun_triggers=[])
        workflow = Workflow(
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


def prepare_inflation_rate(fn: str, currency_to_use: str = "eur") -> pd.Series:
    """
    The function reads-in annual the inflation rates from Eurostat
    https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/dataflow/ESTAT/prc_hicp_aind/1.0?references=descendants&detail=referencepartial&format=sdmx_2.1_generic&compressed=true

    Parameters
    ----------
    fn: str
        file name for the Eurostat inflation rates
    currency_to_use: str
        currency to select for the inflation rate

    Returns
    -------
    pandas.Series
        inflation rates series
    """

    if currency_to_use.casefold() == "usd":
        row_to_use = "United States"
    else:
        row_to_use = "European Union - 27 countries (from 2020)"

    inflation_rate_series = pd.read_excel(
        fn, sheet_name="Sheet 1", index_col=0, header=[8], engine="calamine"
    )
    inflation_rate_series = (inflation_rate_series.loc[row_to_use].dropna()).loc[
        "2001"::
    ]
    inflation_rate_series.rename(
        index=lambda inflation_rate_val: int(inflation_rate_val), inplace=True
    )
    inflation_rate_series = inflation_rate_series.astype(float)
    inflation_rate_series /= 100.0
    return inflation_rate_series


def adjust_for_inflation(
    inflation_rate: pd.Series,
    costs: pd.DataFrame,
    techs: pd.Series,
    eur_year: int,
    col_name: str,
    usa_costs_flag: bool = False,
) -> pd.DataFrame:
    """
    The function adjust the investment costs for the specified techs for inflation.

    Parameters
    ----------
    inflation_rate : pandas.Series
        inflation rates for several years
    costs : pd.DataFrame
        existing cost dataframe
    techs : pd.Series
        technologies
    eur_year : int,
        reference year for which the costs are provided and based on which the inflation adjustment is done
    col_name : str
        column name to which to apply the inflation rate adjustment
    usa_costs_flag: bool
        flag for US specific costs

    Returns
    -------
    pandas.Dataframe
        inflation updated cost dataframe
    """

    def get_factor(inflation_rate_df, ref_year, eur_year_val):
        if (pd.isna(ref_year)) or (ref_year < 1900):
            return np.nan
        if ref_year == eur_year_val:
            return 1
        mean = inflation_rate_df.mean()
        if ref_year < eur_year_val:
            new_index = np.arange(ref_year + 1, eur_year_val + 1)
            df = 1 + inflation_rate_df.reindex(new_index).fillna(mean)
            return df.cumprod().loc[eur_year_val]
        else:
            new_index = np.arange(eur_year_val + 1, ref_year + 1)
            df = 1 + inflation_rate_df.reindex(new_index).fillna(mean)
            return 1 / df.cumprod().loc[ref_year]

    inflation = costs.currency_year.apply(
        lambda x: get_factor(inflation_rate, x, eur_year)
    )

    paras = ["investment", "VOM", "fuel"]

    if usa_costs_flag:
        filter_i = costs.technology.isin(techs) & costs.parameter.isin(paras)
    else:
        filter_i = costs.index.get_level_values(0).isin(
            techs
        ) & costs.index.get_level_values(1).isin(paras)
    costs.loc[filter_i, col_name] = costs.loc[filter_i, col_name].mul(
        inflation.loc[filter_i], axis=0
    )

    return costs


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)

    # Setup a function to handle uncaught exceptions and include them with their stacktrace into logfiles
    def handle_exception(exc_type, exc_value, exc_traceback):
        # Log the exception
        logger = logging.getLogger()
        logger.error(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )

    sys.excepthook = handle_exception
