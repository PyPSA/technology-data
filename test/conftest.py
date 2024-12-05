#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pathlib
import pytest
import shutil
import yaml

_content_temp_file = "content"
_name_temp_file = "hello.txt"
_temp_content_dir = "temp_content_dir"
_sub_temp_content_dir = "sub_temp_content_dir"


@pytest.fixture(scope="function")
def get_temp_file(tmpdir):
    p = tmpdir.join(_name_temp_file)
    p.write(_content_temp_file)
    yield p
    pathlib.Path(p).unlink(missing_ok=True)


@pytest.fixture(scope="function")
def get_temp_folder(tmpdir):
    temp_content_dir = tmpdir.join(_temp_content_dir)
    sub_temp_content_dir = temp_content_dir.join(_sub_temp_content_dir)
    yield sub_temp_content_dir
    shutil.rmtree(str(sub_temp_content_dir))


@pytest.fixture(scope="function")
def config():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def cost_dataframe():
    return pd.DataFrame(
        {
            "technology": ["coal", "coal", "coal", "coal", "coal", "coal", "coal", "coal"],
            "parameter": ["investment", "FOM", "VOM", "fuel", "investment", "discount rate", "co2 intensity", "lifetime"],
            "value": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "unit": ["unit", "unit", "unit", "unit", "unit", "unit", "unit", "unit"],
            "source": ["source", "source", "source", "source", "source", "source", "source", "source"],
            "further description": ["a", "b", "c", "d", "e", "f", "g", "h"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020],
        },
        index=[0, 1, 2, 3, 4, 5, 6, 7],
    )


@pytest.fixture(scope="function")
def atb_cost_dataframe():
    return pd.DataFrame(
        {
            "technology": ["coal", "coal", "coal", "coal", "coal", "coal"],
            "parameter": ["investment", "FOM", "VOM", "fuel", "investment", "discount rate"],
            "value": [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            "unit": ["unit_atb", "unit_atb", "unit_atb", "unit_atb", "unit_atb", "unit_atb"],
            "source": ["source_atb", "source_atb", "source_atb", "source_atb", "source_atb", "source_atb"],
            "further description": ["a", "b", "c", "d", "e", "f"],
            "currency_year": [2020, 2020, 2020, 2020, 2020, 2020],
            "financial_case": ["R&D", "R&D", "R&D", "R&D", "R&D", "R&D"],
            "scenario": ["Moderate", "Moderate", "Moderate", "Moderate", "Moderate", "Moderate"],
            "tax_credit_case": ["ITC", "ITC", "ITC", "ITC", "ITC", "ITC"],
        },
        index=[0, 1, 2, 3, 4, 5],
    )
