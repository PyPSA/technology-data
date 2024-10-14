#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib
import shutil

import pytest
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
def get_config_dict():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict
