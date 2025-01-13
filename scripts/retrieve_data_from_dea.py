# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8
"""
Created on Mon May  4 18:48:11 2020

script retrieves technology data sheets from the danish energy agency

@author: bw0928
"""

import time
import urllib.request

import requests
from bs4 import BeautifulSoup

# %%
prefix = "https://ens.dk"
url = prefix + "/en/our-services/projections-and-models/technology-data/"
path_out = "inputs/"
docu = False

# %%
response = requests.get(url)
soup = BeautifulSoup(response.text, "html.parser")

# %%
search = "https://ens.dk/en/our-services/projections-and-models/technology-data/"
links = soup.findAll("a", {"href": lambda href: href and search in href})

for i in range(len(links)):
    one_a_tag = links[i]
    link_to_site = one_a_tag["href"]
    response2 = requests.get(link_to_site)
    soup2 = BeautifulSoup(response2.text, "html.parser")
    data = soup2.findAll("a", {"href": lambda href: href and ".xlsx" in href})
    docu = soup2.findAll("a", {"href": lambda href: href and ".pdf" in href})

    # get the data
    for j in range(len(data)):
        link_to_data = data[j]["href"]
        download_url = prefix + link_to_data
        urllib.request.urlretrieve(download_url, path_out + link_to_data.split("/")[-1])
        time.sleep(1)

    # get the documentation
    if docu:
        for j in range(len(docu)):
            link_to_docu = docu[j]["href"]

            if prefix not in link_to_docu:
                download_url = prefix + link_to_docu
            else:
                download_url = link_to_docu

            urllib.request.urlretrieve(
                download_url, path_out + "/docu/" + link_to_docu.split("/")[-1]
            )
            time.sleep(1)
