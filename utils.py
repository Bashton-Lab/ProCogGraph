#!/usr/bin/env python

import re
import numpy as np
import pandas as pd
import requests
from urllib.parse import quote
from bs4 import BeautifulSoup
import json 

#this is a core function - should be stored elsewhere.
#a limitation of this function is that it selects only one representative EC number for transferred entries
def get_terminal_record(entry, row, df):
    entry = row.ID
    pattern = r"[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+"
    while row.DE.startswith("Transferred entry: "):
        transfers = re.findall(pattern, row.DE)
        #when updating, if multiple possible transfers, selecting only first.
        row = df.loc[df.ID == transfers[0]].iloc[0]
    return row.ID

def get_csdb_from_glycoct(glycoct):
    if glycoct is np.nan:
        return np.nan
    else:
        url = "http://csdb.glycoscience.ru/database/core/convert_api.php"
        data = {"glycoct":glycoct}
        headers = {'Content-Type': 'application/json'}
        response = requests.get(f"{url}?glycoct={quote(glycoct)}")
        if response.status_code == 200:
            response_data = response.text.replace("<pre>", "")  # Remove the <pre> tag
            lines = response_data.split("\n")  # Split into lines
            csdb_linear = np.nan
            for line in lines:
                if line.startswith("CSDB Linear:"):
                    csdb_linear = line.replace("CSDB Linear:", "").strip()  # Extract the CSDB Linear string
                    break
        else:
            csdb_linear = np.nan
        return csdb_linear
    
def get_glycoct_from_wurcs(wurcs):
    url = "https://api.glycosmos.org/glycanformatconverter/2.8.2/wurcs2glycoct"
    data = {"input":wurcs}
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, headers=headers, data=json.dumps(data))
    
    if response.status_code == 200:
        response_data = response.json()
        if 'message' in response_data and response_data['message'] == 'Returned null.':
            glycoct = np.nan
        else:
            glycoct = response_data['GlycoCT']
    else:
        glycoct = np.nan
    return glycoct

def get_smiles_from_csdb(csdb_linear):
    if csdb_linear is np.nan:
        return np.nan
    else:
        response = requests.get(f"http://csdb.glycoscience.ru/database/core/convert_api.php?csdb={quote(csdb_linear)}&format=smiles")
        mol = np.nan
        smiles = np.nan
        if response.status_code == 200:
            html = response.text
            soup = BeautifulSoup(html, 'html.parser')
            for a in soup.find_all("a"):
                title = a.get('title')
                if title == "Find this structure in ChemSpider":
                    smiles = a.contents[0].strip()
                    break
        else:
            smiles = np.nan   
        return smiles
