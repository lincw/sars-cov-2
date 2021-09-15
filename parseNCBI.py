#!/usr/bin/env python
# Lin Chung-wen
# parse NCBI with ensembl ID

import requests
from bs4 import BeautifulSoup
import pandas as pd

xl_file = pd.read_excel('/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Supplementary_Table_1_20210914LCW.xlsx')

url = "https://www.ncbi.nlm.nih.gov/gene/?term="

for i in xl_file['Ensembl gene ID']:
	response = requests.get(url + i)
	soup = BeautifulSoup(response.text, "html.parser")
	result = soup.find_all("div", class_ = "rprt full-rprt")
	print(i, "\t", end = "")
for name in result:
	symbol = name.select_one("dd", class_ = "noline").getText()
	print(symbol.split("provided by "), "\t", end = "")

	id = name.select_one("em", class_ = "tax").getText()
	print(id)
