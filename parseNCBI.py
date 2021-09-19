#!/usr/bin/env python
# Lin Chung-wen
# parse NCBI with ensembl ID

import requests
from bs4 import BeautifulSoup
import pandas as pd
import sys

sys.stdout = open('/tmp/HuSCI_symbol.txt', 'w')

xl_file = pd.read_excel('/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Supplementary_Table_1_20210914LCW.xlsx')

url = "https://www.ncbi.nlm.nih.gov/gene/?term="

for i in xl_file['Ensembl gene ID']:
	print(i, "\t", end = "")
	response = requests.get(url + i)
	soup = BeautifulSoup(response.text, "html.parser")
	primary = soup.find_all("div", class_ = "rprt-header")
	for gene in primary:
		gene_id = gene.find("span", class_ = "geneid").getText()
		print(gene_id.split(", updated on")[0], "\t", end = "")
	result = soup.find_all("div", class_ = "rprt full-rprt")
	for name in result:
		symbol = name.select_one("dd", class_ = "noline").getText()
		print(symbol.split("provided by ")[0])

sys.stdout.close()