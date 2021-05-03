# Evaluate the gene-disease association (GDA) for HuSCI dataset
# date: 03.05.2021
# author: Lin Chung-wen

## environment
library(gprofiler2)
library(disgenet2r)

#!!!!!!!!!!!!!!!!!!!!!!!!
# Not used not
## custom GMT
# gda <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/DisGeNet_curated_GDA.gmt") # generated 03.05.2021
# gda <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/disgenet.curated.v7.symbols.gmt")

## HuSCI
binary <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126_simple.csv", header = T)
binary_node <- unique(binary[binary$category == "human", "node"])

sars2_gda <- disease_enrichment(entities = binary_node, 
                  vocabulary = "HGNC",
                  database = "CTD_human")
