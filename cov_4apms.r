# date: 2021.04.06
# 2021.09.27
# author: Lin Chung-wen

library(openxlsx)

## loading Human Protein Atlas tissue data ----
huri_tissue <- read.table("~/Documents/INET-work/references/HuRI_binaryPPI/Supplementary_Tables_new/Supplementary Table 22.txt", header = T, sep = "\t")
has_all <- read.xlsx("~/workplace/database/Homo_sapiens/20210409proteinatlas.xlsx")
has_tissue <- has_all[, c(3, 1, 17:18, 20, 83:143)]
has_name <- read.csv("~/workplace/database/Homo_sapiens/RNAname_has_protein_data.csv", header = T)
names(has_tissue) <- c("ensemblID", "Gene", "RNA_tissue_specificity", "RNA_tissue_distribution", "RNA_tissue_NX", has_name$tissue[c(1:61)])
has_organ <- has_name$organ

has_dorwie <- read.delim("~/Documents/INET-work/virus_network/raw/human_split.tsv", header = T, sep = "\t", na.strings = c("", " ", "NA"))

## loading HuSCI ----
husci <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/Paper_VirHostome_CoV2/04_Supplementary\ Information/Supplementary_Table_1.xlsx", sheet = "1b - HuSCI", startRow = 4)
husci <- unique(husci[, c("Ensembl.gene.ID", "Host.protein_symbol")])
names(husci) <- c("Ensembl_ID", "node")
husci_rna <- unique(has_tissue[has_tissue$ensemblID %in% husci$Ensembl_ID,])
row.names(husci_rna) <- husci_rna$Gene

## loading 4 CoV AP-MS data ----
gordon_science <- read.csv("~/Documents/INET-work/virus_network/references/PPIs/Gordon_Science/withEnsemblID.csv", header = T)
gordon_rna <- unique(has_tissue[has_tissue$ensemblID %in% gordon_science$Ensembl_uniprotIDmapping,])
row.names(gordon_rna) <- gordon_rna$Gene

stukalov <- read.csv("~/Documents/INET-work/virus_network/raw/stukalov_APMS.csv", header = T)
stukalov_rna <- unique(has_tissue[has_tissue$ensemblID %in% stukalov$ensemblID,])
row.names(stukalov_rna) <- stukalov_rna$Gene

li <- read.csv("~/Documents/INET-work/virus_network/raw/li_APMS.csv", header = T)
li_rna <- unique(has_tissue[has_tissue$ensemblID %in% li$ensemblID,])
row.names(li_rna) <- li_rna$Gene

nabeel <- read.csv("~/Documents/INET-work/virus_network/raw/nabeel_APMS.csv", header = T)
nabeel_rna <- unique(has_tissue[has_tissue$ensemblID %in% nabeel$Prey_ensembl,])
row.names(nabeel_rna) <- nabeel_rna$Gene

## loading BioID (Laurent et al) ----
bioid1 <- read.csv("~/Documents/INET-work/virus_network/raw/laurent_BioID.csv", header = T) # with identical Gene.name as the human column
bioid_rna <- unique(has_tissue[has_tissue$ensemblID %in% bioid1$ensemblID, ])
row.names(bioid_rna) <- bioid_rna$Gene

## loading BioID (St-Germain et al) ----
bioid_st <- read.csv("~/Documents/INET-work/virus_network/raw/St_BioID.csv", header = T)
bioid_st_rna <- unique(has_tissue[has_tissue$ensemblID %in% bioid_st$ensemblID, ])
row.names(bioid_st_rna) <- bioid_st_rna$PreyGene

## loading BioID (St-Germain et al) ----
bioid_sama <- read.csv("~/Documents/INET-work/virus_network/raw/Samavarchi_BioID.csv", header = T)
bioid_sama_rna <- unique(has_tissue[has_tissue$ensemblID %in% bioid_sama$ensemblID, ])
row.names(bioid_sama_rna) <- bioid_sama_rna$PreyGene

## modified tissue specific data ----
raw_path <- ("~/Documents/INET-work/virus_network/raw")
husci_r <- read.delim(file.path(raw_path, "husci_forR.tsv"), header = T, sep = "\t")
gordon_r <- read.delim(file.path(raw_path, "gordon_forR.tsv"), header = T, sep = "\t")
li_r <- read.delim(file.path(raw_path, "li_forR.tsv"), header = T, sep = "\t")
stukalov_r <- read.delim(file.path(raw_path, "stukalov_forR.tsv"), header = T, sep = "\t")
nabeel_r <- read.delim(file.path(raw_path, "nabeel_forR.tsv"), header = T, sep = "\t")
