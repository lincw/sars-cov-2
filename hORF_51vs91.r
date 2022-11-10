# compare HORFeome 5.1 and 9.1
# 26.09.2022
# Lin Chung-wen

# install ggvenn package if not exist
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
library(openxlsx)

## HORFeome v5.1 downloaded from http://horfdb.dfci.harvard.edu/hv5/
## fasta: http://horfdb.dfci.harvard.edu/hv5/docs/human_orfeome51_fasta.tar.gz
## extract fasta annotation into info.csv
orf51 <- read.csv("~/Documents/INET-work/references/HORFeome/human_orfeome51_info_Ensembl.csv", header = T)
orf91 <- read.table("~/Documents/INET-work/references/HORFeome/HORFeome_all.tab", sep = "\t", header = T)
horf <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/ORF_searchSpace.xlsx", sheet = 1)

# union of orf51 and orf91, based on ENTREZ gene id
uni <- unique(c(orf51$gene_id[!is.na(orf51$gene_id)], orf91$entrez_gene_id[!is.na(orf91$entrez_gene_id)]))
int <- unique(intersect(orf51$gene_id[!is.na(orf51$gene_id)], orf91$entrez_gene_id[!is.na(orf91$entrez_gene_id)]))

overlap <- length(int) / length(uni)

# visualize with venn diagram
orf_list <- list(
    HORFeome_5.1 = orf51$gene_id[!is.na(orf51$gene_id)],
    HORFeome_9.1 = orf91$entrez_gene_id[!is.na(orf91$entrez_gene_id)]
)
ggvenn(orf_list, fill_color = c("blue", "red"), stroke_size = 1, set_name_size = 6, text_size = 6)

## union of orf51 and orf91, based on ORFs with Ensembl gene id
uni_ensembl <- unique(c(orf51[!is.na(orf51$Ensembl_id), "orf_id"], unique(horf$orf_id)))
int_ensembl <- unique(intersect(
    orf51[!is.na(orf51$Ensembl_id), "orf_id"], 
    unique(horf$orf_id)
    ))

overlap_ensembl <- length(int_ensembl) / length(uni_ensembl)

# visualize with venn diagram
orf_list_ensembl <- list(
    HORFeome_5.1 = unique(orf51[!is.na(orf51$Ensembl_id), "orf_id"]),
    HORFeome_9.1 = unique(horf$orf_id)
)
ggvenn(orf_list_ensembl, fill_color = c("blue", "red"), stroke_size = 1, set_name_size = 6, text_size = 6)

## union of orf51 and orf91, based on ORF
uni_orf_level <- unique(c(orf51$orf_id, horf$orf_id))
int_orf_level <- intersect(orf51$orf_id, horf$orf_id)

overlap_orf_level <- length(int_orf_level) / length(uni_orf_level)
