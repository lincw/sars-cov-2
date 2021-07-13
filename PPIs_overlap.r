# PPI overlap calculation
# Lin Chung-wen

######
# packages
library(openxlsx)
library(igraph)
library(dplyr)
######
# data
husci <- read.csv("/Volumes/GoogleDrive/My Drive/VirHostome_CW/network/data/HuSCI_edge.csv", header = T)
ppis <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx"
gordon <- read.xlsx(ppis, sheet = "Gordon")
stukalov <- read.xlsx(ppis, sheet = "Stukalov")
li <- read.xlsx(ppis, sheet = "Li")
nabeel <- read.xlsx(ppis, sheet = "Nabeel")
laurent <- read.xlsx(ppis, sheet = "Laurent")
stgermain <- read.xlsx(ppis, sheet = "St_Germain")
samavarchi <- read.xlsx(ppis, sheet = "Samavarchi")

husci <- mutate_all(husci, .funs = toupper)
gordon <- mutate_all(gordon[, c(2, 4)], .funs = toupper)
stukalov <- mutate_all(stukalov[, c(1, 2)], .funs = toupper)
li <- mutate_all(li[, c(1, 2)], .funs = toupper)
nabeel <- mutate_all(nabeel[, c(1, 3)], .funs = toupper)
laurent <- mutate_all(laurent[, c(1, 2)], .funs = toupper)
stgermain <- mutate_all(stgermain[, c(1, 2)], .funs = toupper)
samavarchi <- mutate_all(samavarchi[, c(1, 2)], .funs = toupper)

######
# generate network graph
husci_g <- graph_from_data_frame(husci, directed = F)
gordon_g <- graph_from_data_frame(gordon, directed = F)
stukalov_g <- graph_from_data_frame(stukalov, directed = F)
li_g <- graph_from_data_frame(li, directed = F)
nabeel_g <- graph_from_data_frame(nabeel, directed = F)
laurent_g <- graph_from_data_frame(laurent, directed = F)
stgermain_g <- graph_from_data_frame(stgermain, directed = F)
samavarchi_g <- graph_from_data_frame(samavarchi, directed = F)

######
# overlap calculation
husci_gordon <- husci_g %s% gordon_g
husci_stukalov <- husci_g %s% stukalov_g
husci_li <- husci_g %s% li_g
husci_nabeel <- husci_g %s% nabeel_g
husci_laurent <- husci_g %s% laurent_g
husci_stgermain <- husci_g %s% stgermain_g
husci_samavarchi <- husci_g %s% samavarchi_g
gordon_stukalov <- gordon_g %s% stukalov_g
gordon_li <- gordon_g %s% li_g
gordon_nabeel <- gordon_g %s% nabeel_g
gordon_laurent <- gordon_g %s% laurent_g
gordon_stgermain <- gordon_g %s% stgermain_g
gordon_samavarchi <- gordon_g %s% samavarchi_g
stukalov_li <- stukalov_g %s% li_g
stukalov_nabeel <- stukalov_g %s% nabeel_g
stukalov_laurent <- stukalov_g %s% laurent_g
stukalov_stgermain <- stukalov_g %s% stgermain_g
stukalov_samavarchi <- stukalov_g %s% samavarchi_g
li_nabeel <- li_g %s% nabeel_g
li_laurent <- li_g %s% laurent_g
li_stgermain <- li_g %s% stgermain_g
li_samavarchi <- li_g %s% samavarchi_g
nabeel_laurent <- nabeel_g %s% laurent_g
nabeel_stgermain <- nabeel_g %s% stgermain_g
nabeel_samavarchi <- nabeel_g %s% stukalov_g
laurent_stgermain <- laurent_g %s% stgermain_g
laurent_samavarchi <- laurent_g %s% samavarchi_g
stgermain_samavarchi <- stgermain_g %s% samavarchi_g

overlap <- list(
    # HuSCI
    husci_gordon = c(ecount(husci_g) - ecount(husci_gordon), ecount(husci_gordon), ecount(gordon_g) - ecount(husci_gordon)),
    husci_stukalov = c(ecount(husci_g) - ecount(husci_stukalov), ecount(husci_stukalov), ecount(stukalov_g) - ecount(husci_stukalov)),
    husci_li = c(ecount(husci_g) - ecount(husci_li), ecount(husci_li), ecount(li_g) - ecount(husci_li)),
    husci_nabeel = c(ecount(husci_g) - ecount(husci_nabeel), ecount(husci_nabeel), ecount(nabeel_g) - ecount(husci_nabeel)),
    husci_laurent = c(ecount(husci_g) - ecount(husci_laurent), ecount(husci_laurent), ecount(laurent_g) - ecount(husci_laurent)),
    husci_stgermain = c(ecount(husci_g) - ecount(husci_stgermain), ecount(husci_stgermain), ecount(stgermain_g) - ecount(husci_stgermain)),
    husci_samavarchi = c(ecount(husci_g) - ecount(husci_samavarchi), ecount(husci_samavarchi), ecount(samavarchi_g) - ecount(husci_samavarchi)),
    # Gordon
    gordon_stukalov = c(ecount(gordon_g) - ecount(gordon_stukalov), ecount(gordon_stukalov), ecount(stukalov_g) - ecount(gordon_stukalov)),
    gordon_li = c(ecount(gordon_g) - ecount(gordon_li), ecount(gordon_li), ecount(li_g) - ecount(gordon_li)),
    gordon_nabeel = c(ecount(gordon_g) - ecount(gordon_nabeel), ecount(gordon_nabeel), ecount(nabeel_g) - ecount(gordon_nabeel)),
    gordon_laurent = c(ecount(gordon_g) - ecount(gordon_laurent), ecount(gordon_laurent), ecount(laurent_g) - ecount(gordon_laurent)),
    gordon_stgermain = c(ecount(gordon_g) - ecount(gordon_stgermain), ecount(gordon_stgermain), ecount(stgermain_g) - ecount(gordon_stgermain)),
    gordon_samavarchi = c(ecount(gordon_g) - ecount(gordon_samavarchi), ecount(gordon_samavarchi), ecount(samavarchi_g) - ecount(gordon_samavarchi)),
    # Stukalov
    stukalov_li = c(ecount(stukalov_g) - ecount(stukalov_li), ecount(stukalov_li), ecount(li_g) - ecount(stukalov_li)),
    stukalov_nabeel = c(ecount(stukalov_g) - ecount(stukalov_nabeel), ecount(stukalov_nabeel), ecount(nabeel_g) - ecount(stukalov_nabeel)),
    stukalov_laurent = c(ecount(stukalov_g) - ecount(stukalov_laurent), ecount(stukalov_laurent), ecount(laurent_g) - ecount(stukalov_laurent)),
    stukalov_stgermain = c(ecount(stukalov_g) - ecount(stukalov_stgermain), ecount(stukalov_stgermain), ecount(stgermain_g) - ecount(stukalov_stgermain)),
    stukalov_samavarchi = c(ecount(stukalov_g) - ecount(stukalov_samavarchi), ecount(stukalov_samavarchi), ecount(samavarchi_g) - ecount(stukalov_samavarchi)),
    # Li
    li_nabeel = c(ecount(li_g) - ecount(li_nabeel), ecount(li_nabeel), ecount(nabeel_g) - ecount(li_nabeel)),
    li_laurent = c(ecount(li_g) - ecount(li_laurent), ecount(li_laurent), ecount(laurent_g) - ecount(li_laurent)),
    li_stgermain = c(ecount(li_g) - ecount(li_stgermain), ecount(li_stgermain), ecount(stgermain_g) - ecount(li_stgermain)),
    li_samavarchi = c(ecount(li_g) - ecount(li_samavarchi), ecount(li_samavarchi), ecount(samavarchi_g) - ecount(li_samavarchi)),
    # Nabeel
    nabeel_laurent = c(ecount(nabeel_g) - ecount(nabeel_laurent), ecount(nabeel_laurent), ecount(laurent_g) - ecount(nabeel_laurent)),
    nabeel_stgermain = c(ecount(nabeel_g) - ecount(nabeel_stgermain), ecount(nabeel_stgermain), ecount(stgermain_g) - ecount(nabeel_stgermain)),
    nabeel_samavarchi = c(ecount(nabeel_g) - ecount(nabeel_samavarchi), ecount(nabeel_samavarchi), ecount(samavarchi_g) - ecount(nabeel_samavarchi)),
    # Laurent
    laurent_stgermain = c(ecount(laurent_g) - ecount(laurent_stgermain), ecount(laurent_stgermain), ecount(stgermain_g) - ecount(laurent_stgermain)),
    laurent_samavarchi = c(ecount(laurent_g) - ecount(laurent_samavarchi), ecount(laurent_samavarchi), ecount(samavarchi_g) - ecount(laurent_samavarchi)),
    # St-Germain
    stgermain_samavarchi = c(ecount(stgermain_g) - ecount(stgermain_samavarchi), ecount(stgermain_samavarchi), ecount(samavarchi_g) - ecount(stgermain_samavarchi))
)
