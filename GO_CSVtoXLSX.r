# export gProfiler CSV to XLSX
# Lin Chung-wen
# 08.07.2021 **12:32pm**--

library(openxlsx)

gordon <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_noIEA_gProfiler.csv", skip = 17, header = T)
stukalov <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_noIEA_gProfiler.csv", skip = 17, header = T)
li <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_noIEA_gProfiler.csv", skip = 17, header = T)
nabeel <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_noIEA_gProfiler.csv", skip = 17, header = T)
laurent <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_noIEA_gProfiler.csv", skip = 17, header = T)
st_germain <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/st_germain_noIEA_gProfiler.csv", skip = 17, header = T)
samavarchi <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_noIEA_gProfiler.csv", skip = 17, header = T)

go <- list(
    Gordon = gordon, Stukalov = stukalov, Li = li, Nabeel = nabeel, Laurent = laurent, "St-Germain" = st_germain, "Samavarchi-Tehrani" = samavarchi)

write.xlsx(go, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/other_PPI_GO.xlsx")
