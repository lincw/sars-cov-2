# GO enrichment analysis for HuRI community, which obtained from OCG algorithm
# Lin Chung-wen
# 29.09.2021

######
# environment
huri_ocg <- readRDS("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_ocg.RDS")
database <- "~/workplace/database"

######
# packages
library(gprofiler2)
library(linkcomm)

######
# functions
source("~/Documents/INET-work/virus_network/src/statCal.r")

funcEnrich <- function(ocg, organism, querylist) {
    go <- list()
    df <- data.frame("name" = names(ocg$clustsizes), "clustsize" = ocg$clustsizes)
    for (i in querylist) {
        query_list <- ocg[[4]][ocg[[4]][, 2] == i, "node"]
        goquery <- gost(query = query_list, organism = organism, correction_method = "bonferroni", evcodes = TRUE)
        goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
        goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
        go[[i]] <- goquery$result
    }
    return(go)
}

goOutput <- function(check_list, go_result) { # get enriched community
        dfa1 <- list()
        for (i in check_list) {
                if (length(go_result[[i]]) != 19) {
                        next
                } else {
                        dfa1[i] <- go_result[i]
                }
        }
        go_out <- do.call(rbind, dfa1)
        return (go_out)
}

######
# GO data
bp_cust <- upload_GMT_file(gmtfile = file.path(database, "hsapien_HuRI_GOBP_EXP.gmt"))
mf_cust <- upload_GMT_file(gmtfile = file.path(database, "hsapien_HuRI_GOMF_EXP.gmt"))
cc_cust <- upload_GMT_file(gmtfile = file.path(database, "hsapien_HuRI_GOCC_EXP.gmt"))

######
# GO enrichment analysis
check_list <- names(huri_ocg$clustsizes[huri_ocg$clustsizes >= 4]) # based on community size only
huri_ocg_gobp <- funcEnrich(huri_ocg, bp_cust, check_list)
huri_ocg_gomf <- funcEnrich(huri_ocg, mf_cust, check_list)
huri_ocg_gocc <- funcEnrich(huri_ocg, cc_cust, check_list)

for (i in 1:length(huri_ocg_gobp)) {
	huri_ocg_gobp[[i]]$community <- names(huri_ocg_gobp[i])
	}
for (i in 1:length(huri_ocg_gomf)) {
	huri_ocg_gomf[[i]]$community <- names(huri_ocg_gomf[i])
	}
for (i in 1:length(huri_ocg_gocc)) {
	huri_ocg_gocc[[i]]$community <- names(huri_ocg_gocc[i])
	}

gobp_out <- goOutput(check_list, huri_ocg_gobp)
gobp_out$source <- "BP"
gomf_out <- goOutput(check_list, huri_ocg_gomf)
gomf_out$source <- "MF"
gocc_out <- goOutput(check_list, huri_ocg_gocc)
gocc_out$source <- "CC"

######
# save result
write.csv(rbind(gobp_out, gomf_out, gocc_out)[, c(19, 3:13, 15:18)], file = "~/Documents/INET-work/virus_network/statistic_results/community/HuRI_GO_Exp_all.csv", row.names = FALSE)
