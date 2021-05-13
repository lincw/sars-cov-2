# title: Display gene cell line specificity profiles, HuSCI + APMS
# date: 2021.04.26

## environment setup ----
library(openxlsx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(pastecs)
library(tidyverse)
library(RVAideMemoire)
source("~/Documents/INET-work/INETscripts/R_functions.R")

subLocaPtoM <- function(x) {
    dt_raw <- str_split(x$Subcellular.main.location, pattern = ", ")
    dt <- data.frame(subcell = tolower(unlist(dt_raw)))
    dt$subcell[is.na(dt$subcell)] <- "not detected"
    dt2 <- merge(sub_major, dt, by.x = "subcellular", by.y = "subcell")
    dt2_t <- prop.table(table(dt2$major))
    return(as.data.frame(dt2_t))
}
subLocaCtoM <- function(x) {
    dt_raw <- str_split(x$Subcellular.main.location, pattern = ", ")
    dt <- data.frame(subcell = tolower(unlist(dt_raw)))
    dt$subcell[is.na(dt$subcell)] <- "not detected"
    dt2 <- merge(sub_major, dt, by.x = "subcellular", by.y = "subcell")
    dt2_t <- table(dt2$major)
    return(as.data.frame(dt2_t))
}

## 1. loading HPA, HuSCI, Gordon, Li, Stukalov, Nabeel and HuRI data ----
# add BioID
source("~/Documents/INET-work/virus_network/src/cov_4apms.r")
has_sub <- has_all[, c(3, 1, 63, 57, 64)]
sub_major <- read.csv("~/Documents/INET-work/references/HPA/HPA_subcell.csv", header = T)
sub_major$major <- tolower(sub_major$major)
sub_major$subcellular <- tolower(sub_major$subcellular)

bioid1 <- read.csv("~/Documents/INET-work/virus_network/raw/laurent_BioID.csv", header = T) # with identical Gene.name as the human column

# St_Germain & Samavarchi already loaded from "cov_4apms.r"

total <- c(length(unique(has_sub$Ensembl)), length(unique(husci$Ensembl_ID)), length(unique(gordon_science$Ensembl_uniprotIDmapping)), length(unique(stukalov$ensemblID)), length(unique(li$ensemblID)), length(unique(nabeel$Prey_ensembl)), length(unique(bioid1$ensemblID)), length(unique(bioid_st$ensemblID)), length(unique(bioid_sama$ensemblID)))

### count ratio of gene subcellular location
path_final <- "~/documents/INET-work/virus_network/"

subset <- list(HPA = has_sub, # 19670
     HuSCI = unique(merge(husci, has_sub, by.x = "Ensembl_ID", by.y = "Ensembl", all.x = TRUE)[, c(1, 4:7)]), # 171
     Gordon = unique(merge(gordon_science, has_sub, by.x = "Ensembl_uniprotIDmapping", by.y = "Ensembl", all.x = TRUE)[, c(1, 5, 11:13)]), # 384
    Stukalov = unique(merge(stukalov, has_sub, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 3, 7:9)]), # 875
    Li = unique(merge(li, has_sub, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1, 3, 9:11)]), # 285
    Nabeel = unique(merge(nabeel, has_sub, by.x = "Prey_ensembl", by.y = "Ensembl", all.x = TRUE)[, c(1, 4, 6:8)]),  # 277
    Laurent = unique(merge(bioid1[, c(4, 2)], has_sub, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1:2, 4:6)]),
    St_Germain = unique(merge(bioid_st[, c(2, 3)], has_sub, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1:2, 4:6)]),
    Samavarchi = unique(merge(bioid_sama[, c(2, 3)], has_sub, by.x = "ensemblID", by.y = "Ensembl", all.x = TRUE)[, c(1:2, 4:6)])
        )
write.xlsx(subset, file = file.path(path_final, "/toSummarize/subcellular_location/data/subcell_raw.xlsx"))

# display 13 major organelles ----
## proportion calculation

# will not be used for the SARS-CoV-2 paper
# subset_prop <- lapply(subset, subLocaP)
# subset_prop_dt <- subset_prop %>% reduce(left_join, by = "dt") %>% setNames(., c("subcell", names(subset_prop)))
# subset_prop_dt$subcell <- tolower(subset_prop_dt$subcell)
#
# subset_p_melt <- melt(subset_prop_dt)
# subset_p_melt$variable <- factor(subset_p_melt$variable, levels = c("HuSCI", "Gordon", "Stukalov", "Li", "Nabeel", "HPA"))

# pdf(file = file.path(path_final, "Y2H_screening/20201104_final/figures/subcellular_location.pdf"), width = 12, height = 10)
# ggplot(subset_p_melt, aes(x = variable, y = value * 100)) +
#     geom_bar(aes(fill = variable), stat = "identity") +
#     labs(x = "", y = "% of genes", title = "proportion of protein subcellular location") +
#     facet_wrap(~ subcell, scales = "free") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         axis.text = element_text(color = "black"),
#         strip.text = element_text(size = 8),
#         legend.position = "none")
# dev.off()
# write.csv(subset_prop_dt, file = file.path(path_final, "statistic_results/5interactome_subcell_location_prop.csv"), row.names = FALSE)

# subset_count <- lapply(subset, subLocaC)
# subset_count_dt <- subset_count %>% reduce(left_join, by = "dt") %>% setNames(., c("subcell", names(subset_count)))
# write.csv(subset_count_dt, file = file.path(path_final, "statistic_results/5interactome_subcell_location_count.csv"), row.names = FALSE)

subset_major_p <- lapply(subset, subLocaPtoM)
subset_major_pt <- subset_major_p %>% reduce(left_join, by = "Var1") %>% setNames(., c("subcell", names(subset_major_p)))
subset_major_pt[is.na(subset_major_pt)] <- 0
subset_m_melt <- melt(subset_major_pt)
subset_m_melt$variable <- factor(subset_m_melt$variable, levels = c("HuSCI", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi",  "HPA"))

subset_major_c <- lapply(subset, subLocaCtoM)
subset_major_ct <- subset_major_c %>% reduce(left_join, by = "Var1") %>% setNames(., c("subcell", names(subset_major_c)))
subset_major_ct[is.na(subset_major_ct)] <- 0

## fisher test ----
f_pvalue <- list()
count <- 1
for (i in 1:14) {
    for (j in 3:10) {
        dt_tmp <- matrix(as.numeric(c(subset_major_ct[i, c(j, 2)], total[c(j - 1, 1)] - subset_major_ct[i, c(j, 2)])), ncol = 2)
        f_pvalue[[count]] <- fisher.test(dt_tmp)$p.value
        count <- count + 1
    }
}
f_pvalue_df <- data.frame(matrix(unlist(f_pvalue), ncol = 8, byrow = TRUE))
names(f_pvalue_df) <- names(subset_major_ct)[c(3:10)]
f_pvalue_adj <- p.adjust(unlist(f_pvalue), method = "fdr")
f_pvalue_adj_df <- data.frame(matrix(unlist(f_pvalue_adj), ncol = 8, byrow = TRUE))
names(f_pvalue_adj_df) <- names(subset_major_ct)[c(3:10)]

write.csv(cbind(subset_major_ct[c(1, 3:10, 2)], subset_major_pt[c(3:10, 2)], f_pvalue_df, f_pvalue_adj_df), file = file.path(path_final, "statistic_results/interactome_subcell_major_location_all.csv"), row.names = FALSE)

colour <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#FF7F00", "#FFFF00", "brown", "grey")

anno <- data.frame(subcell = rep(subset_major_ct$subcell, each = 8),
    xstar = ifelse(unlist(f_pvalue) < 0.05, rep(1:8, 14), NA),
    ystar = ifelse(unlist(f_pvalue) < 0.05, as.vector(t(subset_major_pt[, c(3:10)])), NA),
    lab = ifelse(f_pvalue_adj < 0.05, "*", NA))

ggplot(subset_m_melt, aes(x = variable, y = value * 100)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    scale_fill_manual(values = colour) +
    geom_text(data = anno, aes(x = xstar, y = ystar * 100, label = lab)) +
    labs(x = "", y = "% of genes", title = "proportion of protein major subcellular location") +
    facet_wrap(~ subcell, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 8),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"))
ggsave(file.path(path_final, "/toSummarize/subcellular_location/figure/subcellular_major_location_all.pdf"), width = 12, height = 10)

## output processed data ----
process <- list(raw_count = subset_major_ct,
                percent = subset_major_pt)
write.xlsx(process, file = file.path(path_final, "/toSummarize/subcellular_location/data/subcell_statistics.xlsx"))
